"""Module for deconvolving data from Otto-Bliesner et at. (2016)."""

import datetime
from typing import overload

import cftime
import matplotlib.pyplot as plt
import nc_time_axis  # noqa: F401
import numpy as np
import volcano_base
import xarray as xr

import vdd.get_cesm_data

plt.style.use(
    "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle"
)


@overload
def normalise(arr: np.ndarray) -> np.ndarray:
    ...


@overload
def normalise(arr: xr.DataArray) -> xr.DataArray:
    ...


def normalise(arr: xr.DataArray | np.ndarray) -> xr.DataArray | np.ndarray:
    """Normalise the array to the maximum or minimum value."""
    if abs(np.nanmax(arr)) > abs(np.nanmin(arr)):
        return arr / np.nanmax(arr)
    return arr / np.nanmin(arr)


def check_cesm_output() -> None:
    """Check the CESM2 output data."""
    aod = vdd.get_cesm_data.get_aod_cesm(False)
    rf = vdd.get_cesm_data.get_rf_cesm(False)
    t = vdd.get_cesm_data.get_trefht_cesm(False)
    x = aod.time.data
    if any(x != rf.time.data) or any(x != t.time.data):
        raise ValueError("Times do not match.")
    plt.figure()
    normalise(aod).plot()
    normalise(rf).plot()
    normalise(t).plot()
    # plt.plot(x, normalise(aod))
    # plt.plot(x, normalise(rf))
    # plt.plot(x, normalise(t))
    plt.show()


def find_sampling_rate(dates: np.ndarray) -> float:
    """Find the sampling rate.

    Parameters
    ----------
    dates : np.ndarray
        An array of dates. The dates can be in the following formats:
        - float
        - datetime.datetime
        - cftime._cftime.DatetimeNoLeap

    Returns
    -------
    float
        The sampling rate.

    Raises
    ------
    TypeError
        If the input date format is unknown.
    """
    match dates[0]:
        case float():
            return np.mean(np.diff(dates[3:-3]))
        case datetime.datetime():
            return np.mean(
                np.diff(
                    [
                        eval(datetime.datetime.strftime(i, "%Y+%-m/12+%-d/365"))
                        for i in dates[3:-3]
                    ]
                )
            )
        case cftime._cftime.DatetimeNoLeap():
            cftime2float_ = volcano_base.manipulate.dt2float(dates)
            return np.mean(np.diff(cftime2float_[3:-3]))
        case _:
            raise TypeError(
                f"The input date format is unknown, found {type(dates[0]) = }"
            )


def weighted_monthly_avg(da: xr.DataArray) -> xr.DataArray:
    """Calculate a temporal mean, weighted by days in each month.

    Parameters
    ----------
    da : xr.DataArray
        Input data structure to do temporal average on

    Returns
    -------
    xr.DataArray
        The new data structure with time averaged data

    Notes
    -----
    From
    https://ncar.github.io/esds/posts/2021/yearly-averages-xarray/#wrap-it-up-into-a-function
    """
    # Determine the month length
    month_length = da.time.dt.days_in_month
    # Calculate the weights
    wgts = month_length.groupby("time.month") / month_length.groupby("time.month").sum()
    # Make sure the weights in each year add up to 1
    np.testing.assert_allclose(wgts.groupby("time.month").sum(xr.ALL_DIMS), np.ones(12))
    # Setup our masking for nan values
    cond = da.isnull()
    ones = xr.where(cond, 0.0, 1.0)
    # Calculate the numerator
    obs_sum = (da * wgts).resample(time="MS").sum(dim="time")
    # Calculate the denominator
    ones_out = (ones * wgts).resample(time="MS").sum(dim="time")
    # Return the weighted average
    return obs_sum / ones_out


def pad_before_convolution(
    arr: xr.DataArray, diff_before_zero: float = 0.08493151
) -> xr.DataArray:
    """Pad the array on the left as preparation for a convolution.

    This uses the ``xr.DataArray().pad`` method under the hood, and tries to fix the
    time axis by inserting valid time values. The default behaviour is to set ``np.nan``
    for the padding values, which will hide the padded values in a plot.

    Parameters
    ----------
    arr : xr.DataArray
        The input array to pad.
    diff_before_zero : float, optional
        The time difference between time zero and the time step prior. Since time data
        may be monthly with varying number of days, a constant time step is not always
        wanted. The default is 0.08493151, used for monthly data.

    Returns
    -------
    xr.DataArray
        The padded array.
    """
    n = len(arr.time) - 1
    arr_ = arr.pad(time=(n, 0), mode="constant", constant_values=0)
    match arr.time.data[0]:
        case float():
            time_ = arr.time.data - arr.time.data[0]
        case _:
            time_ = volcano_base.manipulate.dt2float(arr.time.data)
            time_ -= time_[0]
    assert time_[0] == 0
    assert len(arr_) % 2 != 0
    # Assign new time values to the padded array that is twice the length of the
    # original, and that goes from -n*dt to n*dt.
    arr_ = arr_.assign_coords(
        time=np.concatenate((time_[1:] - time_[-1] - diff_before_zero, time_)),
    )
    return arr_


def extend_temp_shape() -> xr.DataArray:
    """Extend the shape of the temperature array."""
    t = vdd.get_cesm_data.get_trefht_cesm()
    t *= -1
    count = 1
    x_base = np.concatenate((np.array([0]), t.time.data))
    assert len(x_base) == int(12 * 20)
    x = x_base + 20 * count
    while x[-1] < 200:
        count += 1
        x = np.append(x, x_base + 20 * count)
    decay = 1.5 * np.exp(-x / 30)
    decay = decay + np.random.default_rng().normal(0, 0.1, len(x))
    signal = np.concatenate((t, decay))
    t_new = xr.DataArray(
        -1 * signal, coords={"time": np.concatenate((t.time.data, x))}, dims=["time"]
    )
    return t_new


def main() -> None:
    """Deconvolve data."""
    # volcano_base.load.get_ob16_outputs()
    temperature_day = volcano_base.load.get_ob16_temperature()
    temperature = weighted_monthly_avg(temperature_day)
    rf_day = volcano_base.load.get_ob16_rf()
    rf = weighted_monthly_avg(rf_day)
    so2_day = volcano_base.load.get_so2_ob16_peak_timeseries(xarray=True) / 3 * 2
    d1 = 190
    so2_day = so2_day.assign_coords(
        time=so2_day.time.data - datetime.timedelta(days=d1)
    )
    so2 = weighted_monthly_avg(so2_day)

    so2, rf_fr, temp = xr.align(so2, rf, temperature)
    # t_so2, err_t_so2 = fppanalysis.deconvolution_methods.RL_gauss_deconvolve(
    #     temperature, so2, 200
    # )
    # t_rf, err_t_rf = fppanalysis.deconvolution_methods.RL_gauss_deconvolve(
    #     temperature, rf, 200
    # )
    # print(t_so2)
    # plt.plot(t_so2)
    # plt.plot(t_rf)
    # plt.show()

    frc = pad_before_convolution(vdd.get_cesm_data.get_rf_cesm())
    tmp = pad_before_convolution(vdd.get_cesm_data.get_trefht_cesm())
    tmp_ext = pad_before_convolution(extend_temp_shape())
    signal = np.convolve(so2, frc, mode="same")
    signal2 = np.convolve(so2, tmp, mode="same")
    signal3 = np.convolve(so2, tmp_ext, mode="same")
    plt.figure()
    normalise(temp).plot()
    plt.plot(so2.time.data, normalise(signal2))
    plt.plot(so2.time.data, normalise(signal3))
    plt.figure()
    normalise(rf_fr).plot()
    plt.plot(so2.time.data, normalise(signal))
    plt.show()


if __name__ == "__main__":
    main()
    # extend_temp_shape()
