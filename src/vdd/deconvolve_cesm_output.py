"""Deconvolve CESM2 output data."""

import fppanalysis
import matplotlib.pyplot as plt
import numpy as np
import volcano_base
import xarray as xr

import vdd.get_cesm_data


def pad_before_convolution(
    arr: xr.DataArray, diff_before_zero: float = 0.08493151
) -> xr.DataArray:
    """Pad the array on the left as preparation for a convolution.

    This uses the ``xr.DataArray().pad`` method under the hood, and tries to fix the
    time axis by inserting valid time values. The default behaviour of
    ``xr.DataArray().pad`` is to set time stamps of ``np.nan`` for the padding values,
    which will hide the padded values in a plot.

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


def normalise(arr: xr.DataArray) -> xr.DataArray:
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


def extend_aod_shape() -> xr.DataArray:
    """Extend the shape of the aerosol optical depth array."""
    aod = vdd.get_cesm_data.get_aod_cesm()
    count = 1
    x_base = np.concatenate((np.array([0]), aod.time.data))
    assert len(x_base) == int(12 * 20)
    x = x_base + 20 * count
    while x[-1] < 200:
        count += 1
        x = np.append(x, x_base + 20 * count)
    decay = np.random.default_rng().normal(0, 0.00003, len(x))
    signal = np.concatenate((aod, decay))
    aod_new = xr.DataArray(
        signal, coords={"time": np.concatenate((aod.time.data, x))}, dims=["time"]
    )
    return aod_new


def extend_rf_shape() -> xr.DataArray:
    """Extend the shape of the radiative forcing array."""
    rf = vdd.get_cesm_data.get_rf_cesm()
    rf *= -1
    count = 1
    x_base = np.concatenate((np.array([0]), rf.time.data))
    assert len(x_base) == int(12 * 20)
    x = x_base + 20 * count
    while x[-1] < 200:
        count += 1
        x = np.append(x, x_base + 20 * count)
    decay = np.random.default_rng().normal(0, 0.5, len(x))
    signal = np.concatenate((rf, decay))
    rf_new = xr.DataArray(
        -1 * signal, coords={"time": np.concatenate((rf.time.data, x))}, dims=["time"]
    )
    return rf_new


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
    decay = decay + np.random.default_rng().normal(0, 0.05, len(x))
    signal = np.concatenate((t, decay))
    t_new = xr.DataArray(
        -1 * signal, coords={"time": np.concatenate((t.time.data, x))}, dims=["time"]
    )
    return t_new


def main():
    """Run the main function."""
    aod = vdd.get_cesm_data.get_aod_cesm()
    rf = vdd.get_cesm_data.get_rf_cesm()
    t = vdd.get_cesm_data.get_trefht_cesm()
    # aod = extend_aod_shape()
    # rf = extend_rf_shape()
    # t = extend_temp_shape()
    x = aod.time.data
    print(len(x))
    if any(x != rf.time.data) or any(x != t.time.data):
        raise ValueError("Times do not match.")
    t_pad = pad_before_convolution(t)
    rf_pad = pad_before_convolution(rf)
    aod_pad = pad_before_convolution(aod)
    # Make all arrays positive
    rf_pad *= -1
    t_pad *= -1
    aod_pad.plot()
    rf_pad.plot()
    t_pad.plot()
    plt.show()
    t_rl_rf, err = fppanalysis.RL_gauss_deconvolve(t_pad, rf_pad, 200)
    t_rl_rf = t_rl_rf.flatten()
    t_rl_aod, err = fppanalysis.RL_gauss_deconvolve(t_pad, aod_pad, 200)
    t_rl_aod = t_rl_aod.flatten()
    # TODO: convolve the found t_rl_rf with RF from OB16, and compare with T from OB16
    plt.figure()
    plt.plot(t_rl_rf, label="RLD(T, RF)")
    plt.legend()
    plt.figure()
    plt.plot(t_rl_aod, label="RLD(T, AOD)")
    plt.legend()
    plt.figure()
    plt.plot(t_pad)
    plt.plot(np.convolve(t_rl_aod, aod_pad, mode="same"))
    plt.plot(np.convolve(t_rl_rf, rf_pad, mode="same"))
    plt.plot(t_pad)

    plt.show()


if __name__ == "__main__":
    main()
    # extend_temp_shape()
    # extend_rf_shape()
    # extend_aod_shape()
