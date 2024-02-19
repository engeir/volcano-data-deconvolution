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


def main():
    """Run the main function."""
    aod = vdd.get_cesm_data.get_aod_cesm(False)
    rf = vdd.get_cesm_data.get_rf_cesm(False)
    t = vdd.get_cesm_data.get_trefht_cesm(False)
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
    plt.plot(t_rl_rf)
    plt.figure()
    plt.plot(t_rl_aod)
    plt.figure()
    plt.plot(t_pad)
    plt.plot(np.convolve(t_rl_aod, aod_pad, mode="same"))
    plt.plot(np.convolve(t_rl_rf, rf_pad, mode="same"))
    plt.plot(t_pad)

    plt.show()


if __name__ == "__main__":
    main()
