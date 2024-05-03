"""Deconvolve CESM2 output data."""

import warnings

import fppanalysis
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import vdd.load
import vdd.utils

warnings.warn(
    "This script is deprecated and may be removed in the future.",
    category=DeprecationWarning,
    stacklevel=1,
)


def check_cesm_output() -> None:
    """Check the CESM2 output data."""
    cesm = vdd.load.CESMData()
    x = cesm.aod.time.data
    if any(x != cesm.rf.time.data) or any(x != cesm.temp.time.data):
        raise ValueError("Times do not match.")
    plt.figure()
    vdd.utils.normalise(cesm.aod).plot()
    vdd.utils.normalise(cesm.rf).plot()
    vdd.utils.normalise(cesm.temp).plot()
    # plt.plot(x, vdd.utils.normalise(cesm.aod))
    # plt.plot(x, vdd.utils.normalise(cesm.rf))
    # plt.plot(x, vdd.utils.normalise(cesm.temp))
    plt.show()


def extend_aod_shape() -> xr.DataArray:
    """Extend the shape of the aerosol optical depth array."""
    cesm = vdd.load.CESMData()
    count = 1
    x_base = np.concatenate((np.array([0]), cesm.aod.time.data))
    assert len(x_base) == int(12 * 20)
    x = x_base + 20 * count
    # Let's 200 years as the maximum time
    end_time = 200
    while x[-1] < end_time:
        count += 1
        x = np.append(x, x_base + 20 * count)
    decay = np.random.default_rng().normal(0, 0.00003, len(x))
    signal: np.ndarray = np.concatenate((cesm.aod.data, decay))
    aod_new = xr.DataArray(
        signal, coords={"time": np.concatenate((cesm.aod.time.data, x))}, dims=["time"]
    )
    return aod_new


def extend_rf_shape() -> xr.DataArray:
    """Extend the shape of the radiative forcing array."""
    cesm = vdd.load.CESMData()
    cesm.rf *= -1
    count = 1
    x_base = np.concatenate((np.array([0]), cesm.rf.time.data))
    assert len(x_base) == int(12 * 20)
    x = x_base + 20 * count
    end_time = 200
    while x[-1] < end_time:
        count += 1
        x = np.append(x, x_base + 20 * count)
    decay = np.random.default_rng().normal(0, 0.5, len(x))
    signal: np.ndarray = np.concatenate((cesm.rf.data, decay))
    rf_new = xr.DataArray(
        -1 * signal,
        coords={"time": np.concatenate((cesm.rf.time.data, x))},
        dims=["time"],
    )
    return rf_new


def extend_temp_shape() -> xr.DataArray:
    """Extend the shape of the temperature array."""
    cesm = vdd.load.CESMData()
    cesm.temp *= -1
    count = 1
    x_base = np.concatenate((np.array([0]), cesm.temp.time.data))
    assert len(x_base) == int(12 * 20)
    x = x_base + 20 * count
    end_time = 200
    while x[-1] < end_time:
        count += 1
        x = np.append(x, x_base + 20 * count)
    decay = 1.5 * np.exp(-x / 30)
    decay += np.random.default_rng().normal(0, 0.05, len(x))
    signal = np.concatenate((cesm.temp, decay))
    t_new = xr.DataArray(
        -1 * signal,
        coords={"time": np.concatenate((cesm.temp.time.data, x))},
        dims=["time"],
    )
    return t_new


def main():
    """Run the main function."""
    cesm = vdd.load.CESMData()
    # aod = extend_aod_shape()
    # rf = extend_rf_shape()
    # t = extend_temp_shape()
    x = cesm.aod.time.data
    print(len(x))
    if any(x != cesm.rf.time.data) or any(x != cesm.temp.time.data):
        raise ValueError("Times do not match.")
    # cesm.temp = pad_before_convolution(cesm.temp)
    # cesm.rf = pad_before_convolution(cesm.rf)
    # cesm.aod = pad_before_convolution(cesm.aod)
    # Make all arrays positive
    rf = cesm.rf * -1
    temp = cesm.temp * -1
    cesm.aod.plot()
    rf.plot()
    temp.plot()
    plt.show()
    t_rl_rf, _err = fppanalysis.RL_gauss_deconvolve(temp, rf, 200)
    t_rl_rf = t_rl_rf.flatten()
    t_rl_aod, _err = fppanalysis.RL_gauss_deconvolve(temp, cesm.aod, 200)
    t_rl_aod = t_rl_aod.flatten()
    # TODO: convolve the found t_rl_rf with RF from OB16, and compare with T from OB16
    plt.figure()
    plt.plot(cesm.aod.time, t_rl_rf, label="RLD(T, RF)")
    plt.legend()
    plt.figure()
    plt.plot(cesm.aod.time, t_rl_aod, label="RLD(T, AOD)")
    plt.legend()
    plt.figure()
    plt.plot(cesm.aod.time, temp)
    plt.plot(cesm.aod.time, np.convolve(t_rl_aod, cesm.aod, mode="same"))
    plt.plot(cesm.aod.time, np.convolve(t_rl_rf, rf, mode="same"))
    plt.plot(cesm.aod.time, temp)

    plt.show()


if __name__ == "__main__":
    main()
    # extend_temp_shape()
    # extend_rf_shape()
    # extend_aod_shape()
