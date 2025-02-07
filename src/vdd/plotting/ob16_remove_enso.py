"""Deconvolve OB16 data, but without ENSO."""

import datetime
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
import volcano_base
import xarray as xr

import vdd.load
import vdd.utils

plt.style.use({"legend.framealpha": 0.8})
_SAVE_DIR = volcano_base.config.SAVE_PATH / "review"
if not _SAVE_DIR.exists():
    _SAVE_DIR.mkdir(parents=False)


_COLORS = plt.rcParams["axes.prop_cycle"].by_key()["color"]


def _combine_times(first: xr.DataArray, second: xr.DataArray) -> xr.DataArray:
    return xr.concat([first, second], dim="time")


def _remove_quadratic_fit(da: xr.DataArray) -> xr.DataArray:
    x, y = volcano_base.manipulate.dt2float(da.time.data), da.data
    out = np.polynomial.polynomial.polyfit(x, y, deg=4)
    fit = out[0] + out[1] * x + out[2] * x**2 + out[3] * x**3 + out[4] * x**4
    da.data -= fit
    return da


def _remove_all_freqs(da: xr.DataArray) -> xr.DataArray:
    plot_freq = False
    base_freq = 0.9811
    vb_remove_season = volcano_base.manipulate.remove_seasonality
    da = vb_remove_season(da, freq=base_freq, plot=plot_freq)
    da = vb_remove_season(da, freq=2 * base_freq, plot=plot_freq)
    da = vb_remove_season(da, freq=3 * base_freq, plot=plot_freq)
    return vb_remove_season(da, freq=4 * base_freq, plot=plot_freq)


def _select_enso_area(da: xr.DataArray) -> xr.DataArray:
    return da.sel(lat=slice(-5, 5), lon=slice(190, 240))


def _regress_a_on_b(a: xr.DataArray, b: xr.DataArray) -> xr.DataArray:
    a_up = a.copy()
    x = sm.add_constant(b.data)  # Adds a constant term to the predictor
    model = sm.OLS(a.data, x).fit()
    predictions = model.predict(x)
    a_up.data -= predictions
    return a_up


def _compute_enso() -> xr.DataArray:
    regex = r"""# Match everything up to "cesm-lme"
        .*cesm-lme/.*g16\.
        # Match for simulation type (group 1)
        (.*)\.00
        # Match for ensemble (group 2)
        (\d)\.cam\.
        # Match for frequency (group 3)
        (..)\.
        # Match for attribute (group 4)
        (TS)\.
        # Match for time segment (group 5)
        (\d{6}-\d{6})"""
    group_names = {
        "Simulation type": "sim",
        "Ensembles": "ensemble",
        "Frequency": "freq",
        "Attribute": "attr",
        "Times": "times",
    }
    reverse_search = (
        pathlib.Path("cesm-lme")
        # b.e11.BLMTRC5CN.f19_g16.VOLC_GRA.001.cam.h0.TS.085001-184912.nc
        / "b.e11.BLMTRC5CN.f19_g16.<sim>.00<ensemble>.cam.<freq>.<attr>.<times><ft>"
    )
    ob16_regex = volcano_base.load.RegexLookup(
        ".nc", group_names, reverse_search, regex
    )
    files = volcano_base.load.FindFiles(regex=ob16_regex)
    # files.avail()
    ctrl: tuple[xr.DataArray, xr.DataArray] = tuple(
        files.copy().find("850forcing").sort("ensemble", "times").load()
    )
    ctrl_arr = _combine_times(*ctrl)
    enso = _select_enso_area(ctrl_arr)
    enso = volcano_base.manipulate.mean_flatten(enso, dims=["lat", "lon"])
    enso = _remove_quadratic_fit(enso)
    return _remove_all_freqs(enso)


def _main() -> None:
    enso = _compute_enso()
    dec_ob16 = vdd.load.DeconvolveOB16("h0")
    temp_before = dec_ob16.temp
    plt.plot(dec_ob16.tau, dec_before := dec_ob16.response_temp_so2)
    # Verify that this is equal to the default
    signal, err = dec_ob16.deconvolve(dec_ob16.temp.data, dec_ob16.so2.data)
    signal = signal.flatten()
    assert np.array_equal(signal, dec_before)  # noqa: S101
    enso = enso.assign_coords(
        time=xr.cftime_range(
            start="0850-01-01", periods=len(enso.data), calendar="noleap", freq="MS"
        )
        + datetime.timedelta(days=15)
    )
    enso, temp_before = xr.align(enso, temp_before)
    enso_response_before, _ = dec_ob16.deconvolve(-enso.data, dec_ob16.so2.data)
    assert enso.shape == temp_before.shape  # noqa: S101
    enso, temp_after = (
        _regress_a_on_b(enso, temp_before),
        _regress_a_on_b(temp_before, enso),
    )
    enso_response_after, _ = dec_ob16.deconvolve(-enso.data, dec_ob16.so2.data)
    dec_after, err = dec_ob16.deconvolve(temp_after.data, dec_ob16.so2.data)
    dec_after = dec_after.flatten()
    plt.plot(dec_ob16.tau, dec_after)
    print(f"Temperatures are equal: {np.array_equal(temp_after, temp_before)}")
    print(
        f"Responses are equal: {np.array_equal(dec_after, dec_before)}. At most "
        f"{abs(dec_after - dec_before).max()}, probably because the OB16 time series is"
        " so long that ENSO is efficiently averaged out by default."
    )
    plt.figure()
    plt.plot(dec_ob16.tau, enso_response_before)
    plt.plot(dec_ob16.tau, enso_response_after)
    plt.xlim((-2, 7))
    # plt.savefig("ob16-enso-response-to-eruptions_negative.pdf")
    plt.show()


if __name__ == "__main__":
    _main()
