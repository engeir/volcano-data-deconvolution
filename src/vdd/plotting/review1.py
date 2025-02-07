"""Implementations used to answers to review 1."""

from __future__ import annotations

import contextlib
from typing import TYPE_CHECKING

import fppanalysis
import matplotlib.pyplot as plt
import numpy as np
import plastik
import volcano_base
import volcano_base.manipulate as vbm
import xarray as xr
from rich.console import Console

import vdd.load
import vdd.utils
from vdd.plotting.deconv_ob16_cesm2 import PlotResponseFunctions

if TYPE_CHECKING:
    from collections.abc import Callable

    from numpy.typing import NDArray

plt.style.use({"legend.framealpha": 0.8})
_SAVE_DIR = volcano_base.config.SAVE_PATH / "review"
if not _SAVE_DIR.exists():
    _SAVE_DIR.mkdir(parents=False)


_COLORS = plt.rcParams["axes.prop_cycle"].by_key()["color"]


def _ob16_data() -> None:  # noqa: C901, PLR0915
    # Load raw data
    ob_dec = vdd.load.DeconvolveOB16(data="h0")
    file = "IVI2LoadingLatHeight501-2000_L18_c20100518.nc"
    if not (fn := volcano_base.config.DATA_PATH / "cesm-lme" / file).exists():
        raise FileNotFoundError
    da = xr.open_dataset(fn).colmass
    very_tiny = 1e-17
    da_nonzero = xr.where(da < very_tiny, np.nan, da)

    # Decide on tropical or not, and create delta forcing arrays for tropical,
    # non-tropical, north and south eruptions.
    tropical = ob_dec.so2.copy()
    north = ob_dec.so2.copy()
    south = ob_dec.so2.copy()
    tropical.data = np.zeros_like(tropical.data)
    north.data = np.zeros_like(north.data)
    south.data = np.zeros_like(south.data)
    nonzero_value_times = ob_dec.so2.where(ob_dec.so2 > 0, drop=True).time.data
    for times, times_dt in zip(
        vbm.dt2float(nonzero_value_times),
        nonzero_value_times,
        strict=True,
    ):
        current_eruption_value = ob_dec.so2.sel(time=times_dt).data
        val = da_nonzero.sel(time=times + 0.5, method="nearest")
        south_ = val.sel(lat=-30, method="nearest").data
        north_ = val.sel(lat=30, method="nearest").data
        match not np.isnan(north_), not np.isnan(south_):
            # We match for _existence_ rather than non-existence on the (Northern, Southern)
            # hemispheres.
            case (True, True):
                tropical.loc[{"time": times_dt}] = current_eruption_value
            case (True, False):
                north.loc[{"time": times_dt}] = current_eruption_value
            case (False, True):
                south.loc[{"time": times_dt}] = current_eruption_value
            case _:
                raise ValueError

    nontropical = north + south
    if not nontropical.equals(ob_dec.so2 - tropical):
        raise ValueError
    # Feeling confident!

    def _plot_found_eruption_arrs() -> None:
        # Plots to verify successful analysis. (Don't know how to test for this!?)
        ydiff = (ymax := ob_dec.so2.data.max()) - (ymin := ob_dec.so2.data.min())
        ylim = (ymin - ydiff * 0.05, ymax + ydiff * 0.05)
        alpha = 0.7
        plt.figure()
        ob_dec.so2.plot()
        plt.ylim(ylim)
        plt.savefig(_SAVE_DIR / "original")
        fig, axs = plastik.figure_grid(rows=1, columns=2)
        plt.sca(axs[0])
        north.plot(label="N")
        south.plot(label="S", alpha=alpha)
        plt.ylim(ylim)
        plt.legend()
        plt.sca(axs[1])
        tropical.plot(label="T")
        nontropical.plot(label="NT", alpha=alpha)
        plt.ylim(ylim)
        plt.legend()
        plt.savefig(_SAVE_DIR / "eruption-arrs")
        plt.show()

    def _plot_pulse_shapes() -> None:
        ob_dec.name = "CESM2 OB16"
        trop_dec = vdd.load.DeconvolveOB16(data="h0")
        trop_dec.so2 = tropical
        trop_dec.name = "CESM2 Trop"
        nontrop_dec = vdd.load.DeconvolveOB16(data="h0")
        nontrop_dec.so2 = nontropical
        nontrop_dec.name = "CESM2 Non-trop"
        north_dec = vdd.load.DeconvolveOB16(data="h0")
        north_dec.so2 = north
        north_dec.name = "CESM2 North"
        south_dec = vdd.load.DeconvolveOB16(data="h0")
        south_dec.so2 = south
        south_dec.name = "CESM2 South"

        def _deconv(
            signal: NDArray[np.float64],
            forcing: NDArray[np.float64],
        ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
            if not len(signal) % 2 or not len(forcing) % 2:
                raise vdd.load.EvenLengthError
            x = np.arange(len(signal)) - len(signal) // 2
            guess = np.heaviside(x, 1)
            # box_width = 20 * 12
            # guess[np.argwhere(x > box_width)] = 0
            n_iters = 2000
            with (
                Console().status("[bold yellow]Deconvolving ...", spinner="point"),
                contextlib.redirect_stdout(None),
            ):
                out, err = fppanalysis.RL_gauss_deconvolve(
                    signal,
                    forcing,
                    initial_guess=guess,
                    iteration_list=n_iters,
                )
            return out, err

        for dec in (ob_dec, trop_dec, nontrop_dec, nontrop_dec, south_dec):
            dec.change_deconvolution_method(_deconv)
        plot = PlotResponseFunctions(
            ob_dec, trop_dec, nontrop_dec, north_dec, south_dec, norm=False
        )
        fig, axs = plastik.figure_grid(rows=2, columns=2)
        for i, ax in enumerate(axs):
            plot.grayscale_plot(ax, "temp", i + 1)
            ax.legend()
            ax.set_xlim((-2, 30))
        fig.savefig(_SAVE_DIR / "select-eruptions")
        plt.show()

    _plot_found_eruption_arrs()
    _plot_pulse_shapes()


# I have no idea how to type this creature.
class Pipe:
    """Simple pipeline/chain implementation."""

    def __init__(self, value, func=None):  # noqa: ANN001, ANN204
        self.value = value
        self.func = func

    def __getitem__(self, func):  # noqa: ANN001, ANN204
        """Implement slice/square bracket functionality."""
        return Pipe(self.value, func)

    def __call__(self, *args, **kwargs):  # noqa: D102, ANN002, ANN003, ANN204
        return Pipe(self.func(self.value, *args, **kwargs))

    def __repr__(self):  # noqa: D105, ANN204
        return f"Pipe({self.value}, {self.func})"

    def __rrshift__(self, *args, **kwargs):  # noqa: ANN204, ANN003, D105, ANN002
        return Pipe(self.func(self.value, *args, **kwargs))


class MetaFunc:
    """Simple meta-function implementing a unix style pipeline.

    From: https://gist.github.com/miyuchina/f463bbd8ad605d8c9a56aef83709ebe2
    """

    def __init__(self, func: Callable) -> None:
        self.func = func

    def __call__[T](self, *args: T) -> T:
        """Define behaviour for calling a class instance."""
        return self.func(*args)

    def __or__(self, other: MetaFunc | Callable) -> MetaFunc:
        """Define behaviour for the pipe operator."""
        return MetaFunc(lambda *args: other(self(*args)))


def _compare_means() -> None:  # noqa: PLR0915, C901
    def finder(strength: str) -> volcano_base.load.FindFiles:
        return (
            volcano_base.load.FindFiles()
            .find("e_BWma1850", "TREFHT", "h0", strength)
            .remove(
                "ens1"
                if strength not in {"tt-2sep", "tt-4sep", "medium-2sep", "medium-4sep"}
                else "ens0",
            )
            # .remove("ens2", "ens4")
            .keep_most_recent()
            .sort("ensemble")
        )

    # Eight eruptions in total.
    control_files = finder("control")
    s26_files = finder("medium")
    s400_files = finder("medium-plus")
    s26_4sep_files = finder("medium-4sep")
    s400_4sep_files = finder("tt-4sep")
    feb_and_aug = 2
    if len(s26_4sep_files) != feb_and_aug and len(s400_4sep_files) != feb_and_aug:
        raise ValueError

    def _2zero(arr: xr.DataArray) -> xr.DataArray:
        """Convert the time to start from zero."""
        if isinstance(arr.time.data[0], float):
            arr = arr.assign_coords(time=arr.time.data - arr.time.data[0])
        else:
            arr = arr.assign_coords(time=vbm.dt2float(arr.time.data) - 1850)
        arr.coords["time"].attrs["long_name"] = "Time-after-eruption"
        arr.coords["time"].attrs["units"] = "yr"
        return arr

    def shift_and_flatten(ff: volcano_base.load.FindFiles) -> list[xr.DataArray]:
        version = "meta"
        if version == "pipe":
            sft = vbm.shift_arrays
            flat = vbm.mean_flatten
            op = vbm.data_array_operation
            dims = ["lat", "lon"]
            return Pipe(ff.load())[flat](dims=dims)[sft](daily=False)[sft](custom=1)[
                op
            ](_2zero).value
        if version == "meta":

            def shift_notday(x: xr.DataArray) -> xr.DataArray:
                return vbm.shift_arrays(x, daily=False)

            def shift_custom(x: xr.DataArray) -> xr.DataArray:
                return vbm.shift_arrays(x, custom=1)

            def data_array_operation(x: xr.DataArray) -> xr.DataArray:
                return vbm.data_array_operation(x, _2zero)

            dims = ["lat", "lon"]

            def flat(x: xr.DataArray) -> xr.DataArray:
                return vbm.mean_flatten(x, dims=dims)

            meta = MetaFunc(flat) | shift_notday | shift_custom | data_array_operation
            return meta(ff.load())
        return op(sft(sft(flat(ff.load(), dims=dims), daily=False), custom=1), _2zero)

    def norm(arr: xr.DataArray) -> xr.DataArray:
        return (arr - arr.mean()) / arr.std()

    end = 4
    control = vbm.mean_flatten(control_files.load()[0], dims=["lat", "lon"])
    control = _2zero(vbm.subtract_climatology(control, control, "time.month")[0])
    control = norm(control.sel(time=slice(None, end)))
    s26_spring, s26_summer, s26_fall, s26_winter = shift_and_flatten(s26_files)
    s26_spring = _2zero(s26_spring.sel(time=slice(None, end)))
    s26_summer = _2zero(s26_summer.sel(time=slice(None, end)))
    s26_fall = _2zero(s26_fall.sel(time=slice(None, end)))
    s26_winter = _2zero(s26_winter.sel(time=slice(None, end)))
    s26_spring = norm(s26_spring)
    s26_summer = norm(s26_summer)
    s26_fall = norm(s26_fall)
    s26_winter = norm(s26_winter)
    s26_all_seasons = norm(s26_winter + s26_summer + s26_spring + s26_fall)
    s26_winter_summer = norm(s26_winter + s26_summer)
    s26_spring_fall = norm(s26_spring + s26_fall)
    s400_spring, s400_summer, s400_fall, s400_winter = shift_and_flatten(s400_files)
    s400_spring = _2zero(s400_spring.sel(time=slice(None, end)))
    s400_summer = _2zero(s400_summer.sel(time=slice(None, end)))
    s400_fall = _2zero(s400_fall.sel(time=slice(None, end)))
    s400_winter = _2zero(s400_winter.sel(time=slice(None, end)))
    s400_spring = norm(s400_spring)
    s400_summer = norm(s400_summer)
    s400_fall = norm(s400_fall)
    s400_winter = norm(s400_winter)
    s400_all_seasons = norm(s400_winter + s400_summer + s400_spring + s400_fall)
    s400_winter_summer = norm(s400_winter + s400_summer)
    s400_spring_fall = norm(s400_spring + s400_fall)
    # Plot
    fig, axs = plastik.figure_grid(rows=1, columns=2)
    plt.sca(axs[0])
    # times, corr = fppanalysis.corr_fun(control, control, 1 / 12)
    # plt.plot(times, corr, label="ctrl")
    times, corr = fppanalysis.corr_fun(s26_all_seasons, s26_all_seasons, 1 / 12)
    plt.plot(times, corr, label="all-all")
    times, corr = fppanalysis.corr_fun(s26_winter_summer, s26_winter_summer, 1 / 12)
    plt.plot(times, corr, label="fa-fa")
    times, corr = fppanalysis.corr_fun(s26_spring_fall, s26_spring_fall, 1 / 12)
    plt.plot(times, corr, label="mn-mn")
    times, corr = fppanalysis.corr_fun(s26_winter_summer, s26_spring_fall, 1 / 12)
    plt.plot(times, corr, label="fa-mn")
    plt.sca(axs[1])
    # times, corr = fppanalysis.corr_fun(control, control, 1 / 12)
    # plt.plot(times, corr, label="ctrl")
    times, corr = fppanalysis.corr_fun(s400_all_seasons, s400_all_seasons, 1 / 12)
    plt.plot(times, corr, label="all-all")
    times, corr = fppanalysis.corr_fun(s400_winter_summer, s400_winter_summer, 1 / 12)
    plt.plot(times, corr, label="fa-fa")
    times, corr = fppanalysis.corr_fun(s400_spring_fall, s400_spring_fall, 1 / 12)
    plt.plot(times, corr, label="mn-mn")
    times, corr = fppanalysis.corr_fun(s400_winter_summer, s400_spring_fall, 1 / 12)
    plt.plot(times, corr, label="fa-mn")
    [ax.set_xlabel("Time lag") for ax in axs]
    [ax.set_ylabel("Correlation") for ax in axs]
    [ax.legend() for ax in axs]
    plt.savefig(_SAVE_DIR / "correlations-wintersummer_springfall")
    plt.show()

    s26_feb, s26_aug = shift_and_flatten(s26_4sep_files)
    s26_feb_1 = s26_feb.sel(time=slice(None, 4))
    s26_feb_2 = _2zero(s26_feb.sel(time=slice(4, 8)))
    s26_aug_1 = s26_aug.sel(time=slice(None, 4))
    s26_aug_2 = _2zero(s26_aug.sel(time=slice(4, 8)))
    s400_feb, s400_aug = shift_and_flatten(s400_4sep_files)
    s400_feb_1 = s400_feb.sel(time=slice(None, 4))
    s400_feb_2 = _2zero(s400_feb.sel(time=slice(4, 8)))
    s400_aug_1 = s400_aug.sel(time=slice(None, 4))
    s400_aug_2 = _2zero(s400_aug.sel(time=slice(4, 8)))
    # Since we expect these eruptions to superpose, normalising should make them all
    # similar!
    s26_feb_1 = norm(s26_feb_1)
    s26_feb_2 = norm(s26_feb_2)
    s26_aug_1 = norm(s26_aug_1)
    s26_aug_2 = norm(s26_aug_2)
    s400_feb_1 = norm(s400_feb_1)
    s400_feb_2 = norm(s400_feb_2)
    s400_aug_1 = norm(s400_aug_1)
    s400_aug_2 = norm(s400_aug_2)
    # Make means from mixing and matching:
    # - s26_feb (1) + s26_aug (1)
    # - s26_feb (1) + s26_aug (2)
    # - s26_feb (2) + s26_aug (1)
    # - s26_feb (2) + s26_aug (2)
    # - s400_feb (1) + s400_aug (1)
    # - s400_feb (1) + s400_aug (2)
    # - s400_feb (2) + s400_aug (1)
    # - s400_feb (2) + s400_aug (2)
    s26_f1_a1 = norm(s26_feb_1 + s26_aug_1)
    s26_f1_a2 = norm(s26_feb_1 + s26_aug_2)
    s26_f2_a1 = norm(s26_feb_2 + s26_aug_1)
    s26_f2_a2 = norm(s26_feb_2 + s26_aug_2)
    s400_f1_a1 = norm(s400_feb_1 + s400_aug_1)
    s400_f1_a2 = norm(s400_feb_1 + s400_aug_2)
    s400_f2_a1 = norm(s400_feb_2 + s400_aug_1)
    s400_f2_a2 = norm(s400_feb_2 + s400_aug_2)
    # s26[0].plot()
    # s26[1].plot()
    # s26_feb.plot()
    # s26_aug.plot()
    # s400_feb.plot()
    # s400_aug.plot()
    # Individually
    fig, axs = plastik.figure_grid(rows=2, columns=2)
    plt.sca(axs[0])
    s26_feb_1.plot(label="Feb_1")
    s26_feb_2.plot(label="Feb_2")
    s26_aug_1.plot(label="Aug_1")
    s26_aug_2.plot(label="Aug_2")
    plt.sca(axs[1])
    s400_feb_1.plot(label="Feb_1")
    s400_feb_2.plot(label="Feb_2")
    s400_aug_1.plot(label="Aug_1")
    s400_aug_2.plot(label="Aug_2")
    plt.sca(axs[2])
    s26_f1_a1.plot(label="f1a1")
    s26_f1_a2.plot(label="f1a2")
    s26_f2_a1.plot(label="f2a1")
    s26_f2_a2.plot(label="f2a2")
    plt.sca(axs[3])
    s400_f1_a1.plot(label="f1a1")
    s400_f1_a2.plot(label="f1a2")
    s400_f2_a1.plot(label="f2a1")
    s400_f2_a2.plot(label="f2a2")
    axs[0].legend()
    axs[1].legend()
    axs[2].legend()
    axs[3].legend()
    plt.savefig(_SAVE_DIR / "individual-eruptions")
    # plt.show()
    # Compute correlation functions.
    fig, axs = plastik.figure_grid(rows=1, columns=2)
    plt.sca(axs[0])
    times, corr = fppanalysis.corr_fun(s26_f1_a1, s26_f1_a1, 1 / 12)
    plt.plot(times, corr, label="f1a1-f1a1")
    times, corr = fppanalysis.corr_fun(s26_f1_a2, s26_f1_a2, 1 / 12)
    plt.plot(times, corr, label="f1a2-f1a2")
    times, corr = fppanalysis.corr_fun(s26_f1_a1, s26_f1_a2, 1 / 12)
    plt.plot(times, corr, label="f1a1-f1a2")
    plt.sca(axs[1])
    times, corr = fppanalysis.corr_fun(s400_f1_a1, s400_f1_a1, 1 / 12)
    plt.plot(times, corr, label="f1a1-f1a1")
    times, corr = fppanalysis.corr_fun(s400_f1_a2, s400_f1_a2, 1 / 12)
    plt.plot(times, corr, label="f1a2-f1a2")
    times, corr = fppanalysis.corr_fun(s400_f1_a1, s400_f1_a2, 1 / 12)
    plt.plot(times, corr, label="f1a1-f1a2")
    [ax.set_xlabel("Time lag") for ax in axs]
    [ax.set_ylabel("Correlation") for ax in axs]
    [ax.legend() for ax in axs]
    plt.savefig(_SAVE_DIR / "correlations")
    plt.show()


def _plot_sea_ice_s26() -> None:
    sims = (
        "control",
        "medium",
        "medium-2sep",
        "medium-4sep",
    )
    new: volcano_base.load.FindFiles = (
        volcano_base.load.FindFiles()
        .find("ICEFRAC", "e_BWma1850", sims, "h0")
        .keep_most_recent()
    )
    fig, axs = plastik.figure_grid(rows=1, columns=2)
    for c, new_ in zip(
        ["grey", _COLORS[0], _COLORS[1], _COLORS[2]], sims, strict=False
    ):
        arrs_: volcano_base.load.FindFiles = new.copy().keep(new_)
        if len(arrs_) > 1:
            large_ens = 4
            arrs_.remove("ens1" if len(arrs_) >= large_ens else "ens0")
        arr_ = arrs_.load()
        arr_ = volcano_base.manipulate.mean_flatten(arr_, dims=["lat", "lon"])
        if new_ == "control":
            ctrl = arr_[0]
            continue
        first = True
        plt.sca(axs[0])
        for a in arr_:
            kwargs = {"label": vdd.utils.name_swap(new_).upper()} if first else {}
            a.plot(c=c, **kwargs)
            first = False
        [a.plot(c=c) for a in arr_]
        plt.sca(axs[1])

        def sub(the_array: xr.DataArray) -> xr.DataArray:
            tmp = (
                volcano_base.manipulate.subtract_climatology(
                    the_array,
                    ctrl,  # noqa: B023
                    groupby="time.month",
                )[0]
                + ctrl.data.mean()  # noqa: B023
            )
            return tmp.assign_attrs(the_array.attrs)

        arr_ = volcano_base.manipulate.data_array_operation(arr_, sub)
        [a.plot(c=c) for a in arr_]
    axs[0].set_xlim((vdd.utils.d2n("1850-01-01"), vdd.utils.d2n("1856-01-01")))
    axs[1].set_xlim((vdd.utils.d2n("1850-01-01"), vdd.utils.d2n("1856-01-01")))
    axs[0].legend()
    plt.savefig(_SAVE_DIR / "sea-ice-verification")
    plt.show()


if __name__ == "__main__":
    # _ob16_data()
    _compare_means()
    # _plot_sea_ice_s26()
