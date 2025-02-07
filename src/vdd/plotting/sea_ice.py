"""Plotting functions for sea ice data."""

import matplotlib.pyplot as plt
import volcano_base
import xarray as xr

from vdd.utils import name_swap as ns

_SAVE_DIR = volcano_base.config.SAVE_PATH / "sea-ice"
if not _SAVE_DIR.exists():
    _SAVE_DIR.mkdir(parents=False)

plt.style.use(
    [
        "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle",
        "vdd.extra",
        "vdd.jgr",
    ],
)


_COLORS = plt.rcParams["axes.prop_cycle"].by_key()["color"]
C = {
    "control": "grey",
    "medium": _COLORS[0],
    "medium-2sep": _COLORS[0],
    "medium-4sep": _COLORS[0],
    "medium-plus": _COLORS[1],
    "tt-2sep": _COLORS[1],
    "tt-4sep": _COLORS[1],
    "strong": _COLORS[2],
    "size5000": _COLORS[3],
}
LS = {
    "control": "-",
    "medium": "-",
    "medium-2sep": "dotted",
    "medium-4sep": "--",
    "medium-plus": "-",
    "tt-2sep": "dotted",
    "tt-4sep": "--",
    "strong": "-",
    "size5000": "-",
}


def plot_sea_ice() -> None:
    """Plot sea ice data."""
    sims = (
        "control",
        "medium",
        "medium-2sep",
        "medium-4sep",
        "medium-plus",
        "tt-2sep",
        "tt-4sep",
        "strong",
        "size5000",
    )
    new: volcano_base.load.FindFiles = (
        volcano_base.load.FindFiles()
        .find("ICEFRAC", "e_BWma1850", sims, "h0")
        .keep_most_recent()
    )
    plt.figure()
    for new_ in sims:
        arrs_: volcano_base.load.FindFiles = new.copy().keep(new_)
        if len(arrs_) > 1:
            large_ens = 4
            arrs_.remove(*["ens1"] if len(arrs_) >= large_ens else "ens0")
        print(arrs_)
        arr_ = arrs_.load()
        arr_ = volcano_base.manipulate.mean_flatten(arr_, dims=["lat", "lon"])
        if new_ == "control":
            ctrl = arr_[0]

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
        arr_ = volcano_base.manipulate.shift_arrays(arr_, daily=False)
        arr = volcano_base.manipulate.get_median(arr_, xarray=True)
        arr = volcano_base.manipulate.weighted_year_avg(arr)[:20]
        arr = arr.assign_coords(
            {
                "time": volcano_base.manipulate.dt2float(arr.time.data) - 1850,
            },
        )
        lab = f"$I_{{\\text{{{ns(new_).upper()}}}}}$"
        if "sep" in lab.lower():
            lab = f"_{lab}"
        if new_ == "control":
            plt.fill_between(
                arr.time.data,
                arr.mean() - arr.std(),
                arr.mean() + arr.std(),
                color=C[new_],
                label=lab,
            )
            plt.fill_between(
                arr.time.data,
                arr.mean() - 2 * arr.std(),
                arr.mean() + 2 * arr.std(),
                color=C[new_],
                label="_hidden",
                alpha=0.6,
            )
        else:
            arr.plot(label=lab, c=C[new_], ls=LS[new_])
    plt.xlabel("Time after first eruption [yr]")
    plt.ylabel("Sea Ice Fraction [1]")
    plt.legend(framealpha=0.4, loc="upper right")
    plt.savefig(_SAVE_DIR / "sea-ice.pdf")
    plt.show()


if __name__ == "__main__":
    plot_sea_ice()
