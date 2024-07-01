"""Plotting functions for sea ice data."""

import matplotlib.pyplot as plt
import volcano_base

from vdd.utils import name_swap as ns

_SAVE_DIR = volcano_base.config.SAVE_PATH / "sea-ice"
if not _SAVE_DIR.exists():
    _SAVE_DIR.mkdir(parents=False)

plt.style.use(
    [
        "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle",
        "vdd.extra",
    ],
)


def plot_sea_ice() -> None:
    """Plot sea ice data."""
    sims = (
        "medium",
        "medium-plus",
        "strong",
        "size5000",
        "medium-2sep",
        "medium-4sep",
        "tt-2sep",
        "tt-4sep",
    )
    new: volcano_base.load.FindFiles = (
        volcano_base.load.FindFiles()
        .find("ICEFRAC", "e_BWma1850", sims, "h0")
        .keep_most_recent()
    )
    plt.figure()
    for new_ in sims:
        arrs_: volcano_base.load.FindFiles = new.copy().keep(new_)
        large_ens = 4
        arrs_.remove("ens1" if len(arrs_) >= large_ens else "ens0")
        print(arrs_)
        arr_ = arrs_.load()
        arr_ = volcano_base.manipulate.shift_arrays(arr_, daily=False)
        arr_ = volcano_base.manipulate.mean_flatten(arr_, dims=["lat", "lon"])
        arr = volcano_base.manipulate.get_median(arr_, xarray=True)
        arr = volcano_base.manipulate.weighted_year_avg(arr)[:20]
        arr = arr.assign_coords(
            {
                "time": volcano_base.manipulate.dt2float(arr.time.data) - 1850,
            },
        )
        lab = f"$I_{{\\text{{{ns(new_).upper()}}}}}$"
        arr.plot(label=lab)
    plt.xlabel("Time after first eruption [yr]")
    plt.ylabel("Sea Ice Fraction [1]")
    plt.legend()
    plt.savefig(_SAVE_DIR / "sea-ice.pdf")
    plt.show()


if __name__ == "__main__":
    plot_sea_ice()
