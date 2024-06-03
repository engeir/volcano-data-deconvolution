"""Plot the simulations ensembles and their volcanoes."""

from collections.abc import Sequence

import matplotlib.pyplot as plt
import volcano_base
from matplotlib.axes import Axes

plt.style.use([
    "vdd.extra",
    "vdd.jgr",
    {
        "text.usetex": True,
        "font.size": 8,
        "axes.labelsize": 8,
        "legend.fontsize": 8,
        "figure.dpi": 300,
        "figure.figsize": (3.37, 2.08277),
        "figure.subplot.left": 0.30,
        "figure.subplot.right": 0.9,
        "figure.subplot.bottom": 0.05,
        "figure.subplot.top": 0.8,
    },
])

_SAVE_DIR = volcano_base.config.SAVE_PATH
if not _SAVE_DIR.exists():
    _SAVE_DIR.mkdir(parents=False)

sim_start = 1 / 12
sim_end = 5 + 1 / 12


def _define_eruption(eruption_time: float, text: str, ax: Axes) -> Axes:
    """Give an eruption date and draw a vertical line."""
    text_y = 18
    axvline_kwargs = {
        "lw": 0.7,
        "ls": "--",
        "color": "grey",
        "alpha": 0.5,
        "zorder": -1,
    }
    ax.text(eruption_time, text_y, text, rotation=90)
    ax.axvline(eruption_time, 0, 1, **axvline_kwargs)
    return ax


def _place_eruption(
    eruption_time: float | Sequence[float],
    vertical_position: float,
    color: str,
    ax: Axes,
) -> Axes:
    arrow_kwargs = {"width": 0, "head_width": 0.5, "head_length": 0.1, "lw": 0.7}
    scatter_kwargs = {"s": 6, "color": "k"}
    ax.arrow(sim_start, vertical_position, sim_end, 0, color=color, **arrow_kwargs)
    vertical_scatter = (
        [vertical_position for _ in eruption_time]
        if isinstance(eruption_time, Sequence)
        else vertical_position
    )
    ax.scatter(eruption_time, vertical_scatter, **scatter_kwargs)  # type: ignore[arg-type]
    return ax


def gantt() -> None:
    """Create a Gantt chart of the simulation timeline."""
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    small_c = colors[0]
    intermediate_c = colors[1]
    strong_c = colors[2]
    extreme_c = colors[3]
    tt2sep_c = colors[4]
    tt4sep_c = colors[5]
    plt.figure()
    ax = plt.gca()
    # Remove all borders
    ax.axis("off")
    # Time line
    ax.arrow(
        sim_start,
        -2,
        sim_end,
        0,
        color="k",
        width=0,
        head_width=0.5,
        head_length=0.1,
        lw=0.5,
    )
    ax.text(sim_start - 0.1, -2, "1850", ha="right", va="center")
    ax.text(sim_end + 0.3, -2, "1871", ha="left", va="center")
    ax.vlines(sim_start, -2 - 0.5, -2 + 0.5, lw=0.5, color="k")
    ax.vlines(sim_end + 0.2, -2 - 0.5, -2 + 0.5, lw=0.5, color="k")
    _define_eruption(2 / 12, "1850-02", ax)
    _define_eruption(5 / 12, "1850-05", ax)
    _define_eruption(8 / 12, "1850-08", ax)
    _define_eruption(11 / 12, "1850-11", ax)
    _define_eruption(1 + 2 / 12, "1851-02", ax)
    _define_eruption(2 + 2 / 12, "1852-02", ax)
    _define_eruption(2 + 8 / 12, "1852-08", ax)
    _define_eruption(4 + 2 / 12, "1854-02", ax)
    _define_eruption(4 + 8 / 12, "1854-08", ax)
    # SMALL
    ax.text(0, 1.5, "SMALL", ha="right", va="center")
    _place_eruption(5 / 12, 0, small_c, ax)
    _place_eruption(8 / 12, 1, small_c, ax)
    _place_eruption(11 / 12, 2, small_c, ax)
    _place_eruption(1 + 2 / 12, 3, small_c, ax)
    # INTERMEDIATE
    ax.text(0, 5.5, "INTERMEDIATE", ha="right", va="center")
    _place_eruption(5 / 12, 4, intermediate_c, ax)
    _place_eruption(8 / 12, 5, intermediate_c, ax)
    _place_eruption(11 / 12, 6, intermediate_c, ax)
    _place_eruption(1 + 2 / 12, 7, intermediate_c, ax)
    # STRONG
    ax.text(0, 9.5, "STRONG", ha="right", va="center")
    _place_eruption(5 / 12, 8, strong_c, ax)
    _place_eruption(8 / 12, 9, strong_c, ax)
    _place_eruption(11 / 12, 10, strong_c, ax)
    _place_eruption(1 + 2 / 12, 11, strong_c, ax)
    # EXTREME
    ax.text(0, 12.5, "EXTREME", ha="right", va="center")
    _place_eruption(5 / 12, 12, extreme_c, ax)
    _place_eruption(11 / 12, 13, strong_c, ax)
    # INT-2SEP
    ax.text(0, 14.5, "INT-2SEP", ha="right", va="center")
    _place_eruption([2 / 12, 2 + 2 / 12], 14, tt2sep_c, ax)
    _place_eruption([8 / 12, 2 + 8 / 12], 15, tt2sep_c, ax)
    # INT-4SEP
    ax.text(0, 16.5, "INT-4SEP", ha="right", va="center")
    _place_eruption([2 / 12, 4 + 2 / 12], 16, tt4sep_c, ax)
    _place_eruption([8 / 12, 4 + 8 / 12], 17, tt4sep_c, ax)
    # ax.set_ylim(-3, 21.5)
    plt.savefig(_SAVE_DIR / "simulation_timeline")


def main() -> None:
    """Run the main function."""
    gantt()
    plt.show()


if __name__ == "__main__":
    main()
