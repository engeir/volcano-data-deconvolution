"""Plot the simulations ensembles and their volcanoes."""

from collections.abc import Sequence

import matplotlib.pyplot as plt
import volcano_base
from matplotlib.axes import Axes

plt.style.use(
    [
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
    ]
)

_SAVE_DIR = volcano_base.config.SAVE_PATH
if not _SAVE_DIR.exists():
    _SAVE_DIR.mkdir(parents=False)

_COLORS = plt.rcParams["axes.prop_cycle"].by_key()["color"]
_SIM_START = 1 / 12


def _define_eruption(
    eruption_time: float, text: str, ax: Axes, *, include_double: bool = True
) -> Axes:
    """Give an eruption date and draw a vertical line."""
    text_y = 22 if include_double else 16
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
    sim_end: float,
    color: str,
    ax: Axes,
) -> Axes:
    arrow_kwargs = {"width": 0, "head_width": 0.5, "head_length": 0.1, "lw": 0.7}
    scatter_kwargs = {"s": 6, "color": "k"}
    ax.arrow(_SIM_START, vertical_position, sim_end, 0, color=color, **arrow_kwargs)
    vertical_scatter = (
        [vertical_position for _ in eruption_time]
        if isinstance(eruption_time, Sequence)
        else vertical_position
    )
    ax.scatter(eruption_time, vertical_scatter, **scatter_kwargs)  # type: ignore[arg-type]
    return ax


def _create_all_eruptions(ax: Axes, *, include_double: bool) -> None:
    _define_eruption(2 / 12, "0-02", ax, include_double=include_double)
    _define_eruption(5 / 12, "0-05", ax, include_double=include_double)
    _define_eruption(8 / 12, "0-08", ax, include_double=include_double)
    _define_eruption(11 / 12, "0-11", ax, include_double=include_double)
    _define_eruption(1 + 2 / 12, "1-02", ax, include_double=include_double)
    if include_double:
        _define_eruption(2 + 2 / 12, "2-02", ax)
        _define_eruption(2 + 8 / 12, "2-08", ax)
        _define_eruption(4 + 2 / 12, "4-02", ax)
        _define_eruption(4 + 8 / 12, "4-08", ax)


def gantt() -> None:  # noqa: PLR0915
    """Create a Gantt chart of the simulation timeline."""
    # sim_end = 2 + 1 / 12  # Used for paper 1.
    sim_end = 5 + 1 / 12  # Used for paper 2.
    long_sims = 4
    small_c = _COLORS[0]
    intermediate_c = _COLORS[1]
    strong_c = _COLORS[2]
    strongn_c = _COLORS[5]
    extreme_c = _COLORS[3] if sim_end > long_sims else _COLORS[6]
    tt2sep_c = _COLORS[4]
    tt4sep_c = _COLORS[5]
    small2sep_c = _COLORS[6]
    small4sep_c = _COLORS[7]
    plt.figure()
    ax = plt.gca()
    # Remove all borders
    ax.axis("off")
    # Time line
    time_line_h = -3
    ax.arrow(
        _SIM_START,
        time_line_h,
        sim_end,
        0,
        color="grey",
        width=0,
        head_width=0.5,
        head_length=0.1,
        lw=0.5,
    )
    ax.text(_SIM_START - 0.1, time_line_h, "0", ha="right", va="center")
    ax.text(sim_end + 0.3, time_line_h, "21", ha="left", va="center")
    ax.vlines(_SIM_START, time_line_h - 0.5, time_line_h + 0.5, lw=0.5, color="grey")
    ax.vlines(sim_end + 0.2, time_line_h - 0.5, time_line_h + 0.5, lw=0.5, color="grey")
    _create_all_eruptions(ax, include_double=sim_end > long_sims)
    # Control
    ax.text(0, -1, "CONTROL", ha="right", va="center")
    ax.arrow(
        _SIM_START,
        -1,
        sim_end,
        0,
        color="k",
        width=0,
        head_width=0.5,
        head_length=0.1,
        lw=0.7,
    )
    # S26
    ax.text(0, 1.5, "S26", ha="right", va="center")
    _place_eruption(5 / 12, 0, sim_end, small_c, ax)
    _place_eruption(8 / 12, 1, sim_end, small_c, ax)
    _place_eruption(11 / 12, 2, sim_end, small_c, ax)
    _place_eruption(1 + 2 / 12, 3, sim_end, small_c, ax)
    # S400
    ax.text(0, 5.5, "S400", ha="right", va="center")
    _place_eruption(5 / 12, 4, sim_end, intermediate_c, ax)
    _place_eruption(8 / 12, 5, sim_end, intermediate_c, ax)
    _place_eruption(11 / 12, 6, sim_end, intermediate_c, ax)
    _place_eruption(1 + 2 / 12, 7, sim_end, intermediate_c, ax)
    # S1629
    ax.text(0, 9.5, "S1629", ha="right", va="center")
    _place_eruption(5 / 12, 8, sim_end, strong_c, ax)
    _place_eruption(8 / 12, 9, sim_end, strong_c, ax)
    _place_eruption(11 / 12, 10, sim_end, strong_c, ax)
    _place_eruption(1 + 2 / 12, 11, sim_end, strong_c, ax)
    if sim_end < long_sims:
        # S1629N
        ax.text(0, 12.5, "S1629N", ha="right", va="center")
        _place_eruption(2 / 12, 12, sim_end, strongn_c, ax)
        _place_eruption(8 / 12, 13, sim_end, strongn_c, ax)
        # S3000
        ax.text(0, 14.5, "S3000", ha="right", va="center")
        _place_eruption(5 / 12, 14, sim_end, extreme_c, ax)
        _place_eruption(11 / 12, 15, sim_end, extreme_c, ax)
        plt.savefig(_SAVE_DIR / "simulation_timeline_mini")
        return
    # S3000
    ax.text(0, 12.5, "S3000", ha="right", va="center")
    _place_eruption(5 / 12, 12, sim_end, extreme_c, ax)
    _place_eruption(11 / 12, 13, sim_end, extreme_c, ax)
    # `S26-2SEP`
    ax.text(0, 14.5, "S26-2SEP", ha="right", va="center")
    _place_eruption([2 / 12, 2 + 2 / 12], 14, sim_end, small2sep_c, ax)
    _place_eruption([8 / 12, 2 + 8 / 12], 15, sim_end, small2sep_c, ax)
    # `S26-4SEP`
    ax.text(0, 16.5, "S26-4SEP", ha="right", va="center")
    _place_eruption([2 / 12, 4 + 2 / 12], 16, sim_end, small4sep_c, ax)
    _place_eruption([8 / 12, 4 + 8 / 12], 17, sim_end, small4sep_c, ax)
    # `S400-2SEP`
    ax.text(0, 18.5, "S400-2SEP", ha="right", va="center")
    _place_eruption([2 / 12, 2 + 2 / 12], 18, sim_end, tt2sep_c, ax)
    _place_eruption([8 / 12, 2 + 8 / 12], 19, sim_end, tt2sep_c, ax)
    # `S400-4SEP`
    ax.text(0, 20.5, "S400-4SEP", ha="right", va="center")
    _place_eruption([2 / 12, 4 + 2 / 12], 20, sim_end, tt4sep_c, ax)
    _place_eruption([8 / 12, 4 + 8 / 12], 21, sim_end, tt4sep_c, ax)
    # `ax.set_ylim(-3, 21.5)`
    plt.savefig(_SAVE_DIR / "simulation_timeline")


def main() -> None:
    """Run the main function."""
    gantt()
    plt.show()


if __name__ == "__main__":
    main()
