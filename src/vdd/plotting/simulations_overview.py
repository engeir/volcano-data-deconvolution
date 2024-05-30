"""Plot the simulations ensembles and their volcanoes."""

import matplotlib.pyplot as plt
import volcano_base

plt.style.use([
    # "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle",
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


def gantt() -> None:
    start_time = 1 / 12
    end_time = 5 + 1 / 12
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
    text_y = 18
    arrow_kwargs = {"width": 0, "head_width": 0.5, "head_length": 0.1, "lw": 0.7}
    axvline_kwargs = {
        "lw": 0.7,
        "ls": "--",
        "color": "grey",
        "alpha": 0.5,
        "zorder": -1,
    }
    scatter_kwargs = {"s": 6, "color": "k"}
    # Time line
    ax.arrow(
        start_time,
        -2,
        end_time,
        0,
        color="k",
        width=0,
        head_width=0.5,
        head_length=0.1,
        lw=0.5,
    )
    ax.text(start_time - 0.1, -2, "1850", ha="right", va="center")
    ax.text(end_time + 0.3, -2, "1871", ha="left", va="center")
    ax.vlines(start_time, -2 - 0.5, -2 + 0.5, lw=0.5, color="k")
    ax.vlines(end_time + 0.2, -2 - 0.5, -2 + 0.5, lw=0.5, color="k")
    ax.text(2 / 12, text_y, "1850-02", rotation=90)
    ax.text(5 / 12, text_y, "1850-05", rotation=90)
    ax.text(8 / 12, text_y, "1850-08", rotation=90)
    ax.text(11 / 12, text_y, "1850-11", rotation=90)
    ax.text(1 + 2 / 12, text_y, "1851-02", rotation=90)
    ax.text(2 + 2 / 12, text_y, "1852-02", rotation=90)
    ax.text(2 + 8 / 12, text_y, "1852-08", rotation=90)
    ax.text(4 + 2 / 12, text_y, "1854-02", rotation=90)
    ax.text(4 + 8 / 12, text_y, "1854-08", rotation=90)
    ax.axvline(2 / 12, **axvline_kwargs)
    ax.axvline(5 / 12, **axvline_kwargs)
    ax.axvline(8 / 12, **axvline_kwargs)
    ax.axvline(11 / 12, **axvline_kwargs)
    ax.axvline(1 + 2 / 12, **axvline_kwargs)
    ax.axvline(2 + 2 / 12, **axvline_kwargs)
    ax.axvline(2 + 8 / 12, **axvline_kwargs)
    ax.axvline(4 + 2 / 12, **axvline_kwargs)
    ax.axvline(4 + 8 / 12, **axvline_kwargs)
    # SMALL
    ax.text(0, 1.5, "SMALL", ha="right", va="center")
    ax.arrow(start_time, 0, end_time, 0, color=small_c, **arrow_kwargs)
    ax.arrow(start_time, 1, end_time, 0, color=small_c, **arrow_kwargs)
    ax.arrow(start_time, 2, end_time, 0, color=small_c, **arrow_kwargs)
    ax.arrow(start_time, 3, end_time, 0, color=small_c, **arrow_kwargs)
    ax.scatter(5 / 12, 0, **scatter_kwargs)
    ax.scatter(8 / 12, 1, **scatter_kwargs)
    ax.scatter(11 / 12, 2, **scatter_kwargs)
    ax.scatter(1 + 2 / 12, 3, **scatter_kwargs)
    # INTERMEDIATE
    ax.text(0, 5.5, "INTERMEDIATE", ha="right", va="center")
    ax.arrow(start_time, 4, end_time, 0, color=intermediate_c, **arrow_kwargs)
    ax.arrow(start_time, 5, end_time, 0, color=intermediate_c, **arrow_kwargs)
    ax.arrow(start_time, 6, end_time, 0, color=intermediate_c, **arrow_kwargs)
    ax.arrow(start_time, 7, end_time, 0, color=intermediate_c, **arrow_kwargs)
    ax.scatter(5 / 12, 4, **scatter_kwargs)
    ax.scatter(8 / 12, 5, **scatter_kwargs)
    ax.scatter(11 / 12, 6, **scatter_kwargs)
    ax.scatter(1 + 2 / 12, 7, **scatter_kwargs)
    # STRONG
    ax.text(0, 9.5, "STRONG", ha="right", va="center")
    ax.arrow(start_time, 8, end_time, 0, color=strong_c, **arrow_kwargs)
    ax.arrow(start_time, 9, end_time, 0, color=strong_c, **arrow_kwargs)
    ax.arrow(start_time, 10, end_time, 0, color=strong_c, **arrow_kwargs)
    ax.arrow(start_time, 11, end_time, 0, color=strong_c, **arrow_kwargs)
    ax.scatter(5 / 12, 8, **scatter_kwargs)
    ax.scatter(8 / 12, 9, **scatter_kwargs)
    ax.scatter(11 / 12, 10, **scatter_kwargs)
    ax.scatter(1 + 2 / 12, 11, **scatter_kwargs)
    # EXTREME
    ax.text(0, 12.5, "EXTREME", ha="right", va="center")
    ax.arrow(start_time, 12, end_time, 0, color=extreme_c, **arrow_kwargs)
    ax.arrow(start_time, 13, end_time, 0, color=extreme_c, **arrow_kwargs)
    ax.scatter(5 / 12, 12, **scatter_kwargs)
    ax.scatter(11 / 12, 13, **scatter_kwargs)
    # INT-2SEP
    ax.text(0, 14.5, "INT-2SEP", ha="right", va="center")
    ax.arrow(start_time, 14, end_time, 0, color=tt2sep_c, **arrow_kwargs)
    ax.arrow(start_time, 15, end_time, 0, color=tt2sep_c, **arrow_kwargs)
    ax.scatter([2 / 12, 2 + 2 / 12], [14, 14], **scatter_kwargs)
    ax.scatter([8 / 12, 2 + 8 / 12], [15, 15], **scatter_kwargs)
    # INT-4SEP
    ax.text(0, 16.5, "INT-4SEP", ha="right", va="center")
    ax.arrow(start_time, 16, end_time, 0, color=tt4sep_c, **arrow_kwargs)
    ax.arrow(start_time, 17, end_time, 0, color=tt4sep_c, **arrow_kwargs)
    ax.scatter([2 / 12, 4 + 2 / 12], [16, 16], **scatter_kwargs)
    ax.scatter([8 / 12, 4 + 8 / 12], [17, 17], **scatter_kwargs)
    # ax.set_ylim(-3, 21.5)
    plt.savefig(_SAVE_DIR / "simulation_timeline")


if __name__ == "__main__":
    gantt()
    plt.show()
