"""Plot the deconvolution comparison between OB16 and CESM2."""

from typing import Self

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import volcano_base
from numpy.typing import NDArray
from rich.console import Console

import vdd.deconvolve_methods
import vdd.load
import vdd.utils
from vdd.utils import name_swap as ns

_SAVE_DIR = volcano_base.config.SAVE_PATH / "deconv_ob16_cesm2"
if not _SAVE_DIR.exists():
    _SAVE_DIR.mkdir(parents=False)

plt.style.use(
    [
        "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle",
        "vdd.extra",
        "vdd.jgr",
        {"text.latex.preamble": r"\usepackage{amsmath,siunitx}"},
    ],
)

DataCESM = vdd.load.CESMData
DecCESM = vdd.load.DeconvolveCESM
# CESM2
console = Console()
_USE_NEW_DECONV = False


def _new_deconv(
    signal: NDArray[np.float64], forcing: NDArray[np.float64]
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    with console.status("Deconvolving with new method"):
        res, err = vdd.deconvolve_methods.alternative_deconv(signal, forcing)
    return res, err


padding = vdd.load.PaddingMethod.NOISE
dec_cesm_int4 = DecCESM(pad_before=padding, cesm=DataCESM(strength="tt-4sep"))
dec_cesm_int2 = DecCESM(pad_before=padding, cesm=DataCESM(strength="tt-2sep"))
dec_cesm_med4 = DecCESM(pad_before=padding, cesm=DataCESM(strength="medium-4sep"))
dec_cesm_med2 = DecCESM(pad_before=padding, cesm=DataCESM(strength="medium-2sep"))
dec_cesm_e = DecCESM(pad_before=padding, cesm=DataCESM(strength="size5000"))
dec_cesm_s = DecCESM(pad_before=padding, cesm=DataCESM(strength="strong"))
dec_cesm_p = DecCESM(pad_before=padding, cesm=DataCESM(strength="medium-plus"))
dec_cesm_m = DecCESM(pad_before=padding, cesm=DataCESM(strength="medium"))
# Original
dec_ob16 = vdd.load.DeconvolveOB16(data="h0", length=int(12 * 1000) + 1)
dec_ob16.name = "OB16"
all_decs = (
    dec_ob16,
    dec_cesm_m,
    dec_cesm_p,
    dec_cesm_s,
    dec_cesm_e,
    dec_cesm_med2,
    dec_cesm_med4,
    dec_cesm_int2,
    dec_cesm_int4,
)
if _USE_NEW_DECONV:
    for dec in all_decs:
        # NOTE: The OB16 time series takes a while to deconvolve with this method.
        # (10-20 minutes.) We skip it by default.
        if dec.name == "OB16":
            continue
        dec.change_deconvolution_method(_new_deconv)


class PlotResponseFunctions:
    """Plot the response functions.

    This takes any number of deconvolution objects as input and creates plots that
    compare the response functions of the deconvolution objects.

    Parameters
    ----------
    *decs : vdd.load.Deconvolve
        The deconvolution objects to compare.
    norm : bool, optional
        Whether to normalise the response functions, by default False.
    """

    def __init__(self: Self, *decs: vdd.load.Deconvolve, norm: bool = False) -> None:
        self.decs = decs
        self.norm = norm

    def plot_rf_so2(
        self: Self,
        fig: mpl.figure.Figure | None = None,
        save_as: str = "rf-so2",  # noqa: ARG002
    ) -> mpl.figure.Figure:
        """Plot the radiative forcing to SO2 response functions."""
        match fig:
            case None:
                return plt.figure()
            case mpl.figure.Figure():
                ax = fig.gca()
                rf_xlim = (-2, 10)
                ax.set_xlim(rf_xlim)
                ax.set_xlabel("Time lag [yr]")
                ax.set_ylabel(
                    f"{"Normalised " if self.norm else ""}RF to {"\n" if self.norm else ""}SO2 response [1]",
                )
                ax.legend()
                # fig.savefig(_SAVE_DIR / f"{save_as}-{"norm" if self.norm else "abs"}")
                return fig
            case _:
                raise ValueError

    def plot_rf_so2_decay(
        self: Self,
        fig: mpl.figure.Figure | None = None,
        save_as: str = "rf-so2_decay",  # noqa: ARG002
    ) -> mpl.figure.Figure:
        """Plot the radiative forcing to SO2 decay response functions."""
        match fig:
            case None:
                return plt.figure()
            case mpl.figure.Figure():
                ax = fig.gca()
                rf_xlim = (-2, 10)
                ax.set_xlim(rf_xlim)
                ax.set_xlabel("Time lag [yr]")
                ax.set_ylabel(
                    f"{"Normalised " if self.norm else ""}RF to {"\n" if self.norm else ""}SO2 burden response [1]",
                )
                ax.legend()
                # fig.savefig(_SAVE_DIR / f"{save_as}-{"norm" if self.norm else "abs"}")
                return fig
            case _:
                raise ValueError

    def plot_temp_so2(
        self: Self,
        fig: mpl.figure.Figure | None = None,
        save_as: str = "temp-so2",  # noqa: ARG002
    ) -> mpl.figure.Figure:
        """Plot the temperature to SO2 response functions."""
        match fig:
            case None:
                return plt.figure()
            case mpl.figure.Figure():
                ax = fig.gca()
                temp_xlim = (-2, 20)
                ax.set_xlim(temp_xlim)
                ax.set_xlabel("Time lag [yr]")
                ax.set_ylabel(
                    f"{"Normalised t" if self.norm else "T"}emperature to {"\n" if self.norm else ""}SO2 response [1]",
                )
                ax.legend()
                # fig.savefig(_SAVE_DIR / f"{save_as}-{"norm" if self.norm else "abs"}")
                return fig
            case _:
                raise ValueError

    def plot_temp_rf(
        self: Self,
        fig: mpl.figure.Figure | None = None,
        save_as: str = "temp-rf",  # noqa: ARG002
    ) -> mpl.figure.Figure:
        """Plot the temperature to radiative forcing response functions."""
        match fig:
            case None:
                return plt.figure()
            case mpl.figure.Figure():
                ax = fig.gca()
                temp_xlim = (-1, 10)
                ax.set_xlim(temp_xlim)
                ax.set_xlabel("Time lag [yr]")
                ax.set_ylabel(
                    f"{"Normalised t" if self.norm else "T"}emperature to {"\n" if self.norm else ""}RF response [1]",
                )
                ax.legend()
                # fig.savefig(_SAVE_DIR / f"{save_as}-{"norm" if self.norm else "abs"}")
                return fig
            case _:
                raise ValueError

    @staticmethod
    def plot_grayscale_highlight(
        fig: mpl.figure.Figure | None = None,
        save_as: str = "temp-so2-gs",
    ) -> mpl.figure.Figure:
        """Plot the temperature to SO2 response functions with a grayscale highlight."""
        rows = 4
        cols = 2
        match fig:
            case None:
                fig, _ = vdd.utils.figure_multiple_rows_columns(
                    rows,
                    cols,
                    share_axes="x",
                )
                return fig
            case mpl.figure.Figure():
                axs = fig.get_axes()
                for ax in axs:
                    temp_xlim = (-2, 20)
                    ax.set_xlim(temp_xlim)
                    ax.set_xlabel("Time lag [yr]")
                    sub = "R" if "rf" in save_as else "T"
                    ax.set_ylabel(f"$\\varphi_{{{sub}}} / \\max\\varphi_{{{sub}}}$")
                    ax.legend()
                new_deconv = "-alternative_deconv" if _USE_NEW_DECONV else ""
                fig.savefig(_SAVE_DIR / f"{save_as}{new_deconv}")
                return fig
            case _:
                raise ValueError

    def _norm_plot(
        self: Self, temp_so2_gs_a: list[mpl.axes.Axes], rf_so2_gs_a: list[mpl.axes.Axes]
    ) -> None:
        for ax_list, res_name in zip(
            [temp_so2_gs_a, rf_so2_gs_a],
            ["temp", "rf"],
            strict=True,
        ):
            for i, ax in enumerate(ax_list):
                i_ = i
                self._grayscale_plot(ax, res_name, i_)

    def _grayscale_plot(self: Self, ax: mpl.axes.Axes, res_name: str, i_: int) -> None:  # noqa: C901
        for dec in self.decs:
            name = ns(dec.name)
            clr = "r"
            name = name.replace("CESM2 ", "").upper()
            response = -1 * getattr(dec, f"response_{res_name}_so2")
            # Effectively a min-max feature scaling
            # https://en.wikipedia.org/wiki/Feature_scaling#Rescaling_(min-max_normalization)
            scale = max(response, key=abs)
            arr = response / scale
            lab = f"$\\varphi_{{{res_name[0].upper()}}}^{{\\text{{{name}}}}}$ ({vdd.utils.s2n(scale)})"
            kwargs = {"c": clr, "zorder": 10, "label": lab, "lw": 1}
            match i_, name:
                case _, "OB16":
                    ax.plot(dec.tau, arr, label=lab, c="k", zorder=5, lw=1)
                case 0, "SMALL":
                    ax.plot(dec.tau, arr, **kwargs)  # type: ignore[arg-type]
                case 1, "INTERMEDIATE":
                    ax.plot(dec.tau, arr, **kwargs)  # type: ignore[arg-type]
                case 2, "STRONG":
                    ax.plot(dec.tau, arr, **kwargs)  # type: ignore[arg-type]
                case 3, "EXTREME":
                    ax.plot(dec.tau, arr, **kwargs)  # type: ignore[arg-type]
                case 4, "SMALL-2SEP":
                    ax.plot(dec.tau, arr, **kwargs)  # type: ignore[arg-type]
                case 5, "SMALL-4SEP":
                    ax.plot(dec.tau, arr, **kwargs)  # type: ignore[arg-type]
                case 6, "INT-2SEP":
                    ax.plot(dec.tau, arr, **kwargs)  # type: ignore[arg-type]
                case 7, "INT-4SEP":
                    ax.plot(dec.tau, arr, **kwargs)  # type: ignore[arg-type]
                case _:
                    ax.plot(dec.tau, arr, label=f"_{lab}", c="gray", lw=0.5)

    def run(self: Self, save_as: list[str] | None = None) -> None:
        """Run the class."""
        # rf_so2 = self.plot_rf_so2()
        # rf_so2_decay = self.plot_rf_so2_decay()
        # temp_so2 = self.plot_temp_so2()
        # temp_rf = self.plot_temp_rf()
        temp_so2_gs = self.plot_grayscale_highlight()
        rf_so2_gs = self.plot_grayscale_highlight()
        # for dec in self.decs:
        #     rf_so2_resp = (
        #         vdd.utils.normalise(dec.response_rf_so2)
        #         if self.norm
        #         else dec.response_rf_so2
        #     )
        #     rf_so2_decay_resp = (
        #         vdd.utils.normalise(dec.response_rf_so2_decay)
        #         if self.norm
        #         else dec.response_rf_so2_decay
        #     )
        #     temp_so2_resp = (
        #         vdd.utils.normalise(dec.response_temp_so2)
        #         if self.norm
        #         else dec.response_temp_so2
        #     )
        #     _temp_rf_resp = dec.deconvolve(
        #         dec.response_temp_so2, dec.response_rf_so2
        #     )[0].flatten()
        #     temp_rf_resp = (
        #         vdd.utils.normalise(_temp_rf_resp) if self.norm else _temp_rf_resp
        #     )
        #     rf_so2.gca().plot(dec.tau, rf_so2_resp, label=ns(dec.name))
        #     rf_so2_decay.gca().plot(dec.tau, rf_so2_decay_resp, label=ns(dec.name))
        #     temp_so2.gca().plot(dec.tau, temp_so2_resp, label=ns(dec.name))
        #     temp_rf.gca().plot(dec.tau, temp_rf_resp, label=ns(dec.name))
        if self.norm:
            self._norm_plot(temp_so2_gs.get_axes(), rf_so2_gs.get_axes())
        save_as = save_as or [
            "rf-so2",
            "rf-so2_decay",
            "temp-so2",
            "temp-rf",
            "temp-so2-gs",
            "rf-so2-gs",
        ]
        # self.plot_rf_so2(rf_so2, save_as[0])
        # self.plot_rf_so2_decay(rf_so2_decay, save_as[1])
        # self.plot_temp_so2(temp_so2, save_as[2])
        # self.plot_temp_rf(temp_rf, save_as[3])
        if self.norm:
            self.plot_grayscale_highlight(temp_so2_gs, save_as[4])
            self.plot_grayscale_highlight(rf_so2_gs, save_as[5])


def main() -> None:
    """Run the main function."""
    PlotResponseFunctions(*all_decs, norm=True).run()
    plt.show()


if __name__ == "__main__":
    main()
