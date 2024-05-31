"""Plot the deconvolution comparison between OB16 and CESM2."""

import cosmoplots
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import volcano_base

import vdd.load
import vdd.utils
from vdd.utils import name_swap as ns

_SAVE_DIR = volcano_base.config.SAVE_PATH / "deconv_ob16_cesm2"
if not _SAVE_DIR.exists():
    _SAVE_DIR.mkdir(parents=False)

plt.style.use([
    "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle",
    "vdd.extra",
    {"text.latex.preamble": r"\usepackage{amsmath,siunitx}"},
])

DataCESM = vdd.load.CESMData
DecCESM = vdd.load.DeconvolveCESM
# CESM2
dec_cesm_4sep = DecCESM(pad_before=True, cesm=DataCESM(strength="tt-4sep"))
dec_cesm_2sep = DecCESM(pad_before=True, cesm=DataCESM(strength="tt-2sep"))
dec_cesm_e = DecCESM(pad_before=True, cesm=DataCESM(strength="size5000"))
dec_cesm_s = DecCESM(pad_before=True, cesm=DataCESM(strength="strong"))
dec_cesm_p = DecCESM(pad_before=True, cesm=DataCESM(strength="medium-plus"))
dec_cesm_m = DecCESM(pad_before=True, cesm=DataCESM(strength="medium"))
# Original
dec_ob16 = vdd.load.DeconvolveOB16(data="h0", length=int(12 * 1000) + 1)
# Experimental (needs vdd.deconvolve_methods.alternative_deconv)
# ob16 = volcano_base.load.OttoBliesner(freq="h0", progress=True)
# dec_ob16 = vdd.load.DeconvolveOB16(data=ob16)
# dec_ob16.change_deconvolution_method(alternative_deconv)
dec_ob16.name = "OB16"
all_decs = (
    dec_ob16,
    dec_cesm_m,
    dec_cesm_p,
    dec_cesm_s,
    dec_cesm_e,
    dec_cesm_2sep,
    dec_cesm_4sep,
)


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

    def __init__(self, *decs: vdd.load.Deconvolve, norm: bool = False):
        self.decs = decs
        self.norm = norm

    def plot_rf_so2(
        self, fig: mpl.figure.Figure | None = None, save_as: str = "rf-so2"
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
                    f"{"Normalised " if self.norm else ""}RF to {"\n" if self.norm else ""}SO2 response [1]"
                )
                ax.legend()
                fig.savefig(_SAVE_DIR / f"{save_as}-{"norm" if self.norm else "abs"}")
                return fig
            case _:
                raise ValueError("rf must be a mpl.figure.Figure or None")

    def plot_rf_so2_decay(
        self, fig: mpl.figure.Figure | None = None, save_as: str = "rf-so2_decay"
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
                    f"{"Normalised " if self.norm else ""}RF to {"\n" if self.norm else ""}SO2 burden response [1]"
                )
                ax.legend()
                fig.savefig(_SAVE_DIR / f"{save_as}-{"norm" if self.norm else "abs"}")
                return fig
            case _:
                raise ValueError("rf must be a mpl.figure.Figure or None")

    def plot_temp_so2(
        self, fig: mpl.figure.Figure | None = None, save_as: str = "temp-so2"
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
                    f"{"Normalised t" if self.norm else "T"}emperature to {"\n" if self.norm else ""}SO2 response [1]"
                )
                ax.legend()
                fig.savefig(_SAVE_DIR / f"{save_as}-{"norm" if self.norm else "abs"}")
                return fig
            case _:
                raise ValueError("temp must be a mpl.figure.Figure or None")

    def plot_temp_rf(
        self, fig: mpl.figure.Figure | None = None, save_as: str = "temp-rf"
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
                    f"{"Normalised t" if self.norm else "T"}emperature to {"\n" if self.norm else ""}RF response [1]"
                )
                ax.legend()
                fig.savefig(_SAVE_DIR / f"{save_as}-{"norm" if self.norm else "abs"}")
                return fig
            case _:
                raise ValueError("temp must be a mpl.figure.Figure or None")

    @staticmethod
    def plot_grayscale_highlight(
        fig: mpl.figure.Figure | None = None, save_as: str = "temp-so2-gs"
    ) -> mpl.figure.Figure:
        """Plot the temperature to SO2 response functions with a grayscale highlight."""
        rows = 3
        cols = 2
        match fig:
            case None:
                fig, _ = cosmoplots.figure_multiple_rows_columns(rows, cols)
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
                fig.savefig(_SAVE_DIR / f"{save_as}")
                return fig
            case _:
                raise ValueError("temp must be a mpl.figure.Figure or None")

    def run(self, save_as: list[str] | None = None) -> None:
        """Run the class."""
        rf_so2 = self.plot_rf_so2()
        rf_so2_decay = self.plot_rf_so2_decay()
        temp_so2 = self.plot_temp_so2()
        temp_rf = self.plot_temp_rf()
        temp_so2_gs = self.plot_grayscale_highlight()
        rf_so2_gs = self.plot_grayscale_highlight()
        rf_so2_a = rf_so2.gca()
        rf_so2_decay_a = rf_so2_decay.gca()
        temp_so2_a = temp_so2.gca()
        temp_rf_a = temp_rf.gca()
        temp_so2_gs_a = temp_so2_gs.get_axes()
        rf_so2_gs_a = rf_so2_gs.get_axes()
        for dec in self.decs:
            rf_so2_resp = (
                vdd.utils.normalise(dec.response_rf_so2)
                if self.norm
                else dec.response_rf_so2
            )
            rf_so2_decay_resp = (
                vdd.utils.normalise(dec.response_rf_so2_decay)
                if self.norm
                else dec.response_rf_so2_decay
            )
            temp_so2_resp = (
                vdd.utils.normalise(dec.response_temp_so2)
                if self.norm
                else dec.response_temp_so2
            )
            _temp_rf_resp = dec._deconv_method(
                dec.response_temp_so2, dec.response_rf_so2
            )[0].flatten()
            temp_rf_resp = (
                vdd.utils.normalise(_temp_rf_resp) if self.norm else _temp_rf_resp
            )
            rf_so2_a.plot(dec.tau, rf_so2_resp, label=ns(dec.name))
            rf_so2_decay_a.plot(dec.tau, rf_so2_decay_resp, label=ns(dec.name))
            temp_so2_a.plot(dec.tau, temp_so2_resp, label=ns(dec.name))
            temp_rf_a.plot(dec.tau, temp_rf_resp, label=ns(dec.name))
        if self.norm:
            for ax_list, res_name in zip(
                [temp_so2_gs_a, rf_so2_gs_a], ["temp", "rf"], strict=True
            ):
                for i, ax in enumerate(ax_list):
                    i_ = i
                    for dec in self.decs:
                        name = ns(dec.name)
                        clr = "r"
                        name = name.replace("CESM2 ", "").upper()
                        response = getattr(dec, f"response_{res_name}_so2")
                        scale = np.nanmax(response)
                        arr = response / scale
                        lab = f"$\\varphi_{{{res_name[0].upper()}}}^{{\\text{{{name}}}}}$ ({vdd.utils.s2n(scale)})"
                        kwargs = {"c": clr, "zorder": 10, "label": lab, "lw": 1}
                        if name == "OB16":
                            ax.plot(dec.tau, arr, label=lab, c="k", zorder=5, lw=1)
                        elif i_ == 0 and name == "SMALL":
                            ax.plot(dec.tau, arr, **kwargs)
                        elif i_ == 1 and name == "INTERMEDIATE":
                            ax.plot(dec.tau, arr, **kwargs)
                        elif i_ == 2 and name == "STRONG":
                            ax.plot(dec.tau, arr, **kwargs)
                        elif i_ == 3 and name == "EXTREME":
                            ax.plot(dec.tau, arr, **kwargs)
                        elif i_ == 4 and name == "INT-2SEP":
                            ax.plot(dec.tau, arr, **kwargs)
                        elif i_ == 5 and name == "INT-4SEP":
                            ax.plot(dec.tau, arr, **kwargs)
                        else:
                            ax.plot(dec.tau, arr, label=f"_{lab}", c="gray", lw=0.5)
        save_as = save_as or [
            "rf-so2",
            "rf-so2_decay",
            "temp-so2",
            "temp-rf",
            "temp-so2-gs",
            "rf-so2-gs",
        ]
        self.plot_rf_so2(rf_so2, save_as[0])
        self.plot_rf_so2_decay(rf_so2_decay, save_as[1])
        self.plot_temp_so2(temp_so2, save_as[2])
        self.plot_temp_rf(temp_rf, save_as[3])
        if self.norm:
            self.plot_grayscale_highlight(temp_so2_gs, save_as[4])
            self.plot_grayscale_highlight(rf_so2_gs, save_as[5])


def _main() -> None:
    save_as = [
        "rf-so2",
        "rf-so2_decay",
        "temp-so2",
        "temp-rf",
        "temp-so2-gs",
        "rf-so2-gs",
    ]
    # for dec in all_decs[1:]:
    #     save_as_ = [
    #         val + "-" + vdd.utils.clean_filename(ns(dec.name)).name for val in save_as
    #     ]
    #     PlotResponseFunctions(dec_ob16, dec, norm=True).run(save_as_)
    #     PlotResponseFunctions(dec_ob16, dec, norm=False).run(save_as_)
    #     plt.close("all")
    # PlotResponseFunctions(*all_decs, norm=False).run()
    PlotResponseFunctions(*all_decs, norm=True).run()
    plt.show()
    return
    files = (
        [_SAVE_DIR / f"rf-so2-{k}.jpg" for k in ["abs", "norm"]],
        [_SAVE_DIR / f"rf-so2_decay-{k}.jpg" for k in ["abs", "norm"]],
        [_SAVE_DIR / f"temp-so2-{k}.jpg" for k in ["abs", "norm"]],
        [_SAVE_DIR / f"temp-rf-{k}.jpg" for k in ["abs", "norm"]],
    )
    for file in files:
        vdd.utils.combine(*file).in_grid(2, 1).save(
            _SAVE_DIR / ns(file[0].name.replace("-abs", ""))
        )
        for f in file:
            f.unlink()
    plt.show()


if __name__ == "__main__":
    _main()
