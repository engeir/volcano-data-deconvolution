"""Plot the deconvolution comparison between OB16 and CESM2."""

import cosmoplots
import matplotlib as mpl
import matplotlib.pyplot as plt
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
    {"text.latex.preamble": r"\usepackage{amsmath}"},
])

DataCESM = vdd.load.CESMData
DecCESM = vdd.load.DeconvolveCESM
# CESM2
dec_cesm_4sep = DecCESM(pad_before=True, cesm=DataCESM(strength="double-overlap"))
dec_cesm_2sep = DecCESM(pad_before=True, cesm=DataCESM(strength="tt-2sep"))
dec_cesm_e = DecCESM(pad_before=True, cesm=DataCESM(strength="size5000"))
dec_cesm_s = DecCESM(pad_before=True, cesm=DataCESM(strength="strong"))
dec_cesm_p = DecCESM(pad_before=True, cesm=DataCESM(strength="medium-plus"))
dec_cesm_m = DecCESM(pad_before=True, cesm=DataCESM(strength="medium"))
# Original
dec_ob16 = vdd.load.DeconvolveOB16(data="h1")
dec_ob16_month = vdd.load.DeconvolveOB16(data="h0")
# Experimental (needs vdd.deconvolve_methods.alternative_deconv)
# ob16_day = volcano_base.load.OttoBliesner(freq="h1", progress=True)
# dec_ob16 = vdd.load.DeconvolveOB16(data=ob16_day)
# dec_ob16.change_deconvolution_method(alternative_deconv)
# ob16_month = volcano_base.load.OttoBliesner(freq="h0", progress=True)
# dec_ob16_month = vdd.load.DeconvolveOB16(data=ob16_month)
# dec_ob16_month.change_deconvolution_method(alternative_deconv)
dec_ob16.name = "OB16"
dec_ob16_month.name = "OB16 month"
all_decs = (
    # dec_ob16,
    dec_ob16_month,
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

    def plot_rf_so2(self, fig: mpl.figure.Figure | None = None) -> mpl.figure.Figure:
        """Plot the radiative forcing to SO2 response functions."""
        match fig:
            case None:
                return plt.figure()
            case mpl.figure.Figure():
                ax = fig.gca()
                rf_xlim = (-2, 10)
                ax.set_xlim(rf_xlim)
                ax.set_xlabel("Time lag ($\\tau$) [yr]")
                ax.set_ylabel(
                    f"{"Normalised " if self.norm else ""}RF to {"\n" if self.norm else ""}SO2 response [1]"
                )
                ax.legend()
                fig.savefig(_SAVE_DIR / f"rf-so2-{"norm" if self.norm else "abs"}.png")
                return fig
            case _:
                raise ValueError("rf must be a mpl.figure.Figure or None")

    def plot_temp_so2(self, fig: mpl.figure.Figure | None = None) -> mpl.figure.Figure:
        """Plot the temperature to SO2 response functions."""
        match fig:
            case None:
                return plt.figure()
            case mpl.figure.Figure():
                ax = fig.gca()
                temp_xlim = (-2, 20)
                ax.set_xlim(temp_xlim)
                ax.set_xlabel("Time lag ($\\tau$) [yr]")
                ax.set_ylabel(
                    f"{"Normalised t" if self.norm else "T"}emperature to {"\n" if self.norm else ""}SO2 response [1]"
                )
                ax.legend()
                fig.savefig(
                    _SAVE_DIR / f"temp-so2-{"norm" if self.norm else "abs"}.png"
                )
                return fig
            case _:
                raise ValueError("temp must be a mpl.figure.Figure or None")

    def plot_temp_rf(self, fig: mpl.figure.Figure | None = None) -> mpl.figure.Figure:
        """Plot the temperature to radiative forcing response functions."""
        match fig:
            case None:
                return plt.figure()
            case mpl.figure.Figure():
                ax = fig.gca()
                temp_xlim = (-2, 20)
                ax.set_xlim(temp_xlim)
                ax.set_xlabel("Time lag ($\\tau$) [yr]")
                ax.set_ylabel(
                    f"{"Normalised t" if self.norm else "T"}emperature to {"\n" if self.norm else ""}RF response [1]"
                )
                ax.legend()
                fig.savefig(_SAVE_DIR / f"temp-rf-{"norm" if self.norm else "abs"}.png")
                return fig
            case _:
                raise ValueError("temp must be a mpl.figure.Figure or None")

    def run(self) -> None:
        """Run the class."""
        rf_so2 = self.plot_rf_so2()
        temp_so2 = self.plot_temp_so2()
        temp_rf = self.plot_temp_rf()
        rf_so2_a = rf_so2.gca()
        temp_so2_a = temp_so2.gca()
        temp_rf_a = temp_rf.gca()
        # so2s = {
        #     "CESM2 strong": 1629,
        #     "CESM2 medium": 26,
        #     "CESM2 medium-plus": 400,
        #     "CESM2 double-overlap": 400,
        #     "CESM2 tt-2sep": 400,
        #     "CESM2 size5000": 3000,
        #     "OB16 month": 200,
        # }
        for dec in self.decs:
            rf_so2_resp = (
                vdd.utils.normalise(dec.response_rf_so2)
                # dec.response_rf_so2 / np.log(1+(400 / so2s[dec.name]))
                # dec.response_rf_so2 / (400 / so2s[dec.name])**(1/1.5)
                # dec.response_rf_so2 / (400 / so2s[dec.name])**(1/2)
                if self.norm
                else dec.response_rf_so2
            )
            temp_so2_resp = (
                vdd.utils.normalise(dec.response_temp_so2)
                # dec.response_rf_so2 / np.log(1+(400 / so2s[dec.name]))
                # dec.response_rf_so2 / (400 / so2s[dec.name])**(1/1.5)
                # dec.response_rf_so2 / (400 / so2s[dec.name])**(1/2)
                if self.norm
                else dec.response_temp_so2
            )
            temp_rf_resp = (
                vdd.utils.normalise(dec.response_temp_rf)
                if self.norm
                else dec.response_temp_rf
            )
            rf_so2_a.plot(dec.tau, rf_so2_resp, label=ns(dec.name))
            temp_so2_a.plot(dec.tau, temp_so2_resp, label=ns(dec.name))
            temp_rf_a.plot(dec.tau, temp_rf_resp, label=ns(dec.name))
        self.plot_rf_so2(rf_so2)
        self.plot_temp_so2(temp_so2)
        self.plot_temp_rf(temp_rf)


if __name__ == "__main__":
    PlotResponseFunctions(*all_decs, norm=True).run()
    PlotResponseFunctions(*all_decs, norm=False).run()
    files = (
        [_SAVE_DIR / f"rf-so2-{k}.png" for k in ["abs", "norm"]],
        [_SAVE_DIR / f"temp-so2-{k}.png" for k in ["abs", "norm"]],
        [_SAVE_DIR / f"temp-rf-{k}.png" for k in ["abs", "norm"]],
    )
    for file in files:
        cosmoplots.combine(*file).in_grid(2, 1).using(fontsize=50).save(
            _SAVE_DIR / ns(file[0].name.replace("-abs", ""))
        )
        for f in file:
            f.unlink()
    plt.show()
