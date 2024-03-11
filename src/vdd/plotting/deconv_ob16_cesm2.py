"""Plot the deconvolution comparison between OB16 and CESM2."""

import matplotlib as mpl
import matplotlib.pyplot as plt

import vdd.load
import vdd.utils

plt.rc("text.latex", preamble=r"\usepackage{amsmath}")
plt.style.use(
    "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle"
)

dec_cesm_s = vdd.load.DeconvolveCESM(pad_before=True)
dec_cesm_s.name = "CESM2 Strong"
dec_cesm_p = vdd.load.DeconvolveCESM(
    pad_before=True, cesm=vdd.load.CESMData(strength="medium-plus")
)
dec_cesm_p.name = "CESM2 Intermediate"
dec_cesm_m = vdd.load.DeconvolveCESM(
    pad_before=True, cesm=vdd.load.CESMData(strength="medium")
)
dec_cesm_m.name = "CESM2 Small"
dec_cesm_e = vdd.load.DeconvolveCESM(
    pad_before=True, cesm=vdd.load.CESMData(strength="size5000")
)
dec_cesm_e.name = "CESM2 Extra Strong"
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
all_decs = (dec_ob16, dec_ob16_month, dec_cesm_m, dec_cesm_p, dec_cesm_s, dec_cesm_e)


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
                fig.savefig(f"rf-so2{"-norm" if self.norm else ""}.png")
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
                fig.savefig(f"temp-so2{"-norm" if self.norm else ""}.png")
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
                fig.savefig(f"temp-rf{"-norm" if self.norm else ""}.png")
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
        for dec in self.decs:
            rf_so2_resp = (
                vdd.utils.normalise(dec.response_rf_so2)
                if self.norm
                else dec.response_rf_so2
            )
            temp_so2_resp = (
                vdd.utils.normalise(dec.response_temp_so2)
                if self.norm
                else dec.response_temp_so2
            )
            temp_rf_resp = (
                vdd.utils.normalise(dec.response_temp_rf)
                if self.norm
                else dec.response_temp_rf
            )
            rf_so2_a.plot(dec.tau, rf_so2_resp, label=dec.name)
            temp_so2_a.plot(dec.tau, temp_so2_resp, label=dec.name)
            temp_rf_a.plot(dec.tau, temp_rf_resp, label=dec.name)
        self.plot_rf_so2(rf_so2)
        self.plot_temp_so2(temp_so2)
        self.plot_temp_rf(temp_rf)


if __name__ == "__main__":
    PlotResponseFunctions(*all_decs, norm=True).run()
    plt.show()
    PlotResponseFunctions(*all_decs, norm=False).run()
    plt.show()
