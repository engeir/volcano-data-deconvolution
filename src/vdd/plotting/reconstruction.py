"""Plotting functions for the reconstruction of the data."""

# Steps:
#
# 1. [x] Load OB16 monthly data. (Since deconvolution of T against RF is problematic.)
# 2. [x] Deconvolve T against SO2 and T against RF.
# 3. [x] Reconstruct T by convolving with SO2 and RF.
# 4. [x] Get hold of OB16 control run data for T.
# 5. [x] Compute residuals.
# 6. [x] Compute the spectrum of the residuals and the control run data.
# 7. [x] Look at the correlation between the residuals and the reconstructed T.
# 8. [x] Compute the difference between the peaks of the residuals and the reconstructed
#        T. Is it symmetric?
# 9. [x] Compare against CESM2 responses as well, but proper normalisation must be done.
#        This could be to use OB16 response amplitude as the true amplitude, to let the
#        OB16 response sum (layman's integral) be the true sum, or to normalise every
#        array we come across. Try all. (Amplitude seems to work best.)

import pathlib
import warnings
from functools import cached_property
from typing import Literal

import cosmoplots
import fppanalysis
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as ssi
import scipy.stats
import volcano_base
import xarray as xr
from rich import print as rprint

import vdd.load
import vdd.utils

_SAVE_DIR = volcano_base.config.SAVE_PATH / "reconstruction"
if not _SAVE_DIR.exists():
    _SAVE_DIR.mkdir(parents=False)

plt.rc("text.latex", preamble=r"\usepackage{amsmath}")
plt.style.use(
    "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle"
)


def _setup() -> (
    tuple[volcano_base.load.OttoBliesner, tuple[vdd.load.Reconstructor, ...]]
):
    # Used as a reference for the reconstruction.
    ob16_month = volcano_base.load.OttoBliesner(freq="h0", progress=True)
    dec_ob16 = vdd.load.DeconvolveOB16(data=ob16_month)
    dec_ob16.name = "OB16 month"
    recs: tuple[vdd.load.Reconstructor, ...] = ()
    # OB16 month -----------------------------------------------------------------------
    rec_ob16 = vdd.load.DeconvolveOB16(normalise=False, data=ob16_month)
    rec_ob16.name = "OB16 month"
    recs += (rec_ob16.dump_reconstructor(),)
    # OB16 cut offs --------------------------------------------------------------------
    co_ob16_temp_so2 = vdd.load.CutOff(dec_ob16, ("temp", "so2"))
    co_ob16_temp_rf = vdd.load.CutOff(dec_ob16, ("temp", "rf"))
    meta_years = {12 * i for i in [2, 4, 8, 16]}
    for co in (co_ob16_temp_so2, co_ob16_temp_rf):
        co.cut_off(meta_years)
    for meta_year_ in meta_years:
        meta_year = int(meta_year_)
        rec_co = co_ob16_temp_so2.dump_reconstructor(
            cut=meta_year, temp_rf=co_ob16_temp_rf.cuts[str(meta_year)].response.data
        )
        num = str(meta_year // 12)
        num = "0" * (2 - len(num)) + num
        rec_co.name = f"cut off {num} years"
        recs += (rec_co,)
    # CESM2 ----------------------------------------------------------------------------
    for cesm_ in (
        "tt-2sep",
        "medium",
        "medium-plus",
        "size5000",
        "double-overlap",
        "strong",
    ):
        cesm = vdd.load.CESMData(strength=cesm_)  # type: ignore
        rec_cesm = vdd.load.DeconvolveCESM(normalise=False, pad_before=True, cesm=cesm)
        recs += (rec_cesm.dump_reconstructor(),)
    return ob16_month, recs


class PlotReconstruction:
    """Plot the reconstruction of the data.

    Parameters
    ----------
    ob16 : volcano_base.load.OttoBliesner
        OB16 monthly data. This is what we compare against and what we try to
        reconstruct.
    reconstruction : vdd.load.Reconstructor
        The reconstruction object. We use the response functions from this object to
        reconstruct the OB16 data via convolution of response functions.

    Attributes
    ----------
    ob16 : volcano_base.load.OttoBliesner
        OB16 monthly data.
    normalise : bool
        Whether to normalise the data or not.
    dec_ob16 : vdd.load.DeconvolveOB16
        Deconvolved OB16 data.
    reconstruction : vdd.load.Reconstructor
        The reconstruction object.
    sim_name : pathlib.Path
        Name of the simulation extracted from the reconstruction object, used for saving
        figures.
    """

    def __init__(
        self,
        ob16: volcano_base.load.OttoBliesner,
        reconstruction: vdd.load.Reconstructor,
    ):
        self.ob16 = ob16
        self.normalise = reconstruction.normalise
        dec_ob16 = vdd.load.DeconvolveOB16(normalise=self.normalise, data=ob16)
        dec_ob16.name = "OB16 month"
        self.dec_ob16 = dec_ob16
        self.reconstruction = reconstruction
        self.sim_name = vdd.utils.name_swap(
            vdd.utils.clean_filename(reconstruction.name)
        )

    @cached_property
    def temp_control(self) -> xr.DataArray:
        """Temperature control."""
        temp_control = xr.align(self.ob16.temperature_control, self.dec_ob16.temp)[0]
        if self.normalise:
            temp_control = (temp_control - temp_control.mean()) / temp_control.std()
        return temp_control

    @cached_property
    def rec_temp_so2(self) -> np.ndarray:
        """Reconstructed temperature from SO2."""
        response = self.reconstruction.response_temp_so2
        response[: len(response) // 2] = 0
        if not self.normalise:
            # response = response / np.sum(response) * np.sum(self.dec_ob16.response_temp_so2)
            response = response / response.max() * self.dec_ob16.response_temp_so2.max()
        convolved = np.convolve(self.dec_ob16.so2.data, response, "same")
        if self.normalise:
            # convolved = vdd.utils.normalise(convolved)
            convolved = (convolved - convolved.mean()) / convolved.std()
        return convolved

    @cached_property
    def rec_temp_rf(self) -> np.ndarray:
        """Reconstructed temperature from radiative forcing."""
        response = self.reconstruction.response_temp_rf
        response[: len(response) // 2] = 0
        if not self.normalise:
            # response = response / np.sum(response) * np.sum(self.dec_ob16.response_temp_rf)
            response = response / response.max() * self.dec_ob16.response_temp_rf.max()
        convolved = np.convolve(self.dec_ob16.rf.data, response, "same")
        if self.normalise:
            # convolved = vdd.utils.normalise(convolved)
            convolved = (convolved - convolved.mean()) / convolved.std()
        return convolved

    @cached_property
    def residual_so2(self) -> np.ndarray:
        """Compute the residual of the SO2.

        Returns
        -------
        np.ndarray
            Residual of the SO2.
        """
        return self.dec_ob16.temp.data - self.rec_temp_so2

    @cached_property
    def residual_rf(self) -> np.ndarray:
        """Compute the residual of the radiative forcing.

        Returns
        -------
        np.ndarray
            Residual of the radiative forcing.
        """
        return self.dec_ob16.temp.data - self.rec_temp_rf

    @cached_property
    def _peaks_tup(self) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Compute the peaks of the residual of the SO2."""
        return self._find_peaks()

    @property
    def peaks_time(self) -> np.ndarray:
        """Time array for the peaks.

        Returns
        -------
        np.ndarray
            Time array for the peaks.
        """
        return self._peaks_tup[0]

    @property
    def peaks_original(self) -> np.ndarray:
        """Peaks of the original temperature.

        Returns
        -------
        np.ndarray
            Peaks of the original temperature.
        """
        return self._peaks_tup[1]

    @property
    def peaks_so2(self) -> np.ndarray:
        """Peaks of the temperature from SO2.

        Returns
        -------
        np.ndarray
            Peaks of the temperature from SO2.
        """
        return self._peaks_tup[2]

    @property
    def peaks_rf(self) -> np.ndarray:
        """Peaks of the temperature from radiative forcing.

        Returns
        -------
        np.ndarray
            Peaks of the temperature from radiative forcing.
        """
        return self._peaks_tup[3]

    def _find_peaks(
        self, view: bool = False
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Compute the peaks of the OB16 temperature."""
        aligned_arrays = self.ob16.aligned_arrays
        so2_temp = aligned_arrays["so2-temperature"]
        _idx_temp = np.argwhere(so2_temp.data > 0)
        peaks_original = self.dec_ob16.temp.data[_idx_temp].flatten()
        peaks_time = self.dec_ob16.temp.time.data[_idx_temp].flatten()
        peaks_so2 = self.rec_temp_so2[_idx_temp].flatten()
        peaks_rf = self.rec_temp_rf[_idx_temp].flatten()
        if view:
            aligned_arrays["so2-start"].plot()
            so2_temp.plot()
            self.dec_ob16.temp.plot()
            plt.plot(self.dec_ob16.temp.time, self.rec_temp_so2)
            plt.plot(self.dec_ob16.temp.time, self.rec_temp_rf)
            plt.scatter(peaks_time, peaks_original)
            plt.scatter(peaks_time, peaks_so2)
            plt.scatter(peaks_time, peaks_rf)
            plt.show()
        return peaks_time, peaks_original, peaks_so2, peaks_rf

    def plot_reconstruction_temp(self, fig: mpl.figure.Figure | None = None) -> None:
        """Plot the reconstruction of the data."""
        abso = plt.figure()
        abso_a = abso.gca()
        abso_a.set_xlabel("Time [yr]")
        abso_a.set_ylabel("Absolute")
        # norm = plt.figure()
        # norm_a = norm.gca()
        # norm_a.set_ylabel("Normalised")
        # nor2 = plt.figure()
        # nor2_a = nor2.gca()
        # nor2_a.set_ylabel("Normalised to OB16 response amplitude")
        time_ = self.dec_ob16.temp.time
        temp = self.dec_ob16.temp
        # rt2 = self.dec_ob16.response_temp_so2.max()
        # rtr = self.dec_ob16.response_temp_rf.max()
        abso_a.plot(time_, self.temp_control, label="OB16 control")
        abso_a.plot(time_, temp.data, label="OB16 month")
        # norm_a.plot(time_, vdd.utils.normalise(temp.data), label="OB16 month")
        # nor2_a.plot(time_, temp.data, label="OB16 month")
        rec = self.reconstruction
        # conv_norm_temp_so2 = np.convolve(
        #     so2.data,
        #     rec.response_temp_so2 / rec.response_temp_so2.max() * rt2,
        #     "same",
        # )
        # conv_norm_temp_rf = np.convolve(
        #     rf.data, rec.response_temp_rf / rec.response_temp_rf.max() * rtr, "same"
        # )
        lso2 = f"SO2 ({rec.name})"
        lrf = f"RF ({rec.name})"
        abso_a.plot(time_, self.rec_temp_so2, label=lso2)
        abso_a.plot(time_, self.rec_temp_rf, label=lrf)
        # norm_a.plot(time_, vdd.utils.normalise(self.rec_temp_so2), label=lso2)
        # norm_a.plot(time_, vdd.utils.normalise(self.rec_temp_rf), label=lrf)
        # nor2_a.plot(time_, conv_norm_temp_so2, label=lso2)
        # nor2_a.plot(time_, conv_norm_temp_rf, label=lrf)
        abso_a.set_xlim((-790 * 365, -650 * 365))
        abso_a.legend(framealpha=0.5)
        # norm_a.legend()
        # nor2_a.legend()
        plt.savefig(_SAVE_DIR / f"{self.sim_name}-temp-reconstructed.png")

    def correlation(self) -> None:
        """Compute the correlation between the residuals and temperature."""
        corr_so2_time, corr_so2 = fppanalysis.corr_fun(
            self.residual_so2, self.dec_ob16.temp.data, 1 / 12
        )
        corr_rf_time, corr_rf = fppanalysis.corr_fun(
            self.residual_rf, self.dec_ob16.temp.data, 1 / 12
        )
        plt.figure()
        plt.plot(corr_so2_time, corr_so2, label="SO2", alpha=0.7)
        plt.plot(corr_rf_time, corr_rf, label="RF", alpha=0.7)
        plt.xlabel("Time lag ($\\tau$) [yr]")
        plt.ylabel("Correlation between residual \nand original temperature")
        plt.legend()
        plt.savefig(
            _SAVE_DIR / f"{self.sim_name}-correlation-residual-reconstructed.png"
        )

    @staticmethod
    def _spectrum_1d(signal: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Calculate the one sided spectrum of the signal.

        Parameters
        ----------
        signal : np.ndarray
            Signal to calculate the spectrum of.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            Frequency and power of the signal.
        """
        signal = (signal - signal.mean()) / signal.std()
        sample_frequency = 12
        frequency, power = ssi.welch(
            signal, sample_frequency, nperseg=2**11, return_onesided=False
        )
        frequency_plus = frequency[frequency > 0]
        power_plus = power[frequency > 0]
        return np.asarray(frequency_plus[1:]), np.asarray(power_plus[1:])

    def spectrum(self) -> None:
        """Compare the spectrum of the residuals and the control temperature."""
        f_so2, p_so2 = self._spectrum_1d(self.residual_so2)
        f_rf, p_rf = self._spectrum_1d(self.residual_rf)
        f_control, p_control = self._spectrum_1d(self.temp_control.data)
        f_orig, p_orig = self._spectrum_1d(self.dec_ob16.temp.data)
        plt.figure()
        plt.plot(f_so2, p_so2, label="SO2", alpha=0.5)
        plt.plot(f_rf, p_rf, label="RF", alpha=0.5)
        plt.plot(f_control, p_control, label="Control", alpha=0.5)
        plt.plot(f_orig, p_orig, label="OB16", alpha=0.5)
        # Suppress the warning
        warnings.filterwarnings("ignore")
        cosmoplots.change_log_axis_base(plt.gca(), "both")
        warnings.resetwarnings()
        plt.xlabel("Frequency")
        plt.ylabel("Power spectral density")
        plt.legend()
        plt.savefig(_SAVE_DIR / f"{self.sim_name}-spectrum-residual-control_temp.png")

    def spectrum_parts(self) -> None:
        """View the spectrum of the response functions and the input data."""
        f_so2, p_so2 = self._spectrum_1d(self.reconstruction.response_temp_so2)
        f_rf, p_rf = self._spectrum_1d(self.reconstruction.response_temp_rf)
        f_control, p_control = self._spectrum_1d(self.temp_control.data)
        f_orig_so2, p_orig_so2 = self._spectrum_1d(self.dec_ob16.so2.data)
        f_orig_rf, p_orig_rf = self._spectrum_1d(self.dec_ob16.rf.data)
        plt.figure()
        plt.plot(f_so2, p_so2, label="$\\phi_{SO2}$", alpha=0.5)
        plt.plot(f_rf, p_rf, label="$\\phi_{RF}$", alpha=0.5)
        plt.plot(f_control, p_control, label="Control", alpha=0.5)
        plt.plot(f_orig_so2, p_orig_so2, label="SO2", alpha=0.5)
        plt.plot(f_orig_rf, p_orig_rf, label="RF", alpha=0.5)
        # Suppress the warning
        warnings.filterwarnings("ignore")
        cosmoplots.change_log_axis_base(plt.gca(), "both")
        warnings.resetwarnings()
        plt.xlabel("Frequency")
        plt.ylabel("Power spectral density")
        plt.legend()
        plt.savefig(_SAVE_DIR / f"{self.sim_name}-spectrum-response-input.png")

    @staticmethod
    def _peak_difference_ttest(
        so2_basis: np.ndarray, rf_basis: np.ndarray
    ) -> tuple[str, str]:
        # Specify the value to test for symmetry
        test_value = 0
        # Perform a one-sample t-test to check for symmetry around the test value
        result_so2 = scipy.stats.ttest_1samp(so2_basis, popmean=test_value)
        t_statistic_so2, p_value_so2 = result_so2.statistic, result_so2.pvalue
        result_rf = scipy.stats.ttest_1samp(rf_basis, popmean=test_value)
        t_statistic_rf, p_value_rf = result_rf.statistic, result_rf.pvalue

        def info(name, p_value) -> None:
            rprint(
                f"[blue][bold]{name}[/bold]: I can with [/blue][red]"
                f"{(1 - p_value) * 100:.4f}% confidence[/red][blue] say that the "
                f"distribution does not have a mean of {test_value}[/blue]"
            )

        # Check if the p-value is less than a significance level (e.g., 0.05) to
        # determine symmetry (meaning a confidence level of 95%)
        sl = 0.01
        reject = f"The distribution does not have a mean of {test_value} (confidence of {int((1 - sl) * 100)}%)"
        reject_no = f"I cannot with at least {int((1 - sl) * 100)}% confidence say that the distribution does not have a mean of {test_value}"
        if p_value_so2 < sl:
            print(reject)
        else:
            print(reject_no)
        print(t_statistic_so2, p_value_so2)
        info("SO2", p_value_so2)
        if p_value_rf < sl:
            print(reject)
        else:
            print(reject_no)
        print(t_statistic_rf, p_value_rf)
        info("RF", p_value_rf)
        return p_value_so2, p_value_rf

    def peak_difference_analysis(self) -> None:  # noqa:PLR0914
        """Plot the difference between the reconstructed and the original peaks."""
        so2_basis = self.peaks_so2 - self.peaks_original
        rf_basis = self.peaks_rf - self.peaks_original
        so2_conf, rf_conf = self._peak_difference_ttest(so2_basis, rf_basis)
        pdf_so2, cdf_so2, bin_centers_so2 = fppanalysis.distribution(
            so2_basis, 30, ccdf=False
        )
        pdf_rf, cdf_rf, bin_centers_rf = fppanalysis.distribution(
            rf_basis, 30, ccdf=False
        )
        stats_so2 = scipy.stats.describe(so2_basis)
        stats_rf = scipy.stats.describe(rf_basis)
        fit_so2 = scipy.stats.norm.fit(so2_basis)
        fit_rf = scipy.stats.norm.fit(rf_basis)
        dist_so2 = scipy.stats.skewnorm(
            a=stats_so2.skewness, loc=stats_so2.mean, scale=np.sqrt(stats_so2.variance)
        )
        dist_rf = scipy.stats.skewnorm(
            a=stats_rf.skewness, loc=stats_rf.mean, scale=np.sqrt(stats_rf.variance)
        )
        self._peak_difference_plot(
            (bin_centers_so2, bin_centers_rf, pdf_so2, pdf_rf),
            ((fit_so2, fit_rf), (dist_so2, dist_rf)),
            "pdf",
            txt=(so2_conf, rf_conf),
        )
        self._peak_difference_plot(
            (bin_centers_so2, bin_centers_rf, cdf_so2, cdf_rf),
            ((fit_so2, fit_rf), (dist_so2, dist_rf)),
            "cdf",
            txt=(so2_conf, rf_conf),
        )

    def _peak_difference_plot(  # noqa:PLR0914
        self,
        fpp_out: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
        fits: tuple[tuple[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]],
        dist: Literal["pdf", "cdf"],
        txt: tuple[str, str],
    ) -> None:
        b_so2, b_rf, f_s, f_r = fpp_out
        norm_fit, skewnorm_fit = fits
        prop_cycle = plt.rcParams["axes.prop_cycle"]
        colors = prop_cycle.by_key()["color"]
        norm_so2 = getattr(scipy.stats.norm, dist)(b_so2, *norm_fit[0])
        norm_rf = getattr(scipy.stats.norm, dist)(b_rf, *norm_fit[1])
        skewnorm_so2 = getattr(skewnorm_fit[0], dist)(b_so2)
        skewnorm_rf = getattr(skewnorm_fit[1], dist)(b_rf)
        plt.figure()
        ax = plt.gca()
        ax.bar(b_so2, f_s, width=0.01, label=f"SO2 (p-value: {txt[0]:.4f})", alpha=0.5)
        ax.bar(b_rf, f_r, width=0.01, label=f"RF (p-value: {txt[1]:.4f})", alpha=0.5)
        bar_hand, bar_lab = ax.get_legend_handles_labels()
        # Norm
        (norm,) = ax.plot(b_so2, norm_so2, c=colors[0], label="_Norm SO2", alpha=0.5)
        ax.plot(b_rf, norm_rf, c=colors[1], label="_Norm RF", alpha=0.5)
        # Skewnorm
        (skewnorm,) = plt.plot(
            b_so2, skewnorm_so2, "--", c=colors[0], label="_Skewnorm SO2", alpha=0.5
        )
        ax.plot(b_rf, skewnorm_rf, "--", c=colors[1], label="_Skewnorm RF", alpha=0.5)
        bar_legend = ax.legend(bar_hand, bar_lab, loc="upper left", framealpha=0.5)
        norm_loc = "center left" if dist == "cdf" else "upper right"
        norm_legend = ax.legend(
            [norm, skewnorm], ["Norm", "Skewnorm"], loc=norm_loc, framealpha=0.5
        )
        norm_legend.legend_handles[0].set_color("black")  # type: ignore
        norm_legend.legend_handles[1].set_color("black")  # type: ignore
        ax.add_artist(bar_legend)
        ax.add_artist(norm_legend)
        # Make the plot symmetric around 0
        xlim = np.abs(plt.gca().get_xlim()).max()
        plt.xlim((-xlim, xlim))
        plt.ylabel(dist.upper())
        plt.xlabel("Difference between the peaks")
        plt.savefig(_SAVE_DIR / f"{self.sim_name}-peak-difference-{dist}.png")


class PlotManyReconstructions:
    """Wrapper class that creates any number of reconstructions."""

    def __init__(
        self, ob16: volcano_base.load.OttoBliesner, *recs: vdd.load.Reconstructor
    ):
        self.ob16 = ob16
        self.recs = recs

    def run(self) -> None:
        """Run the reconstructions."""
        for rec in self.recs:
            rec_class = PlotReconstruction(self.ob16, rec)
            rec_class.plot_reconstruction_temp()
            rec_class.peak_difference_analysis()
            rec_class.correlation()
            rec_class.spectrum()
            rec_class.spectrum_parts()
            plt.close("all")


if __name__ == "__main__":
    ob16, recs = _setup()
    plot_recs = PlotManyReconstructions(ob16, *recs)
    plot_recs.run()
    # We take extra care of the cut-off data.
    files = pathlib.Path(_SAVE_DIR).glob("cut-off*")
    temp = []
    pdf = []
    cdf = []
    power = []
    correlation = []
    for f in files:
        match f.name:
            case _ if "temp-reconstructed" in f.name:
                temp.append(f)
            case _ if "peak-difference-pdf" in f.name:
                pdf.append(f)
            case _ if "peak-difference-cdf" in f.name:
                cdf.append(f)
            case _ if "spectrum-residual-control_temp" in f.name:
                power.append(f)
            case _ if "correlation-residual-reconstructed" in f.name:
                correlation.append(f)
            case _:
                pass
    temp.sort()
    pdf.sort()
    cdf.sort()
    power.sort()
    correlation.sort()
    for list_, name in zip(
        (temp, pdf, cdf, power, correlation),
        ("temp", "pdf", "cdf", "power", "correlation"),
        strict=True,
    ):
        cosmoplots.combine(*list_).in_grid(2, len(list_) // 2).using(fontsize=50).save(
            _SAVE_DIR / f"cut-off-{name}-combined.png"
        )
        for f in list_:
            f.unlink()
