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

import contextlib
import datetime
import warnings
from functools import cached_property
from typing import Literal, Self

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
from rich.console import Console
from rich.table import Table

import vdd.load
import vdd.utils
from vdd.load import Normalise, PaddingMethod
from vdd.utils import name_swap as ns

_SAVE_DIR = volcano_base.config.SAVE_PATH / "reconstruction"
if not _SAVE_DIR.exists():
    _SAVE_DIR.mkdir(parents=False)

plt.style.use(
    [
        "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle",
        "vdd.extra",
        "vdd.jgr",
    ],
)

COLORS = plt.rcParams["axes.prop_cycle"].by_key()["color"]

DataCESM = vdd.load.CESMData
DecCESM = vdd.load.DeconvolveCESM
# CESM2
padding = PaddingMethod.NOISE
dec_cesm_m = DecCESM(
    pad_before=padding, cesm=DataCESM(dims=["lat", "lon"], strength="medium")
)
# OB16
dec_ob16_month = vdd.load.DeconvolveOB16(data="h0")
dec_ob16_month.name = "OB16 month"
all_decs = (dec_cesm_m, dec_ob16_month)


class ReconstructOB16:
    """Class that reconstructs the temperature of OB16 from CESM2 simulations."""

    def __init__(self: Self, *decs: vdd.load.Deconvolve) -> None:
        self.ob16 = vdd.load.DeconvolveOB16(data="h0")
        self.ob16.name = "OB16 month"
        self.decs = decs

    def plot_temperature(self: Self) -> None:
        """Plot the reconstructed temperatures."""
        xlim = (
            vdd.utils.d2n(datetime.datetime(1250, 1, 1, 0, 0, tzinfo=datetime.UTC)),
            vdd.utils.d2n(datetime.datetime(1310, 1, 1, 0, 0, tzinfo=datetime.UTC)),
        )
        all_f = plt.figure()
        all_a = all_f.gca()
        all_a.plot(self.ob16.temp.time, self.ob16.temp, label=self.ob16.name)
        all_zoom_f = plt.figure()
        all_zoom_a = all_zoom_f.gca()
        all_zoom_a.plot(self.ob16.temp.time, self.ob16.temp, label=self.ob16.name)
        all_zoom_a.set_xlim(xlim)
        res: list[tuple[str, str, str]] = []
        for dec in self.decs:
            res = self._plot_temperature_single(dec, res, (all_a, all_zoom_a))
        table = Table(
            title="Difference between reconstructed temperature from OB16 and other simulations",
        )
        table.add_column("Simulation name", justify="left", style="cyan", no_wrap=True)
        table.add_column("Raw response", justify="center", style="magenta")
        table.add_column("Scaled response", justify="center", style="magenta")
        for r_ in res:
            table.add_row(*r_)
        console = Console()
        console.print(table)
        all_a.legend()
        all_f.savefig(_SAVE_DIR / "reconstruct_from_all.jpg")
        all_zoom_a.legend()
        all_zoom_f.savefig(_SAVE_DIR / "reconstruct_from_all_zoom.jpg")

    def _plot_temperature_single(
        self: Self,
        dec: vdd.load.Deconvolve,
        res: list[tuple[str, str, str]],
        axs: tuple[mpl.axes.Axes, mpl.axes.Axes],
    ) -> list[tuple[str, str, str]]:
        """Plot the reconstructed temperature for a single simulation."""
        xlim = (
            vdd.utils.d2n(datetime.datetime(1250, 1, 1, 0, 0, tzinfo=datetime.UTC)),
            vdd.utils.d2n(datetime.datetime(1310, 1, 1, 0, 0, tzinfo=datetime.UTC)),
        )
        fn = ns(vdd.utils.clean_filename(dec.name))
        inv_f = plt.figure()
        inv_a = inv_f.gca()
        inv_zoom_f = plt.figure()
        inv_zoom_a = inv_zoom_f.gca()
        inv_zoom_a.set_xlim(xlim)
        response = dec.response_temp_so2
        response_scaled = response / response.max() * self.ob16.response_temp_so2.max()
        new_temp = np.convolve(self.ob16.so2, response, mode="same")
        new_temp_scaled = np.convolve(self.ob16.so2, response_scaled, mode="same")
        axs[0].plot(self.ob16.temp.time, new_temp_scaled, label=ns(dec.name))
        axs[1].plot(self.ob16.temp.time, new_temp_scaled, label=ns(dec.name))
        inv_a.plot(self.ob16.temp.time, self.ob16.temp, label="OB16 temperature")
        inv_a.plot(self.ob16.temp.time, new_temp_scaled, label="Scaled response")
        inv_a.legend()
        inv_f.savefig(_SAVE_DIR / f"reconstruct_from_{fn}.jpg")
        inv_zoom_a.plot(self.ob16.temp.time, self.ob16.temp, label="OB16 temperature")
        inv_zoom_a.plot(self.ob16.temp.time, new_temp_scaled, label="Scaled response")
        inv_zoom_a.legend()
        inv_zoom_f.savefig(_SAVE_DIR / f"reconstruct_from_{fn}_zoom.jpg")
        # Print the distance away from the reconstructed
        ob16_temp = np.convolve(self.ob16.so2, self.ob16.response_temp_so2, mode="same")
        ob16_diff = np.abs(ob16_temp - new_temp).sum()
        ob16_diff_scaled = np.abs(ob16_temp - new_temp_scaled).sum()
        res.append((ns(dec.name), f"{ob16_diff:.2f}", f"{ob16_diff_scaled:.2f}"))
        return res


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
    """

    def __init__(
        self: Self,
        ob16: volcano_base.load.OttoBliesner,
        reconstruction: vdd.load.Reconstructor,
    ) -> None:
        self.ob16 = ob16
        self.normalise = reconstruction.normalise
        self.dec_ob16 = vdd.load.DeconvolveOB16(
            normalise=self.normalise,
            data=ob16,
            length=12001,
        )
        self.dec_ob16.name = "OB16 month"
        self.reconstruction = reconstruction
        self.reconstruction.name = vdd.utils.name_swap(self.reconstruction.name)
        self.sim_name = vdd.utils.name_swap(
            vdd.utils.clean_filename(reconstruction.name),
        )

    @cached_property
    def temp_control(self: Self) -> xr.DataArray:
        """Temperature control."""
        temp_control = xr.align(self.ob16.temperature_control, self.dec_ob16.temp)[0]
        match self.normalise:
            case Normalise.MEAN_STD:
                temp_control = (temp_control - temp_control.mean()) / temp_control.std()
        return temp_control

    @cached_property
    def rec_temp_so2(self: Self) -> np.ndarray:
        """Reconstructed temperature from SO2."""
        response = self.reconstruction.response_temp_so2
        response[: len(response) // 2] = 0
        match self.normalise:
            case Normalise.NO:
                response = (
                    response / response.max() * self.dec_ob16.response_temp_so2.max()
                )
        convolved = np.convolve(self.dec_ob16.so2.data, response, "same")
        match self.normalise:
            case Normalise.MEAN_STD:
                convolved = (convolved - convolved.mean()) / convolved.std()
        return convolved

    @cached_property
    def rec_temp_rf(self: Self) -> np.ndarray:
        """Reconstructed temperature from radiative forcing."""
        response = self.reconstruction.response_temp_rf
        response[: len(response) // 2] = 0
        match self.normalise:
            case Normalise.NO:
                response = (
                    response / response.max() * self.dec_ob16.response_temp_rf.max()
                )
        convolved = np.convolve(self.dec_ob16.rf.data, response, "same")
        match self.normalise:
            case Normalise.MEAN_STD:
                convolved = (convolved - convolved.mean()) / convolved.std()
        return convolved

    @cached_property
    def residual_so2(self: Self) -> np.ndarray:
        """Compute the residual of the SO2.

        Returns
        -------
        np.ndarray
            Residual of the SO2.
        """
        return -1 * self.dec_ob16.temp.data - -1 * self.rec_temp_so2

    @cached_property
    def residual_rf(self: Self) -> np.ndarray:
        """Compute the residual of the radiative forcing.

        Returns
        -------
        np.ndarray
            Residual of the radiative forcing.
        """
        return -1 * self.dec_ob16.temp.data - -1 * self.rec_temp_rf

    @cached_property
    def _peaks_tup(self: Self) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Compute the peaks of the residual of the SO2."""
        return self._find_peaks()

    @property
    def peaks_time(self: Self) -> np.ndarray:
        """Time array for the peaks.

        Returns
        -------
        np.ndarray
            Time array for the peaks.
        """
        return self._peaks_tup[0]

    @property
    def peaks_original(self: Self) -> np.ndarray:
        """Peaks of the original temperature.

        Returns
        -------
        np.ndarray
            Peaks of the original temperature.
        """
        return self._peaks_tup[1]

    @property
    def peaks_so2(self: Self) -> np.ndarray:
        """Peaks of the temperature from SO2.

        Returns
        -------
        np.ndarray
            Peaks of the temperature from SO2.
        """
        return self._peaks_tup[2]

    @property
    def peaks_rf(self: Self) -> np.ndarray:
        """Peaks of the temperature from radiative forcing.

        Returns
        -------
        np.ndarray
            Peaks of the temperature from radiative forcing.
        """
        return self._peaks_tup[3]

    def _find_peaks(
        self: Self,
        *,
        view: bool = False,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Compute the peaks of the OB16 temperature."""
        aligned_arrays = self.ob16.aligned_arrays
        s, e = self.dec_ob16.start_pt, self.dec_ob16.end_pt
        so2_temp = aligned_arrays["so2-temperature"][s:e]
        _idx_temp = np.argwhere(so2_temp.data > 0)
        peaks_original = self.dec_ob16.temp.data[_idx_temp].flatten()
        peaks_time = self.dec_ob16.temp.time.data[_idx_temp].flatten()
        # We could add 'peak_shift' to '_idx_temp' to account for the later peak of the
        # CESM2 pulse functions, but since we already scale the pulse function, this
        # will clearly yield a result very similar to OB16, and be of little interest.
        peaks_so2 = self.rec_temp_so2[_idx_temp].flatten()
        peaks_rf = self.rec_temp_rf[_idx_temp].flatten()
        if view:
            aligned_arrays["so2-start"][s:e].plot()
            so2_temp.plot()
            self.dec_ob16.temp.plot()  # type: ignore[call-arg]
            plt.plot(self.dec_ob16.temp.time, self.rec_temp_so2)
            plt.plot(self.dec_ob16.temp.time, self.rec_temp_rf)
            plt.scatter(peaks_time, peaks_original)
            plt.scatter(peaks_time, peaks_so2)
            plt.scatter(peaks_time, peaks_rf)
            plt.show()
        return peaks_time, peaks_original, peaks_so2, peaks_rf

    def _get_name(self: Self) -> str:
        match self.reconstruction.name:
            case name_ if "26" in name_.lower():
                name = "S26"
            case name_ if "400" in name_.lower():
                name = "S400"
            case name_ if "ob16" in name_.lower():
                name = "OB16"
            case _:
                name = "OTHER"
        return name

    def plot_reconstruction_temp(self: Self, ax: mpl.axes.Axes) -> mpl.axes.Axes:
        """Plot the reconstruction of the data."""
        ax.set_xlabel("Time [yr]")
        ax.set_ylabel("$T$ [K]")
        time_ = self.dec_ob16.temp.time
        temp = self.dec_ob16.temp
        aligned_arrays = self.ob16.aligned_arrays
        s, e = self.dec_ob16.start_pt, self.dec_ob16.end_pt
        so2_temp = aligned_arrays["so2-temperature"][s:e]
        _idx_temp = np.argwhere(so2_temp.data > 0)
        peaks_time = self.dec_ob16.temp.time.data[_idx_temp].flatten()
        ax.axvline(
            peaks_time[0],
            lw=0.5,
            c="k",
            label="$s_{k,\\text{OB16}}+10\\,\\mathrm{months}$",
        )
        [ax.axvline(x_, c="k", lw=0.5) for x_ in peaks_time[1:]]
        ax.plot(
            time_,
            -1 * self.temp_control,
            c=COLORS[1],
            label="$T_{\\text{CONTROL}}$",
        )
        ax.plot(time_, -1 * temp.data, c=COLORS[0], label="$T_{\\text{OB16}}$")
        l_so2 = (
            f"$\\varphi_T^{{\\text{{{self._get_name()}}}}}\\ast S_{{\\text{{OB16}}}}$"
        )
        ax.plot(time_, -1 * self.rec_temp_so2, c=COLORS[2], label=l_so2)
        xlim = (
            vdd.utils.d2n(datetime.datetime(1250, 1, 1, 0, 0, tzinfo=datetime.UTC)),
            vdd.utils.d2n(datetime.datetime(1310, 1, 1, 0, 0, tzinfo=datetime.UTC)),
        )
        ax.set_xlim(xlim)
        ax.set_ylim((-2.1, 0.6))
        ax.legend(framealpha=1.0)
        return ax

    def correlation(self: Self, ax: mpl.axes.Axes) -> mpl.axes.Axes:
        """Compute the correlation between the residuals and temperature."""
        corr_self_time, corr_self = fppanalysis.corr_fun(
            -1 * self.dec_ob16.temp.data,
            -1 * self.dec_ob16.temp.data,
            1 / 12,
        )
        corr_so2_time, corr_so2 = fppanalysis.corr_fun(
            -1 * self.dec_ob16.temp.data,
            -1 * self.rec_temp_so2,
            1 / 12,
            # self.residual_so2, self.dec_ob16.temp.data, 1 / 12
        )
        corr_ctrl_time, corr_ctrl = fppanalysis.corr_fun(
            -1 * self.dec_ob16.temp.data,
            -1 * self.dec_ob16.temp_control,
            1 / 12,
        )
        rprint(f"[bold]Lag 0 correlation[/bold]: {np.max(corr_self) =:.4f}")
        rprint(f"[bold]Lag 0 correlation[/bold]: {np.max(corr_so2) =:.4f}")
        rprint(f"[bold]Lag 0 correlation[/bold]: {np.max(corr_ctrl) =:.4f}")
        ax.plot(corr_ctrl_time, corr_ctrl, c=COLORS[1], label="$T_{\\text{CONTROL}}$")
        ax.plot(corr_self_time, corr_self, c=COLORS[0], label="$T_{\\text{OB16}}$")
        l_so2 = (
            f"$\\varphi_T^{{\\text{{{self._get_name()}}}}}\\ast S_{{\\text{{OB16}}}}$"
        )
        ax.plot(corr_so2_time, corr_so2, c=COLORS[2], label=l_so2)
        ax.set_xlim((-100, 100))
        ax.set_xlabel("Time lag [yr]")
        ax.set_ylabel("Correlation with \n$T_{\\text{OB16}}$")
        ax.legend()
        return ax

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
            signal,
            sample_frequency,
            nperseg=2**11,
            return_onesided=False,
        )
        frequency_plus = frequency[frequency > 0]
        power_plus = power[frequency > 0]
        return np.asarray(frequency_plus[1:]), np.asarray(power_plus[1:])

    def spectrum(self: Self, ax: mpl.axes.Axes) -> mpl.axes.Axes:
        """Compare the spectrum of the residuals and the control temperature."""
        f_so2_res, p_so2_res = self._spectrum_1d(self.residual_so2)
        f_so2_rec, p_so2_rec = self._spectrum_1d(-1 * self.rec_temp_so2)
        f_control, p_control = self._spectrum_1d(-1 * self.temp_control.data)
        f_orig, p_orig = self._spectrum_1d(-1 * self.dec_ob16.temp.data)
        f_so2, p_so2 = self._spectrum_1d(self.reconstruction.response_temp_so2)
        ax.plot(
            f_control,
            p_control,
            c=COLORS[1],
            label="$T_{\\text{CONTROL}}$",
        )
        ax.plot(f_orig, p_orig, c=COLORS[0], label="$T_{\\text{OB16}}$")
        lrec = (
            f"$\\varphi_T^{{\\text{{{self._get_name()}}}}}\\ast S_{{\\text{{OB16}}}}$"
        )
        ax.plot(f_so2_rec, p_so2_rec, c=COLORS[2], label=lrec)
        lresidual = f"$T_{{\\text{{OB16}}}}-\\varphi_T^{{\\text{{{self._get_name()}}}}}\\ast S_{{\\text{{OB16}}}}$"
        ax.plot(f_so2_res, p_so2_res, c=COLORS[4], label=lresidual)
        ax.plot(
            f_so2,
            p_so2,
            c=COLORS[3],
            label=f"$\\varphi_{{T}}^{{\\text{{{self._get_name()}}}}}$",
        )
        return self._plot_power(ax, "Frequency [yr$^{-1}$]")

    def spectrum_parts(self: Self, ax: mpl.axes.Axes) -> mpl.axes.Axes:
        """View the spectrum of the response functions and the input data."""
        f_so2, p_so2 = self._spectrum_1d(self.reconstruction.response_temp_so2)
        f_control, p_control = self._spectrum_1d(self.temp_control.data)
        f_orig_so2, p_orig_so2 = self._spectrum_1d(self.dec_ob16.so2.data)
        ax.plot(
            f_so2,
            p_so2,
            label=f"$\\varphi_{{T}}^{{\\text{{{self._get_name()}}}}}$",
            alpha=0.5,
        )
        ax.plot(f_control, p_control, label="$T_{\\text{CONTROL}}$", alpha=0.5)
        ax.plot(f_orig_so2, p_orig_so2, label="SO2 TS", alpha=0.5)
        return self._plot_power(ax, "Frequency")

    @staticmethod
    def _plot_power(ax: mpl.axes.Axes, xlab: str) -> mpl.axes.Axes:
        warnings.filterwarnings("ignore")
        cosmoplots.change_log_axis_base(ax, "both")
        warnings.resetwarnings()
        ax.set_xlabel(xlab)
        ax.set_ylabel("Power spectral density")
        ax.legend()
        return ax

    def _peak_difference_ttest(
        self: Self,
        so2_basis: np.ndarray,
        rf_basis: np.ndarray,
    ) -> tuple[str, str]:
        # Specify the value to test for symmetry
        test_value = 0
        # Perform a one-sample t-test to check for symmetry around the test value
        result_so2 = scipy.stats.ttest_1samp(so2_basis, popmean=test_value)
        t_statistic_so2, p_value_so2 = result_so2.statistic, result_so2.pvalue
        result_rf = scipy.stats.ttest_1samp(rf_basis, popmean=test_value)
        t_statistic_rf, p_value_rf = result_rf.statistic, result_rf.pvalue

        sl = 0.05  # Significance level

        def info(name: str, p_value: float) -> None:
            color = "red" if p_value < sl else "green"
            rprint(
                f"[blue][bold]{self.sim_name}/{name}[/bold]: I can with [/blue][{color}]"
                f"{(1 - p_value) * 100:.4f}% confidence[/{color}][blue] say that the "
                f"distribution does not have a mean of {test_value}[/blue]",
            )

        # Check if the p-value is less than a significance level (e.g., 0.05) to
        # determine symmetry (meaning a confidence level of 95%)
        reject = f"The distribution does not have a mean of {test_value} (confidence of {int((1 - sl) * 100)}%)"
        reject_no = f"I cannot with at least {int((1 - sl) * 100)}% confidence say that the distribution does not have a mean of {test_value}"
        if p_value_so2 < sl:
            rprint(reject)
        else:
            rprint(reject_no)
        rprint(t_statistic_so2, p_value_so2)
        info("SO2", p_value_so2)
        if p_value_rf < sl:
            rprint(reject)
        else:
            rprint(reject_no)
        rprint(t_statistic_rf, p_value_rf)
        info("CTRL", p_value_rf)
        return p_value_so2, p_value_rf

    def peak_difference_analysis(
        self: Self,
        ax1: mpl.axes.Axes,
        ax2: mpl.axes.Axes,
    ) -> tuple[mpl.axes.Axes, mpl.axes.Axes]:
        """Plot the difference between the reconstructed and the original peaks."""
        self._diff_of_max_peak()
        so2_basis = -1 * (self.peaks_so2 - self.peaks_original)
        ctrl_basis = -1 * self.temp_control
        so2_conf, ctrl_conf = self._peak_difference_ttest(so2_basis, ctrl_basis)  # type: ignore[arg-type]
        pdf_so2, cdf_so2, bin_centers_so2 = fppanalysis.distribution(
            so2_basis,
            30,
            ccdf=False,
        )
        pdf_ctrl, cdf_ctrl, bin_centers_ctrl = fppanalysis.distribution(
            ctrl_basis,
            30,
            ccdf=False,
        )
        stats_so2 = scipy.stats.describe(so2_basis)
        stats_ctrl = scipy.stats.describe(ctrl_basis)
        fit_so2 = scipy.stats.norm.fit(so2_basis)
        fit_ctrl = scipy.stats.norm.fit(ctrl_basis)
        dist_so2 = scipy.stats.skewnorm(
            a=stats_so2.skewness,
            loc=stats_so2.mean,
            scale=np.sqrt(stats_so2.variance),
        )
        dist_ctrl = scipy.stats.skewnorm(
            a=stats_ctrl.skewness,
            loc=stats_ctrl.mean,
            scale=np.sqrt(stats_ctrl.variance),
        )
        axpdf = self._peak_difference_plot(
            ax1,
            (bin_centers_so2, bin_centers_ctrl, pdf_so2, pdf_ctrl),
            ((fit_so2, fit_ctrl), (dist_so2, dist_ctrl)),
            "pdf",
            txt=(so2_conf, ctrl_conf),
        )
        axcdf = self._peak_difference_plot(
            ax2,
            (bin_centers_so2, bin_centers_ctrl, cdf_so2, cdf_ctrl),
            ((fit_so2, fit_ctrl), (dist_so2, dist_ctrl)),
            "cdf",
            txt=(so2_conf, ctrl_conf),
        )
        return axpdf, axcdf

    def _peak_difference_plot(
        self: Self,
        ax: mpl.axes.Axes,
        fpp_out: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
        fits: tuple[tuple[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]],
        dist: Literal["pdf", "cdf"],
        txt: tuple[str, str],
    ) -> mpl.axes.Axes:
        b_so2, b_rf, f_s, f_r = fpp_out
        norm_fit, _ = fits
        norm_so2 = getattr(scipy.stats.norm, dist)(b_so2, *norm_fit[0])
        norm_rf = getattr(scipy.stats.norm, dist)(b_rf, *norm_fit[1])
        kw = {"width": 0.01, "alpha": 1.0}
        ax3 = plt.figure().gca()
        ax3.bar(
            b_rf,
            f_r,
            color=COLORS[1],
            label=f"$T_{{\\text{{CONTROL}}}}$ (p-value: {txt[1]:.4f})",
            **kw,  # type: ignore[arg-type]
        )
        lresidual = f"$T_{{\\text{{OB16}}}}-\\varphi_T^{{\\text{{{self._get_name()}}}}}\\ast S_{{\\text{{OB16}}}}$"
        ax3.bar(
            b_so2,
            f_s,
            color=COLORS[4],
            label=f"{lresidual} (p-value: {txt[0]:.4f})",
            **kw,  # type: ignore[arg-type]
        )
        bar_hand, bar_lab = ax3.get_legend_handles_labels()
        # Norm
        ax3.plot(b_rf, norm_rf, c=COLORS[1], label="_Norm CTRL")
        ax3.plot(b_so2, norm_so2, c=COLORS[4], label="_Norm SO2")
        ax3.legend(bar_hand, bar_lab, loc="upper left", framealpha=0.5)
        # Make the plot symmetric around 0
        xlim = np.abs(ax3.get_xlim()).max()
        ax3.set_xlim((-xlim, xlim))
        ax3.set_ylabel(dist.upper())
        ax3.set_xlabel("$T$ [K]")
        lspread = f"$T_{{\\text{{OB16}}}}(\\varphi_T^{{\\text{{{self._get_name()}}}}}\\ast S_{{\\text{{OB16}}}})$"
        norm = 10
        shift_x = 0.15
        lims = (-1.55, 0.6)
        std = self.temp_control.data.std() / 2
        ax.annotate(
            "",
            xy=(lims[1], lims[1]),
            xytext=(lims[0], lims[0]),
            arrowprops={
                "shrinkA": 0,
                "shrinkB": 0,
                "arrowstyle": "->",
                "lw": 0.7,
                "color": "grey",
            },
        )
        ax.fill(
            [lims[0] - std, lims[1] - std, lims[1] + std, lims[0] + std],
            [lims[0] + std, lims[1] + std, lims[1] - std, lims[0] - std],
            "grey",
            alpha=0.3,
        )
        ax.fill(
            [
                lims[0] - 2 * std,
                lims[1] - 2 * std,
                lims[1] + 2 * std,
                lims[0] + 2 * std,
            ],
            [
                lims[0] + 2 * std,
                lims[1] + 2 * std,
                lims[1] - 2 * std,
                lims[0] - 2 * std,
            ],
            "grey",
            alpha=0.3,
        )
        ax.annotate(
            r"$\Delta T$ [K]",
            xy=(0.4 + shift_x, -0.4 + shift_x),
            xytext=(-0.45 + shift_x, 0.45 + shift_x),
            arrowprops={
                "shrinkA": 0,
                "shrinkB": 0,
                "lw": 0.7,
                "relpos": (1, 0),
                "color": "grey",
                "arrowstyle": "<-",
            },
            c="grey",
            ha="right",
            va="bottom",
            size=6,
        )
        bc_s = b_rf / 2
        bc_sx, bc_sy = bc_s + shift_x, bc_s - shift_x
        bs_s = b_so2 / 2
        bs_sx, bs_sy = bs_s + shift_x, bs_s - shift_x
        ax.scatter(bc_sx + f_r / norm, f_r / norm - bc_sy, s=3, c=COLORS[1], zorder=5)
        ax.scatter(bs_sx + f_s / norm, f_s / norm - bs_sy, s=3, c=COLORS[4], zorder=5)
        ax.plot(bc_sx + norm_rf / norm, norm_rf / norm - bc_sy, c=COLORS[1])
        ax.plot(bs_sx + norm_so2 / norm, norm_so2 / norm - bs_sy, c=COLORS[4])
        ax.scatter(
            -self.peaks_so2,
            -self.peaks_original,
            marker="*",
            c=COLORS[4],
            s=12,
            zorder=10,
            label=lspread,
        )
        ax.set_ylabel("$T_{{\\text{{OB16}}}}$ peak [K]")
        ax.set_xlabel(
            f"$\\varphi_T^{{\\text{{{self._get_name()}}}}}\\ast S_{{\\text{{OB16}}}}$ peak [K]",
        )
        ax.legend(loc="lower right")
        return ax

    def _diff_of_max_peak(self: Self) -> None:
        """Print the difference between the peaks of the original and the reconstructed."""
        peaks_so2 = scipy.signal.savgol_filter(self.rec_temp_so2, 12, 3)
        peaks_orig = scipy.signal.savgol_filter(self.dec_ob16.temp.data, 12, 3)
        reasonable_diff = peaks_so2.max() - peaks_orig.max()
        idx = np.argmax(self.peaks_so2)
        diff_max = (self.peaks_so2 - self.peaks_original)[idx]
        time_max = self.peaks_time[idx]
        rprint(
            f"[bold]{self._get_name()}: Difference at the largest volcanic eruption at {time_max}[/bold]: "
            f"{diff_max:.4f} K, but the more reasonable difference is {reasonable_diff:.4f} K",
        )


class PlotManyReconstructions:
    """Wrapper class that creates any number of reconstructions."""

    def __init__(
        self: Self,
        ob16: volcano_base.load.OttoBliesner,
        *recs: vdd.load.Reconstructor,
    ) -> None:
        self.ob16 = ob16
        self.recs = recs

    def run(self: Self) -> None:
        """Run the reconstructions."""
        for rec in self.recs:
            rec_class = PlotReconstruction(self.ob16, rec)
            rec_class.peak_difference_analysis(plt.figure().gca(), plt.figure().gca())
            rec_class.correlation(plt.figure().gca())
            rec_class.spectrum(plt.figure().gca())
            rec_class.spectrum_parts(plt.figure().gca())
            rec_class.plot_reconstruction_temp(plt.figure().gca())
            plt.show()
            plt.close("all")


def _plot_reconstructed_temperature() -> None:
    # Makes most sense when comparing all CESM2 simulations against each other, with
    # OB16 as the control.
    rec = ReconstructOB16(*all_decs)
    rec.plot_temperature()
    plt.show()


def _plot_many_reconstructions() -> None:
    ob16 = dec_ob16_month.data
    rec_ob16_ = dec_ob16_month.dump_reconstructor()
    rec_small_ = dec_cesm_m.dump_reconstructor()
    zero_like = 1e-4
    valid_until = 7
    rec_ob16_.response_temp_so2[
        np.intersect1d(
            np.argwhere(rec_ob16_.response_temp_so2 < zero_like).flatten(),
            np.argwhere(rec_ob16_.tau > valid_until).flatten(),
        ).flatten()[0] :
    ] = 0
    with contextlib.suppress(Exception):
        # In case we never reach zero, for example when running with the large eruption
        # simulation output
        rec_small_.response_temp_so2[
            np.intersect1d(
                np.argwhere(rec_small_.response_temp_so2 < zero_like).flatten(),
                np.argwhere(rec_small_.tau > valid_until).flatten(),
            ).flatten()[0] :
        ] = 0
    rec_ob16 = PlotReconstruction(ob16, rec_ob16_)
    rec_small = PlotReconstruction(ob16, rec_small_)
    _plot_individual(rec_ob16, rec_small)
    # _plot_all(rec_ob16, rec_small)


def _plot_individual(
    rec_ob16: PlotReconstruction,
    rec_small: PlotReconstruction,
) -> None:
    figpdf, axspdf = cosmoplots.figure_multiple_rows_columns(1, 2)
    figcdf, axscdf = cosmoplots.figure_multiple_rows_columns(1, 2)
    rec_ob16.peak_difference_analysis(*[axspdf[0], axscdf[0]])
    rec_small.peak_difference_analysis(*[axspdf[1], axscdf[1]])
    figpdf.savefig(_SAVE_DIR / "compare-historical-size-peak-difference-pdf")
    figcdf.savefig(_SAVE_DIR / "compare-historical-size-peak-difference-cdf")
    figcorr, axcorr = cosmoplots.figure_multiple_rows_columns(1, 2)
    rec_ob16.correlation(axcorr[0])
    rec_small.correlation(axcorr[1])
    figcorr.savefig(
        _SAVE_DIR / "compare-historical-size-correlation-residual-reconstructed",
    )
    figsp, axsp = cosmoplots.figure_multiple_rows_columns(1, 2)
    rec_ob16.spectrum(axsp[0])
    rec_small.spectrum(axsp[1])
    figsp.savefig(_SAVE_DIR / "compare-historical-size-spectrum-residual-control_temp")
    figrec, axrec = cosmoplots.figure_multiple_rows_columns(1, 2)
    rec_ob16.plot_reconstruction_temp(axrec[0])
    rec_small.plot_reconstruction_temp(axrec[1])
    figrec.savefig(_SAVE_DIR / "compare-historical-size-temp-reconstructed")


def _plot_all(rec_ob16: PlotReconstruction, rec_small: PlotReconstruction) -> None:
    # All in one
    figtot, axtot = cosmoplots.figure_multiple_rows_columns(4, 2)
    rec_ob16.plot_reconstruction_temp(axtot[0])
    rec_small.plot_reconstruction_temp(axtot[1])
    rec_ob16.peak_difference_analysis(*[axtot[2], plt.figure().gca()])
    rec_small.peak_difference_analysis(*[axtot[3], plt.figure().gca()])
    rec_ob16.spectrum(axtot[4])
    rec_small.spectrum(axtot[5])
    rec_ob16.correlation(axtot[6])
    rec_small.correlation(axtot[7])
    figtot.savefig(_SAVE_DIR / "compare-historical-size-all-in-one")


def _main() -> None:
    # _plot_reconstructed_temperature()
    _plot_many_reconstructions()
    plt.show()


if __name__ == "__main__":
    _main()
