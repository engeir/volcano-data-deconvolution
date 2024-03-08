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
# 8. [x] Compute the difference between the peaks of the residuals and the reconstructed T.
#        Is it symmetric?

from functools import cached_property

import cosmoplots
import fppanalysis
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as ssi
import scipy.stats
import volcano_base
import xarray as xr

import vdd.load
import vdd.utils

plt.rc("text.latex", preamble=r"\usepackage{amsmath}")
plt.style.use(
    "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle"
)

ob16_month = volcano_base.load.OttoBliesner(freq="h0", progress=True)
dec_ob16 = vdd.load.DeconvolveOB16(data=ob16_month)
dec_ob16.name = "OB16 month"


class PlotReconstruction:
    """Plot the reconstruction of the data."""

    def __init__(
        self, ob16: volcano_base.load.OttoBliesner, reconstruction: vdd.load.Deconvolve
    ):
        self.ob16 = ob16
        dec_ob16 = vdd.load.DeconvolveOB16(data=ob16)
        dec_ob16.name = "OB16 month"
        self.dec_ob16 = dec_ob16
        self.reconstruction = reconstruction

    @cached_property
    def temp_control(self) -> xr.DataArray:
        """Temperature control."""
        temp_control = xr.align(self.ob16.temperature_control, self.dec_ob16.temp)[0]
        return temp_control

    @cached_property
    def rec_temp_so2(self) -> np.ndarray:
        """Reconstructed temperature from SO2."""
        return np.convolve(
            self.dec_ob16.so2.data, self.reconstruction.response_temp_so2, "same"
        )

    @cached_property
    def rec_temp_rf(self) -> np.ndarray:
        """Reconstructed temperature from radiative forcing."""
        return np.convolve(
            self.dec_ob16.rf.data, self.reconstruction.response_temp_rf, "same"
        )

    @cached_property
    def residual_so2(self) -> np.ndarray:
        """Compute the residual of the SO2."""
        return self.dec_ob16.temp.data - self.rec_temp_so2

    @cached_property
    def residual_rf(self) -> np.ndarray:
        """Compute the residual of the radiative forcing."""
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
        temp = aligned_arrays["temperature"]
        _idx_temp = np.argwhere(so2_temp.data > 0)
        peaks_original = temp.data[_idx_temp].flatten()
        peaks_time = temp.time.data[_idx_temp].flatten()
        peaks_so2 = self.rec_temp_so2[_idx_temp].flatten()
        peaks_rf = self.rec_temp_rf[_idx_temp].flatten()
        if view:
            aligned_arrays["so2-start"].plot()
            so2_temp.plot()
            temp.plot()
            plt.plot(temp.time, self.rec_temp_so2)
            plt.plot(temp.time, self.rec_temp_rf)
            plt.scatter(peaks_time, peaks_original)
            plt.scatter(peaks_time, peaks_so2)
            plt.scatter(peaks_time, peaks_rf)
            plt.show()
        return peaks_time, peaks_original, peaks_so2, peaks_rf

    def plot_reconstruction(self, fig: mpl.figure.Figure | None = None) -> None:
        """Plot the reconstruction of the data."""
        abso = plt.figure()
        abso_a = abso.gca()
        abso_a.set_ylabel("Absolute")
        norm = plt.figure()
        norm_a = norm.gca()
        norm_a.set_ylabel("Normalised")
        nor2 = plt.figure()
        nor2_a = nor2.gca()
        nor2_a.set_ylabel("Normalised to OB16 response amplitude")
        time_ = self.dec_ob16.temp.time
        so2 = self.dec_ob16.so2
        rf = self.dec_ob16.rf
        temp = self.dec_ob16.temp
        rt2 = self.dec_ob16.response_temp_so2.max()
        rtr = self.dec_ob16.response_temp_rf.max()
        abso_a.plot(time_, self.temp_control, label="OB16 control")
        abso_a.plot(time_, temp.data, label="OB16 month")
        norm_a.plot(time_, vdd.utils.normalise(temp.data), label="OB16 month")
        nor2_a.plot(time_, temp.data, label="OB16 month")
        rec = self.reconstruction
        conv_temp_so2 = np.convolve(so2.data, rec.response_temp_so2, "same")
        conv_temp_rf = np.convolve(rf.data, rec.response_temp_rf, "same")
        conv_norm_temp_so2 = np.convolve(
            so2.data,
            rec.response_temp_so2 / rec.response_temp_so2.max() * rt2,
            "same",
        )
        conv_norm_temp_rf = np.convolve(
            rf.data, rec.response_temp_rf / rec.response_temp_rf.max() * rtr, "same"
        )
        lso2 = f"From SO2 ({rec.name})"
        lrf = f"From RF ({rec.name})"
        abso_a.plot(time_, conv_temp_so2, label=lso2)
        abso_a.plot(time_, conv_temp_rf, label=lrf)
        norm_a.plot(time_, vdd.utils.normalise(conv_temp_so2), label=lso2)
        norm_a.plot(time_, vdd.utils.normalise(conv_temp_rf), label=lrf)
        nor2_a.plot(time_, conv_norm_temp_so2, label=lso2)
        nor2_a.plot(time_, conv_norm_temp_rf, label=lrf)
        abso_a.legend()
        norm_a.legend()
        nor2_a.legend()

    def correlation(self) -> None:
        """Compute the correlation between the residuals and the reconstructed T."""
        corr_so2_time, corr_so2 = fppanalysis.corr_fun(
            self.residual_so2, self.rec_temp_so2, 1 / 12
        )
        corr_rf_time, corr_rf = fppanalysis.corr_fun(
            self.residual_rf, self.rec_temp_rf, 1 / 12
        )
        plt.figure()
        plt.plot(corr_so2_time, corr_so2, label="SO2")
        plt.plot(corr_rf_time, corr_rf, label="RF")
        plt.legend()

    def _spectrum_1d(self, signal: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
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
            signal, sample_frequency, nperseg=2**13, return_onesided=False
        )
        frequency_plus = frequency[frequency > 0]
        power_plus = power[frequency > 0]
        return np.asarray(frequency_plus[1:]), np.asarray(power_plus[1:])

    def spectrum(self) -> None:
        """Compare the spectrum of the residuals and the control temperature."""
        f_so2, p_so2 = self._spectrum_1d(self.residual_so2)
        f_rf, p_rf = self._spectrum_1d(self.residual_rf)
        f_control, p_control = self._spectrum_1d(self.temp_control.data)
        plt.figure()
        plt.loglog(f_so2, p_so2, label="SO2")
        plt.loglog(f_rf, p_rf, label="RF")
        plt.loglog(f_control, p_control, label="Control")
        cosmoplots.change_log_axis_base(plt.gca(), "both")
        plt.xlabel("Frequency")
        plt.ylabel("Power spectral density")
        plt.legend()

    def peak_difference_analysis(self) -> None:
        """Plot the difference between the reconstructed and the original peaks."""
        # Generate a sample array of values
        # data = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        so2_basis = self.peaks_original - self.peaks_so2
        rf_basis = self.peaks_original - self.peaks_rf
        # Specify the value to test for symmetry
        test_value = 0
        # Perform a one-sample t-test to check for symmetry around the test value
        result_so2 = scipy.stats.ttest_1samp(so2_basis, test_value)
        t_statistic_so2, p_value_so2 = result_so2.statistic, result_so2.pvalue
        result_rf = scipy.stats.ttest_1samp(rf_basis, test_value)
        t_statistic_rf, p_value_rf = result_rf.statistic, result_rf.pvalue
        # Check if the p-value is less than a significance level (e.g., 0.05) to
        # determine symmetry (meaning a confidence level of 95%)
        sl = 0.05
        statement = f"symmetrically distributed around {test_value} (confidence of {int((1-sl)*100)}%)"
        if p_value_so2 < sl:
            print(f"Values are not {statement}")
        else:
            print(f"Values are {statement}")
        print(t_statistic_so2, p_value_so2)
        if p_value_rf < sl:
            print(f"Values are not {statement}")
        else:
            print(f"Values are {statement}")
        print(t_statistic_rf, p_value_rf)
        plt.figure()
        plt.xlabel("Difference between the peaks")
        plt.scatter(so2_basis, np.zeros_like(so2_basis), label="SO2")
        plt.scatter(rf_basis, np.zeros_like(rf_basis), label="RF")
        plt.figure()
        plt.ylabel("PDF")
        plt.xlabel("Difference between the peaks")
        pdf_so2, cdf_so2, bin_centers_so2 = fppanalysis.distribution(
            so2_basis, 30, ccdf=False
        )
        pdf, cdf, bin_centers = fppanalysis.distribution(rf_basis, 30, ccdf=False)
        plt.bar(bin_centers_so2, pdf_so2, width=0.01, label="SO2", alpha=0.5)
        plt.bar(bin_centers, pdf, width=0.01, label="RF", alpha=0.5)
        plt.legend()
        plt.figure()
        plt.ylabel("CDF")
        plt.xlabel("Difference between the peaks")
        plt.bar(bin_centers_so2, cdf_so2, width=0.01, label="SO2", alpha=0.5)
        plt.bar(bin_centers, cdf, width=0.01, label="RF", alpha=0.5)
        plt.legend()


if __name__ == "__main__":
    rec_class = PlotReconstruction(ob16_month, dec_ob16)
    rec_class.plot_reconstruction()
    # rec_class.peak_difference_analysis()
    rec_class.correlation()
    rec_class.spectrum()
    plt.show()
