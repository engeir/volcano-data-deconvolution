"""Functions that fetches data from files and returns it as an xarray.DataArray."""

from collections.abc import Callable
from functools import cached_property
from typing import Any, Literal

import cftime
import fppanalysis
import matplotlib.pyplot as plt
import numpy as np
import volcano_base
import xarray as xr

import vdd.utils


def _convert_time_start_zero(arr: xr.DataArray) -> xr.DataArray:
    """Convert the time to start from zero."""
    arr = arr.assign_coords(time=volcano_base.manipulate.dt2float(arr.time.data) - 1850)
    arr.coords["time"].attrs["long_name"] = "Time-after-eruption"
    arr.coords["time"].attrs["units"] = "yr"
    return arr


class CESMData:
    """Class to store the CESM2 data."""

    def __init__(
        self,
        strength: Literal["strong", "medium", "medium-plus", "size5000"] = "strong",
    ) -> None:
        self.strength = strength

    @cached_property
    def so2(self) -> xr.DataArray:
        """Get the CESM2 SO2 data.

        Returns
        -------
        xr.DataArray
            The SO2 data.
        """
        x = self.rf.time.data
        y = np.zeros_like(x)
        match self.strength:
            case "strong":
                y[0] = 1629
            case "medium-plus":
                y[0] = 400
            case "medium":
                y[0] = 26
            case "size5000":
                y[0] = 3000
            case _:
                vdd.utils.never_called(self.strength)
        return xr.DataArray(y, coords={"time": x}, dims=["time"])

    @cached_property
    def aod(self) -> xr.DataArray:
        """Get the CESM2 AOD data.

        Returns
        -------
        xr.DataArray
            The AOD data.
        """
        return self._get_aod_cesm()

    @cached_property
    def rf(self) -> xr.DataArray:
        """Get the CESM2 RF data.

        Returns
        -------
        xr.DataArray
            The RF data.
        """
        return self._get_rf_cesm()

    @cached_property
    def temp(self) -> xr.DataArray:
        """Get the CESM2 temperature data.

        Returns
        -------
        xr.DataArray
            The temperature data.
        """
        return self._get_trefht_cesm()

    def _get_aod_cesm(self, plot_example: bool = False) -> xr.DataArray:
        """Get the CESM2 AOD data."""
        data = (
            volcano_base.load.FindFiles()
            .find("e_fSST1850", "AODVISstdn", "h0", self.strength)
            .remove("ens1")
            .keep_most_recent()
        )
        files = data.load()
        files = volcano_base.manipulate.mean_flatten(files, dims=["lat", "lon"])
        files = volcano_base.manipulate.shift_arrays(files, daily=False)
        files = volcano_base.manipulate.shift_arrays(files, custom=1)
        files = volcano_base.manipulate.subtract_mean_of_tail(files)
        files = volcano_base.manipulate.data_array_operation(
            files, _convert_time_start_zero
        )
        files = list(xr.align(*files))
        if plot_example:
            plt.figure()
            [f.plot() for f in files]
            plt.show()
        mean_array = volcano_base.manipulate.get_median(files, xarray=True)
        if plot_example:
            plt.figure()
            mean_array.plot()
            plt.show()
        return mean_array.dropna("time")

    def _get_rf_cesm(self, plot_example: bool = False) -> xr.DataArray:
        """Get the CESM2 RF data."""
        control: volcano_base.load.FindFiles = (
            volcano_base.load.FindFiles()
            .find("e_fSST1850", "h0", "control", {"FLNT", "FSNT"}, "ens1")
            .sort("ensemble", "attr")
        )
        data: volcano_base.load.FindFiles = (
            volcano_base.load.FindFiles()
            .find("e_fSST1850", {"FLNT", "FSNT"}, "h0", self.strength)
            .remove("ens1")
            .keep_most_recent()
            .sort("ensemble", "attr")
        )
        c_size = 1
        match (control.copy().keep("FLNT"), control.copy().keep("FSNT")):
            case (c_flnt, c_fsnt) if len(c_flnt) == c_size and len(c_fsnt) == c_size:
                c_flnt_xr = c_flnt.load()[0]
                c_fsnt_xr = c_fsnt.load()[0]
            case _:
                raise ValueError("Control data not found.")
        d_size = 2 if self.strength == "size5000" else 4
        match (data.copy().keep("FLNT"), data.copy().keep("FSNT")):
            case (flnt, fsnt) if len(flnt) == d_size and len(fsnt) == d_size:
                flnt_xr = flnt.load()
                fsnt_xr = fsnt.load()
            case _:
                raise ValueError("Data not found.")
        c_fsnt_xr, c_flnt_xr = volcano_base.manipulate.mean_flatten(
            [c_fsnt_xr, c_flnt_xr], dims=["lat", "lon"]
        )
        fsnt_xr = volcano_base.manipulate.mean_flatten(fsnt_xr, dims=["lat", "lon"])
        flnt_xr = volcano_base.manipulate.mean_flatten(flnt_xr, dims=["lat", "lon"])

        def remove_control(rf: xr.DataArray) -> xr.DataArray:
            rf, c_s, c_l = xr.align(rf, c_fsnt_xr, c_flnt_xr)
            rf.data = rf.data - (c_s.data - c_l.data)
            return rf.assign_attrs(attr="RF")

        def difference(
            short: list[xr.DataArray], long: list[xr.DataArray]
        ) -> list[xr.DataArray]:
            rf_: list[xr.DataArray] = []
            rf_ap = rf_.append
            for s, ell in zip(short, long, strict=True):
                s_, long_ = xr.align(s, ell)
                s_.data -= long_.data
                rf_ap(s_)
            return volcano_base.manipulate.data_array_operation(rf_, remove_control)

        rf = difference(fsnt_xr, flnt_xr)
        rf = volcano_base.manipulate.shift_arrays(rf, daily=False)
        rf = volcano_base.manipulate.shift_arrays(rf, custom=1)
        rf = volcano_base.manipulate.subtract_mean_of_tail(rf)
        rf = volcano_base.manipulate.data_array_operation(rf, _convert_time_start_zero)
        rf = list(xr.align(*rf))
        if plot_example:
            plt.figure()
            [f.plot() for f in rf]
            plt.show()
        mean_array = volcano_base.manipulate.get_median(rf, xarray=True)
        if plot_example:
            plt.figure()
            plt.plot(mean_array)
            plt.show()
        return mean_array.dropna("time")

    def _get_trefht_cesm(self, plot_example: bool = False) -> xr.DataArray:
        """Get the CESM2 temperature data."""
        data = (
            volcano_base.load.FindFiles()
            .find("e_BWma1850", "TREFHT", "h0", self.strength)
            .remove("ens1")
            .keep_most_recent()
        )
        files = data.load()
        files = volcano_base.manipulate.mean_flatten(files, dims=["lat", "lon"])
        files = volcano_base.manipulate.shift_arrays(files, daily=False)
        files = volcano_base.manipulate.shift_arrays(files, custom=1)

        def subtract_control_mean(arr: xr.DataArray) -> xr.DataArray:
            return arr - volcano_base.config.MEANS["TREFHT"]

        files = volcano_base.manipulate.data_array_operation(
            files, subtract_control_mean
        )
        files = volcano_base.manipulate.data_array_operation(
            files, _convert_time_start_zero
        )
        files = list(xr.align(*files))
        if plot_example:
            plt.figure()
            [f.plot() for f in files]
            plt.show()
        mean_array = volcano_base.manipulate.get_median(files, xarray=True)
        if plot_example:
            plt.figure()
            plt.plot(mean_array)
            plt.show()
        return mean_array.dropna("time")


class Deconvolve:
    """Class for deconvolving data."""

    tau: np.ndarray
    so2: xr.DataArray
    rf: xr.DataArray
    temp: xr.DataArray
    name: str = "Deconvolve"

    def __init__(self) -> None:
        kwargs = {"iteration_list": 200}
        self._deconv_method: Callable[
            [np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]
        ] = lambda signal, forcing: fppanalysis.RL_gauss_deconvolve(
            signal, forcing, **kwargs
        )

    def change_deconvolution_method(
        self,
        method: Callable[[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]],
        *args: Any,
        **kwargs: Any,
    ) -> None:
        """Change the deconvolution method.

        Parameters
        ----------
        method : Callable[[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]]
            The new deconvolution method.
        *args : Any
            The positional arguments to pass to the method, if any. Not including the
            signal and forcing.
        **kwargs : Any
            The keyword arguments to pass to the method, if any.
        """

        def method_(
            signal: np.ndarray, forcing: np.ndarray
        ) -> tuple[np.ndarray, np.ndarray]:
            return method(signal, forcing, *args, **kwargs)

        self._deconv_method = method_

    @cached_property
    def _response_rf_so2_tup(self) -> tuple[np.ndarray, np.ndarray]:
        """Deconvolve the RF signal with the SO2 signal."""
        return self._deconv_method(self.rf.data, self.so2.data)

    @property
    def response_rf_so2(self) -> np.ndarray:
        """Deconvolve the RF signal with the SO2 signal.

        Returns
        -------
        np.ndarray
            The deconvolved RF to SO2 signal.
        """
        return self._response_rf_so2_tup[0].flatten()

    @property
    def response_rf_so2_err(self) -> np.ndarray:
        """Deconvolve the RF signal with the SO2 signal."""
        return self._response_rf_so2_tup[1]

    @cached_property
    def _response_temp_so2_tup(self) -> tuple[np.ndarray, np.ndarray]:
        """Deconvolve the temperature signal with the SO2 signal."""
        return self._deconv_method(self.temp.data, self.so2.data)

    @property
    def response_temp_so2(self) -> np.ndarray:
        """Deconvolve the temperature signal with the SO2 signal.

        Returns
        -------
        np.ndarray
            The deconvolved temperature to SO2 signal.
        """
        return self._response_temp_so2_tup[0].flatten()

    @property
    def response_temp_so2_err(self) -> np.ndarray:
        """Deconvolve the temperature signal with the SO2 signal."""
        return self._response_temp_so2_tup[1]

    @cached_property
    def _response_temp_rf_tup(self) -> tuple[np.ndarray, np.ndarray]:
        """Deconvolve the temperature signal with the RF signal."""
        return self._deconv_method(self.temp.data, self.rf.data)

    @property
    def response_temp_rf(self) -> np.ndarray:
        """Deconvolve the temperature signal with the RF signal.

        Returns
        -------
        np.ndarray
            The deconvolved temperature to RF signal.
        """
        return self._response_temp_rf_tup[0].flatten()

    @property
    def response_temp_rf_err(self) -> np.ndarray:
        """Deconvolve the temperature signal with the RF signal."""
        return self._response_temp_rf_tup[1]

    def plot_dec_rf_with_so2(self) -> None:
        """Deconvolve the RF signal with the SO2 signal."""
        # Quick comparison
        plt.figure()
        vdd.utils.normalise(self.rf).plot()
        vdd.utils.normalise(self.so2).plot()
        plt.figure()
        plt.plot(self.tau, self.response_rf_so2)
        plt.figure()
        plt.semilogy(self.response_rf_so2_err)

    def plot_dec_temp_with_so2(self) -> None:
        """Deconvolve the temperature signal with the SO2 signal."""
        # Quick comparison
        plt.figure()
        vdd.utils.normalise(self.temp).plot()
        vdd.utils.normalise(self.so2).plot()
        plt.figure()
        plt.plot(self.tau, self.response_temp_so2)
        plt.figure()
        plt.semilogy(self.response_temp_so2_err)


class DeconvolveCESM(Deconvolve):
    """Class for deconvolving data from CESM2.

    Parameters
    ----------
    pad_before : bool, optional
        Whether to pad with zeros all time series before the convolution, by default
        False
    cesm : CESMData | None, optional
        The CESM data to use. If not given, the "strong" strength data will be used.
    """

    def __init__(self, pad_before: bool = False, cesm: CESMData | None = None) -> None:
        super().__init__()
        cesm_ = CESMData() if cesm is None else cesm
        self.so2 = (
            vdd.utils.pad_before_convolution(cesm_.so2) if pad_before else cesm_.so2
        )
        self.tau = (
            self.so2.time.data
            if pad_before
            else self.so2.time.data - self.so2.time.data[len(self.so2.time.data) // 2]
        )
        self.aod = (
            vdd.utils.pad_before_convolution(cesm_.aod) if pad_before else cesm_.aod
        )
        self.rf = (
            vdd.utils.pad_before_convolution(cesm_.rf * -1)
            if pad_before
            else cesm_.rf * -1
        )
        self.temp = (
            vdd.utils.pad_before_convolution(cesm_.temp * -1)
            if pad_before
            else cesm_.temp * -1
        )

    @cached_property
    def _response_rf_aod_tup(self) -> tuple[np.ndarray, np.ndarray]:
        """Deconvolve the RF signal with the AOD signal."""
        return self._deconv_method(self.rf.data, self.aod.data)

    @property
    def response_rf_aod(self) -> np.ndarray:
        """Deconvolve the RF signal with the AOD signal.

        Returns
        -------
        np.ndarray
            The deconvolved RF to AOD signal.
        """
        return self._response_rf_aod_tup[0].flatten()

    @property
    def response_rf_aod_err(self) -> np.ndarray:
        """Deconvolve the RF signal with the AOD signal."""
        return self._response_rf_aod_tup[1]

    @cached_property
    def _response_temp_aod_tup(self) -> tuple[np.ndarray, np.ndarray]:
        """Deconvolve the temperature signal with the AOD signal."""
        return self._deconv_method(self.temp.data, self.aod.data)

    @property
    def response_temp_aod(self) -> np.ndarray:
        """Deconvolve the temperature signal with the AOD signal.

        Returns
        -------
        np.ndarray
            The deconvolved temperature to AOD signal.
        """
        return self._response_temp_aod_tup[0].flatten()

    @property
    def response_temp_aod_err(self) -> np.ndarray:
        """Deconvolve the temperature signal with the AOD signal."""
        return self._response_temp_aod_tup[1]


class DeconvolveOB16(Deconvolve):
    """Class for deconvolving data from Otto-Bliesner et at. (2016).

    Parameters
    ----------
    data : volcano_base.load.OttoBliesner | Literal["h0", "h1"], optional
        The OB16 data class to use. If not given, the h1 (daily) data will be loaded.

    Raises
    ------
    ValueError
        If the data is invalid.
    """

    def __init__(
        self, data: volcano_base.load.OttoBliesner | Literal["h0", "h1"] = "h1"
    ) -> None:
        super().__init__()
        match data:
            case "h0":
                ob16 = volcano_base.load.OttoBliesner(freq="h0", progress=True)
            case "h1":
                ob16 = volcano_base.load.OttoBliesner(freq="h1", progress=True)
            case volcano_base.load.OttoBliesner():
                ob16 = data
            case _:
                raise ValueError(f"Invalid data: {data}")
        self.so2 = ob16.aligned_arrays["so2-start"]
        tau = self.so2.time.data - (
            self.so2.time.data[len(self.so2.time.data) // 2]
            - cftime.DatetimeNoLeap(0, 1, 1, has_year_zero=True, calendar="noleap")
        )
        self.tau = volcano_base.manipulate.dt2float(tau)
        self.rf = ob16.aligned_arrays["rf"]
        self.temp = ob16.aligned_arrays["temperature"]
