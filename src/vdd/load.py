"""Functions that fetches data from files and returns it as an xarray.DataArray."""

import warnings
from abc import ABC, abstractmethod, abstractproperty
from collections.abc import Callable, Iterable
from functools import cached_property
from typing import Any, Literal, Self

import cftime
import cosmoplots
import fppanalysis
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy
import scipy.signal as ssi
import volcano_base
import xarray as xr
from pydantic import BaseModel, Field
from rich import print as rprint

import vdd.utils

type T_RF = tuple[Literal["temp"], Literal["rf"]]  # type: ignore
type T_SO2 = tuple[Literal["temp"], Literal["so2"]]  # type: ignore
type RF_SO2 = tuple[Literal["rf"], Literal["so2"]]  # type: ignore
type T_Strengths = Literal[  # type: ignore
    "strong", "medium", "medium-plus", "size5000", "tt-2sep", "double-overlap"
]


def _convert_time_start_zero(arr: xr.DataArray) -> xr.DataArray:
    """Convert the time to start from zero."""
    arr = arr.assign_coords(time=volcano_base.manipulate.dt2float(arr.time.data) - 1850)
    arr.coords["time"].attrs["long_name"] = "Time-after-eruption"
    arr.coords["time"].attrs["units"] = "yr"
    return arr


class Reconstructor(BaseModel):
    """Class that keeps a static set of arrays.

    The arrays represent time series that have already been deconvolved completely, thus
    we assume no further deconvolution analysis is needed.

    Parameters
    ----------
    tau : np.ndarray
        The time axis for the deconvolution.
    time : xr.DataArray
        The time axis for the response functions.
    response_temp_so2 : np.ndarray
        The temperature to SO2 response function.
    response_temp_rf : np.ndarray
        The temperature to RF response function.
    forcing : xr.DataArray
        The forcing time series.
    signal : xr.DataArray
        The signal time series.
    name : str
        The name of the reconstructor.
    normalise : bool
        Whether the data is normalised.
    """

    class Config:
        """Configuration for the Reconstructor BaseModel object."""

        validate_assignment = True
        strict = True
        arbitrary_types_allowed = True

    tau: np.ndarray = Field(..., alias="tau", frozen=True)
    time: xr.DataArray = Field(..., alias="time", frozen=True)
    response_temp_so2: np.ndarray = Field(..., alias="response_temp_so2", frozen=True)
    response_temp_rf: np.ndarray = Field(..., alias="response_temp_rf", frozen=True)
    forcing: xr.DataArray = Field(..., alias="forcing", frozen=True)
    signal: xr.DataArray = Field(..., alias="signal", frozen=True)
    name: str = Field(..., frozen=False)
    normalise: bool = Field(..., frozen=True)


class CESMData(BaseModel):
    """Class to store the CESM2 data.

    Parameters
    ----------
    strength : Literal["strong", "medium", "medium-plus", "size5000", "tt-2sep", "double-overlap"], optional
        The strength of the eruption, by default "strong".
    dims : list[str], optional
        The dimensions of the data to average over, by default ["lat", "lon"]. An empty
        list keeps all dimensions.
    """

    strength: T_Strengths = Field(default="strong", frozen=True)
    dims: list[str] = Field(default=["lat", "lon"], frozen=True)

    class Config:
        """Configuration for the CESMData BaseModel object."""

        validate_assignment = True
        frozen = True
        extra = "forbid"
        strict = True

    @cached_property
    def so2(self) -> xr.DataArray:
        """Get the CESM2 SO2 data.

        Returns
        -------
        xr.DataArray
            The SO2 data.
        """
        x = self._get_aod_cesm().time.data
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
            case "tt-2sep":
                y[0] = 400
                y[24] = 400
            case "double-overlap":
                y[0] = 400
                y[48] = 400
            case _:
                vdd.utils.never_called(self.strength)
        out = xr.DataArray(y, coords={"time": x}, dims=["time"])
        return self._align_arrays("so2", out)

    @cached_property
    def tmso2(self) -> xr.DataArray:
        """Get the CESM2 TMSO2 data.

        Returns
        -------
        xr.DataArray
            The TMSO2 data.
        """
        out = self._get_tmso2_cesm()
        return self._align_arrays("tmso2", out)

    @cached_property
    def aod(self) -> xr.DataArray:
        """Get the CESM2 AOD data.

        Returns
        -------
        xr.DataArray
            The AOD data.
        """
        out = self._get_aod_cesm()
        return self._align_arrays("aod", out)

    @cached_property
    def rf(self) -> xr.DataArray:
        """Get the CESM2 RF data.

        Returns
        -------
        xr.DataArray
            The RF data.
        """
        out = self._get_rf_cesm()
        return self._align_arrays("rf", out)

    @cached_property
    def rf_control(self) -> xr.DataArray:
        """Get the CESM2 control RF data.

        Returns
        -------
        xr.DataArray
            The control RF data.
        """
        if not hasattr(self, "_rf_control"):
            _ = self.rf
        out = self._rf_control
        return self._align_arrays("rf_control", out)

    @cached_property
    def temp(self) -> xr.DataArray:
        """Get the CESM2 temperature data.

        Returns
        -------
        xr.DataArray
            The temperature data.
        """
        out = self._get_trefht_cesm()
        return self._align_arrays("temp", out)

    @cached_property
    def temperature_control(self) -> xr.DataArray:
        """Get the CESM2 control temperature data.

        Returns
        -------
        xr.DataArray
            The control temperature data.
        """
        if not hasattr(self, "_temperature_control"):
            _ = self.temp
        out = self._temperature_control
        return self._align_arrays("temperature_control", out)

    def initialise_data(self) -> None:
        """Initialise the data, ensuring that it is loaded and aligned."""
        _ = (
            self.so2,
            self.tmso2,
            self.aod,
            self.rf,
            self.temp,
            self.temperature_control,
            self.rf_control,
        )

    def _align_arrays(self, new: str, new_obj: xr.DataArray) -> xr.DataArray:
        out = list(
            set(self.__dict__.keys())
            & {"temperature_control", "rf_control", "temp", "tmso2", "so2", "rf", "aod"}
        )
        if out:
            aligned = xr.align(new_obj, *[getattr(self, o) for o in out])
            self.__dict__.update({o: a for o, a in zip(out, aligned[1:], strict=True)})
        else:
            aligned = (new_obj,)
        return aligned[0]

    def _get_tmso2_cesm(self, plot_example: bool = False) -> xr.DataArray:
        """Get the CESM2 TMSO2 data."""
        data = (
            volcano_base.load.FindFiles()
            .find("e_fSST1850", "TMSO2", "h0", self.strength)
            .remove("ens1" if self.strength != "tt-2sep" else "ens0")
            .keep_most_recent()
        )
        files = data.load()
        shift = 35 if self.strength == "double-overlap" else None
        files = volcano_base.manipulate.mean_flatten(files, dims=self.dims)
        files = volcano_base.manipulate.shift_arrays(files, custom=shift, daily=False)
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

    def _get_aod_cesm(self, plot_example: bool = False) -> xr.DataArray:
        """Get the CESM2 AOD data."""
        data = (
            volcano_base.load.FindFiles()
            .find("e_fSST1850", "AODVISstdn", "h0", self.strength)
            .remove("ens1" if self.strength != "tt-2sep" else "ens0")
            .keep_most_recent()
        )
        files = data.load()
        shift = 35 if self.strength == "double-overlap" else None
        files = volcano_base.manipulate.mean_flatten(files, dims=self.dims)
        files = volcano_base.manipulate.shift_arrays(files, custom=shift, daily=False)
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

    def _load_rf_cesm(
        self,
    ) -> tuple[xr.DataArray, xr.DataArray, list[xr.DataArray], list[xr.DataArray]]:
        control_data: volcano_base.load.FindFiles = (
            volcano_base.load.FindFiles()
            .find("e_fSST1850", "h0", "control", {"FLNT", "FSNT"}, "ens1")
            .sort("ensemble", "attr")
        )
        data: volcano_base.load.FindFiles = (
            volcano_base.load.FindFiles()
            .find("e_fSST1850", {"FLNT", "FSNT"}, "h0", self.strength)
            .remove("ens1" if self.strength != "tt-2sep" else "ens0")
            .keep_most_recent()
            .sort("ensemble", "attr")
        )
        c_size = 1
        match (control_data.copy().keep("FLNT"), control_data.copy().keep("FSNT")):
            case (c_flnt, c_fsnt) if len(c_flnt) == c_size and len(c_fsnt) == c_size:
                c_flnt_xr = c_flnt.load()[0]
                c_fsnt_xr = c_fsnt.load()[0]
            case _:
                raise ValueError("Control data not found.")
        match self.strength:
            case "size5000" | "tt-2sep":
                d_size = 2
            case "double-overlap":
                d_size = 1
            case _:
                d_size = 4
        match (data.copy().keep("FLNT"), data.copy().keep("FSNT")):
            case (flnt, fsnt) if len(flnt) == d_size and len(fsnt) == d_size:
                flnt_xr = flnt.load()
                fsnt_xr = fsnt.load()
            case _:
                raise ValueError("Data not found.")
        c_fsnt_xr, c_flnt_xr = volcano_base.manipulate.mean_flatten(
            [c_fsnt_xr, c_flnt_xr], dims=self.dims
        )
        fsnt_xr = volcano_base.manipulate.mean_flatten(fsnt_xr, dims=self.dims)
        flnt_xr = volcano_base.manipulate.mean_flatten(flnt_xr, dims=self.dims)
        return c_fsnt_xr, c_flnt_xr, fsnt_xr, flnt_xr

    def _get_rf_cesm(self, plot_example: bool = False) -> xr.DataArray:
        """Get the CESM2 RF data."""
        c_fsnt_xr, c_flnt_xr, fsnt_xr, flnt_xr = self._load_rf_cesm()

        def remove_control(rf: xr.DataArray) -> xr.DataArray:
            rf, c_s, c_l = xr.align(rf, c_fsnt_xr, c_flnt_xr)
            rf.data = rf.data - (c_s.data - c_l.data)
            # Combine attributes, overwrite if needed
            original_attrs = rf.attrs
            new_attrs = {
                "attr": "RF",
                "long_name": "Net radiative flux (residual) at top-of-model",
            }
            combined_attrs = {**original_attrs, **new_attrs}
            return rf.assign_attrs(**combined_attrs).rename("RESTOM")

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
        shift = 35 if self.strength == "double-overlap" else None
        rf = volcano_base.manipulate.shift_arrays(rf, custom=shift, daily=False)
        rf = volcano_base.manipulate.shift_arrays(rf, custom=1)
        rf = volcano_base.manipulate.subtract_mean_of_tail(rf)
        rf = volcano_base.manipulate.data_array_operation(rf, _convert_time_start_zero)
        rf = list(xr.align(*rf))
        # Create control run time series
        c_total = c_fsnt_xr.copy()
        c_total.data -= c_flnt_xr.data
        c_total = c_total.assign_attrs(attr="RF")
        c_total = volcano_base.manipulate.subtract_mean_of_tail([c_total])[0]
        c_total = volcano_base.manipulate.data_array_operation(
            [c_total], _convert_time_start_zero
        )[0]
        c_total.data = np.asarray(c_total.data)
        # Fourier
        # self._rf_control = volcano_base.manipulate.remove_seasonality(
        #     c_total.dropna("time"), freq=1.09, radius=0.01
        # )
        # Climatology
        c_total_float = c_total.dropna("time")
        c_total_dt = c_total_float.assign_coords(
            time=volcano_base.manipulate.float2dt(c_total_float.time, freq="MS")
        )
        c_total_dt = volcano_base.manipulate.subtract_climatology(
            c_total_dt, c_total_dt, "time.month"
        )[0]
        self._rf_control = c_total_dt.assign_coords(time=c_total_float.time)
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
        control_data: volcano_base.load.FindFiles = (
            volcano_base.load.FindFiles()
            .find("e_BWma1850", "h0", "control", "TREFHT", "ens0")
            .sort("ensemble", "attr")
        )
        if len(control_data) != 1:
            raise ValueError("Control data not found.")
        control_l = control_data.load()
        control_l = volcano_base.manipulate.mean_flatten(control_l, dims=self.dims)
        data = (
            volcano_base.load.FindFiles()
            .find("e_BWma1850", "TREFHT", "h0", self.strength)
            .remove("ens1" if self.strength != "tt-2sep" else "ens0")
            .keep_most_recent()
        )
        files = data.load()
        orig_attrs = files[0].attrs
        shift = 35 if self.strength == "double-overlap" else None
        files = volcano_base.manipulate.mean_flatten(files, dims=self.dims)
        control = control_l[0]
        files = volcano_base.manipulate.shift_arrays(files, custom=shift, daily=False)
        files = volcano_base.manipulate.shift_arrays(files, custom=1)

        def subtract_control_array(arr: xr.DataArray) -> xr.DataArray:
            arr, c = xr.align(arr, control)
            arr.data = arr.data - c.data
            return arr

        def subtract_control_mean(arr: xr.DataArray) -> xr.DataArray:
            return arr - volcano_base.config.MEANS["TREFHT"]

        control_l = volcano_base.manipulate.data_array_operation(
            control_l, _convert_time_start_zero
        )
        # Fourier
        # self._temperature_control = volcano_base.manipulate.remove_seasonality(
        #     volcano_base.manipulate.get_median(control_l, xarray=True).dropna("time")
        #     - volcano_base.config.MEANS["TREFHT"],
        #     freq=1,
        #     radius=0.1,
        # )
        # Climatology
        control_float = volcano_base.manipulate.get_median(
            control_l, xarray=True
        ).dropna("time")
        control_dt = control_float.assign_coords(
            time=volcano_base.manipulate.float2dt(control_float.time, freq="MS")
        )
        control_dt = volcano_base.manipulate.subtract_climatology(
            control_dt, control_dt, "time.month"
        )[0]
        self._temperature_control = control_dt.assign_coords(time=control_float.time)
        subtract_control = (
            subtract_control_array if len(data) == 1 else subtract_control_mean
        )
        files = volcano_base.manipulate.data_array_operation(files, subtract_control)
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
        return mean_array.dropna("time").assign_attrs(**orig_attrs)


class _PostInitCaller(ABC, type):
    def __call__(cls, *args, **kwargs):
        obj = type.__call__(cls, *args, **kwargs)
        obj.__post_init__()
        return obj


class Deconvolve(metaclass=_PostInitCaller):
    """Class for deconvolving data."""

    name: str = "Deconvolve"

    def __init__(self, normalise: bool = False) -> None:
        def _deconv(signal, forcing) -> tuple[np.ndarray, np.ndarray]:
            guess = np.heaviside(np.arange(len(signal)) - len(signal) // 2, 1)
            kwargs = {"initial_guess": guess, "iteration_list": 1000}
            return fppanalysis.RL_gauss_deconvolve(signal, forcing, **kwargs)

        self._deconv_method: Callable[
            [np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]
        ] = _deconv
        self.normalise = normalise

    def __post_init__(self) -> None:
        """Run the normalisation if needed."""
        if self.normalise:
            self._update_if_normalise()

    @abstractmethod
    def _update_if_normalise(self) -> None: ...

    @abstractproperty
    def tau(self) -> np.ndarray:
        """Time axis for the deconvolution."""
        raise NotImplementedError

    @abstractproperty
    def so2(self) -> xr.DataArray:
        """SO2 time series data."""
        raise NotImplementedError

    @abstractproperty
    def so2_decay(self) -> xr.DataArray:
        """SO2 time series data."""
        raise NotImplementedError

    @abstractproperty
    def rf(self) -> xr.DataArray:
        """Radiative forcing time series data."""
        raise NotImplementedError

    @abstractproperty
    def temp(self) -> xr.DataArray:
        """Temperature time series data."""
        raise NotImplementedError

    @abstractproperty
    def temp_control(self) -> xr.DataArray:
        """Temperature control time series data."""
        raise NotImplementedError

    def dump_reconstructor(self) -> Reconstructor:
        """Dump the current deconvolution object into a reconstructor object."""
        return Reconstructor(
            tau=self.tau,
            time=self.temp.time,
            response_temp_so2=self.response_temp_so2,
            response_temp_rf=self.response_temp_rf,
            forcing=self.rf,
            signal=self.temp,
            name=self.name,
            normalise=self.normalise,
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

        Examples
        --------
        The default deconvolution method could be re-set as follows:
        >>> def default_deconv(signal, forcing) -> tuple[np.ndarray, np.ndarray]:
        ...     guess = np.heaviside(np.arange(len(signal)) - len(signal) // 2, 1)
        ...     kwargs = {"initial_guess": guess, "iteration_list": 1000}
        ...     return fppanalysis.RL_gauss_deconvolve(signal, forcing, **kwargs)
        >>> dec = vdd.load.Deconvolve()
        >>> dec.change_deconvolution_method(default_deconv)

        Or if the heaviside function is not needed, one could use the following:
        >>> kwargs = {"iteration_list": 1000}
        >>> dec.change_deconvolution_method(fppanalysis.RL_gauss_deconvolve, **kwargs)

        Or even
        >>> dec.change_deconvolution_method(fppanalysis.RL_gauss_deconvolve, iteration_list=1000)
        """

        def method_(
            signal: np.ndarray, forcing: np.ndarray
        ) -> tuple[np.ndarray, np.ndarray]:
            return method(signal, forcing, *args, **kwargs)

        self._deconv_method = method_

    @cached_property
    def _response_rf_so2_tup(self) -> tuple[np.ndarray, np.ndarray]:
        """Deconvolve the RF signal with the SO2 signal."""
        signal, err = self._deconv_method(self.rf.data, self.so2.data)
        if self.normalise:
            # signal = vdd.utils.normalise(signal)
            signal = (signal - signal.mean()) / signal.std()
        return signal, err

    @property
    def response_rf_so2(self) -> np.ndarray:
        """Deconvolve the RF signal with the SO2 delta signal.

        Returns
        -------
        np.ndarray
            The deconvolved RF to SO2 signal.
        """
        return self._response_rf_so2_tup[0].flatten()

    @property
    def _response_rf_so2_err(self) -> np.ndarray:
        """Deconvolve the RF signal with the SO2 signal."""
        return self._response_rf_so2_tup[1]

    @cached_property
    def _response_temp_so2_tup(self) -> tuple[np.ndarray, np.ndarray]:
        """Deconvolve the temperature signal with the SO2 signal."""
        signal, err = self._deconv_method(self.temp.data, self.so2.data)
        if self.normalise:
            # signal = vdd.utils.normalise(signal)
            signal = (signal - signal.mean()) / signal.std()
        return signal, err

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
    def _response_temp_so2_err(self) -> np.ndarray:
        """Deconvolve the temperature signal with the SO2 signal."""
        return self._response_temp_so2_tup[1]

    @cached_property
    def _response_rf_so2_decay_tup(self) -> tuple[np.ndarray, np.ndarray]:
        """Deconvolve the RF signal with the SO2 decay signal."""
        signal, err = self._deconv_method(self.rf.data, self.so2_decay.data)
        if self.normalise:
            # signal = vdd.utils.normalise(signal)
            signal = (signal - signal.mean()) / signal.std()
        return signal, err

    @property
    def response_rf_so2_decay(self) -> np.ndarray:
        """Deconvolve the RF signal with the SO2 decay signal.

        Returns
        -------
        np.ndarray
            The deconvolved RF to SO2 signal.
        """
        return self._response_rf_so2_decay_tup[0].flatten()

    @property
    def _response_rf_so2_decay_err(self) -> np.ndarray:
        """Deconvolve the RF signal with the SO2 decay signal."""
        return self._response_rf_so2_decay_tup[1]

    @cached_property
    def _response_temp_so2_decay_tup(self) -> tuple[np.ndarray, np.ndarray]:
        """Deconvolve the temperature signal with the SO2 decay signal."""
        signal, err = self._deconv_method(self.temp.data, self.so2_decay.data)
        if self.normalise:
            # signal = vdd.utils.normalise(signal)
            signal = (signal - signal.mean()) / signal.std()
        return signal, err

    @property
    def response_temp_so2_decay(self) -> np.ndarray:
        """Deconvolve the temperature signal with the SO2 decay signal.

        Returns
        -------
        np.ndarray
            The deconvolved temperature to SO2 decay signal.
        """
        return self._response_temp_so2_decay_tup[0].flatten()

    @property
    def _response_temp_so2_decay_err(self) -> np.ndarray:
        """Deconvolve the temperature signal with the SO2 decay signal."""
        return self._response_temp_so2_decay_tup[1]

    @cached_property
    def _response_temp_rf_tup(self) -> tuple[np.ndarray, np.ndarray]:
        """Deconvolve the temperature signal with the RF signal."""
        signal, err = self._deconv_method(self.temp.data, self.rf.data)
        if self.normalise:
            # signal = vdd.utils.normalise(signal)
            signal = (signal - signal.mean()) / signal.std()
        return signal, err

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
    def _response_temp_rf_err(self) -> np.ndarray:
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
        plt.semilogy(self._response_rf_so2_err)

    def plot_dec_temp_with_so2(self) -> None:
        """Deconvolve the temperature signal with the SO2 signal."""
        # Quick comparison
        plt.figure()
        vdd.utils.normalise(self.temp).plot()
        vdd.utils.normalise(self.so2).plot()
        plt.figure()
        plt.plot(self.tau, self.response_temp_so2)
        plt.figure()
        plt.semilogy(self._response_temp_so2_err)

    def plot_dec_temp_with_rf(self) -> None:
        """Deconvolve the temperature signal with the SO2 signal."""
        # Quick comparison
        plt.figure()
        vdd.utils.normalise(self.temp).plot()
        vdd.utils.normalise(self.rf).plot()
        plt.figure()
        plt.plot(self.tau, self.response_temp_rf)
        plt.figure()
        plt.semilogy(self._response_temp_rf_err)


class DeconvolveCESM(Deconvolve):
    """Class for deconvolving data from CESM2.

    Parameters
    ----------
    normalise : bool, optional
        Whether to normalise the data, by default False.
    pad_before : bool, optional
        Whether to pad with zeros all time series before the convolution, by default
        False.
    cesm : CESMData | None, optional
        The CESM data to use. If not given, the "strong" strength data will be used.
    """

    def __init__(
        self,
        normalise: bool = False,
        pad_before: bool = False,
        cesm: CESMData | None = None,
    ) -> None:
        super().__init__(normalise)
        self._data = CESMData() if cesm is None else cesm
        # Since the time series can sometimes be of different length, we make sure to
        # call all of them before using them here. That way, they will be updated within
        # the CESMData class before assigning them here.
        self.pad_before = pad_before
        self.name = f"CESM2 {self._data.strength}"

    def _update_if_normalise(self) -> None:
        self.so2 = (self.so2 - self.so2.mean()) / self.so2.std()
        self.tmso2 = (self.tmso2 - self.tmso2.mean()) / self.tmso2.std()
        self.aod = (self.aod - self.aod.mean()) / self.aod.std()
        self.rf = (self.rf - self.rf.mean()) / self.rf.std()
        self.temp = (self.temp - self.temp.mean()) / self.temp.std()
        self.temp_control = (
            self.temp_control - self.temp_control.mean()
        ) / self.temp_control.std()
        # self.so2 = vdd.utils.normalise(self.so2)
        # self.aod = vdd.utils.normalise(self.aod)
        # self.rf = vdd.utils.normalise(self.rf)
        # self.temp = vdd.utils.normalise(self.temp)

    @cached_property
    def so2(self) -> xr.DataArray:
        """SO2 time series data."""
        self._data.initialise_data()
        return (
            vdd.utils.pad_before_convolution(self._data.so2)
            if self.pad_before
            else self._data.so2
        )

    @cached_property
    def so2_decay(self) -> xr.DataArray:
        """SO2 time series data."""
        self._data.initialise_data()
        return (
            vdd.utils.pad_before_convolution(self._data.tmso2)
            if self.pad_before
            else self._data.tmso2
        )

    @cached_property
    def tau(self) -> np.ndarray:
        """Time axis for the deconvolution."""
        self._data.initialise_data()
        return (
            self.so2.time.data
            if self.pad_before
            else self.so2.time.data - self.so2.time.data[len(self.so2.time.data) // 2]
        )

    @cached_property
    def tmso2(self) -> xr.DataArray:
        """SO2 column burden."""
        self._data.initialise_data()
        return (
            vdd.utils.pad_before_convolution(self._data.tmso2)
            if self.pad_before
            else self._data.tmso2
        )

    @cached_property
    def aod(self) -> xr.DataArray:
        """Aerosol optical depth time series data."""
        self._data.initialise_data()
        return (
            vdd.utils.pad_before_convolution(self._data.aod)
            if self.pad_before
            else self._data.aod
        )

    @cached_property
    def rf(self) -> xr.DataArray:
        """Radiative forcing time series data."""
        self._data.initialise_data()
        attrs = self._data.rf.attrs
        return (
            vdd.utils.pad_before_convolution(self._data.rf * -1)
            if self.pad_before
            else self._data.rf * -1
        ).assign_attrs(**attrs)

    @cached_property
    def rf_control(self) -> xr.DataArray:
        """RF time series data."""
        self._data.initialise_data()
        return (
            vdd.utils.pad_before_convolution(self._data.rf_control * -1)
            if self.pad_before
            else self._data.rf_control * -1
        )

    @cached_property
    def temp(self) -> xr.DataArray:
        """Temperature time series data."""
        self._data.initialise_data()
        attrs = self._data.temp.attrs
        return (
            vdd.utils.pad_before_convolution(self._data.temp * -1)
            if self.pad_before
            else self._data.temp * -1
        ).assign_attrs(**attrs)

    @cached_property
    def temp_control(self) -> xr.DataArray:
        """Temperature control time series data."""
        self._data.initialise_data()
        return (
            vdd.utils.pad_before_convolution(self._data.temperature_control * -1)
            if self.pad_before
            else self._data.temperature_control * -1
        )

    @cached_property
    def _response_aod_so2_tup(self) -> tuple[np.ndarray, np.ndarray]:
        """Deconvolve the AOD signal with the SO2 signal."""
        signal, err = self._deconv_method(self.aod.data, self.so2.data)
        if self.normalise:
            # signal = vdd.utils.normalise(signal)
            signal = (signal - signal.mean()) / signal.std()
        return signal, err

    @property
    def response_aod_so2(self) -> np.ndarray:
        """Deconvolve the AOD signal with the SO2 signal.

        Returns
        -------
        np.ndarray
            The deconvolved AOD to SO2 signal.
        """
        return self._response_aod_so2_tup[0].flatten()

    @property
    def _response_aod_so2_err(self) -> np.ndarray:
        """Deconvolve the AOD signal with the SO2 signal."""
        return self._response_aod_so2_tup[1]

    @cached_property
    def _response_rf_aod_tup(self) -> tuple[np.ndarray, np.ndarray]:
        """Deconvolve the RF signal with the AOD signal."""
        signal, err = self._deconv_method(self.rf.data, self.aod.data)
        if self.normalise:
            # signal = vdd.utils.normalise(signal)
            signal = (signal - signal.mean()) / signal.std()
        return signal, err

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
    def _response_rf_aod_err(self) -> np.ndarray:
        """Deconvolve the RF signal with the AOD signal."""
        return self._response_rf_aod_tup[1]

    @cached_property
    def _response_temp_aod_tup(self) -> tuple[np.ndarray, np.ndarray]:
        """Deconvolve the temperature signal with the AOD signal."""
        signal, err = self._deconv_method(self.temp.data, self.aod.data)
        if self.normalise:
            # signal = vdd.utils.normalise(signal)
            signal = (signal - signal.mean()) / signal.std()
        return signal, err

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
    def _response_temp_aod_err(self) -> np.ndarray:
        """Deconvolve the temperature signal with the AOD signal."""
        return self._response_temp_aod_tup[1]


class DeconvolveOB16(Deconvolve):
    """Class for deconvolving data from Otto-Bliesner et at. (2016).

    Parameters
    ----------
    data : volcano_base.load.OttoBliesner | Literal["h0", "h1"], optional
        The OB16 data class to use. If not given, the h1 (daily) data will be loaded.
    normalise : bool, optional
        Whether to normalise the data, by default False.
    length : int, optional
        After 1850, the SO2 dataset mismatches with the simulation output. This
        parameter specifies how many items should be included in all the arrays. Default
        is 0, which means all the data will be included.

    Raises
    ------
    ValueError
        If the data is invalid.
    """

    def __init__(
        self,
        data: volcano_base.load.OttoBliesner | Literal["h0", "h1"] = "h1",
        normalise: bool = False,
        length: int = 0,
    ) -> None:
        super().__init__(normalise)
        match data:
            case "h0":
                self.data = volcano_base.load.OttoBliesner(freq="h0", progress=True)
            case "h1":
                self.data = volcano_base.load.OttoBliesner(freq="h1", progress=True)
            case volcano_base.load.OttoBliesner():
                self.data = data
            case _:
                raise ValueError(f"Invalid data: {data}")
        self._start_pt = 0
        if length:
            self._end_pt = length if length - self._start_pt % 2 else length + 1
        else:
            self._end_pt = -1

    def _update_if_normalise(self) -> None:
        self.so2 = (self.so2 - self.so2.mean()) / self.so2.std()
        self.so2_decay = (self.so2_decay - self.so2_decay.mean()) / self.so2_decay.std()
        self.rf = (self.rf - self.rf.mean()) / self.rf.std()
        self.rf_control = (
            self.rf_control - self.rf_control.mean()
        ) / self.rf_control.std()
        self.temp = (self.temp - self.temp.mean()) / self.temp.std()
        self.temp_control = (
            self.temp_control - self.temp_control.mean()
        ) / self.temp_control.std()
        # self.so2 = vdd.utils.normalise(self.so2)
        # self.rf = vdd.utils.normalise(self.rf)
        # self.temp = vdd.utils.normalise(self.temp)

    @cached_property
    def so2(self) -> xr.DataArray:
        """SO2 time series data."""
        return self.data.aligned_arrays["so2-start"][self._start_pt : self._end_pt]

    @cached_property
    def so2_decay(self) -> xr.DataArray:
        """SO2 time series data."""
        return (
            self.data.aligned_arrays["so2-decay-start"][self._start_pt : self._end_pt]
            / 510e3
        )

    @cached_property
    def tau(self) -> np.ndarray:
        """Time axis for the deconvolution."""
        tau = self.so2.time.data - (
            self.so2.time.data[len(self.so2.time.data) // 2]
            - cftime.DatetimeNoLeap(0, 1, 1, has_year_zero=True, calendar="noleap")
        )
        return np.asarray(volcano_base.manipulate.dt2float(tau))

    @cached_property
    def rf(self) -> xr.DataArray:
        """Radiative forcing time series data."""
        return self.data.aligned_arrays["rf"][self._start_pt : self._end_pt]

    @cached_property
    def rf_control(self) -> xr.DataArray:
        """RF time series data."""
        return xr.align(self.data.rf_control, self.rf)[0][self._start_pt : self._end_pt]

    @cached_property
    def temp(self) -> xr.DataArray:
        """Temperature time series data."""
        return self.data.aligned_arrays["temperature"][self._start_pt : self._end_pt]

    @cached_property
    def temp_control(self) -> xr.DataArray:
        """Temperature time series data."""
        return xr.align(self.data.temperature_control, self.temp)[0][
            self._start_pt : self._end_pt
        ]


class CutOff:
    """Cut off the response functions of a deconvolution object."""

    def __init__(self, dec: Deconvolve, arrays: T_RF | T_SO2 | RF_SO2) -> None:
        self.dec = dec
        self.ts_specifier: T_RF | T_SO2 | RF_SO2 = arrays
        self.cuts: dict[str, xr.Dataset] = {}
        self.ensembles: dict[str, xr.Dataset] = {}

    def dump_reconstructor(
        self,
        cut: int,
        temp_so2: np.ndarray | None = None,
        temp_rf: np.ndarray | None = None,
    ) -> Reconstructor:
        """Dump the current deconvolution object into a reconstructor object.

        Parameters
        ----------
        cut : int
            The cut-off time lag to use.
        temp_so2 : np.ndarray | None, optional
            The response function of the SO2 signal, by default None.
        temp_rf : np.ndarray | None, optional
            The response function of the RF signal, by default None.

        Returns
        -------
        Reconstructor
            The reconstructor object.

        Raises
        ------
        ValueError
            If not exactly one of temp_so2 or temp_rf is given.
        """
        kwargs = {
            "tau": self.dec.tau,
            "time": self.dec.rf.time,
            "forcing": self.forcing,
            "signal": self.output,
            "name": self.dec.name,
            "normalise": self.dec.normalise,
        }
        match temp_so2, temp_rf:
            case (None, None) | (np.ndarray(), np.ndarray()):
                raise ValueError("Exactly one of temp_so2 or temp_rf must be given.")
            case np.ndarray(), None:
                kwargs["response_temp_so2"] = temp_so2
                kwargs["response_temp_rf"] = self.cuts[str(cut)].response.data
            case None, np.ndarray():
                kwargs["response_temp_so2"] = self.cuts[str(cut)].response.data
                kwargs["response_temp_rf"] = temp_rf
            case _:
                raise ValueError(
                    f"I do not recognise the types of the response functions. Got: {temp_so2}, {temp_rf}"
                )
        return Reconstructor(**kwargs)

    @cached_property
    def response(self) -> np.ndarray:
        """The response function in the convolution."""
        out = getattr(
            self.dec, f"response_{self.ts_specifier[0]}_{self.ts_specifier[1]}"
        )
        out[self.dec.tau <= 0] = 0
        return out

    @cached_property
    def forcing(self) -> xr.DataArray:
        """The forcing time series in the convolution."""
        return getattr(self.dec, self.ts_specifier[1])

    @cached_property
    def output(self) -> xr.DataArray:
        """The final output time series of the convolution."""
        return getattr(self.dec, self.ts_specifier[0])

    @cached_property
    def control(self) -> xr.DataArray:
        """The control time series in the convolution."""
        return getattr(self.dec, f"{self.ts_specifier[0]}_control")

    def cut_off(self, cutoff: int | Iterable[int]) -> Self:
        """Cut off the response function at a given time lag."""
        match cutoff:
            case int():
                self._single_cut_off(cutoff)
            case Iterable():
                for c in set(cutoff):
                    if not isinstance(c, int):
                        raise ValueError(
                            "cutoff must be an integer or a sequence of integers."
                        )
                    self._single_cut_off(c)
            case _:
                raise ValueError("cutoff must be an integer or a sequence of integers.")
        return self

    def _single_cut_off(self, cutoff: int) -> None:
        if str(cutoff) in self.cuts:
            return
        r_cut = self.response.copy()
        r_cut[len(r_cut) // 2 + cutoff :] = 0
        tau = self.dec.tau
        time = self.output.time
        temp_r = np.convolve(self.forcing, r_cut, "same")
        ds = xr.Dataset(
            {
                "response": ("tau", r_cut, {"label": f"cut {cutoff}"}),
                "temp_rec": ("time", temp_r, {"label": f"temp rec {cutoff}"}),
            },
            coords={"tau": tau, "time": time},
        )
        self.cuts[str(cutoff)] = ds

    def generate_ensembles(self, n: int) -> None:
        """Generate an ensemble of response function estimates."""
        if not self.cuts:
            raise ValueError("No cuts have been made.")
        iters = np.arange(200)
        for k, v in self.cuts.items():
            if k in self.ensembles:
                continue
            arrays: dict[str, tuple] = {}
            for i in range(n):
                temp_rec = v.temp_rec.copy()
                temp_random = fppanalysis.signal_rand_phase(self.control.data)
                temp_rec += temp_random
                res_rec, err = fppanalysis.RL_gauss_deconvolve(
                    temp_rec, self.forcing, len(iters) - 1
                )
                r_cut_rec = res_rec.flatten()
                r_cut_rec[self.dec.tau <= 0] = 0
                arrays[f"response_{i}"] = ("tau", r_cut_rec, {"label": f"response {i}"})
                arrays[f"iters_{i}"] = ("iters", err.flatten(), {"label": f"err {i}"})
                # arrays[f"temp_{i}"] = rec
            self.ensembles[k] = xr.Dataset(
                arrays,
                coords={"tau": self.dec.tau, "time": self.output.time, "iters": iters},
            )


class TSComparison:
    """Easily compare the statistics of two time series."""

    def __init__(
        self,
        original: xr.DataArray,
        reconstructed: np.ndarray,
        peaks: np.ndarray,
    ) -> None:
        self.orig = original
        self.rec = reconstructed
        assert original.data.shape == reconstructed.shape
        self.peaks = peaks

    def _find_peaks(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Compute the peaks of the time series."""
        _idx = np.argwhere(self.peaks.data > 0)
        peak_times = self.orig.time.data[_idx].flatten()
        peaks_ts1 = self.orig.data[_idx].flatten()
        peaks_ts2 = self.rec[_idx].flatten()
        return peak_times, peaks_ts1, peaks_ts2

    @cached_property
    def _peaks_tup(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Compute the peaks of the residual of the SO2."""
        return self._find_peaks()

    @property
    def peaks_orig(self) -> np.ndarray:
        """Peaks of the temperature from SO2.

        Returns
        -------
        np.ndarray
            Peaks of the temperature from SO2.
        """
        return self._peaks_tup[1]

    @property
    def peaks_rec(self) -> np.ndarray:
        """Peaks of the temperature from radiative forcing.

        Returns
        -------
        np.ndarray
            Peaks of the temperature from radiative forcing.
        """
        return self._peaks_tup[2]

    @cached_property
    def residual(self) -> np.ndarray:
        """Compute the residual of the SO2.

        Returns
        -------
        np.ndarray
            Residual of the SO2.
        """
        return self.orig.data - self.rec

    def correlation(self) -> None:
        """Compute the correlation between the residuals and temperature."""
        corr_time, corr = fppanalysis.corr_fun(self.residual, self.orig.data, 1 / 12)
        # corr_ts2_time, corr_ts2 = fppanalysis.corr_fun(
        #     self.residual, self.rec, 1 / 12
        # )
        plt.figure()
        plt.plot(corr_time, corr, label="Residual / TS1", alpha=0.7)
        # plt.plot(corr_ts2_time, corr_ts2, label="Residual / TS2", alpha=0.7)
        plt.xlabel("Time lag ($\\tau$) [yr]")
        plt.ylabel("Correlation between residual \nand original temperature")
        plt.legend()
        # plt.savefig(
        #     _SAVE_DIR / f"{self.sim_name}-correlation-residual-reconstructed.png"
        # )

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
        f_res, p_res = self._spectrum_1d(self.residual)
        f_orig, p_orig = self._spectrum_1d(self.orig.data)
        f_rec, p_rec = self._spectrum_1d(self.rec)
        plt.figure()
        plt.plot(f_res, p_res, label="Residual", alpha=0.5)
        plt.plot(f_orig, p_orig, label="Original", alpha=0.5)
        plt.plot(f_rec, p_rec, label="Reconstructed", alpha=0.5)
        # Suppress the warning
        warnings.filterwarnings("ignore")
        cosmoplots.change_log_axis_base(plt.gca(), "both")
        warnings.resetwarnings()
        plt.xlabel("Frequency")
        plt.ylabel("Power spectral density")
        plt.legend()
        # plt.savefig(_SAVE_DIR / f"{self.sim_name}-spectrum-residual-control_temp.png")

    @staticmethod
    def _peak_difference_ttest(basis: np.ndarray) -> str:
        # Specify the value to test for symmetry
        test_value = 0
        # Perform a one-sample t-test to check for symmetry around the test value
        result = scipy.stats.ttest_1samp(basis, popmean=test_value)
        t_statistic, p_value = result.statistic, result.pvalue

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
        if p_value < sl:
            print(reject)
        else:
            print(reject_no)
        print(t_statistic, p_value)
        info("Reconstructed", p_value)
        return p_value

    def peak_difference_analysis(self) -> None:  # noqa:PLR0914
        """Plot the difference between the reconstructed and the original peaks."""
        basis = self.peaks_orig.data - self.peaks_rec
        ttest_res = self._peak_difference_ttest(basis)
        pdf, cdf, bin_centers = fppanalysis.distribution(basis, 30, ccdf=False)
        stats = scipy.stats.describe(basis)
        fit = scipy.stats.norm.fit(basis)
        dist = scipy.stats.skewnorm(
            a=stats.skewness, loc=stats.mean, scale=np.sqrt(stats.variance)
        )
        self._peak_difference_plot(
            (bin_centers, pdf),
            (fit, dist),
            "pdf",
            txt=ttest_res,
        )
        self._peak_difference_plot(
            (bin_centers, cdf),
            (fit, dist),
            "cdf",
            txt=ttest_res,
        )

    @staticmethod
    def _peak_difference_plot(
        dist_data: tuple[np.ndarray, np.ndarray],
        fits: tuple[np.ndarray, np.ndarray],
        dist: Literal["pdf", "cdf"],
        txt: str,
    ) -> None:
        norm_fit, skewnorm_fit = fits
        prop_cycle = plt.rcParams["axes.prop_cycle"]
        colors = prop_cycle.by_key()["color"]
        norm_so2 = getattr(scipy.stats.norm, dist)(dist_data[0], *norm_fit)
        skewnorm_so2 = getattr(skewnorm_fit, dist)(dist_data[0])
        plt.figure()
        ax = plt.gca()
        ax.bar(
            dist_data[0],
            dist_data[1],
            width=0.01,
            label=f"SO2 (p-value: {txt:.4f})",
            alpha=0.5,
        )
        bar_hand, bar_lab = ax.get_legend_handles_labels()
        # Norm
        (norm,) = ax.plot(
            dist_data[0], norm_so2, c=colors[0], label="_Norm SO2", alpha=0.5
        )
        # Skewnorm
        (skewnorm,) = plt.plot(
            dist_data[0],
            skewnorm_so2,
            "--",
            c=colors[0],
            label="_Skewnorm SO2",
            alpha=0.5,
        )
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
        # plt.savefig(_SAVE_DIR / f"{self.sim_name}-peak-difference-{dist}.png")

    def plot_reconstructions(self, fig: mpl.figure.Figure | None = None) -> None:
        """Plot the reconstruction of the data."""
        if fig is None:
            fig = plt.figure()
        ax = fig.gca()
        ax.set_xlabel("Time [yr]")
        ax.set_ylabel("Absolute")
        time_ = self.orig.time
        ax.plot(time_, self.orig.data, label="Original")
        ax.plot(time_, self.rec, label="Reconstriction")
        # ax.set_xlim((-790 * 365, -650 * 365))
        ax.legend(framealpha=0.5)
        # plt.savefig(_SAVE_DIR / f"{self.sim_name}-temp-reconstructed.png")
