"""Functions that fetches data from files and returns it as an xarray.DataArray."""

from abc import ABC, abstractmethod, abstractproperty
from collections.abc import Callable
from functools import cached_property
from typing import Any, Literal

import cftime
import fppanalysis
import matplotlib.pyplot as plt
import numpy as np
import volcano_base
import xarray as xr
from pydantic import BaseModel, Field

import vdd.utils


def _convert_time_start_zero(arr: xr.DataArray) -> xr.DataArray:
    """Convert the time to start from zero."""
    arr = arr.assign_coords(time=volcano_base.manipulate.dt2float(arr.time.data) - 1850)
    arr.coords["time"].attrs["long_name"] = "Time-after-eruption"
    arr.coords["time"].attrs["units"] = "yr"
    return arr


class CESMData(BaseModel):
    """Class to store the CESM2 data.

    Parameters
    ----------
    strength : Literal["strong", "medium", "medium-plus", "size5000", "tt-2sep", "double-overlap"], optional
        The strength of the eruption, by default "strong".
    """

    strength: Literal[
        "strong", "medium", "medium-plus", "size5000", "tt-2sep", "double-overlap"
    ] = Field(
        default="strong",
        frozen=True,
    )

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
            self._get_trefht_cesm()
        out = self._temperature_control
        return self._align_arrays("temperature_control", out)

    def initialise_data(self) -> None:
        """Initialise the data, ensuring that it is loaded and aligned."""
        _ = self.so2, self.aod, self.rf, self.temp, self.temperature_control

    def _align_arrays(self, new: str, new_obj: xr.DataArray) -> xr.DataArray:
        out = list(
            set(self.__dict__.keys())
            & {"temperature_control", "temp", "so2", "rf", "aod"}
        )
        if out:
            aligned = xr.align(new_obj, *[getattr(self, o) for o in out])
            self.__dict__.update({o: a for o, a in zip(out, aligned[1:], strict=True)})
        else:
            aligned = (new_obj,)
        return aligned[0]

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
        files = volcano_base.manipulate.mean_flatten(files, dims=["lat", "lon"])
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
            .remove("ens1" if self.strength != "tt-2sep" else "ens0")
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
        shift = 35 if self.strength == "double-overlap" else None
        rf = volcano_base.manipulate.shift_arrays(rf, custom=shift, daily=False)
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
        control_data: volcano_base.load.FindFiles = (
            volcano_base.load.FindFiles()
            .find("e_BWma1850", "h0", "control", "TREFHT", "ens0")
            .sort("ensemble", "attr")
        )
        if len(control_data) != 1:
            raise ValueError("Control data not found.")
        control_l = control_data.load()
        control_l = volcano_base.manipulate.mean_flatten(control_l, dims=["lat", "lon"])
        data = (
            volcano_base.load.FindFiles()
            .find("e_BWma1850", "TREFHT", "h0", self.strength)
            .remove("ens1" if self.strength != "tt-2sep" else "ens0")
            .keep_most_recent()
        )
        files = data.load()
        shift = 35 if self.strength == "double-overlap" else None
        files = volcano_base.manipulate.mean_flatten(files, dims=["lat", "lon"])
        control_l = volcano_base.manipulate.shift_arrays(
            control_l, custom=shift, daily=False
        )
        control = volcano_base.manipulate.shift_arrays(control_l, custom=1)[0]
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
        self._temperature_control = volcano_base.manipulate.remove_seasonality(
            volcano_base.manipulate.get_median(control_l, xarray=True).dropna("time")
            - volcano_base.config.MEANS["TREFHT"],
            freq=1,
            radius=0.1,
        )
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
        return mean_array.dropna("time")


class _PostInitCaller(ABC, type):
    def __call__(cls, *args, **kwargs):
        obj = type.__call__(cls, *args, **kwargs)
        obj.__post_init__()
        return obj


class Deconvolve(metaclass=_PostInitCaller):
    """Class for deconvolving data."""

    name: str = "Deconvolve"

    def __init__(self, normalise: bool = False) -> None:
        kwargs = {"iteration_list": 200}
        self._deconv_method: Callable[
            [np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]
        ] = lambda signal, forcing: fppanalysis.RL_gauss_deconvolve(
            signal, forcing, **kwargs
        )
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
        signal, err = self._deconv_method(self.rf.data, self.so2.data)
        if self.normalise:
            # signal = vdd.utils.normalise(signal)
            signal = (signal - signal.mean()) / signal.std()
        return signal, err

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
    def response_temp_so2_err(self) -> np.ndarray:
        """Deconvolve the temperature signal with the SO2 signal."""
        return self._response_temp_so2_tup[1]

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
    def tau(self) -> np.ndarray:
        """Time axis for the deconvolution."""
        self._data.initialise_data()
        return (
            self.so2.time.data
            if self.pad_before
            else self.so2.time.data - self.so2.time.data[len(self.so2.time.data) // 2]
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
        return (
            vdd.utils.pad_before_convolution(self._data.rf * -1)
            if self.pad_before
            else self._data.rf * -1
        )

    @cached_property
    def temp(self) -> xr.DataArray:
        """Temperature time series data."""
        self._data.initialise_data()
        return (
            vdd.utils.pad_before_convolution(self._data.temp * -1)
            if self.pad_before
            else self._data.temp * -1
        )

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
    def response_aod_so2_err(self) -> np.ndarray:
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
    def response_rf_aod_err(self) -> np.ndarray:
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
    def response_temp_aod_err(self) -> np.ndarray:
        """Deconvolve the temperature signal with the AOD signal."""
        return self._response_temp_aod_tup[1]


class DeconvolveOB16(Deconvolve):
    """Class for deconvolving data from Otto-Bliesner et at. (2016).

    Parameters
    ----------
    normalise : bool, optional
        Whether to normalise the data, by default False.
    data : volcano_base.load.OttoBliesner | Literal["h0", "h1"], optional
        The OB16 data class to use. If not given, the h1 (daily) data will be loaded.

    Raises
    ------
    ValueError
        If the data is invalid.
    """

    def __init__(
        self,
        normalise: bool = False,
        data: volcano_base.load.OttoBliesner | Literal["h0", "h1"] = "h1",
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

    def _update_if_normalise(self) -> None:
        self.so2 = (self.so2 - self.so2.mean()) / self.so2.std()
        self.rf = (self.rf - self.rf.mean()) / self.rf.std()
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
        return self.data.aligned_arrays["so2-start"]

    @cached_property
    def tau(self) -> np.ndarray:
        """Time axis for the deconvolution."""
        tau = self.so2.time.data - (
            self.so2.time.data[len(self.so2.time.data) // 2]
            - cftime.DatetimeNoLeap(0, 1, 1, has_year_zero=True, calendar="noleap")
        )
        return volcano_base.manipulate.dt2float(tau)

    @cached_property
    def rf(self) -> xr.DataArray:
        """Radiative forcing time series data."""
        return self.data.aligned_arrays["rf"]

    @cached_property
    def temp(self) -> xr.DataArray:
        """Temperature time series data."""
        return self.data.aligned_arrays["temperature"]

    @cached_property
    def temp_control(self) -> xr.DataArray:
        """Temperature time series data."""
        return xr.align(self.data.temperature_control, self.temp)[0]
