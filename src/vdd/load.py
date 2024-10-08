"""Functions that fetches data from files and returns it as an xarray.DataArray."""

import contextlib
import datetime
import warnings
from abc import ABC, abstractmethod
from collections.abc import Callable, Iterable
from enum import Enum
from functools import cached_property
from typing import Literal, Self

import cftime
import cosmoplots
import fppanalysis
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import rich.console
import rich.live
import rich.panel
import rich.progress
import scipy as sp
import scipy.signal as ssi
import volcano_base
import xarray as xr
from numpy.typing import NDArray
from pydantic import BaseModel, Field
from rich import print as rprint
from rich.console import Console
from rich.table import Table

import vdd.utils
from vdd.utils import name_swap as ns

type T_RF = tuple[Literal["temp"], Literal["rf"]]  # type: ignore[valid-type]
type T_SO2 = tuple[Literal["temp"], Literal["so2"]]  # type: ignore[valid-type]
type RF_SO2 = tuple[Literal["rf"], Literal["so2"]]  # type: ignore[valid-type]
type T_Strengths = Literal[  # type: ignore[valid-type]
    "strong",
    "medium",
    "medium-2sep",
    "medium-4sep",
    "medium-plus",
    "size5000",
    "tt-2sep",
    "tt-4sep",
]


console = Console()


def _convert_time_start_zero(arr: xr.DataArray) -> xr.DataArray:
    """Convert the time to start from zero."""
    arr = arr.assign_coords(time=volcano_base.manipulate.dt2float(arr.time.data) - 1850)
    arr.coords["time"].attrs["long_name"] = "Time-after-eruption"
    arr.coords["time"].attrs["units"] = "yr"
    return arr


class DataAttributeNotFoundError(FileNotFoundError):
    """Data is not found in the file."""

    def __init__(self: Self) -> None:
        super().__init__("Control data not found.")


class PaddingMethod(Enum):
    """Specify how the data should be padded.

    Attributes
    ----------
    NO : bool
        False, the data should not be padded.
    ZEROS : bool
        True, the data should be zero padded.
    NOISE : str
        The data should be padded using it control counterpart.
    """

    NO: bool = False
    ZEROS: bool = True
    NOISE: str = "noise"


class Normalise(Enum):
    """Specify how the data should be normalised.

    Attributes
    ----------
    AMPLITUDE : str
        Normalise the data to within 0 and 1.
    MEAN_STD : str
        Normalise the statistics of the data to have zero-mean and unit standard
        deviation.
    NO : bool
        Do not normalise or otherwise alter the data.
    """

    AMPLITUDE: str = "minmax"
    MEAN_STD: str = "mean-std"
    NO: bool = False


class Reconstructor(BaseModel):
    """Class that keeps a static set of arrays.

    The arrays represent time series that have already been deconvolved completely, thus
    we assume no further deconvolution analysis is needed.

    Attributes
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
    normalise : Normalise
        Whether the data is normalised.
    """

    class Config:
        """Configuration for the Reconstructor BaseModel object.

        Attributes
        ----------
        validate_assignment : bool
            Validate all assignments of the base class.
        strict : bool
            Be strict.
        arbitrary_types_allowed : bool
            Any types can be used.
        """

        validate_assignment: bool = True
        strict: bool = True
        arbitrary_types_allowed: bool = True

    tau: np.ndarray = Field(..., alias="tau", frozen=True)
    time: xr.DataArray = Field(..., alias="time", frozen=True)
    response_temp_so2: np.ndarray = Field(..., alias="response_temp_so2", frozen=True)
    response_temp_rf: np.ndarray = Field(..., alias="response_temp_rf", frozen=True)
    forcing: xr.DataArray = Field(..., alias="forcing", frozen=True)
    signal: xr.DataArray = Field(..., alias="signal", frozen=True)
    name: str = Field(..., frozen=False)
    normalise: Normalise = Field(..., frozen=True)


class CESMData(BaseModel):
    """Class to store the CESM2 data.

    Attributes
    ----------
    strength : T_Strengths
        The strength of the eruption, by default "strong".
    dims : list[str], optional
        The dimensions of the data to average over, by default ["lat", "lon"]. An empty
        list keeps all dimensions.
    """

    strength: T_Strengths = Field(default="strong", frozen=True)
    dims: list[str] = Field(default=["lat", "lon"], frozen=True)

    class Config:
        """Configuration for the CESMData BaseModel object.

        Attributes
        ----------
        validate_assignment : bool
            Validate all assignments of the base class.
        frozen : bool
            Do not allow attributes to change.
        extra : str
            Do not enable extra features.
        strict : bool
            Be strict.
        """

        validate_assignment: bool = True
        frozen: bool = True
        extra: str = "forbid"
        strict: bool = True

    @cached_property
    def so2(self: Self) -> xr.DataArray:  # noqa: C901
        """Get the CESM2 SO2 data.

        Returns
        -------
        xr.DataArray
            The SO2 data.

        Raises
        ------
        AttributeError
            If no data is found.
        """
        # This can probably be avoided by using the returns package, so that a
        # match/case statement could be used. But hey, it works. And this project is
        # very soon finished.(???)
        try:
            x = self._get_aod_cesm().time.data
        except AttributeError:
            try:
                x = self._get_trefht_cesm().time.data
            except AttributeError as e:
                msg = "No data found."
                raise AttributeError(msg) from e
        y = np.zeros_like(x)
        match self.strength:
            case "strong":
                y[0] = 1629
            case "medium-plus":
                y[0] = 400
            case "medium-2sep":
                y[0] = 26
                y[24] = 26
            case "medium-4sep":
                y[0] = 26
                y[48] = 26
            case "medium":
                y[0] = 26
            case "size5000":
                y[0] = 3000
            case "tt-2sep":
                y[0] = 400
                y[24] = 400
            case "tt-4sep":
                y[0] = 400
                y[48] = 400
            case _:
                vdd.utils.never_called(self.strength)
        out = xr.DataArray(y, coords={"time": x}, dims=["time"])
        out = out if len(out) % 2 else out[:-1]
        return self._align_arrays("so2", out)

    @cached_property
    def tmso2(self: Self) -> xr.DataArray:
        """Get the CESM2 TMSO2 data.

        Returns
        -------
        xr.DataArray
            The TMSO2 data.
        """
        out = self._get_tmso2_cesm()
        out = out if len(out) % 2 else out[:-1]
        return self._align_arrays("tmso2", out)

    @cached_property
    def tmso2_control(self: Self) -> xr.DataArray:
        """Get the CESM2 control TMSO2 data.

        Returns
        -------
        xr.DataArray
            The control TMSO2 data.
        """
        if not hasattr(self, "_tmso2_control"):
            _ = self.tmso2
        out = self._tmso2_control
        out = out if len(out) % 2 else out[:-1]
        return self._align_arrays("tmso2_control", out)

    @cached_property
    def aod(self: Self) -> xr.DataArray:
        """Get the CESM2 AOD data.

        Returns
        -------
        xr.DataArray
            The AOD data.
        """
        out = self._get_aod_cesm()
        out = out if len(out) % 2 else out[:-1]
        return self._align_arrays("aod", out)

    @cached_property
    def aod_control(self: Self) -> xr.DataArray:
        """Get the CESM2 control AOD data.

        Returns
        -------
        xr.DataArray
            The control AOD data.
        """
        if not hasattr(self, "_aod_control"):
            _ = self.aod
        out = self._aod_control
        out = out if len(out) % 2 else out[:-1]
        return self._align_arrays("aod_control", out)

    @cached_property
    def rf(self: Self) -> xr.DataArray:
        """Get the CESM2 RF data.

        Returns
        -------
        xr.DataArray
            The RF data.
        """
        out = self._get_rf_cesm()
        out = out if len(out) % 2 else out[:-1]
        return self._align_arrays("rf", out)

    @cached_property
    def rf_control(self: Self) -> xr.DataArray:
        """Get the CESM2 control RF data.

        Returns
        -------
        xr.DataArray
            The control RF data.
        """
        if not hasattr(self, "_rf_control"):
            _ = self.rf
        out = self._rf_control
        out = out if len(out) % 2 else out[:-1]
        return self._align_arrays("rf_control", out)

    @cached_property
    def temp(self: Self) -> xr.DataArray:
        """Get the CESM2 temperature data.

        Returns
        -------
        xr.DataArray
            The temperature data.
        """
        out = self._get_trefht_cesm()
        out = out if len(out) % 2 else out[:-1]
        return self._align_arrays("temp", out)

    @cached_property
    def temperature_control(self: Self) -> xr.DataArray:
        """Get the CESM2 control temperature data.

        Returns
        -------
        xr.DataArray
            The control temperature data.
        """
        if not hasattr(self, "_temperature_control"):
            _ = self.temp
        out = self._temperature_control
        out = out if len(out) % 2 else out[:-1]
        return self._align_arrays("temperature_control", out)

    def initialise_data(self: Self) -> None:
        """Initialise the data, ensuring that it is loaded and aligned."""
        _ = (
            self.so2,
            self.tmso2,
            self.aod,
            self.rf,
            self.temp,
            self.temperature_control,
            self.rf_control,
            self.aod_control,
            self.tmso2_control,
        )

    def _align_arrays(self: Self, new: str, new_obj: xr.DataArray) -> xr.DataArray:  # noqa: ARG002
        if out := list(
            set(self.__dict__.keys())
            & {
                "temperature_control",
                "rf_control",
                "aod_control",
                "tmso2_control",
                "temp",
                "tmso2",
                "so2",
                "rf",
                "aod",
            },
        ):
            aligned = xr.align(new_obj, *[getattr(self, o) for o in out])
            self.__dict__.update(dict(zip(out, aligned[1:], strict=True)))
        else:
            aligned = (new_obj,)
        return aligned[0]

    def _get_tmso2_cesm(self: Self, *, plot_example: bool = False) -> xr.DataArray:
        """Get the CESM2 TMSO2 data."""
        data = (
            volcano_base.load.FindFiles()
            .find("e_fSST1850", "TMSO2", "h0", self.strength)
            .remove(
                "ens1"
                if self.strength
                not in {"tt-2sep", "tt-4sep", "medium-2sep", "medium-4sep"}
                else "ens0",
            )
            .keep_most_recent()
        )
        files = data.load()
        files = volcano_base.manipulate.mean_flatten(files, dims=self.dims)
        files = volcano_base.manipulate.shift_arrays(files, daily=False)
        files = volcano_base.manipulate.shift_arrays(files, custom=1)
        files = volcano_base.manipulate.subtract_mean_of_tail(files)
        files = volcano_base.manipulate.data_array_operation(
            files,
            _convert_time_start_zero,
        )
        files = list(xr.align(*files))
        # TMSO2 control data
        control_data: volcano_base.load.FindFiles = (
            volcano_base.load.FindFiles()
            .find("e_fSST1850", "h0", "control", "TMSO2", "ens1")
            .sort("ensemble", "attr")
        )
        if len(control_data) != 1:
            raise DataAttributeNotFoundError
        control_l = control_data.load()
        control_l = volcano_base.manipulate.mean_flatten(control_l, dims=self.dims)
        control_l = volcano_base.manipulate.data_array_operation(
            control_l,
            _convert_time_start_zero,
        )
        # Climatology
        control_float = volcano_base.manipulate.get_median(
            control_l,
            xarray=True,
        ).dropna("time")
        control_dt = control_float.assign_coords(
            time=volcano_base.manipulate.float2dt(control_float.time, freq="MS"),
        )
        control_dt = volcano_base.manipulate.subtract_climatology(
            control_dt,
            control_dt,
            "time.month",
        )[0]
        self._tmso2_control = control_dt.assign_coords(time=control_float.time)
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

    def _get_aod_cesm(self: Self, *, plot_example: bool = False) -> xr.DataArray:
        """Get the CESM2 AOD data."""
        data = (
            volcano_base.load.FindFiles()
            .find("e_fSST1850", "AODVISstdn", "h0", self.strength)
            .remove(
                "ens1"
                if self.strength
                not in {"tt-2sep", "tt-4sep", "medium-2sep", "medium-4sep"}
                else "ens0",
            )
            .keep_most_recent()
        )
        files = data.load()
        files = volcano_base.manipulate.mean_flatten(files, dims=self.dims)
        files = volcano_base.manipulate.shift_arrays(files, daily=False)
        files = volcano_base.manipulate.shift_arrays(files, custom=1)
        files = volcano_base.manipulate.subtract_mean_of_tail(files)
        files = volcano_base.manipulate.data_array_operation(
            files,
            _convert_time_start_zero,
        )
        files = list(xr.align(*files))
        # AOD control data
        control_data: volcano_base.load.FindFiles = (
            volcano_base.load.FindFiles()
            .find("e_fSST1850", "h0", "control", "AODVISstdn", "ens1")
            .sort("ensemble", "attr")
        )
        if len(control_data) != 1:
            raise DataAttributeNotFoundError
        control_l = control_data.load()
        control_l = volcano_base.manipulate.mean_flatten(control_l, dims=self.dims)
        control_l = volcano_base.manipulate.data_array_operation(
            control_l,
            _convert_time_start_zero,
        )
        # Climatology
        control_float = volcano_base.manipulate.get_median(
            control_l,
            xarray=True,
        ).dropna("time")
        control_dt = control_float.assign_coords(
            time=volcano_base.manipulate.float2dt(control_float.time, freq="MS"),
        )
        control_dt = volcano_base.manipulate.subtract_climatology(
            control_dt,
            control_dt,
            "time.month",
        )[0]
        self._aod_control = control_dt.assign_coords(time=control_float.time)
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
        self: Self,
    ) -> tuple[xr.DataArray, xr.DataArray, list[xr.DataArray], list[xr.DataArray]]:
        control_data: volcano_base.load.FindFiles = (
            volcano_base.load.FindFiles()
            .find("e_fSST1850", "h0", "control", {"FLNT", "FSNT"}, "ens1")
            .sort("ensemble", "attr")
        )
        data: volcano_base.load.FindFiles = (
            volcano_base.load.FindFiles()
            .find("e_fSST1850", {"FLNT", "FSNT"}, "h0", self.strength)
            .remove(
                "ens1"
                if self.strength
                not in {"tt-2sep", "tt-4sep", "medium-2sep", "medium-4sep"}
                else "ens0",
            )
            .remove(
                *["ens2", "ens4"]
                if self.strength in {"tt-2sep", "tt-4sep", "medium-2sep", "medium-4sep"}
                else "ens6",
            )
            .keep_most_recent()
            .sort("ensemble", "attr")
        )
        c_size = 1
        match (control_data.copy().keep("FLNT"), control_data.copy().keep("FSNT")):
            case (c_flnt, c_fsnt) if len(c_flnt) == c_size and len(c_fsnt) == c_size:
                c_flnt_xr = c_flnt.load()[0]
                c_fsnt_xr = c_fsnt.load()[0]
            case _:
                raise DataAttributeNotFoundError
        match self.strength:
            case "size5000" | "tt-2sep" | "tt-4sep" | "medium-2sep" | "medium-4sep":
                d_size = 2
            case _:
                d_size = 4
        match (data.copy().keep("FLNT"), data.copy().keep("FSNT")):
            case (flnt, fsnt) if len(flnt) == d_size and len(fsnt) == d_size:
                flnt_xr = flnt.load()
                fsnt_xr = fsnt.load()
            case _:
                msg = f"Data for {self.strength} not found."
                raise ValueError(msg)
        c_fsnt_xr, c_flnt_xr = volcano_base.manipulate.mean_flatten(
            [c_fsnt_xr, c_flnt_xr],
            dims=self.dims,
        )
        fsnt_xr = volcano_base.manipulate.mean_flatten(fsnt_xr, dims=self.dims)
        flnt_xr = volcano_base.manipulate.mean_flatten(flnt_xr, dims=self.dims)
        return c_fsnt_xr, c_flnt_xr, fsnt_xr, flnt_xr

    def _get_rf_cesm(self: Self, *, plot_example: bool = False) -> xr.DataArray:
        """Get the CESM2 RF data."""
        c_fsnt_xr, c_flnt_xr, fsnt_xr, flnt_xr = self._load_rf_cesm()

        def remove_control(rf: xr.DataArray) -> xr.DataArray:
            rf, c_s, c_l = xr.align(rf, c_fsnt_xr, c_flnt_xr)
            rf.data -= c_s.data - c_l.data
            # Combine attributes, overwrite if needed
            original_attrs = rf.attrs
            new_attrs = {
                "attr": "RF",
                "long_name": "Net radiative flux (residual) at top-of-model",
            }
            combined_attrs = {**original_attrs, **new_attrs}
            return rf.assign_attrs(**combined_attrs).rename("RESTOM")

        def difference(
            short: list[xr.DataArray],
            long: list[xr.DataArray],
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
        # Create control run time series
        c_total = c_fsnt_xr.copy()
        c_total.data -= c_flnt_xr.data
        c_total = c_total.assign_attrs(attr="RF")
        c_total = volcano_base.manipulate.subtract_mean_of_tail([c_total])[0]
        c_total = volcano_base.manipulate.data_array_operation(
            [c_total],
            _convert_time_start_zero,
        )[0]
        c_total.data = np.asarray(c_total.data)
        # Fourier
        # self._rf_control = volcano_base.manipulate.remove_seasonality(
        #     c_total.dropna("time"), freq=1.09, radius=0.01
        # )
        # Climatology
        c_total_float = c_total.dropna("time")
        c_total_dt = c_total_float.assign_coords(
            time=volcano_base.manipulate.float2dt(c_total_float.time, freq="MS"),
        )
        c_total_dt = volcano_base.manipulate.subtract_climatology(
            c_total_dt,
            c_total_dt,
            "time.month",
        )[0]
        self._rf_control = c_total_dt.assign_coords(time=c_total_float.time)
        if plot_example:
            plt.figure()
            [f.plot() for f in rf]  # type: ignore[call-arg]
            plt.show()
        mean_array = volcano_base.manipulate.get_median(rf, xarray=True)
        if plot_example:
            plt.figure()
            plt.plot(mean_array)
            plt.show()
        return mean_array.dropna("time")

    def _get_trefht_cesm(self: Self, *, plot_example: bool = False) -> xr.DataArray:
        """Get the CESM2 temperature data."""
        control_data: volcano_base.load.FindFiles = (
            volcano_base.load.FindFiles()
            .find("e_BWma1850", "h0", "control", "TREFHT", "ens0")
            .sort("ensemble", "attr")
        )
        if len(control_data) != 1:
            raise DataAttributeNotFoundError
        control_l = control_data.load()
        control_l = volcano_base.manipulate.mean_flatten(control_l, dims=self.dims)
        data = (
            volcano_base.load.FindFiles()
            .find("e_BWma1850", "TREFHT", "h0", self.strength)
            .remove(
                "ens1"
                if self.strength
                not in {"tt-2sep", "tt-4sep", "medium-2sep", "medium-4sep"}
                else "ens0",
            )
            .keep_most_recent()
        )
        files = data.load()
        orig_attrs = files[0].attrs
        files = volcano_base.manipulate.mean_flatten(files, dims=self.dims)
        control = control_l[0]
        files = volcano_base.manipulate.shift_arrays(files, daily=False)
        files = volcano_base.manipulate.shift_arrays(files, custom=1)

        def subtract_control_array(arr: xr.DataArray) -> xr.DataArray:
            arr, c = xr.align(arr, control)
            arr.data -= c.data
            return arr

        def subtract_control_mean(arr: xr.DataArray) -> xr.DataArray:
            return arr - volcano_base.config.MEANS["TREFHT"]

        control_l = volcano_base.manipulate.data_array_operation(
            control_l,
            _convert_time_start_zero,
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
            control_l,
            xarray=True,
        ).dropna("time")
        control_dt = control_float.assign_coords(
            time=volcano_base.manipulate.float2dt(control_float.time, freq="MS"),
        )
        control_dt = volcano_base.manipulate.subtract_climatology(
            control_dt,
            control_dt,
            "time.month",
        )[0]
        self._temperature_control = control_dt.assign_coords(time=control_float.time)
        subtract_control = (
            subtract_control_array if len(data) == 1 else subtract_control_mean
        )
        files = volcano_base.manipulate.data_array_operation(files, subtract_control)
        files = volcano_base.manipulate.data_array_operation(
            files,
            _convert_time_start_zero,
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


class _PostInitCaller[T](ABC, type):  # type: ignore[valid-type]
    def __call__(cls, *args, **kwargs) -> T:  # type: ignore[name-defined] # noqa: ANN002,ANN003
        obj: T = type.__call__(cls, *args, **kwargs)  # type: ignore[name-defined]
        obj.__post_init__()
        return obj


class EvenLengthError(Exception):
    """Exception raised for even length arrays."""

    def __init__(
        self: Self, message: str = "The arrays must have an odd length."
    ) -> None:
        self.message = message
        super().__init__(self.message)


class Deconvolve(metaclass=_PostInitCaller):
    """Class for deconvolving data.

    Parameters
    ----------
    normalise : Normalise
        A Normalise enum specifying how normalisation should be done.

    Attributes
    ----------
    name : str
        The name of the class.
    progress : bool
        Show progress bars.
    """

    name: str = "Deconvolve"
    progress: bool = True

    def __init__(self: Self, normalise: Normalise = Normalise.NO) -> None:
        def _deconv(
            signal: NDArray[np.float64],
            forcing: NDArray[np.float64],
        ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
            if not len(signal) % 2 or not len(forcing) % 2:
                raise EvenLengthError
            guess = np.heaviside(np.arange(len(signal)) - len(signal) // 2, 1)
            n_iters = 2000
            if self.progress:
                with (
                    console.status("[bold yellow]Deconvolving ...", spinner="point"),
                    contextlib.redirect_stdout(None),
                ):
                    out, err = fppanalysis.RL_gauss_deconvolve(
                        signal,
                        forcing,
                        initial_guess=guess,
                        iteration_list=n_iters,
                    )
            else:
                out, err = fppanalysis.RL_gauss_deconvolve(
                    signal,
                    forcing,
                    initial_guess=guess,
                    iteration_list=n_iters,
                )
            return out, err

        self._deconvolve: Callable[
            [NDArray[np.float64], NDArray[np.float64]],
            tuple[NDArray[np.float64], NDArray[np.float64]],
        ] = _deconv
        self.normalise = normalise

    def __post_init__(self: Self) -> None:
        """Run the normalisation if needed."""
        if self.normalise:
            self._update_if_normalise()

    @abstractmethod
    def _update_if_normalise(self: Self) -> None: ...

    @property
    @abstractmethod
    def tau(self: Self) -> NDArray[np.float64]:
        """Time axis for the deconvolution."""
        raise NotImplementedError

    @property
    @abstractmethod
    def so2(self: Self) -> xr.DataArray:
        """SO2 time series data."""
        raise NotImplementedError

    @property
    @abstractmethod
    def so2_decay(self: Self) -> xr.DataArray:
        """SO2 time series data."""
        raise NotImplementedError

    @property
    @abstractmethod
    def rf(self: Self) -> xr.DataArray:
        """Radiative forcing time series data."""
        raise NotImplementedError

    @property
    @abstractmethod
    def temp(self: Self) -> xr.DataArray:
        """Temperature time series data."""
        raise NotImplementedError

    @property
    @abstractmethod
    def temp_control(self: Self) -> xr.DataArray:
        """Temperature control time series data."""
        raise NotImplementedError

    def deconvolve(
        self: Self, signal: NDArray[np.float64], forcing: NDArray[np.float64]
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Deconvolve the signal with the forcing."""
        return self._deconvolve(signal, forcing)

    def dump_reconstructor(self: Self) -> Reconstructor:
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
        self: Self,
        method: Callable[[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]],
        method_args: list | None = None,
        method_kwargs: dict | None = None,
    ) -> None:
        """Change the deconvolution method.

        Parameters
        ----------
        method : Callable[[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]]
            The new deconvolution method.
        method_args : list | None
            The positional arguments to pass to the method, if any. Not including the
            signal and forcing.
        method_kwargs : dict | None
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
        >>> dec.change_deconvolution_method(fppanalysis.RL_gauss_deconvolve, method_kwargs=kwargs)

        Or even
        >>> dec.change_deconvolution_method(
        ...     fppanalysis.RL_gauss_deconvolve, {"iteration_list": 1000}
        ... )
        """

        def method_(
            signal: np.ndarray,
            forcing: np.ndarray,
        ) -> tuple[np.ndarray, np.ndarray]:
            match method_args, method_kwargs:
                case None, None:
                    return method(signal, forcing)
                case _, None:
                    return method(signal, forcing, *method_args)
                case None, _:
                    return method(signal, forcing, **method_kwargs)
                case _:
                    raise ValueError

        self._deconvolve = method_

    @cached_property
    def _response_rf_so2_tup(self: Self) -> tuple[np.ndarray, np.ndarray]:
        """Deconvolve the RF signal with the SO2 signal."""
        _ = self.tau
        signal, err = self.deconvolve(self.rf.data, self.so2.data)
        match self.normalise:
            case Normalise.AMPLITUDE:
                signal = vdd.utils.normalise(signal)
            case Normalise.MEAN_STD:
                signal = (signal - signal.mean()) / signal.std()
        return signal, err

    @property
    def response_rf_so2(self: Self) -> np.ndarray:
        """Deconvolve the RF signal with the SO2 delta signal.

        Returns
        -------
        np.ndarray
            The deconvolved RF to SO2 signal.
        """
        return self._response_rf_so2_tup[0].flatten()

    @property
    def response_rf_so2_err(self: Self) -> np.ndarray:
        """Deconvolve the RF signal with the SO2 signal."""
        return self._response_rf_so2_tup[1]

    @cached_property
    def _response_temp_so2_tup(self: Self) -> tuple[np.ndarray, np.ndarray]:
        """Deconvolve the temperature signal with the SO2 signal."""
        _ = self.tau
        signal, err = self.deconvolve(self.temp.data, self.so2.data)
        match self.normalise:
            case Normalise.AMPLITUDE:
                signal = vdd.utils.normalise(signal)
            case Normalise.MEAN_STD:
                signal = (signal - signal.mean()) / signal.std()
        return signal, err

    @property
    def response_temp_so2(self: Self) -> np.ndarray:
        """Deconvolve the temperature signal with the SO2 signal.

        Returns
        -------
        np.ndarray
            The deconvolved temperature to SO2 signal.
        """
        return self._response_temp_so2_tup[0].flatten()

    @property
    def response_temp_so2_err(self: Self) -> np.ndarray:
        """Deconvolve the temperature signal with the SO2 signal."""
        return self._response_temp_so2_tup[1]

    @cached_property
    def _response_rf_so2_decay_tup(self: Self) -> tuple[np.ndarray, np.ndarray]:
        """Deconvolve the RF signal with the SO2 decay signal."""
        _ = self.tau
        signal, err = self.deconvolve(self.rf.data, self.so2_decay.data)
        match self.normalise:
            case Normalise.AMPLITUDE:
                signal = vdd.utils.normalise(signal)
            case Normalise.MEAN_STD:
                signal = (signal - signal.mean()) / signal.std()
        return signal, err

    @property
    def response_rf_so2_decay(self: Self) -> np.ndarray:
        """Deconvolve the RF signal with the SO2 decay signal.

        Returns
        -------
        np.ndarray
            The deconvolved RF to SO2 signal.
        """
        return self._response_rf_so2_decay_tup[0].flatten()

    @property
    def response_rf_so2_decay_err(self: Self) -> np.ndarray:
        """Deconvolve the RF signal with the SO2 decay signal."""
        return self._response_rf_so2_decay_tup[1]

    @cached_property
    def _response_temp_so2_decay_tup(self: Self) -> tuple[np.ndarray, np.ndarray]:
        """Deconvolve the temperature signal with the SO2 decay signal."""
        _ = self.tau
        signal, err = self.deconvolve(self.temp.data, self.so2_decay.data)
        match self.normalise:
            case Normalise.AMPLITUDE:
                signal = vdd.utils.normalise(signal)
            case Normalise.MEAN_STD:
                signal = (signal - signal.mean()) / signal.std()
        return signal, err

    @property
    def response_temp_so2_decay(self: Self) -> np.ndarray:
        """Deconvolve the temperature signal with the SO2 decay signal.

        Returns
        -------
        np.ndarray
            The deconvolved temperature to SO2 decay signal.
        """
        return self._response_temp_so2_decay_tup[0].flatten()

    @property
    def response_temp_so2_decay_err(self: Self) -> np.ndarray:
        """Deconvolve the temperature signal with the SO2 decay signal."""
        return self._response_temp_so2_decay_tup[1]

    @cached_property
    def _response_temp_rf_tup(self: Self) -> tuple[np.ndarray, np.ndarray]:
        """Deconvolve the temperature signal with the RF signal."""
        _ = self.tau
        signal, err = self.deconvolve(self.temp.data, self.rf.data)
        match self.normalise:
            case Normalise.AMPLITUDE:
                signal = vdd.utils.normalise(signal)
            case Normalise.MEAN_STD:
                signal = (signal - signal.mean()) / signal.std()
        return signal, err

    @property
    def response_temp_rf(self: Self) -> np.ndarray:
        """Deconvolve the temperature signal with the RF signal.

        Returns
        -------
        np.ndarray
            The deconvolved temperature to RF signal.
        """
        return self._response_temp_rf_tup[0].flatten()

    @property
    def response_temp_rf_err(self: Self) -> np.ndarray:
        """Deconvolve the temperature signal with the RF signal."""
        return self._response_temp_rf_tup[1]

    def plot_dec_rf_with_so2(self: Self) -> None:
        """Deconvolve the RF signal with the SO2 signal."""
        # Quick comparison
        plt.figure()
        vdd.utils.normalise(self.rf).plot()  # type: ignore[call-arg]
        vdd.utils.normalise(self.so2).plot()  # type: ignore[call-arg]
        plt.figure()
        plt.plot(self.tau, self.response_rf_so2)
        plt.figure()
        plt.semilogy(self.response_rf_so2_err)

    def plot_dec_temp_with_so2(self: Self) -> None:
        """Deconvolve the temperature signal with the SO2 signal."""
        # Quick comparison
        plt.figure()
        vdd.utils.normalise(self.temp).plot()  # type: ignore[call-arg]
        vdd.utils.normalise(self.so2).plot()  # type: ignore[call-arg]
        plt.figure()
        plt.plot(self.tau, self.response_temp_so2)
        plt.figure()
        plt.semilogy(self.response_temp_so2_err)

    def plot_dec_temp_with_rf(self: Self) -> None:
        """Deconvolve the temperature signal with the SO2 signal."""
        # Quick comparison
        plt.figure()
        vdd.utils.normalise(self.temp).plot()  # type: ignore[call-arg]
        vdd.utils.normalise(self.rf).plot()  # type: ignore[call-arg]
        plt.figure()
        plt.plot(self.tau, self.response_temp_rf)
        plt.figure()
        plt.semilogy(self.response_temp_rf_err)


class DeconvolveCESM(Deconvolve):
    """Class for deconvolving data from CESM2.

    Parameters
    ----------
    normalise : Normalise, optional
        Whether to normalise the data, by default false.
    pad_before : PaddingMethod, optional
        Whether to pad with zeros all time series before the convolution, by default
        False.
    cesm : CESMData | None, optional
        The CESM data to use. If not given, the "strong" strength data will be used.
    """

    def __init__(
        self: Self,
        normalise: Normalise = Normalise.NO,
        pad_before: PaddingMethod = PaddingMethod.NO,
        cesm: CESMData | None = None,
    ) -> None:
        super().__init__(normalise)
        self._data = CESMData() if cesm is None else cesm
        # Since the time series can sometimes be of different length, we make sure to
        # call all of them before using them here. That way, they will be updated within
        # the CESMData class before assigning them here.
        self.pad_before = pad_before
        self.name = f"CESM2 {self._data.strength}"

    def _update_if_normalise(self: Self) -> None:
        match self.normalise:
            case Normalise.MEAN_STD:
                self._normalise_mean_std()
            case Normalise.AMPLITUDE:
                self._normalise_amplitude()

    def _normalise_amplitude(self: Self) -> None:
        self.so2 = vdd.utils.normalise(self.so2)
        self.aod = vdd.utils.normalise(self.aod)
        self.rf = vdd.utils.normalise(self.rf)
        self.temp = vdd.utils.normalise(self.temp)

    def _normalise_mean_std(self: Self) -> None:
        self.so2 = (self.so2 - self.so2.mean()) / self.so2.std()
        self.tmso2 = (self.tmso2 - self.tmso2.mean()) / self.tmso2.std()
        self.aod = (self.aod - self.aod.mean()) / self.aod.std()
        self.rf = (self.rf - self.rf.mean()) / self.rf.std()
        self.temp = (self.temp - self.temp.mean()) / self.temp.std()
        self.temp_control = (
            self.temp_control - self.temp_control.mean()
        ) / self.temp_control.std()

    @cached_property
    def so2(self: Self) -> xr.DataArray:
        """SO2 time series data."""
        return (
            vdd.utils.pad_before_convolution(self._data.so2)
            if self.pad_before
            else self._data.so2
        )

    @cached_property
    def so2_decay(self: Self) -> xr.DataArray:
        """SO2 time series data."""
        return (
            vdd.utils.pad_before_convolution(self._data.tmso2)
            if self.pad_before
            else self._data.tmso2
        )

    @cached_property
    def tau(self: Self) -> np.ndarray:
        """Time axis for the deconvolution."""
        return (
            self.so2.time.data
            if self.pad_before
            else self.so2.time.data - self.so2.time.data[len(self.so2.time.data) // 2]
        )

    @cached_property
    def tmso2(self: Self) -> xr.DataArray:
        """SO2 column burden."""
        pad_values = (
            self._data.tmso2_control.data[:-1]
            if self.pad_before == PaddingMethod.NOISE
            else None
        )
        return (
            vdd.utils.pad_before_convolution(self._data.tmso2, pad_values=pad_values)
            if self.pad_before
            else self._data.tmso2
        )

    @cached_property
    def tmso2_control(self: Self) -> xr.DataArray:
        """TMSO2 time series data."""
        pad_values = (
            self._data.tmso2_control.data[:-1]
            if self.pad_before == PaddingMethod.NOISE
            else None
        )
        return (
            vdd.utils.pad_before_convolution(
                self._data.tmso2_control,
                pad_values=pad_values,
            )
            if self.pad_before
            else self._data.tmso2_control
        )

    @cached_property
    def aod(self: Self) -> xr.DataArray:
        """Aerosol optical depth time series data."""
        pad_values = (
            self._data.aod_control.data[:-1]
            if self.pad_before == PaddingMethod.NOISE
            else None
        )
        return (
            vdd.utils.pad_before_convolution(self._data.aod, pad_values=pad_values)
            if self.pad_before
            else self._data.aod
        )

    @cached_property
    def aod_control(self: Self) -> xr.DataArray:
        """AOD time series data."""
        pad_values = (
            self._data.aod_control.data[:-1]
            if self.pad_before == PaddingMethod.NOISE
            else None
        )
        return (
            vdd.utils.pad_before_convolution(
                self._data.aod_control,
                pad_values=pad_values,
            )
            if self.pad_before
            else self._data.aod_control
        )

    @cached_property
    def rf(self: Self) -> xr.DataArray:
        """Radiative forcing time series data."""
        pad_values = (
            self._data.rf_control.data[:-1]
            if self.pad_before == PaddingMethod.NOISE
            else None
        )
        attrs = self._data.rf.attrs
        return (
            vdd.utils.pad_before_convolution(self._data.rf * -1, pad_values=pad_values)
            if self.pad_before
            else self._data.rf * -1
        ).assign_attrs(**attrs)

    @cached_property
    def rf_control(self: Self) -> xr.DataArray:
        """RF time series data."""
        pad_values = (
            self._data.rf_control.data[:-1]
            if self.pad_before == PaddingMethod.NOISE
            else None
        )
        attrs = self._data.rf.attrs
        return (
            vdd.utils.pad_before_convolution(
                self._data.rf_control * -1,
                pad_values=pad_values,
            )
            if self.pad_before
            else self._data.rf_control * -1
        ).assign_attrs(**attrs)

    @cached_property
    def temp(self: Self) -> xr.DataArray:
        """Temperature time series data."""
        pad_values = (
            self._data.temperature_control.data[:-1]
            if self.pad_before == PaddingMethod.NOISE
            else None
        )
        attrs = self._data.temp.attrs
        return (
            vdd.utils.pad_before_convolution(
                self._data.temp * -1,
                pad_values=pad_values,
            )
            if self.pad_before
            else self._data.temp * -1
        ).assign_attrs(**attrs)

    @cached_property
    def temp_control(self: Self) -> xr.DataArray:
        """Temperature control time series data."""
        pad_values = (
            self._data.temperature_control.data[:-1]
            if self.pad_before == PaddingMethod.NOISE
            else None
        )
        attrs = self._data.temp.attrs
        return (
            vdd.utils.pad_before_convolution(
                self._data.temperature_control * -1,
                pad_values=pad_values,
            )
            if self.pad_before
            else self._data.temperature_control * -1
        ).assign_attrs(**attrs)

    @cached_property
    def _response_aod_so2_tup(self: Self) -> tuple[np.ndarray, np.ndarray]:
        """Deconvolve the AOD signal with the SO2 signal."""
        _ = self.tau
        signal, err = self.deconvolve(self.aod.data, self.so2.data)
        match self.normalise:
            case Normalise.AMPLITUDE:
                signal = vdd.utils.normalise(signal)
            case Normalise.MEAN_STD:
                signal = (signal - signal.mean()) / signal.std()
        return signal, err

    @property
    def response_aod_so2(self: Self) -> np.ndarray:
        """Deconvolve the AOD signal with the SO2 signal.

        Returns
        -------
        np.ndarray
            The deconvolved AOD to SO2 signal.
        """
        return self._response_aod_so2_tup[0].flatten()

    @property
    def response_aod_so2_err(self: Self) -> np.ndarray:
        """Deconvolve the AOD signal with the SO2 signal."""
        return self._response_aod_so2_tup[1]

    @cached_property
    def _response_rf_aod_tup(self: Self) -> tuple[np.ndarray, np.ndarray]:
        """Deconvolve the RF signal with the AOD signal."""
        _ = self.tau
        signal, err = self.deconvolve(self.rf.data, self.aod.data)
        match self.normalise:
            case Normalise.AMPLITUDE:
                signal = vdd.utils.normalise(signal)
            case Normalise.MEAN_STD:
                signal = (signal - signal.mean()) / signal.std()
        return signal, err

    @property
    def response_rf_aod(self: Self) -> np.ndarray:
        """Deconvolve the RF signal with the AOD signal.

        Returns
        -------
        np.ndarray
            The deconvolved RF to AOD signal.
        """
        return self._response_rf_aod_tup[0].flatten()

    @property
    def response_rf_aod_err(self: Self) -> np.ndarray:
        """Deconvolve the RF signal with the AOD signal."""
        return self._response_rf_aod_tup[1]

    @cached_property
    def _response_temp_aod_tup(self: Self) -> tuple[np.ndarray, np.ndarray]:
        """Deconvolve the temperature signal with the AOD signal."""
        _ = self.tau
        signal, err = self.deconvolve(self.temp.data, self.aod.data)
        match self.normalise:
            case Normalise.AMPLITUDE:
                signal = vdd.utils.normalise(signal)
            case Normalise.MEAN_STD:
                signal = (signal - signal.mean()) / signal.std()
        return signal, err

    @property
    def response_temp_aod(self: Self) -> np.ndarray:
        """Deconvolve the temperature signal with the AOD signal.

        Returns
        -------
        np.ndarray
            The deconvolved temperature to AOD signal.
        """
        return self._response_temp_aod_tup[0].flatten()

    @property
    def response_temp_aod_err(self: Self) -> np.ndarray:
        """Deconvolve the temperature signal with the AOD signal."""
        return self._response_temp_aod_tup[1]


class DeconvolveOB16(Deconvolve):
    """Class for deconvolving data from Otto-Bliesner et at. (2016).

    Parameters
    ----------
    data : volcano_base.load.OttoBliesner | Literal["h0", "h1"], optional
        The OB16 data class to use. If not given, the h1 (daily) data will be loaded.
    normalise : Normalise, optional
        Whether to normalise the data, by default False.
    length : int | None, optional
        After 1850, the SO2 dataset mismatches with the simulation output. This
        parameter specifies how many items should be included in all the arrays. Default
        is 12001, which means data up to 1000 years will be used.
    """

    def __init__(
        self: Self,
        data: volcano_base.load.OttoBliesner | Literal["h0", "h1"] = "h1",
        normalise: Normalise = Normalise.NO,
        length: int | None = None,
    ) -> None:
        super().__init__(normalise)
        match data:
            case "h0":
                self.data = volcano_base.load.OttoBliesner(freq="h0", progress=True)
                length = 12001 if length is None else length
            case "h1":
                self.data = volcano_base.load.OttoBliesner(freq="h1", progress=True)
                length = 365001 if length is None else length
            case volcano_base.load.OttoBliesner():
                self.data = data
            case _:
                vdd.utils.never_called(data)
        self.start_pt = 0
        self.end_pt: int | None = None
        if length:
            self.end_pt = length if (length - self.start_pt) % 2 else length + 1

    @staticmethod
    def _find_good_endpoints(arr: xr.DataArray) -> xr.DataArray:
        range_ = arr.data.max() - arr.data.min()
        # Make the start and end be equal to within 10% of the range, and make the next
        # elements have equals difference to the previous within 10%.
        start = arr.data[0]
        i_ = 0
        for i, val in enumerate(arr.data[::-1]):
            i_ = i
            if abs(val - start) < range_ * 0.1:
                break
        return arr[: len(arr) - i_]

    def _update_if_normalise(self: Self) -> None:
        match self.normalise:
            case Normalise.MEAN_STD:
                self._normalise_mean_std()
            case Normalise.AMPLITUDE:
                self._normalise_amplitude()

    def _normalise_amplitude(self: Self) -> None:
        self.so2 = vdd.utils.normalise(self.so2)
        self.rf = vdd.utils.normalise(self.rf)
        self.temp = vdd.utils.normalise(self.temp)

    def _normalise_mean_std(self: Self) -> None:
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

    @cached_property
    def so2(self: Self) -> xr.DataArray:
        """SO2 time series data."""
        return self.data.aligned_arrays["so2-start"][self.start_pt : self.end_pt]

    @cached_property
    def so2_decay(self: Self) -> xr.DataArray:
        """SO2 time series data."""
        return (
            self.data.aligned_arrays["so2-decay-start"][self.start_pt : self.end_pt]
            / 510e3
        )

    @cached_property
    def tau(self: Self) -> np.ndarray:
        """Time axis for the deconvolution."""
        tau = self.so2.time.data - (
            self.so2.time.data[len(self.so2.time.data) // 2]
            - cftime.DatetimeNoLeap(0, 1, 1, has_year_zero=True, calendar="noleap")
        )
        return np.asarray(volcano_base.manipulate.dt2float(tau))

    @cached_property
    def rf(self: Self) -> xr.DataArray:
        """Radiative forcing time series data."""
        return self.data.aligned_arrays["rf"][self.start_pt : self.end_pt]

    @cached_property
    def rf_control(self: Self) -> xr.DataArray:
        """RF time series data."""
        return xr.align(self.data.rf_control, self.rf)[0][self.start_pt : self.end_pt]

    @cached_property
    def temp(self: Self) -> xr.DataArray:
        """Temperature time series data."""
        return self.data.aligned_arrays["temperature"][self.start_pt : self.end_pt]

    @cached_property
    def temp_control(self: Self) -> xr.DataArray:
        """Temperature time series data."""
        return xr.align(self.data.temperature_control, self.temp)[0][
            self.start_pt : self.end_pt
        ]


class CutOff:
    """Cut off the response functions of a deconvolution object."""

    def __init__(self: Self, dec: Deconvolve, arrays: T_RF | T_SO2 | RF_SO2) -> None:
        self.dec = dec
        self.dec.progress = False
        self.ts_specifier: T_RF | T_SO2 | RF_SO2 = arrays
        self.cuts: dict[str, xr.Dataset] = {}
        self.ensembles: dict[str, xr.Dataset] = {}

    def dump_reconstructor(
        self: Self,
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
                msg = "Exactly one of temp_so2 or temp_rf must be given."
                raise ValueError(msg)
            case np.ndarray(), None:
                kwargs["response_temp_so2"] = temp_so2
                kwargs["response_temp_rf"] = self.cuts[str(cut)].response.data
            case None, np.ndarray():
                kwargs["response_temp_so2"] = self.cuts[str(cut)].response.data
                kwargs["response_temp_rf"] = temp_rf
            case _:
                raise ValueError
        return Reconstructor(**kwargs)

    @cached_property
    def response(self: Self) -> np.ndarray:
        """The response function in the convolution."""
        out = getattr(
            self.dec,
            f"response_{self.ts_specifier[0]}_{self.ts_specifier[1]}",
        )
        out[self.dec.tau <= 0] = 0
        return out

    @cached_property
    def forcing(self: Self) -> xr.DataArray:
        """The forcing time series in the convolution."""
        return getattr(self.dec, self.ts_specifier[1])

    @cached_property
    def output(self: Self) -> xr.DataArray:
        """The final output time series of the convolution."""
        return getattr(self.dec, self.ts_specifier[0])

    @cached_property
    def control(self: Self) -> xr.DataArray:
        """The control time series in the convolution."""
        return getattr(self.dec, f"{self.ts_specifier[0]}_control")

    def cut_off(self: Self, cutoff: int | Iterable[int]) -> Self:
        """Cut off the response function at a given time lag."""
        match cutoff:
            case int():
                self._single_cut_off(cutoff)
            case Iterable():
                self._check_for_duplicates(cutoff)
                for c in cutoff:
                    if not isinstance(c, int):
                        raise TypeError
                    self._single_cut_off(c)
            case _:
                raise TypeError
        return self

    @staticmethod
    def _check_for_duplicates(cutoff: Iterable[int]) -> None:
        try:
            for _ in zip(cutoff, set(cutoff), strict=True):
                continue
        except ValueError:
            warnings.warn("There are duplicates in the cut-off sequence.", stacklevel=1)

    def _single_cut_off(self: Self, cutoff: int) -> None:
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

    @staticmethod
    def _setup_progress_bar() -> (
        tuple[rich.progress.Progress, rich.progress.Progress, rich.console.Group]
    ):
        p1 = rich.progress.Progress(
            rich.progress.TextColumn("[progress.description]{task.description}"),
            rich.progress.SpinnerColumn(),
            rich.progress.BarColumn(),
            rich.progress.TaskProgressColumn(),
            rich.progress.MofNCompleteColumn(),
            rich.progress.TimeRemainingColumn(elapsed_when_finished=True),
            expand=True,
            redirect_stdout=False,
            # transient=True,
        )
        p2 = rich.progress.Progress(
            rich.progress.TextColumn("[progress.description]{task.description}"),
            rich.progress.SpinnerColumn(),
            rich.progress.BarColumn(),
            rich.progress.TaskProgressColumn(),
            rich.progress.MofNCompleteColumn(),
            rich.progress.TimeRemainingColumn(elapsed_when_finished=True),
            expand=True,
            redirect_stdout=False,
            # transient=True,
        )
        return (
            p1,
            p2,
            rich.console.Group(
                rich.panel.Panel(
                    rich.console.Group(
                        p1,
                        p2,
                    )
                )
            ),
        )

    def generate_ensembles(self: Self, n: int) -> None:
        """Generate an ensemble of response function estimates."""
        if not self.cuts:
            msg = "No cuts have been made."
            raise ValueError(msg)
        cut_panel, ens_panel, group = self._setup_progress_bar()
        cut_panel_id = cut_panel.add_task("", total=len(self.cuts))
        with rich.live.Live(group, transient=True):
            for k, v in self.cuts.items():
                cut_panel.update(
                    cut_panel_id,
                    description=f"Working on cut-off at {k}",
                )
                if k in self.ensembles:
                    continue
                arrays: dict[str, tuple] = {}
                ens_panel_id = ens_panel.add_task("", total=n)
                for i in range(n):
                    ens_panel.update(ens_panel_id, description=f"Creating ensemble {i}")
                    temp_rec = v.temp_rec.copy()
                    temp_random = fppanalysis.signal_rand_phase(self.control.data)
                    temp_rec += temp_random
                    res_rec, err = self.dec.deconvolve(temp_rec, self.forcing.data)
                    r_cut_rec = res_rec.flatten()
                    r_cut_rec[self.dec.tau <= 0] = 0
                    arrays[f"response_{i}"] = (
                        "tau",
                        r_cut_rec,
                        {"label": f"response {i}"},
                    )
                    arrays[f"iters_{i}"] = (
                        "iters",
                        err.flatten(),
                        {"label": f"err {i}"},
                    )
                    ens_panel.advance(ens_panel_id)
                ens_panel.update(ens_panel_id, visible=False)
                self.ensembles[k] = xr.Dataset(
                    arrays,
                    coords={
                        "tau": self.dec.tau,
                        "time": self.output.time,
                        "iters": np.arange(len(err)),
                    },
                )
                cut_panel.advance(cut_panel_id)


class ReconstructOB16:
    """Class that reconstructs the temperature of OB16 from CESM2 simulations."""

    def __init__(
        self: Self, *decs: Deconvolve, base_forcing: Literal["so2", "rf"]
    ) -> None:
        self.ob16 = DeconvolveOB16(data="h0")
        self.ob16.name = "OB16 month"
        self.decs = decs
        self.base = base_forcing

    def plot_temperature(self: Self) -> tuple[mpl.figure.Figure, mpl.figure.Figure]:
        """Plot the reconstructed temperatures.

        Returns
        -------
        tuple[mpl.figure.Figure, mpl.figure.Figure]
            The full and zoomed in figures.
        """
        xlim = (
            vdd.utils.d2n(
                datetime.datetime(1250, 1, 1, 0, 0, tzinfo=datetime.UTC),
            ),
            vdd.utils.d2n(
                datetime.datetime(1350, 1, 1, 0, 0, tzinfo=datetime.UTC),
            ),
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
        console.print(table)
        all_a.legend()
        # all_f.savefig(_SAVE_DIR / "reconstruct_from_all.jpg")
        all_zoom_a.legend()
        # all_zoom_f.savefig(_SAVE_DIR / "reconstruct_from_all_zoom.jpg")
        return all_f, all_zoom_f

    def _plot_temperature_single(
        self: Self,
        dec: Deconvolve,
        res: list[tuple[str, str, str]],
        axs: tuple[mpl.axes.Axes, mpl.axes.Axes],
    ) -> list[tuple[str, str, str]]:
        """Plot the reconstructed temperature for a single simulation."""
        xlim = (
            vdd.utils.d2n(
                datetime.datetime(1250, 1, 1, 0, 0, tzinfo=datetime.UTC),
            ),
            vdd.utils.d2n(
                datetime.datetime(1350, 1, 1, 0, 0, tzinfo=datetime.UTC),
            ),
        )
        # fn = ns(vdd.utils.clean_filename(dec.name))
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
        inv_a.plot(self.ob16.temp.time, new_temp, label="Raw response")
        inv_a.plot(self.ob16.temp.time, new_temp_scaled, label="Scaled response")
        inv_a.legend()
        # inv_f.savefig(_SAVE_DIR / f"reconstruct_from_{fn}.jpg")
        inv_zoom_a.plot(self.ob16.temp.time, self.ob16.temp, label="OB16 temperature")
        inv_zoom_a.plot(self.ob16.temp.time, new_temp, label="Raw response")
        inv_zoom_a.plot(self.ob16.temp.time, new_temp_scaled, label="Scaled response")
        inv_zoom_a.legend()
        # inv_zoom_f.savefig(_SAVE_DIR / f"reconstruct_from_{fn}_zoom.jpg")
        # Print the distance away from the reconstructed
        rob16 = self.ob16.response_temp_so2
        ob16_temp = np.convolve(self.ob16.so2, rob16, mode="same")
        ob16_diff = np.abs(ob16_temp - new_temp).sum()
        ob16_diff_scaled = np.abs(ob16_temp - new_temp_scaled).sum()
        res.append((ns(dec.name), f"{ob16_diff:.2f}", f"{ob16_diff_scaled:.2f}"))
        return res


class TSComparison:
    """Easily compare the statistics of two time series."""

    def __init__(
        self: Self,
        original: xr.DataArray,
        reconstructed: np.ndarray,
        peaks: np.ndarray,
    ) -> None:
        self.orig = original
        self.rec = reconstructed
        if original.data.shape != reconstructed.shape:
            raise EvenLengthError
        self.peaks = peaks

    def _find_peaks(self: Self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Compute the peaks of the time series."""
        _idx = np.argwhere(self.peaks > 0)
        peak_times = self.orig.time.data[_idx].flatten()
        peaks_ts1 = self.orig.data[_idx].flatten()
        peaks_ts2 = self.rec[_idx].flatten()
        return peak_times, peaks_ts1, peaks_ts2

    @cached_property
    def _peaks_tup(self: Self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Compute the peaks of the residual of the SO2."""
        return self._find_peaks()

    @property
    def peaks_orig(self: Self) -> np.ndarray:
        """Peaks of the temperature from SO2.

        Returns
        -------
        np.ndarray
            Peaks of the temperature from SO2.
        """
        return self._peaks_tup[1]

    @property
    def peaks_rec(self: Self) -> np.ndarray:
        """Peaks of the temperature from radiative forcing.

        Returns
        -------
        np.ndarray
            Peaks of the temperature from radiative forcing.
        """
        return self._peaks_tup[2]

    @cached_property
    def residual(self: Self) -> np.ndarray:
        """Compute the residual of the SO2.

        Returns
        -------
        np.ndarray
            Residual of the SO2.
        """
        return self.orig.data - self.rec

    def correlation(self: Self) -> None:
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
            signal,
            sample_frequency,
            nperseg=2**11,
            return_onesided=False,
        )
        frequency_plus = frequency[frequency > 0]
        power_plus = power[frequency > 0]
        return np.asarray(frequency_plus[1:]), np.asarray(power_plus[1:])

    def spectrum(self: Self) -> None:
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
        result = sp.stats.ttest_1samp(basis, popmean=test_value)
        t_statistic, p_value = result.statistic, result.pvalue

        def info(name: str, p_value: float) -> None:
            rprint(
                f"[blue][bold]{name}[/bold]: I can with [/blue][red]"
                f"{(1 - p_value) * 100:.4f}% confidence[/red][blue] say that the "
                f"distribution does not have a mean of {test_value}[/blue]",
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

    def peak_difference_analysis(self: Self) -> None:
        """Plot the difference between the reconstructed and the original peaks."""
        basis = self.peaks_orig.data - self.peaks_rec
        ttest_res = self._peak_difference_ttest(basis)
        pdf, cdf, bin_centers = fppanalysis.distribution(basis, 30, ccdf=False)
        stats = sp.stats.describe(basis)
        fit = sp.stats.norm.fit(basis)
        dist = sp.stats.skewnorm(
            a=stats.skewness,
            loc=stats.mean,
            scale=np.sqrt(stats.variance),
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
        norm_so2 = getattr(sp.stats.norm, dist)(dist_data[0], *norm_fit)
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
            dist_data[0],
            norm_so2,
            c=colors[0],
            label="_Norm SO2",
            alpha=0.5,
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
            [norm, skewnorm],
            ["Norm", "Skewnorm"],
            loc=norm_loc,
            framealpha=0.5,
        )
        norm_legend.legend_handles[0].set_color("black")  # type: ignore[union-attr]
        norm_legend.legend_handles[1].set_color("black")  # type: ignore[union-attr]
        ax.add_artist(bar_legend)
        ax.add_artist(norm_legend)
        # Make the plot symmetric around 0
        xlim = np.abs(plt.gca().get_xlim()).max()
        plt.xlim((-xlim, xlim))
        plt.ylabel(dist.upper())
        plt.xlabel("Difference between the peaks")
        # plt.savefig(_SAVE_DIR / f"{self.sim_name}-peak-difference-{dist}.png")

    def plot_reconstructions(self: Self, fig: mpl.figure.Figure | None = None) -> None:
        """Plot the reconstruction of the data."""
        if fig is None:
            fig = plt.figure()
        ax = fig.gca()
        ax.set_xlabel("Time after first eruption [yr]")
        ax.set_ylabel("Absolute")
        time_ = self.orig.time
        ax.plot(time_, self.orig.data, label="Original")
        ax.plot(time_, self.rec, label="Reconstriction")
        # ax.set_xlim((-790 * 365, -650 * 365))
        ax.legend(framealpha=0.5)
        # plt.savefig(_SAVE_DIR / f"{self.sim_name}-temp-reconstructed.png")
