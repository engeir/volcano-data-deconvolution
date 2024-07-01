import matplotlib as mpl
import numpy as np
import volcano_base
import xarray as xr
from _typeshed import Incomplete as Incomplete
from abc import ABC, abstractmethod
from collections.abc import Callable as Callable, Iterable
from enum import Enum
from functools import cached_property as cached_property
from numpy.typing import NDArray as NDArray
from pydantic import BaseModel
from typing import Literal, Self

T_RF: Incomplete
T_SO2: Incomplete
RF_SO2: Incomplete
T_Strengths: Incomplete
console: Incomplete

class DataAttributeNotFoundError(FileNotFoundError):
    def __init__(self) -> None: ...

class PaddingMethod(Enum):
    NO: bool
    ZEROS: bool
    NOISE: str

class Normalise(Enum):
    AMPLITUDE: str
    MEAN_STD: str
    NO: bool

class Reconstructor(BaseModel):
    class Config:
        validate_assignment: bool
        strict: bool
        arbitrary_types_allowed: bool
    tau: np.ndarray
    time: xr.DataArray
    response_temp_so2: np.ndarray
    response_temp_rf: np.ndarray
    forcing: xr.DataArray
    signal: xr.DataArray
    name: str
    normalise: Normalise

class CESMData(BaseModel):
    strength: T_Strengths
    dims: list[str]
    class Config:
        validate_assignment: bool
        frozen: bool
        extra: str
        strict: bool
    @cached_property
    def so2(self) -> xr.DataArray: ...
    @cached_property
    def tmso2(self) -> xr.DataArray: ...
    @cached_property
    def tmso2_control(self) -> xr.DataArray: ...
    @cached_property
    def aod(self) -> xr.DataArray: ...
    @cached_property
    def aod_control(self) -> xr.DataArray: ...
    @cached_property
    def rf(self) -> xr.DataArray: ...
    @cached_property
    def rf_control(self) -> xr.DataArray: ...
    @cached_property
    def temp(self) -> xr.DataArray: ...
    @cached_property
    def temperature_control(self) -> xr.DataArray: ...
    def initialise_data(self) -> None: ...

class _PostInitCaller(ABC, type):
    def __call__(cls, *args, **kwargs) -> T: ...

class EvenLengthError(Exception):
    message: Incomplete
    def __init__(self, message: str = 'The arrays must have an odd length.') -> None: ...

class Deconvolve(metaclass=_PostInitCaller):
    name: str
    normalise: Incomplete
    def __init__(self, normalise: Normalise = ...) -> None: ...
    def __post_init__(self) -> None: ...
    @property
    @abstractmethod
    def tau(self) -> NDArray[np.float64]: ...
    @property
    @abstractmethod
    def so2(self) -> xr.DataArray: ...
    @property
    @abstractmethod
    def so2_decay(self) -> xr.DataArray: ...
    @property
    @abstractmethod
    def rf(self) -> xr.DataArray: ...
    @property
    @abstractmethod
    def temp(self) -> xr.DataArray: ...
    @property
    @abstractmethod
    def temp_control(self) -> xr.DataArray: ...
    def deconvolve(self, signal: NDArray[np.float64], forcing: NDArray[np.float64]) -> tuple[NDArray[np.float64], NDArray[np.float64]]: ...
    def dump_reconstructor(self) -> Reconstructor: ...
    def change_deconvolution_method(self, method: Callable[[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]], method_args: list | None = None, method_kwargs: dict | None = None) -> None: ...
    @property
    def response_rf_so2(self) -> np.ndarray: ...
    @property
    def response_temp_so2(self) -> np.ndarray: ...
    @property
    def response_rf_so2_decay(self) -> np.ndarray: ...
    @property
    def response_temp_so2_decay(self) -> np.ndarray: ...
    @property
    def response_temp_rf(self) -> np.ndarray: ...
    def plot_dec_rf_with_so2(self) -> None: ...
    def plot_dec_temp_with_so2(self) -> None: ...
    def plot_dec_temp_with_rf(self) -> None: ...

class DeconvolveCESM(Deconvolve):
    pad_before: Incomplete
    name: Incomplete
    def __init__(self, normalise: Normalise = ..., pad_before: PaddingMethod = ..., cesm: CESMData | None = None) -> None: ...
    @cached_property
    def so2(self) -> xr.DataArray: ...
    @cached_property
    def so2_decay(self) -> xr.DataArray: ...
    @cached_property
    def tau(self) -> np.ndarray: ...
    @cached_property
    def tmso2(self) -> xr.DataArray: ...
    @cached_property
    def tmso2_control(self) -> xr.DataArray: ...
    @cached_property
    def aod(self) -> xr.DataArray: ...
    @cached_property
    def aod_control(self) -> xr.DataArray: ...
    @cached_property
    def rf(self) -> xr.DataArray: ...
    @cached_property
    def rf_control(self) -> xr.DataArray: ...
    @cached_property
    def temp(self) -> xr.DataArray: ...
    @cached_property
    def temp_control(self) -> xr.DataArray: ...
    @property
    def response_aod_so2(self) -> np.ndarray: ...
    @property
    def response_rf_aod(self) -> np.ndarray: ...
    @property
    def response_temp_aod(self) -> np.ndarray: ...

class DeconvolveOB16(Deconvolve):
    data: Incomplete
    start_pt: int
    end_pt: Incomplete
    def __init__(self, data: volcano_base.load.OttoBliesner | Literal['h0', 'h1'] = 'h1', normalise: Normalise = ..., length: int | None = None) -> None: ...
    @cached_property
    def so2(self) -> xr.DataArray: ...
    @cached_property
    def so2_decay(self) -> xr.DataArray: ...
    @cached_property
    def tau(self) -> np.ndarray: ...
    @cached_property
    def rf(self) -> xr.DataArray: ...
    @cached_property
    def rf_control(self) -> xr.DataArray: ...
    @cached_property
    def temp(self) -> xr.DataArray: ...
    @cached_property
    def temp_control(self) -> xr.DataArray: ...

class CutOff:
    dec: Incomplete
    ts_specifier: Incomplete
    cuts: Incomplete
    ensembles: Incomplete
    def __init__(self, dec: Deconvolve, arrays: T_RF | T_SO2 | RF_SO2) -> None: ...
    def dump_reconstructor(self, cut: int, temp_so2: np.ndarray | None = None, temp_rf: np.ndarray | None = None) -> Reconstructor: ...
    @cached_property
    def response(self) -> np.ndarray: ...
    @cached_property
    def forcing(self) -> xr.DataArray: ...
    @cached_property
    def output(self) -> xr.DataArray: ...
    @cached_property
    def control(self) -> xr.DataArray: ...
    def cut_off(self, cutoff: int | Iterable[int]) -> Self: ...
    def generate_ensembles(self, n: int) -> None: ...

class ReconstructOB16:
    ob16: Incomplete
    decs: Incomplete
    base: Incomplete
    def __init__(self, *decs: Deconvolve, base_forcing: Literal['so2', 'rf']) -> None: ...
    def plot_temperature(self) -> tuple[mpl.figure.Figure, mpl.figure.Figure]: ...

class TSComparison:
    orig: Incomplete
    rec: Incomplete
    peaks: Incomplete
    def __init__(self, original: xr.DataArray, reconstructed: np.ndarray, peaks: np.ndarray) -> None: ...
    @property
    def peaks_orig(self) -> np.ndarray: ...
    @property
    def peaks_rec(self) -> np.ndarray: ...
    @cached_property
    def residual(self) -> np.ndarray: ...
    def correlation(self) -> None: ...
    def spectrum(self) -> None: ...
    def peak_difference_analysis(self) -> None: ...
    def plot_reconstructions(self, fig: mpl.figure.Figure | None = None) -> None: ...