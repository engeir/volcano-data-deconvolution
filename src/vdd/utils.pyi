import cosmoplots
import datetime
import matplotlib as mpl
import numpy as np
import pathlib
import re
import xarray as xr
from _typeshed import Incomplete as Incomplete
from numpy.typing import NDArray as NDArray
from typing import Literal, Never, NoReturn, Self, overload

class LiteralError(ValueError):
    def __init__(self, message: str) -> None: ...

class FigureGrid:
    rows: Incomplete
    columns: Incomplete
    labels: Incomplete
    def __init__(self, rows: int, columns: int, labels: list[str] | None = None) -> None: ...
    def using(self, *, pos: tuple[float, float] | None = None, share_axes: Literal['x', 'y', 'both'] | None = None, columns_first: bool | None = None) -> Self: ...
    def get_grid(self, **kwargs: dict) -> tuple[mpl.figure.Figure, list[mpl.axes.Axes]]: ...

def figure_multiple_rows_columns(rows: int, columns: int, share_axes: Literal['x', 'y', 'both'] | None = None, *, columns_first: bool = False) -> tuple[mpl.figure.Figure, list[mpl.axes.Axes]]: ...
def combine(*files: str | pathlib.Path) -> cosmoplots.Combine: ...
def n2sci(num: float, decimal: int = 2) -> str: ...
def d2n(date: datetime.datetime | str) -> float: ...
def name_translator(name: re.Match) -> str: ...
@overload
def name_swap(name: pathlib.Path) -> pathlib.Path: ...
@overload
def name_swap(name: str) -> str: ...
def never_called(value: Never) -> NoReturn: ...
def clean_filename(filename: str) -> pathlib.Path: ...
@overload
def normalise(arr: NDArray[np.float64]) -> NDArray[np.float64]: ...
@overload
def normalise(arr: xr.DataArray) -> xr.DataArray: ...
def weighted_monthly_avg(da: xr.DataArray) -> xr.DataArray: ...
def pad_before_convolution(arr: xr.DataArray, diff_before_zero: float = 0.08493151, pad_values: np.ndarray | None = None) -> xr.DataArray: ...
