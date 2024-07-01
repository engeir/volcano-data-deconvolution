"""Utility functions for the volcano-data-deconvolution package."""

import datetime
import pathlib
import re
from typing import Literal, Never, NoReturn, Self, overload

import cftime
import cosmoplots
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import volcano_base
import xarray as xr
from numpy.typing import NDArray

plt.style.use(
    [
        "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle",
        "vdd.jgr",
        "vdd.extra",
    ],
)


class LiteralError(ValueError):
    """Unknown literal value."""

    def __init__(self: Self, message: str) -> None:
        msg = f"Incorrect literal provided. Got: {message}"
        super().__init__(msg)


class FigureGrid:
    """Return a figure with axes which is appropriate for (rows, columns) subfigures.

    Parameters
    ----------
    rows : int
        The number of rows in the figure
    columns : int
        The number of columns in the figure
    labels : list[str] | None
        The labels to be applied to each subfigure. Defaults to (a), (b), (c), ...
    """

    def __init__(
        self: Self,
        rows: int,
        columns: int,
        labels: list[str] | None = None,
    ) -> None:
        self.rows = rows
        self.columns = columns
        self.labels = labels
        self._pos: tuple[float, float] = (-0.2, 0.95)
        self._share_axes: Literal["x", "y", "both"] | None = None
        self._columns_first: bool = False

    def _calculate_figsize(self: Self) -> tuple[float, float]:
        """Calculate the figure size based on the number of rows and columns."""
        full_cols = 3.37 * self.columns
        squash_cols = 3.37 * self.columns - (self.columns - 1) * 3.37 * 0.25
        full_rows = 2.08277 * self.rows
        squash_rows = 2.08277 * self.rows - (self.rows - 1) * 2.08277 * 0.25
        match self._share_axes:
            case None:
                return full_rows, full_cols
            case "x":
                return squash_rows, full_cols
            case "y":
                return full_rows, squash_cols
            case "both":
                return squash_rows, squash_cols
            case _:
                raise LiteralError(self._share_axes)

    def _update_labels(self: Self) -> list[str]:
        labels = (
            [
                rf"$\mathrm{{({chr(97 + ell)})}}$"
                for ell in range(self.rows * self.columns)
            ]
            if not self.labels or len(self.labels) != int(self.rows * self.columns)
            else self.labels.copy()
        )
        if self._columns_first:
            labels = [
                labels[i * self.rows + j]
                for j in range(self.rows)
                for i in range(self.columns)
            ]
        return labels

    def using(
        self: Self,
        *,
        pos: tuple[float, float] | None = None,
        share_axes: Literal["x", "y", "both"] | None = None,
        columns_first: bool | None = None,
    ) -> Self:
        """Set text properties.

        The properties must be given as keyword arguments to take effect.

        Parameters
        ----------
        pos : tuple[float, float] | None
            The position in the subfigure relative to `gravity`. Default is `(10.0,
            10.0)`.
        share_axes : Literal["x", "y", "both"] | None
            Use a shared axis for the given direction. Default is `None`.
        columns_first : bool | None
            If the labels should be placed in a columns-first order. Default is `False`,
            meaning rows are numbered first.

        Returns
        -------
        Self
            The object itself.
        """
        self._pos = pos or self._pos
        self._share_axes = share_axes or self._share_axes
        self._columns_first = columns_first or self._columns_first
        return self

    def get_grid(
        self: Self,
        **kwargs: dict,
    ) -> tuple[mpl.figure.Figure, list[mpl.axes.Axes]]:
        """Return a figure with axes appropriate for (rows, columns) subfigures.

        Parameters
        ----------
        **kwargs : dict
            Additional keyword arguments to be passed to Axes.text.

        Returns
        -------
        mpl.figure.Figure
            The figure object
        list[mpl.axes.Axes]
            A list of all the axes objects owned by the figure
        """
        full_height, full_width = self._calculate_figsize()
        fig = plt.figure(figsize=(full_width, full_height))
        axes = []
        labels = self._update_labels()
        for r in range(self.rows):
            if self._share_axes in {"x", "both"}:
                rel_height = 0.75 + 0.25 / self.rows
                height = 0.75 / self.rows / rel_height
                bottom_pad = 0.2 / self.rows / rel_height
                bottom = bottom_pad + height * (self.rows - 1 - r)
            else:
                bottom_pad = 0.2 / self.rows
                height = 0.75 / self.rows
                bottom = bottom_pad + (self.rows - 1 - r) / self.rows
            for c in range(self.columns):
                if self._share_axes in {"y", "both"}:
                    rel_width = 0.75 + 0.25 / self.columns
                    width = 0.75 / self.columns / rel_width
                    left_pad = 0.2 / self.columns / rel_width
                    left = left_pad + width * c
                else:
                    left_pad = 0.2 / self.columns
                    width = 0.75 / self.columns
                    left = left_pad + c / self.columns
                axes.append(fig.add_axes((left, bottom, width, height)))
                if self._share_axes in {"x", "both"} and r != self.rows - 1:
                    axes[-1].set_xticklabels([])
                if self._share_axes in {"y", "both"} and c != 0:
                    axes[-1].set_yticklabels([])
                axes[-1].text(
                    self._pos[0],
                    self._pos[1],
                    labels[self.columns * r + c],
                    transform=axes[-1].transAxes,
                    **kwargs,
                )
        return fig, axes


def figure_multiple_rows_columns(
    rows: int,
    columns: int,
    share_axes: Literal["x", "y", "both"],
    *,
    columns_first: bool = False,
) -> tuple[mpl.figure.Figure, list[mpl.axes.Axes]]:
    """Return a figure with axes which is appropriate for (rows, columns) subfigures.

    Parameters
    ----------
    rows : int
        The number of rows in the figure
    columns : int
        The number of columns in the figure
    share_axes : Literal["x", "y", "both"]
        Share the axes in the figure. Defaults to not sharing.
    columns_first : bool
        If the labels should be placed in a columns-first order. Default is `False`.

    Returns
    -------
    mpl.figure.Figure
        The figure object
    list[mpl.axes.Axes]
        A list of all the axes objects owned by the figure

    Note
    ----
    This functions is deprecated, and works only as a wrapper for the `FigureGrid`
    class.
    """
    return (
        FigureGrid(rows, columns)
        .using(share_axes=share_axes, columns_first=columns_first)
        .get_grid()
    )


def combine(*files: str | pathlib.Path) -> cosmoplots.Combine:
    """Give all files that should be combined.

    Parameters
    ----------
    *files : str | pathlib.Path
        A file path that can be read by pathlib.Path.

    Returns
    -------
    cosmoplots.Combine
        An instance of the Combine class.

    Examples
    --------
    Load the files and subsequently call the methods that updates the properties.

    >>> combine(
    ...     "file1.png", "file2.png", "file3.png", "file4.png", "file5.png", "file6.png"
    ... ).using(fontsize=120).in_grid(w=2, h=3).with_labels(
    ...     "(a)", "(b)", "(c)", "(d)", "(e)", "(f)"
    ... ).save()

    All (global) methods except from `save` and `help` return the object itself, and can
    be chained together.
    """
    return cosmoplots.combine(*files).using(fontsize=8)


def s2n(num: float, decimal: int = 2) -> str:
    """Convert a number to scientific notation."""
    return f"\\num{{{num:.{decimal}e}}}"


def d2n(date: datetime.datetime) -> float:
    """Convert a datetime to a number, using 2000-01-01 as the reference date."""
    unit = "days since 2000-01-01"
    return cftime.date2num(date, units=unit, calendar="noleap", has_year_zero=True)


def name_translator(name: re.Match) -> str:
    """Translate the first match group from the old naming convention to the new one."""
    match name[0]:
        case "medium":
            return "small"
        case "medium-plus":
            return "intermediate"
        case "size5000":
            return "extreme"
        case "strong":
            return name[0]
        case "tt-2sep":
            return "int-2sep"
        case "double-overlap" | "tt-4sep":
            return "int-4sep"
        case _:
            msg = f"Unknown name: {name}"
            raise ValueError(msg)


@overload
def name_swap(name: pathlib.Path) -> pathlib.Path: ...


@overload
def name_swap(name: str) -> str: ...


def name_swap(name: str | pathlib.Path) -> str | pathlib.Path:
    """Replace any occurrence of the old name with the new name."""
    regex = r"strong|medium-plus|medium|size5000|tt-2sep|tt-4sep|double-overlap"
    match name:
        case str():
            return re.sub(
                regex,
                name_translator,
                name,
            )
        case pathlib.Path():
            return pathlib.Path(
                re.sub(
                    regex,
                    name_translator,
                    str(name),
                ),
            )


def never_called(value: Never) -> NoReturn:
    """Raise an error if a value is passed to a function that should never be called."""
    # The function is useful when running mypy. If, in a series of if/elif or
    # match/case, a variable is not fully handled, mypy will complain and say that the
    # variable is of the wrong type when this function is called in the final `else`
    # clause.
    msg = "Code is unreachable."
    raise AssertionError(msg, value)


def clean_filename(filename: str) -> pathlib.Path:
    """Replace non-alphanumeric characters with a hyphen to create a clean filename."""
    # Replace all non-alphanumeric characters, whitespace, and certain special
    # characters with "-"
    cleaned_filename = re.sub(r"[^\w\s.-]", "-", filename)
    # Replace multiple whitespace characters with a single "-"
    cleaned_filename = re.sub(r"\s+", "-", cleaned_filename)
    # Remove consecutive hyphens
    cleaned_filename = re.sub(r"-+", "-", cleaned_filename)
    return pathlib.Path(cleaned_filename.lower())


@overload
def normalise(arr: NDArray[np.float64]) -> NDArray[np.float64]: ...


@overload
def normalise(arr: xr.DataArray) -> xr.DataArray: ...


def normalise(
    arr: xr.DataArray | NDArray[np.float64],
) -> xr.DataArray | NDArray[np.float64]:
    """Normalise the array to the maximum or minimum value."""
    if abs(np.nanmax(arr)) > abs(np.nanmin(arr)):
        return arr / np.nanmax(arr)
    return arr / np.nanmin(arr)


def weighted_monthly_avg(da: xr.DataArray) -> xr.DataArray:
    """Calculate a temporal mean, weighted by days in each month.

    Parameters
    ----------
    da : xr.DataArray
        Input data structure to do temporal average on

    Returns
    -------
    xr.DataArray
        The new data structure with time averaged data

    Notes
    -----
    From
    https://ncar.github.io/esds/posts/2021/yearly-averages-xarray/#wrap-it-up-into-a-function
    """
    # Determine the month length
    month_length = da.time.dt.days_in_month
    # Calculate the weights
    wgts = month_length.groupby("time.month") / month_length.groupby("time.month").sum()
    # Make sure the weights in each year add up to 1
    np.testing.assert_allclose(wgts.groupby("time.month").sum(xr.ALL_DIMS), np.ones(12))
    # Setup our masking for nan values
    cond = da.isna()
    ones = xr.where(cond, 0.0, 1.0)
    # Calculate the numerator
    obs_sum = (da * wgts).resample(time="MS").sum(dim="time")
    # Calculate the denominator
    ones_out = (ones * wgts).resample(time="MS").sum(dim="time")
    # Return the weighted average
    return obs_sum / ones_out


def pad_before_convolution(
    arr: xr.DataArray,
    diff_before_zero: float = 0.08493151,
    pad_values: np.ndarray | None = None,
) -> xr.DataArray:
    """Pad the array on the left with zeros as preparation for a convolution.

    This uses the ``xr.DataArray().pad`` method under the hood, and tries to fix the
    time axis by inserting valid time values. The default behaviour of
    ``xr.DataArray().pad`` is to set ``np.nan`` for the padding values, which will hide
    the padded values in a plot.

    Parameters
    ----------
    arr : xr.DataArray
        The input array to pad.
    diff_before_zero : float, optional
        The time difference between time zero and the time step prior. Since time data
        may be monthly with varying number of days, a constant time step is not always
        wanted. The default is 0.08493151, used for monthly data.
    pad_values : np.ndarray | None
        The values to pad the array with. Defaults to zeros.

    Returns
    -------
    xr.DataArray
        The padded array.

    Raises
    ------
    ValueError
        If the new time does not start with a zero or the length is not divisible by 2.
    IndexError
        If the padding length is not equal to the difference in length between the new
        and original array.
    """
    n = len(arr.time) - 1
    arr_ = arr.pad(time=(n, 0), mode="constant", constant_values=0)
    if pad_values is None:
        pad_values_ = np.zeros(len(arr_) - len(arr))
    elif len(pad_values) == len(arr_) - len(arr):
        pad_values_ = pad_values
    else:
        raise IndexError
    arr_.data[: len(pad_values_)] = pad_values_
    match arr.time.data[0]:
        case float():
            time_ = arr.time.data - arr.time.data[0]
        case _:
            time_ = volcano_base.manipulate.dt2float(arr.time.data)
            time_ -= time_[0]
    if time_[0] != 0 or len(arr_) % 2 == 0:
        raise ValueError
    # Assign new time values to the padded array that is twice the length of the
    # original, and that goes from -n*dt to n*dt.
    return arr_.assign_coords(
        time=np.concatenate((time_[1:] - time_[-1] - diff_before_zero, time_)),
    )
