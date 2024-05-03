"""Utility functions for the volcano-data-deconvolution package."""

import datetime
import pathlib
import re
from typing import Never, NoReturn, overload

import cftime
import cosmoplots
import matplotlib.pyplot as plt
import numpy as np
import volcano_base
import xarray as xr

plt.style.use([
    "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle",
    "vdd.jgr",
    "vdd.extra",
])


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


def d2n(date: datetime.datetime) -> float:
    """Convert a datetime to a number, using 2000-01-01 as the reference date."""
    unit = "days since 2000-01-01"
    return cftime.date2num(date, units=unit, calendar="noleap", has_year_zero=True)


def name_translator(name: re.Match) -> str:
    """Translate the first match group from the old naming convention to the new one."""
    match name.group(0):
        case "medium":
            return "small"
        case "medium-plus":
            return "intermediate"
        case "size5000":
            return "extreme"
        case "tt-2sep" | "tt-4sep" | "strong":
            return name.group(0)
        case "double-overlap":
            return "tt-4sep"
        case _:
            raise ValueError(f"Unknown name: {name}")


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
                )
            )


def never_called(value: Never) -> NoReturn:
    """Raise an error if a value is passed to a function that should never be called."""
    # The function is useful when running mypy. If, in a series of if/elif or
    # match/case, a variable is not fully handled, mypy will complain and say that the
    # variable is of the wrong type when this function is called in the final `else`
    # clause.
    raise AssertionError("Code is unreachable.")


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
def normalise(arr: np.ndarray) -> np.ndarray: ...


@overload
def normalise(arr: xr.DataArray) -> xr.DataArray: ...


def normalise(arr: xr.DataArray | np.ndarray) -> xr.DataArray | np.ndarray:
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
    cond = da.isnull()
    ones = xr.where(cond, 0.0, 1.0)
    # Calculate the numerator
    obs_sum = (da * wgts).resample(time="MS").sum(dim="time")
    # Calculate the denominator
    ones_out = (ones * wgts).resample(time="MS").sum(dim="time")
    # Return the weighted average
    return obs_sum / ones_out


def pad_before_convolution(
    arr: xr.DataArray, diff_before_zero: float = 0.08493151
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

    Returns
    -------
    xr.DataArray
        The padded array.
    """
    n = len(arr.time) - 1
    arr_ = arr.pad(time=(n, 0), mode="constant", constant_values=0)
    match arr.time.data[0]:
        case float():
            time_ = arr.time.data - arr.time.data[0]
        case _:
            time_ = volcano_base.manipulate.dt2float(arr.time.data)
            time_ -= time_[0]
    assert time_[0] == 0
    assert len(arr_) % 2 != 0
    # Assign new time values to the padded array that is twice the length of the
    # original, and that goes from -n*dt to n*dt.
    arr_ = arr_.assign_coords(
        time=np.concatenate((time_[1:] - time_[-1] - diff_before_zero, time_)),
    )
    return arr_
