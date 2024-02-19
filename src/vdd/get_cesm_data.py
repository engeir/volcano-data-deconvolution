"""Functions that fetches data from files and returns it as an xarray.DataArray."""

import matplotlib.pyplot as plt
import volcano_base
import xarray as xr


def convert_time_start_zero(arr: xr.DataArray) -> xr.DataArray:
    """Convert the time to start from zero."""
    arr = arr.assign_coords(time=volcano_base.manipulate.dt2float(arr.time.data) - 1850)
    arr.coords["time"].attrs["long_name"] = "Time-after-eruption"
    arr.coords["time"].attrs["units"] = "yr"
    return arr


def get_aod_cesm(plot_example: bool = False) -> xr.DataArray:
    """Get the CESM2 AOD data."""
    data = (
        volcano_base.load.FindFiles()
        .find("e_fSST1850", "AODVISstdn", "h0", "strong")
        .remove("ens1")
        .keep_most_recent()
    )
    files = data.load()
    files = volcano_base.manipulate.mean_flatten(files, dims=["lat", "lon"])
    files = volcano_base.manipulate.shift_arrays(files, daily=False)
    files = volcano_base.manipulate.shift_arrays(files, custom=1)
    files = volcano_base.manipulate.subtract_mean_of_tail(files)
    files = volcano_base.manipulate.data_array_operation(files, convert_time_start_zero)
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


def get_rf_cesm(plot_example: bool = False) -> xr.DataArray:
    """Get the CESM2 RF data."""
    control: volcano_base.load.FindFiles = (
        volcano_base.load.FindFiles()
        .find("e_fSST1850", "h0", "control", {"FLNT", "FSNT"}, "ens1")
        .sort("ensemble", "attr")
    )
    data: volcano_base.load.FindFiles = (
        volcano_base.load.FindFiles()
        .find("e_fSST1850", {"FLNT", "FSNT"}, "h0", "strong")
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
    rf = volcano_base.manipulate.shift_arrays(rf, daily=False)
    rf = volcano_base.manipulate.shift_arrays(rf, custom=1)
    rf = volcano_base.manipulate.subtract_mean_of_tail(rf)
    rf = volcano_base.manipulate.data_array_operation(rf, convert_time_start_zero)
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


def get_trefht_cesm(plot_example: bool = False) -> xr.DataArray:
    """Get the CESM2 temperature data."""
    data = (
        volcano_base.load.FindFiles()
        .find("e_BWma1850", "TREFHT", "h0", "strong")
        .remove("ens1")
        .keep_most_recent()
    )
    files = data.load()
    files = volcano_base.manipulate.mean_flatten(files, dims=["lat", "lon"])
    files = volcano_base.manipulate.shift_arrays(files, daily=False)
    files = volcano_base.manipulate.shift_arrays(files, custom=1)

    def subtract_control_mean(arr: xr.DataArray) -> xr.DataArray:
        return arr - volcano_base.config.MEANS["TREFHT"]

    files = volcano_base.manipulate.data_array_operation(files, subtract_control_mean)
    files = volcano_base.manipulate.data_array_operation(files, convert_time_start_zero)
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
