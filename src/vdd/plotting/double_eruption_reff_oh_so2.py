"""Create plots of some extra key eruption parameters."""

from typing import ClassVar, NamedTuple

import matplotlib.pyplot as plt
import volcano_base
import volcano_base.manipulate as vbm
import xarray as xr
from volcano_base.load import FindFiles

import vdd.reff

_SAVE_DIR = volcano_base.config.SAVE_PATH / "extra-attrs"
if not _SAVE_DIR.exists():
    _SAVE_DIR.mkdir(parents=False)

plt.style.use(
    [
        "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle",
        "vdd.extra",
        "vdd.jgr",
    ],
)
COLORS = plt.rcParams["axes.prop_cycle"].by_key()["color"]


class PlotArgs(NamedTuple):
    """Definition of plotting parameters."""

    ls: str
    c: str
    lab: str


class ReffPlot:
    """Plot the aerosol effective radius."""

    _SHOW = True
    FINDER = (
        FindFiles()
        .find({"ens1", "ens3"}, {"tt-2sep", "tt-4sep", "medium-2sep", "medium-4sep"})
        .sort("ensemble", "sim")
        .copy
    )
    sad_ = FINDER().keep("SAD_AERO")
    reff_ = FINDER().keep("REFF_AERO")
    temp_ = FINDER().keep("T")
    _PLOT_DICT: ClassVar[dict[str, PlotArgs]] = {
        "S26, 2sep": PlotArgs(ls="-", c=COLORS[0], lab="S26-2SEP"),
        "S26, 4sep": PlotArgs(ls="-", c=COLORS[1], lab="S26-4SEP"),
        "S400, 2sep": PlotArgs(ls="-", c=COLORS[2], lab="S400-2SEP"),
        "S400, 4sep": PlotArgs(ls="-", c=COLORS[3], lab="S400-4SEP"),
    }

    def print(self) -> None:
        """Print all data that is being used."""
        print(self.sad_)
        print(self.reff_)
        print(self.temp_)

    def ens2median(self, arr: list[xr.DataArray]) -> xr.DataArray:
        """Combine a list of arrays from different (known) ensembles into a median."""
        arr = vbm.shift_arrays(arr, daily=False)
        arr = vbm.shift_arrays(arr, custom=1, daily=False)
        arr_ = vbm.get_median(arr, xarray=True)
        arr_ = arr_[: int(12 * 16)]
        return arr_.assign_coords(time=vbm.dt2float(arr_.time.data) - 1850)

    def compute(self) -> xr.Dataset:
        """Compute the effective radius for all simulations."""
        sad = self.sad_.load()
        reff = self.reff_.load()
        temp = self.temp_.load()
        reff_s21 = vdd.reff.Reff(reff[0], temp[0], sad[0]).calculate_reff()
        reff_s41 = vdd.reff.Reff(reff[1], temp[1], sad[1]).calculate_reff()
        reff_m21 = vdd.reff.Reff(reff[2], temp[2], sad[2]).calculate_reff()
        reff_m41 = vdd.reff.Reff(reff[3], temp[3], sad[3]).calculate_reff()
        reff_s23 = vdd.reff.Reff(reff[4], temp[4], sad[4]).calculate_reff()
        reff_s43 = vdd.reff.Reff(reff[5], temp[5], sad[5]).calculate_reff()
        reff_m23 = vdd.reff.Reff(reff[6], temp[6], sad[6]).calculate_reff()
        reff_m43 = vdd.reff.Reff(reff[7], temp[7], sad[7]).calculate_reff()
        e2m = self.ens2median
        names = list(self._PLOT_DICT.keys())
        reff_s2 = e2m([reff_s21, reff_s23]).rename(names[0])
        reff_s4 = e2m([reff_s41, reff_s43]).rename(names[1])
        reff_m2 = e2m([reff_m21, reff_m23]).rename(names[2])
        reff_m4 = e2m([reff_m41, reff_m43]).rename(names[3])
        return xr.merge([reff_s2, reff_s4, reff_m2, reff_m4])

    def load(self) -> xr.Dataset:
        """Get or generate the data as a xarray data set."""
        file = _SAVE_DIR / "reff.nc"
        if file.exists():
            return xr.load_dataset(file)
        self._SHOW = False
        ds = self.compute().assign_attrs(
            {"description": "Global stratospheric mean aerosol effective radius."}
        )
        ds.to_netcdf(file)
        return ds

    def plot(self) -> None:
        """Plot the computed arrays."""
        plt.figure()
        plt.semilogy()
        data = self.load()
        for name, r in data.data_vars.items():
            plot_args = self._PLOT_DICT[name]
            r.plot(label=plot_args.lab, c=plot_args.c, ls=plot_args.ls)  # type: ignore[call-arg]
        plt.xlabel("Time after first eruption [yr]")
        plt.ylabel(r"$R_{\text{EFF}}$ $[\mathrm{\mu m}]$")
        plt.legend()
        plt.savefig(_SAVE_DIR / "reff")
        if self._SHOW:
            plt.show()
        plt.close("all")


class OHPlot:
    """Create plots of the OH CESM output field.

    Attributes
    ----------
    oh_c : FindFiles
        Object holding the file keys to the control simulation.
    oh_m : FindFiles
        Object holding the file keys to the smallest eruption simulation.
    oh_p : FindFiles
        Object holding the file keys to the intermediate eruption simulation.
    oh_s : FindFiles
        Object holding the file keys to the large eruption simulation.
    oh_e : FindFiles
        Object holding the file keys to the extreme eruption simulation.
    oh_m2 : FindFiles
        Object holding the file keys to the smallest 2-year double eruption simulation.
    oh_m4 : FindFiles
        Object holding the file keys to the smallest 4-year double eruption simulation.
    oh_p2 : FindFiles
        Object holding the file keys to the intermediate 2-year double eruption simulation.
    oh_p4 : FindFiles
        Object holding the file keys to the intermediate 4-year double eruption simulation.
    """

    _SHOW = True
    _FINDER: FindFiles = (
        FindFiles().find("OH", "e_fSST1850").sort("sim", "ensemble").copy
    )
    oh_c: FindFiles = _FINDER().keep("control")
    # Only ens5 start in 1850 in the following three experiments. The rest were saved
    # from 1859 onwards.
    oh_m: FindFiles = _FINDER().keep("medium", "ens5")
    oh_p: FindFiles = _FINDER().keep("medium-plus", "ens5")
    oh_s: FindFiles = _FINDER().keep("strong", "ens5")
    oh_e: FindFiles = _FINDER().keep("size5000")
    oh_m2: FindFiles = _FINDER().keep("medium-2sep")
    oh_m4: FindFiles = _FINDER().keep("medium-4sep")
    oh_p2: FindFiles = _FINDER().keep("tt-2sep", {"ens1", "ens3"})
    oh_p4: FindFiles = _FINDER().keep("tt-4sep", {"ens1", "ens3"})
    _PLOT_DICT: ClassVar[dict[str, PlotArgs]] = {
        "CONTROL": PlotArgs(ls="-", c=COLORS[0], lab="CONTROL"),
        "S26": PlotArgs(ls="-", c=COLORS[1], lab="S26"),
        "S400": PlotArgs(ls="-", c=COLORS[2], lab="S400"),
        "S1629": PlotArgs(ls="-", c=COLORS[3], lab="S1629"),
        "S3000": PlotArgs(ls="-", c=COLORS[4], lab="S300"),
        "S26, 2sep": PlotArgs(ls=":", c=COLORS[1], lab="_S26-2SEP"),
        "S26, 4sep": PlotArgs(ls="--", c=COLORS[1], lab="_S26-4SEP"),
        "S400, 2sep": PlotArgs(ls=":", c=COLORS[2], lab="_S400-2SEP"),
        "S400, 4sep": PlotArgs(ls="--", c=COLORS[2], lab="_S400-4SEP"),
    }

    @staticmethod
    def _remove_lev(arr: xr.DataArray) -> xr.DataArray:
        return arr.sum(dim="lev")

    def ens2median(self, arr: list[xr.DataArray]) -> xr.DataArray:
        """Combine a list of arrays from different (known) ensembles into a median."""
        arr = vbm.shift_arrays(arr, daily=False)
        arr = vbm.shift_arrays(arr, custom=1, daily=False)
        arr = vbm.mean_flatten(arr, dims=["lat", "lon"])
        arr = vbm.data_array_operation(arr, self._remove_lev)
        arr_ = vbm.get_median(arr, xarray=True)
        arr_ = arr_[: int(12 * 16)]
        return arr_.assign_coords(time=vbm.dt2float(arr_.time.data) - 1850)

    def print_available(self) -> None:
        """Print all available data."""
        print(self._FINDER())

    def print(self) -> None:
        """Print all data that is being used."""
        print(self.oh_c)
        print(self.oh_m)
        print(self.oh_p)
        print(self.oh_s)
        print(self.oh_e)
        print(self.oh_m2)
        print(self.oh_m4)
        print(self.oh_p2)
        print(self.oh_p4)

    def compute(self) -> xr.Dataset:
        """Compute the global stratospheric mean OH for all simulations."""
        e2m = self.ens2median
        names = list(self._PLOT_DICT.keys())
        oh_c_xr = e2m(self.oh_c.load()).rename(names[0])
        oh_m_xr = e2m(self.oh_m.load()).rename(names[1])
        oh_p_xr = e2m(self.oh_p.load()).rename(names[2])
        oh_s_xr = e2m(self.oh_s.load()).rename(names[3])
        oh_e_xr = e2m(self.oh_e.load()).rename(names[4])
        oh_m2_xr = e2m(self.oh_m2.load()).rename(names[5])
        oh_m4_xr = e2m(self.oh_m4.load()).rename(names[6])
        oh_p2_xr = e2m(self.oh_p2.load()).rename(names[7])
        oh_p4_xr = e2m(self.oh_p4.load()).rename(names[8])
        return xr.merge(
            [
                oh_c_xr,
                oh_m_xr,
                oh_p_xr,
                oh_s_xr,
                oh_e_xr,
                oh_m2_xr,
                oh_m4_xr,
                oh_p2_xr,
                oh_p4_xr,
            ]
        )

    def load(self) -> xr.Dataset:
        """Get or generate the data as a xarray data set."""
        file = _SAVE_DIR / "oh.nc"
        if file.exists():
            return xr.load_dataset(file)
        self._SHOW = False
        ds = self.compute().assign_attrs(
            {"description": "Global stratospheric mean OH concentration."}
        )
        ds.to_netcdf(file)
        return ds

    def plot(self) -> None:
        """Plot the computed arrays."""
        plt.figure()
        plt.semilogy()
        data = self.load()
        for name, o in data.data_vars.items():
            plot_args = self._PLOT_DICT[name]
            o.plot(label=plot_args.lab, color=plot_args.c, ls=plot_args.ls)  # type: ignore[call-arg]
        plt.legend()
        plt.xlabel("Time after first eruption [yr]")
        plt.ylabel("OH concentration")
        plt.savefig(_SAVE_DIR / "oh")
        if self._SHOW:
            plt.show()
        plt.close("all")


class SO2BurdenPlot:
    """Plot the SO2 column burden."""

    _SHOW = True
    FINDER = (
        FindFiles()
        .find("TMSO2", "e_fSST1850", "h0")
        .keep_most_recent()
        .sort("sim", "ensemble")
        .copy
    )
    so2_c = FINDER().keep("control", "ens1")
    so2_m = FINDER().keep("medium", {f"ens{i}" for i in [2, 3, 4, 5]})
    so2_m2 = FINDER().keep("medium-2sep")
    so2_m4 = FINDER().keep("medium-4sep")
    so2_p = FINDER().keep("medium-plus", {f"ens{i}" for i in [2, 3, 4, 5]})
    so2_s = FINDER().keep("strong", {f"ens{i}" for i in [2, 3, 4, 5]})
    so2_e = FINDER().keep("size5000")
    so2_p2 = FINDER().keep("tt-2sep", {"ens1", "ens3"})
    so2_p4 = FINDER().keep("tt-4sep", {"ens1", "ens3"})
    _PLOT_DICT: ClassVar[dict[str, PlotArgs]] = {
        "CONTROL": PlotArgs(ls="-", c=COLORS[0], lab="CONTROL"),
        "S26": PlotArgs(ls="-", c=COLORS[1], lab="S26"),
        "S400": PlotArgs(ls="-", c=COLORS[2], lab="S400"),
        "S1629": PlotArgs(ls="-", c=COLORS[3], lab="S1629"),
        "S3000": PlotArgs(ls="-", c=COLORS[4], lab="S300"),
        "S26, 2sep": PlotArgs(ls=":", c=COLORS[1], lab="_S26-2SEP"),
        "S26, 4sep": PlotArgs(ls="--", c=COLORS[1], lab="_S26-4SEP"),
        "S400, 2sep": PlotArgs(ls=":", c=COLORS[2], lab="_S400-2SEP"),
        "S400, 4sep": PlotArgs(ls="--", c=COLORS[2], lab="_S400-4SEP"),
    }

    @staticmethod
    def ens2median(arr: list[xr.DataArray]) -> xr.DataArray:
        """Combine a list of arrays from different (known) ensembles into a median."""
        arr = vbm.shift_arrays(arr, daily=False)
        arr = vbm.shift_arrays(arr, custom=1, daily=False)
        arr = vbm.mean_flatten(arr, dims=["lat", "lon"])
        arr_ = vbm.get_median(arr, xarray=True)
        arr_ = arr_[: int(12 * 10)]
        return arr_.assign_coords(time=vbm.dt2float(arr_.time.data) - 1850)

    def print(self) -> None:
        """Print all data that is being used."""
        print(self.so2_c)
        print(self.so2_m)
        print(self.so2_p)
        print(self.so2_s)
        print(self.so2_e)
        print(self.so2_m2)
        print(self.so2_m4)
        print(self.so2_p2)
        print(self.so2_p4)

    def compute(self) -> xr.Dataset:
        """Compute the global mean TMSO2 for all simulations."""
        e2m = self.ens2median
        names = list(self._PLOT_DICT.keys())
        so2_c_xr = e2m(self.so2_c.load()).rename(names[0])
        so2_m_xr = e2m(self.so2_m.load()).rename(names[1])
        so2_p_xr = e2m(self.so2_p.load()).rename(names[2])
        so2_s_xr = e2m(self.so2_s.load()).rename(names[3])
        so2_e_xr = e2m(self.so2_e.load()).rename(names[4])
        so2_m2_xr = e2m(self.so2_m2.load()).rename(names[5])
        so2_m4_xr = e2m(self.so2_m4.load()).rename(names[6])
        so2_p2_xr = e2m(self.so2_p2.load()).rename(names[7])
        so2_p4_xr = e2m(self.so2_p4.load()).rename(names[8])
        return xr.merge(
            [
                so2_c_xr,
                so2_m_xr,
                so2_p_xr,
                so2_s_xr,
                so2_e_xr,
                so2_m2_xr,
                so2_m4_xr,
                so2_p2_xr,
                so2_p4_xr,
            ],
        )

    def load(self) -> xr.Dataset:
        """Get or generate the data as a xarray data set."""
        file = _SAVE_DIR / "so2-burden.nc"
        if file.exists():
            return xr.load_dataset(file)
        self._SHOW = False
        ds = self.compute().assign_attrs({"description": "Global mean SO2 burden."})
        ds.to_netcdf(file)
        return ds

    def plot(self) -> None:
        """Plot the computed arrays."""
        plt.figure()
        plt.semilogy()
        data = self.load()
        for name, s in data.data_vars.items():
            plot_args = self._PLOT_DICT[name]
            s.plot(label=plot_args.lab, color=plot_args.c, ls=plot_args.ls)  # type: ignore[call-arg]
        plt.legend()
        plt.xlabel("Time after first eruption [yr]")
        plt.savefig(_SAVE_DIR / "so2-burden")
        if self._SHOW:
            plt.show()
        plt.close("all")


def _main() -> None:
    ReffPlot().plot()
    SO2BurdenPlot().plot()
    OHPlot().plot()


if __name__ == "__main__":
    _main()
