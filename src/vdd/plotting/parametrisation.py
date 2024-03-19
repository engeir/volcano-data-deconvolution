"""Parametrise the temperature from SO2, via RF."""

# TODO: Obtain an estimate of the T to RF response function.
#
# - [x] Get response functions representing T to SO2 and RF to SO2.
# - [x] Try deconvolving right away to obtain the T to RF response function. Both using
#       T and RF, and using the response functions.
# - [ ] Verify that the response functions in the previous step are zero at negative
#       time lags.
# - [ ] Try smoothing the response functions to obtain a more stable result.

import datetime

import cftime
import fppanalysis
import matplotlib.pyplot as plt
import numpy as np
import volcano_base
import xarray as xr
from rich.console import Console
from rich.table import Table

import vdd.load

plt.rc("text.latex", preamble=r"\usepackage{amsmath}")
plt.style.use(
    "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle"
)

_SAVE_DIR = volcano_base.config.SAVE_PATH / "parametrisation"
if not _SAVE_DIR.exists():
    _SAVE_DIR.mkdir(parents=False)

DataCESM = vdd.load.CESMData
DecCESM = vdd.load.DeconvolveCESM
# CESM2
dec_cesm_4sep = DecCESM(pad_before=True, cesm=DataCESM(strength="double-overlap"))
dec_cesm_2sep = DecCESM(pad_before=True, cesm=DataCESM(strength="tt-2sep"))
dec_cesm_e = DecCESM(pad_before=True, cesm=DataCESM(strength="size5000"))
dec_cesm_s = DecCESM(pad_before=True, cesm=DataCESM(strength="strong"))
dec_cesm_p = DecCESM(pad_before=True, cesm=DataCESM(strength="medium-plus"))
dec_cesm_m = DecCESM(pad_before=True, cesm=DataCESM(strength="medium"))
# OB16
dec_ob16_month = vdd.load.DeconvolveOB16(data="h0")
dec_ob16_month.name = "OB16 month"
all_decs = (
    dec_cesm_4sep,
    dec_cesm_2sep,
    dec_cesm_e,
    dec_cesm_s,
    dec_cesm_p,
    # dec_cesm_m,
    dec_ob16_month,
)


class ReconstructOB16:
    """Class that reconstructs the temperature of OB16 from CESM2 simulations."""

    def __init__(self, *decs: vdd.load.Deconvolve) -> None:
        self.ob16 = vdd.load.DeconvolveOB16(data="h0")
        self.ob16.name = "OB16 month"
        self.decs = decs

    def plot_temperature(self) -> None:
        """Plot the reconstructed temperatures."""

        def d2n(date: datetime.datetime) -> float:
            unit = "days since 2000-01-01"
            return cftime.date2num(
                date, units=unit, calendar="noleap", has_year_zero=True
            )

        xlim = (
            d2n(datetime.datetime(1250, 1, 1, 0, 0)),
            d2n(datetime.datetime(1350, 1, 1, 0, 0)),
        )
        all_f = plt.figure()
        all_a = all_f.gca()
        all_a.plot(self.ob16.temp.time, self.ob16.temp, label=self.ob16.name)
        all_zoom_f = plt.figure()
        all_zoom_a = all_zoom_f.gca()
        all_zoom_a.plot(self.ob16.temp.time, self.ob16.temp, label=self.ob16.name)
        all_zoom_a.set_xlim(xlim)
        res = []
        for dec in self.decs:
            fn = vdd.utils.clean_filename(dec.name)
            inv_f = plt.figure()
            inv_a = inv_f.gca()
            inv_zoom_f = plt.figure()
            inv_zoom_a = inv_zoom_f.gca()
            inv_zoom_a.set_xlim(xlim)
            response = dec.response_temp_rf
            response_scaled = (
                response / response.max() * self.ob16.response_temp_rf.max()
            )
            new_temp = np.convolve(self.ob16.rf, response, mode="same")
            new_temp_scaled = np.convolve(self.ob16.rf, response_scaled, mode="same")
            all_a.plot(self.ob16.temp.time, new_temp_scaled, label=dec.name)
            all_zoom_a.plot(self.ob16.temp.time, new_temp_scaled, label=dec.name)
            inv_a.plot(self.ob16.temp.time, self.ob16.temp, label="OB16 temperature")
            inv_a.plot(self.ob16.temp.time, new_temp_scaled, label="Scaled response")
            inv_a.plot(self.ob16.temp.time, new_temp, label="Raw response")
            inv_a.legend()
            inv_f.savefig(_SAVE_DIR / f"reconstruct_from_{fn}.png")
            inv_zoom_a.plot(
                self.ob16.temp.time, self.ob16.temp, label="OB16 temperature"
            )
            inv_zoom_a.plot(
                self.ob16.temp.time, new_temp_scaled, label="Scaled response"
            )
            inv_zoom_a.plot(self.ob16.temp.time, new_temp, label="Raw response")
            inv_zoom_a.legend()
            inv_zoom_f.savefig(_SAVE_DIR / f"reconstruct_from_{fn}_zoom.png")
            # Print the distance away from the reconstructed
            rob16 = self.ob16.response_temp_rf
            ob16_temp = np.convolve(self.ob16.rf, rob16, mode="same")
            ob16_diff = np.abs(ob16_temp - new_temp).sum()
            ob16_diff_scaled = np.abs(ob16_temp - new_temp_scaled).sum()
            res.append((dec.name, f"{ob16_diff:.2f}", f"{ob16_diff_scaled:.2f}"))
        table = Table(
            title="Difference between reconstructed temperature from OB16 and other simulations"
        )
        table.add_column("Simulation name", justify="left", style="cyan", no_wrap=True)
        table.add_column("Raw response", justify="center", style="magenta")
        table.add_column("Scaled response", justify="center", style="magenta")
        for r_ in res:
            table.add_row(*r_)
        console = Console()
        console.print(table)
        all_a.legend()
        all_f.savefig(_SAVE_DIR / "reconstruct_from_all.png")
        all_zoom_a.legend()
        all_zoom_f.savefig(_SAVE_DIR / "reconstruct_from_all_zoom.png")


class PlotParametrisation:
    """Parametrise the functional relationship between temperature and RF."""

    def __init__(self, *decs: vdd.load.Deconvolve):
        self.decs = decs
        self.results: dict[str, xr.Dataset] = {}

    @staticmethod
    def _strategy2(dec: vdd.load.Deconvolve) -> np.ndarray:
        """Obtain the response function from deconvolving `res_t_so2` with `res_rf_so2`."""
        signal = dec.response_temp_so2
        kern = dec.response_rf_so2
        return fppanalysis.RL_gauss_deconvolve(signal, kern, 200)[0].flatten()

    @staticmethod
    def _strategy3(dec: vdd.load.Deconvolve) -> np.ndarray:
        """Obtain the corrected response function from deconvolving `res_t_so2` with `res_rf_so2`."""
        signal = dec.response_temp_so2
        kern = dec.response_rf_so2
        signal[: len(signal) // 2] = 0
        kern[: len(kern) // 2] = 0
        return fppanalysis.RL_gauss_deconvolve(signal, kern, 200)[0].flatten()

    @staticmethod
    def _strategy4(dec: vdd.load.Deconvolve) -> np.ndarray:
        """Compute the same as strategy 3, but we smooth the kernel."""
        signal = dec.response_temp_so2
        kern = dec.response_rf_so2
        signal[: len(signal) // 2] = 0
        pad = 6
        kern = np.pad(kern, pad)
        kern = fppanalysis.run_mean(kern, pad)
        kern[: len(kern) // 2] = 0
        return fppanalysis.RL_gauss_deconvolve(signal, kern, 200)[0].flatten()

    def strategy(self) -> None:
        """Obtain the temperature response function to RF."""
        for dec in self.decs:
            orig = dec.response_temp_rf
            s2 = self._strategy2(dec)
            s3 = self._strategy3(dec)
            # s4 = self._strategy4(dec)
            dataset = xr.Dataset(
                {
                    "s1": (
                        "tau",
                        orig,
                        {
                            "long_name": "Response function",
                            "label": "dec(T, RF)",
                            "ls": "-",
                        },
                    ),
                    "s2": (
                        "tau",
                        s2,
                        {
                            "long_name": "Response function",
                            "label": "dec(dec(T, SO2), dec(RF, SO2))",
                            "ls": "--",
                        },
                    ),
                    "s3": (
                        "tau",
                        s3,
                        {
                            "long_name": "Response function",
                            "label": "dec(dec(T, SO2), dec(RF, SO2)) corrected",
                            "ls": ":",
                        },
                    ),
                    # "s4": ("tau", s4),
                },
                coords={"tau": dec.tau},
            )
            self.results[dec.name] = dataset

    def plot_simulations(self) -> None:
        """Plot the responses and their means from each simulation."""
        ens = plt.figure()
        ens_a = ens.gca()
        for key, value in self.results.items():
            # Individual plot
            ind = plt.figure()
            ind_a = ind.gca()
            ind_a.set_xlim((-2, 21))
            for v in value.data_vars.values():
                ind_a.plot(v.tau, v, label=v.attrs["label"], ls=v.attrs["ls"])
            ind_a.legend()
            ind.savefig(
                _SAVE_DIR / f"parametrisation_{vdd.utils.clean_filename(key)}.png"
            )
            # Ensemble plot
            mean = np.mean([v for v in value.data_vars.values()], axis=0)
            ens_a.plot(value.tau, mean, label=f"{key} mean")
        ens_a.set_xlim((-2, 21))
        ens_a.legend()
        ens.savefig(_SAVE_DIR / "parametrisation_ensemble.png")

    def plot_method_comparisons(self) -> None:
        """Plot the mean for each response strategy and the full mean."""
        final = plt.figure()
        final_a = final.gca()
        method = plt.figure()
        method_a = method.gca()
        final_xr: list[tuple[xr.DataArray, ...]] = []
        final_xr_app = final_xr.append
        final_mean: list[np.ndarray] = []
        final_mean_app = final_mean.append
        for i in range(len(list(self.results.values())[0].data_vars.keys())):
            xr_aligned: tuple[xr.DataArray, ...] = xr.align(
                *volcano_base.manipulate.approx_align(
                    *[value[f"s{i+1}"] for value in self.results.values()]
                )
            )
            final_xr_app(xr_aligned)
            xr_mean = np.mean(xr_aligned, axis=0)
            final_mean_app(xr_mean)
            # label = xr_aligned[0].attrs["label"]
            label = f"s{i+1} mean"
            ls = xr_aligned[0].attrs["ls"]
            method_a.plot(xr_aligned[0].tau, xr_mean, label=label, ls=ls)
        final_plot = np.mean(final_mean, axis=0)
        final_a.plot(final_xr[0][0].tau, final_plot, label="Full mean")
        final_a.set_xlim((-2, 21))
        final_a.legend()
        final.savefig(_SAVE_DIR / "parametrisation_final.png")
        method_a.set_xlim((-2, 21))
        method_a.legend()
        method.savefig(_SAVE_DIR / "parametrisation_method.png")


def _plot_reconstructed_temperature() -> None:
    rec = ReconstructOB16(*all_decs)
    rec.plot_temperature()
    plt.show()


def _plot_response_functions() -> None:
    pp = PlotParametrisation(*all_decs)
    pp.strategy()
    pp.plot_simulations()
    plt.show()
    pp.plot_method_comparisons()
    plt.show()


def main():
    """Run the main script."""
    _plot_reconstructed_temperature()
    _plot_response_functions()


if __name__ == "__main__":
    main()
