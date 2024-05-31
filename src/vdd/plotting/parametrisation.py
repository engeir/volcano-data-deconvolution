"""Explore ways of estimating the temperature response pulse function to RF.

Parametrise the temperature from SO2, via RF.
"""

# INFO: Obtain an estimate of the T to RF response function.
# - [x] Get response functions representing T to SO2 and RF to SO2.
# - [x] Try deconvolving right away to obtain the T to RF response function. Both using
#       T and RF, and using the response functions.
# - [x] Verify that the response functions in the previous step are zero at negative
#       time lags.
# - [x] Try smoothing the response functions to obtain a more stable result. (Bad.)

import fppanalysis
import matplotlib.pyplot as plt
import numpy as np
import volcano_base
import xarray as xr

import vdd.load
from vdd.utils import name_swap as ns

plt.style.use([
    "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle",
    "vdd.extra",
])

_SAVE_DIR = volcano_base.config.SAVE_PATH / "parametrisation"
if not _SAVE_DIR.exists():
    _SAVE_DIR.mkdir(parents=False)

DataCESM = vdd.load.CESMData
DecCESM = vdd.load.DeconvolveCESM
# CESM2
dec_cesm_4sep = DecCESM(pad_before=True, cesm=DataCESM(strength="tt-4sep"))
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
    dec_cesm_m,
    dec_ob16_month,
)


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
        return dec._deconv_method(signal, kern)[0].flatten()

    @staticmethod
    def _strategy3(dec: vdd.load.Deconvolve) -> np.ndarray:
        """Obtain the corrected response function from deconvolving `res_t_so2` with `res_rf_so2`."""
        signal = dec.response_temp_so2
        kern = dec.response_rf_so2
        signal[: len(signal) // 2] = 0
        kern[: len(kern) // 2] = 0
        return dec._deconv_method(signal, kern)[0].flatten()

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
        return dec._deconv_method(signal, kern)[0].flatten()

    def strategy(self) -> None:
        """Obtain the temperature response function to RF."""
        for dec in self.decs:
            orig = dec.response_temp_rf
            s2 = self._strategy2(dec)
            # s3 = self._strategy3(dec)
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
                    # "s3": (
                    #     "tau",
                    #     s3,
                    #     {
                    #         "long_name": "Response function",
                    #         "label": "dec(dec(T, SO2), dec(RF, SO2)) corrected",
                    #         "ls": ":",
                    #     },
                    # ),
                    # "s4": ("tau", s4),
                },
                coords={"tau": dec.tau},
            )
            self.results[ns(dec.name)] = dataset

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
                _SAVE_DIR / ns(f"parametrisation_{vdd.utils.clean_filename(key)}.jpg")
            )
            # Ensemble plot
            mean = np.mean([v for v in value.data_vars.values()], axis=0)
            ens_a.plot(value.tau, mean, label=f"{key} mean")
        ens_a.set_xlim((-2, 21))
        ens_a.legend()
        ens.savefig(_SAVE_DIR / "parametrisation_ensemble.jpg")

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
                *volcano_base.manipulate.approx_align(*[
                    value[f"s{i + 1}"] for value in self.results.values()
                ])
            )
            final_xr_app(xr_aligned)
            xr_mean = np.mean(xr_aligned, axis=0)
            final_mean_app(xr_mean)
            # label = xr_aligned[0].attrs["label"]
            label = f"s{i + 1} mean"
            ls = xr_aligned[0].attrs["ls"]
            method_a.plot(xr_aligned[0].tau, xr_mean, label=label, ls=ls)
        final_plot = np.mean(final_mean, axis=0)
        final_a.plot(final_xr[0][0].tau, final_plot, label="Full mean")
        final_a.set_xlim((-2, 21))
        final_a.legend()
        final.savefig(_SAVE_DIR / "parametrisation_final.jpg")
        method_a.set_xlim((-2, 21))
        method_a.legend()
        method.savefig(_SAVE_DIR / "parametrisation_method.jpg")


def _plot_response_functions() -> None:
    pp = PlotParametrisation(*all_decs)
    pp.strategy()
    pp.plot_simulations()
    plt.show()
    # Combine the response functions I'm interested in
    vdd.utils.combine(*[
        _SAVE_DIR / f"parametrisation_{ns(file)}.jpg"
        for file in [
            "ob16-month",
            "cesm2-medium-plus",
            "cesm2-strong",
            "cesm2-size5000",
            "cesm2-tt-2sep",
            "cesm2-tt-4sep",
        ]
    ]).in_grid(2, 3).save(_SAVE_DIR / "parametrisation_combined.jpg")
    pp.plot_method_comparisons()
    plt.show()


def main():
    """Run the main script."""
    _plot_response_functions()


if __name__ == "__main__":
    main()
