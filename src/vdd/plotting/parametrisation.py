"""Parametrise the temperature from SO2, via RF."""

# TODO: Obtain an estimate of the T to RF response function.
#
# - [x] Get response functions representing T to SO2 and RF to SO2.
# - [x] Try deconvolving right away to obtain the T to RF response function. Both using
#       T and RF, and using the response functions.
# - [ ] Verify that the response functions in the previous step are zero at negative
#       time lags.
# - [ ] Try smoothing the response functions to obtain a more stable result.

import fppanalysis
import matplotlib.pyplot as plt
import numpy as np
import volcano_base

import vdd.load

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
dec_ob16 = vdd.load.DeconvolveOB16(data="h1")
dec_ob16_month = vdd.load.DeconvolveOB16(data="h0")
dec_ob16.name = "OB16"
dec_ob16_month.name = "OB16 month"
all_decs = (
    dec_cesm_4sep,
    dec_cesm_2sep,
    dec_cesm_e,
    dec_cesm_s,
    dec_cesm_p,
    dec_cesm_m,
    dec_ob16,
    dec_ob16_month,
)


class PlotParametrisation:
    """Parametrise the functional relationship between temperature and RF."""

    def __init__(self, *decs: vdd.load.Deconvolve):
        self.decs = decs
        self.results: dict[str, dict[str, np.ndarray]] = {}

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

    def strategy(self) -> None:
        """Obtain the temperature response function to RF."""
        for dec in self.decs:
            orig = dec.response_temp_rf
            s2 = self._strategy2(dec)
            s3 = self._strategy3(dec)
            self.results[dec.name] = {"tau": dec.tau, "s1": orig, "s2": s2, "s3": s3}


def main():
    """Run the main script."""
    pp = PlotParametrisation(*all_decs)
    pp.strategy()

    for key, value in pp.results.items():
        plt.figure()
        plt.suptitle(key)
        plt.xlim((-2, 21))
        x = value["tau"]
        s1 = value["s1"]
        s2 = value["s2"]
        s3 = value["s3"]
        plt.plot(x, s1, label="dec(T, RF)")
        plt.plot(x, s2, "--", label="dec(dec(T, SO2), dec(RF, SO2))")
        plt.plot(x, s3, ":", label="dec(dec(T, SO2), dec(RF, SO2)) corrected")
        plt.legend()
        plt.savefig(_SAVE_DIR / f"parametrisation_{key}.png")


if __name__ == "__main__":
    main()
