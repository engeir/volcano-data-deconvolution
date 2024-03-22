"""Analysis of how cutting off the response functions affect the result."""

# TODO: Using the ideas from `volcano-deconv`, implement a simple CutOff class that will
# repeatedly cut off estimated response functions, recreate, add noise (for example
# phase randomised input signal), and from this create an ensemble of response
# estimates.
#
# 1. Implement the CutOff class.
#    - [ ] Owns a Deconvolve object.
#    - [ ] Uses only one response function and corresponding signal (RF or temperature).

from collections.abc import Iterable
from functools import cached_property
from typing import Literal, Self

import fppanalysis
import matplotlib.pyplot as plt
import numpy as np
import plastik
import volcano_base
import xarray as xr

import vdd.load
import vdd.utils

plt.rc("text.latex", preamble=r"\usepackage{amsmath}")
plt.style.use(
    "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle"
)
_SAVE_DIR = volcano_base.config.SAVE_PATH / "cut_off"

if not _SAVE_DIR.exists():
    _SAVE_DIR.mkdir(parents=False)

type T_RF = tuple[Literal["temp"], Literal["rf"]]  # type: ignore
type T_SO2 = tuple[Literal["temp"], Literal["so2"]]  # type: ignore
# Only got temp ctrl from Otto-Bliesner dataset
type RF_SO2 = tuple[Literal["rf"], Literal["so2"]]  # type: ignore

DataCESM = vdd.load.CESMData
DecCESM = vdd.load.DeconvolveCESM
# CESM2
dec_cesm_4sep = DecCESM(pad_before=True, cesm=DataCESM(strength="double-overlap"))
dec_cesm_2sep = DecCESM(pad_before=True, cesm=DataCESM(strength="tt-2sep"))
dec_cesm_s = DecCESM(pad_before=True, cesm=DataCESM(strength="strong"))
dec_ob16_month = vdd.load.DeconvolveOB16(data="h0")
dec_ob16_month.name = "OB16 month"
all_decs = (dec_cesm_4sep, dec_cesm_2sep, dec_cesm_s, dec_ob16_month)


class CutOff:
    """Cut off the response functions of a deconvolution object."""

    def __init__(self, dec: vdd.load.Deconvolve, arrays: T_RF | T_SO2 | RF_SO2) -> None:
        self.dec = dec
        self.ts_specifier: T_RF | T_SO2 | RF_SO2 = arrays
        self.cuts: dict[str, xr.Dataset] = {}
        self.ensembles: dict[str, xr.Dataset] = {}

    @cached_property
    def response(self) -> np.ndarray:
        """The response function in the convolution."""
        out = getattr(
            self.dec, f"response_{self.ts_specifier[0]}_{self.ts_specifier[1]}"
        )
        out[self.dec.tau <= 0] = 0
        return out

    @cached_property
    def forcing(self) -> xr.DataArray:
        """The forcing time series in the convolution."""
        return getattr(self.dec, self.ts_specifier[1])

    @cached_property
    def output(self) -> xr.DataArray:
        """The final output time series of the convolution."""
        return getattr(self.dec, self.ts_specifier[0])

    @cached_property
    def control(self) -> xr.DataArray:
        """The control time series in the convolution."""
        return getattr(self.dec, f"{self.ts_specifier[0]}_control")

    def cut_off(self, cutoff: int | Iterable[int]) -> Self:
        """Cut off the response function at a given time lag."""
        match cutoff:
            case int():
                self._single_cut_off(cutoff)
            case Iterable():
                for c in set(cutoff):
                    if not isinstance(c, int):
                        raise ValueError(
                            "cutoff must be an integer or a sequence of integers."
                        )
                    self._single_cut_off(c)
            case _:
                raise ValueError("cutoff must be an integer or a sequence of integers.")
        return self

    def _single_cut_off(self, cutoff: int) -> None:
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

    def generate_ensembles(self, n: int) -> None:
        """Generate an ensemble of response function estimates."""
        if not self.cuts:
            raise ValueError("No cuts have been made.")
        iters = np.arange(200)
        for k, v in self.cuts.items():
            if k in self.ensembles:
                continue
            arrays: dict[str, tuple] = {}
            for i in range(n):
                temp_rec = v.temp_rec.copy()
                temp_random = fppanalysis.signal_rand_phase(self.control.data)
                temp_rec += temp_random
                res_rec, err = fppanalysis.RL_gauss_deconvolve(
                    temp_rec, self.forcing, len(iters) - 1
                )
                r_cut_rec = res_rec.flatten()
                r_cut_rec[self.dec.tau <= 0] = 0
                arrays[f"response_{i}"] = ("tau", r_cut_rec, {"label": f"response {i}"})
                arrays[f"iters_{i}"] = ("iters", err.flatten(), {"label": f"err {i}"})
                # arrays[f"temp_{i}"] = rec
            self.ensembles[k] = xr.Dataset(
                arrays,
                coords={"tau": self.dec.tau, "time": self.output.time, "iters": iters},
            )


class PlotCutOff:
    """Plot the results of the CutOff class for any deconvolution object."""

    def __init__(self, *cut_offs: CutOff) -> None:
        self.cut_offs = cut_offs

    def call_cut_offs(self, method: str, *args, **kwargs) -> None:
        """Call a method on all CutOff objects."""
        match method, args, kwargs:
            case method, (arg1,), _:
                for co in self.cut_offs:
                    getattr(co, method)(arg1)

    def plot(self) -> None:
        """Plot the results of the CutOff class."""
        if not self.cut_offs:
            raise ValueError(
                "No CutOff objects have been created. Run `call_cut_offs`."
            )
        for co in self.cut_offs:
            if not co.cuts:
                raise ValueError(
                    "No cuts have been made. Run `call_cut_offs('cut_off', ...)`."
                )
            self._plot_single(co)
            plt.close("all")

    @staticmethod
    def _plot_single(co: CutOff) -> None:
        """Plot the results of the CutOff class."""
        for k, v in co.cuts.items():
            resp_f = plt.figure()
            resp_a = resp_f.gca()
            temp_f = plt.figure()
            temp_a = temp_f.gca()
            resp_a.axvline(int(k) / 12, c="k", ls="--")
            resp_a.plot(co.dec.tau, co.response, label="Response")
            resp_a.plot(v.tau, v.response, label=v.response.attrs["label"])
            temp_a.plot(co.output.time, co.output, label="Temp")
            temp_a.plot(v.time, v.temp_rec, label=v.temp_rec.attrs["label"])
            percs = []
            for i, j in co.ensembles[k].items():
                if "response" in str(i):
                    percs.append(j)
                    # label = j.attrs["label"]
                    # resp_a.plot(j.tau, j, label=f"_{label}", c="grey", alpha=0.3)
            percs_np = np.asarray(percs)
            resp_a = plastik.percentiles(
                co.dec.tau, percs_np, plot_median=False, ax=resp_a
            )
            resp_a.set_xlabel("Time lag ($\\tau$) [yr]")
            resp_a.set_ylabel("Response [1]")
            resp_a.set_xlim((-2, 20))
            ymax = co.response.max()
            resp_a.set_ylim((ymax * (-0.05), ymax * 1.05))
            resp_a.legend(loc="upper right", framealpha=0.9)
            num = "0" * (3 - len(k)) + k
            name = vdd.utils.clean_filename(co.dec.name)
            ts = vdd.utils.clean_filename("-".join(co.ts_specifier))
            resp_f.savefig(_SAVE_DIR / f"{name}_resp_{ts}_{num}.png")
            temp_f.savefig(_SAVE_DIR / f"{name}_temp_{ts}_{num}.png")


def _use_cut_off() -> None:
    # OB16
    co_ob16_rf_so2 = CutOff(dec_ob16_month, ("rf", "so2"))
    co_ob16_temp_so2 = CutOff(dec_ob16_month, ("temp", "so2"))
    co_ob16_temp_rf = CutOff(dec_ob16_month, ("temp", "rf"))
    # CESM2 strong
    co_cesm_s_rf_so2 = CutOff(dec_cesm_s, ("rf", "so2"))
    co_cesm_s_temp_so2 = CutOff(dec_cesm_s, ("temp", "so2"))
    co_cesm_s_temp_rf = CutOff(dec_cesm_s, ("temp", "rf"))
    # CESM2 2sep
    co_2sep_rf_so2 = CutOff(dec_cesm_2sep, ("rf", "so2"))
    co_2sep_temp_so2 = CutOff(dec_cesm_2sep, ("temp", "so2"))
    co_2sep_temp_rf = CutOff(dec_cesm_2sep, ("temp", "rf"))
    # CESM2 4sep
    co_4sep_rf_so2 = CutOff(dec_cesm_4sep, ("rf", "so2"))
    co_4sep_temp_so2 = CutOff(dec_cesm_4sep, ("temp", "so2"))
    co_4sep_temp_rf = CutOff(dec_cesm_4sep, ("temp", "rf"))
    pco = PlotCutOff(
        co_ob16_rf_so2,
        co_ob16_temp_so2,
        co_ob16_temp_rf,
        co_cesm_s_rf_so2,
        co_cesm_s_temp_so2,
        co_cesm_s_temp_rf,
        co_2sep_rf_so2,
        co_2sep_temp_so2,
        co_2sep_temp_rf,
        co_4sep_rf_so2,
        co_4sep_temp_so2,
        co_4sep_temp_rf,
    )
    pco.call_cut_offs("cut_off", {12 * i for i in [2, 4, 8, 16]})
    pco.call_cut_offs("generate_ensembles", 50)
    pco.plot()


def _main() -> None:
    _use_cut_off()
    plt.show()


if __name__ == "__main__":
    _main()
