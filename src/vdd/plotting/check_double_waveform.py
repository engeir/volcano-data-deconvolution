"""Check how well response functions are estimated from double waveforms.

We also check how well we are able to recreate the double waveform time series from the
response functions of other simulations.
"""

from typing import Literal

import cosmoplots
import matplotlib.pyplot as plt
import numpy as np
import volcano_base
import xarray as xr
from matplotlib import patheffects

import vdd.load

plt.rc("text.latex", preamble=r"\usepackage{amsmath}")
plt.style.use(
    "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle"
)
_COLORS = plt.rcParams["axes.prop_cycle"].by_key()["color"]

_SAVE_DIR = volcano_base.config.SAVE_PATH / "waveform"
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


def check_waveform_responses(*decs: vdd.load.DeconvolveCESM) -> None:
    """Check how well response functions are estimated from double waveforms."""
    aod_f = plt.figure()
    aod_a = aod_f.gca()
    rf_f = plt.figure()
    rf_a = rf_f.gca()
    temp_f = plt.figure()
    temp_a = temp_f.gca()
    for dec in decs:
        # Check how well we are able to recreate the double waveform time series from the
        # response functions of other simulations.
        tau = dec.tau
        name = "2sep" if "2sep" in dec.name else "4sep"
        aod_a.plot(tau, dec.response_aod_so2, label=name)
        rf_a.plot(tau, dec.response_rf_so2, label=name)
        temp_a.plot(tau, dec.response_temp_so2, label=name)
    aod_a.set_xlabel("Time lag ($\\tau$) [yr]")
    aod_a.set_ylabel("Aerosol optical depth response to SO2")
    rf_a.set_xlabel("Time lag ($\\tau$) [yr]")
    rf_a.set_ylabel("Radiative forcing response to SO2")
    temp_a.set_xlabel("Time lag ($\\tau$) [yr]")
    temp_a.set_ylabel("Temperature response to SO2")
    aod_a.set_xlim((-2, 15))
    rf_a.set_xlim((-2, 15))
    temp_a.set_xlim((-2, 15))
    aod_a.legend()
    rf_a.legend()
    temp_a.legend()
    files = [_SAVE_DIR / f"responses_{k}.png" for k in ["aod", "rf", "temp"]]
    aod_f.savefig(files[0])
    rf_f.savefig(files[1])
    temp_f.savefig(files[2])
    plt.show()


def curve_fit_aod(aod, alpha, beta):
    """Fit values using AOD as input."""
    return alpha * np.log(beta * aod + 1)


class CheckRecreatedWaveforms:
    """Check how well we are able to recreate the double waveform time series."""

    def __init__(
        self,
        *decs: vdd.load.DeconvolveCESM,
        scale_by_aod: Literal["log", "log-inside", "root"] | bool = False,
    ):
        self.decs = decs
        self.scale_by_aod = scale_by_aod
        self.figs = {"aod": plt.figure(), "rf": plt.figure(), "temp": plt.figure()}
        self.axs = {k: v.gca() for k, v in self.figs.items()}
        self.response_rf = dec_cesm_s.response_rf_so2

    def run_loop(self) -> None:
        """Run the main loop."""
        for i, dec in enumerate(self.decs):
            so2_new = self._get_so2_new(dec)
            for attr in ["aod", "rf", "temp"]:
                self._run_attr_loop(dec, attr, so2_new, i)
        [self.axs[k].set_xlabel("Time [yr]") for k in self.axs]
        [self.axs[k].set_xlim((-2, 25)) for k in self.axs]
        [self.axs[k].legend(framealpha=0.5) for k in self.axs]
        self.axs["aod"].set_ylabel("Aerosol optical depth [1]")
        self.axs["rf"].set_ylabel("Radiative Forcing [W/m$^2$]")
        self.axs["temp"].set_ylabel("Temperature [K]")
        corrected = f"-aod-{self.scale_by_aod}-corrected" if self.scale_by_aod else ""
        files = [
            _SAVE_DIR / f"recreated_waveforms_{k}{corrected}.png"
            for k in ["aod", "rf", "temp"]
        ]
        self.figs["aod"].savefig(files[0])
        self.figs["rf"].savefig(files[1])
        self.figs["temp"].savefig(files[2])
        cosmoplots.combine(*files).in_grid(1, 3).using(fontsize=50).save(
            _SAVE_DIR / f"responses_combined{corrected}.png"
        )
        for f in files:
            f.unlink()
        plt.show()

    def _get_so2_new(self, dec: vdd.load.DeconvolveCESM) -> xr.DataArray:
        so2 = dec.so2
        so2_new = so2.copy()
        idx = so2_new > 0
        extra_scale_arr = dec.aod[idx] / dec.aod.max()
        match self.scale_by_aod:
            case "log":
                so2_new[idx] = (
                    so2_new[idx] * np.log(1 + 1 - extra_scale_arr) / np.log(2)
                    if self.scale_by_aod
                    else so2_new[idx]
                )
            case "log-inside":
                so2_new[idx] = (
                    so2_new[idx] * 1 - np.log(1 + extra_scale_arr) / np.log(2)
                    if self.scale_by_aod
                    else so2_new[idx]
                )
            case "root":
                so2_new[idx] = (
                    so2_new[idx] * (1 - extra_scale_arr) ** 0.5
                    if self.scale_by_aod
                    else so2_new[idx]
                )
        return so2_new

    def _run_attr_loop(
        self, dec: vdd.load.DeconvolveCESM, attr: str, so2_new: xr.DataArray, i: int
    ) -> None:
        name = "2sep" if "2sep" in dec.name else "4sep"
        so2 = dec.so2
        diff_len = len(self.response_rf) - len(so2)
        pe = [
            patheffects.Stroke(linewidth=1, foreground=_COLORS[i]),
            patheffects.Normal(),
        ]
        arr = getattr(dec, attr)
        dec_resp = getattr(dec, f"response_{attr}_so2")
        resp_arr = getattr(dec_cesm_s, f"response_{attr}_so2")[
            diff_len // 2 : -diff_len // 2
        ]
        # Scale r_arr
        scale_arr = np.max(dec_resp) / np.max(resp_arr)
        resp_arr = resp_arr * scale_arr
        # Scale so2 for new array
        rec_same = np.convolve(dec_resp, so2, "same")
        rec_new = np.convolve(resp_arr, so2_new, "same")
        plot = self.axs[attr].plot
        if attr == "aod":
            # Area of the Earth is 5.1e14 m2
            plot(dec.tau, dec.tmso2 * 5.1e3, label="TMSO2")
        elif attr == "rf":
            # params, _ = curve_fit(curve_fit_aod, dec.aod, dec.rf)
            # print(params)
            # [53.04185578  0.24951242]  # 2sep
            # [59.44007657  0.23123325]  # 4sep
            params = [22.34317318, 0.97876114]  # Single peaks
            plot(
                dec.tau,
                curve_fit_aod(dec.aod, *params),
                c="r",
                lw=0.5,
                label=f"{name} AOD",
            )
        elif attr == "temp":
            # params, _ = curve_fit(curve_fit_aod, dec.aod, dec.temp)
            # print(params)
            # [ 1.18746129 28.08689995]  # 2sep
            # [ 1.03776542 65.61337388]  # 4sep
            # params = [3.02039794, 1.19424641]  # Single peaks
            params = [22.34317318, 0.97876114]  # Single peaks, RF/AOD
            # resp_temp_rf = dec_cesm_s.response_temp_rf[
            #     diff_len // 2 : -diff_len // 2
            # ]
            resp_temp_rf = dec.response_temp_rf
            temp = np.convolve(
                resp_temp_rf, curve_fit_aod(dec.aod, *params), mode="same"
            )
            plot(dec.tau, temp, c="r", lw=0.5, label=f"{name} AOD")
        plot(dec.tau, arr, c="k", lw=0.5, label=f"{name} original")
        kwargs = {"ls": "--", "c": "k", "lw": 0.5, "path_effects": pe}
        plot(dec.tau, rec_same, label=f"{name} reconstruct self", **kwargs)  # type: ignore
        plot(dec.tau, rec_new, c=_COLORS[i], label=f"{name} reconstruc other")


def main() -> None:
    """Run the main script."""
    # check_waveform_responses(dec_cesm_2sep, dec_cesm_4sep)
    CheckRecreatedWaveforms(dec_cesm_2sep, dec_cesm_4sep, scale_by_aod="log").run_loop()
    CheckRecreatedWaveforms(
        dec_cesm_2sep, dec_cesm_4sep, scale_by_aod="log-inside"
    ).run_loop()
    CheckRecreatedWaveforms(
        dec_cesm_2sep, dec_cesm_4sep, scale_by_aod="root"
    ).run_loop()
    CheckRecreatedWaveforms(dec_cesm_2sep, dec_cesm_4sep).run_loop()


if __name__ == "__main__":
    main()
