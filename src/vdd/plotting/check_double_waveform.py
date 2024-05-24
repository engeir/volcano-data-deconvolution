"""Check how well response functions are estimated from double waveforms.

We also check how well we are able to recreate the double waveform time series from the
response functions of other simulations.
"""

from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
import volcano_base
import xarray as xr

import vdd.load
import vdd.utils

plt.style.use([
    "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle",
    "vdd.extra",
])
_COLORS = plt.rcParams["axes.prop_cycle"].by_key()["color"]

_SAVE_DIR = volcano_base.config.SAVE_PATH / "waveform"
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
    aod_a.set_xlim((-2, 21))
    rf_a.set_xlim((-2, 21))
    temp_a.set_xlim((-2, 21))
    aod_a.legend()
    rf_a.legend()
    temp_a.legend()
    files = [_SAVE_DIR / f"responses_{k}.jpg" for k in ["aod", "rf", "temp"]]
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
        single_waveform: vdd.load.DeconvolveCESM = dec_cesm_p,
        scale_by_aod: Literal["log", "log-inside", "root"] | bool = False,
    ):
        self.single_waveform = single_waveform
        self.decs = decs
        self.scale_by_aod = scale_by_aod
        self.figs, self.axs = vdd.utils.figure_multiple_rows_columns(
            3, 2, share_axes="x"
        )
        self.keys = {"aod": (0, 1), "rf": (2, 3), "temp": (4, 5)}

    def run_loop(self) -> None:
        """Run the main loop."""
        for i, dec in enumerate(self.decs):
            so2_new = self._get_so2_new(dec)
            for attr in ["aod", "rf", "temp"]:
                self._run_attr_loop(dec, attr, so2_new, i)
        [ax.set_xlabel("Time [yr]") for ax in self.axs]
        [ax.set_xlim((-1, 21)) for ax in self.axs]
        [ax.legend(framealpha=0.5) for ax in self.axs]
        [self.axs[i].set_ylabel("Aerosol optical depth [1]") for i in self.keys["aod"]]
        [self.axs[i].set_ylabel("Radiative forcing [W/m$^2$]") for i in self.keys["rf"]]
        [self.axs[i].set_ylabel("Temperature anomaly [K]") for i in self.keys["temp"]]
        corrected = f"-aod-{self.scale_by_aod}-corrected" if self.scale_by_aod else ""
        self.figs.savefig(_SAVE_DIR / f"responses_combined{corrected}")
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
        idx = 0 if name == "2sep" else 1
        c_idx = 4 if name == "2sep" else 5
        so2 = dec.so2
        diff_len = len(self.single_waveform.response_rf_so2) - len(so2)
        arr = getattr(dec, attr)
        if diff_len < 0:
            diff_len = abs(diff_len)
            dec_resp = getattr(dec, f"response_{attr}_so2")[
                diff_len // 2 : -diff_len // 2
            ]
            resp_arr = getattr(self.single_waveform, f"response_{attr}_so2")
        elif diff_len > 0:
            dec_resp = getattr(dec, f"response_{attr}_so2")
            resp_arr = getattr(self.single_waveform, f"response_{attr}_so2")[
                diff_len // 2 : -diff_len // 2
            ]
        else:
            dec_resp = getattr(dec, f"response_{attr}_so2")
            resp_arr = getattr(self.single_waveform, f"response_{attr}_so2")
        rec_same = np.convolve(dec_resp, so2, "same")
        rec_new = np.convolve(resp_arr, so2_new, "same")
        plot = self.axs[self.keys[attr][idx]].plot
        plot(
            dec.tau,
            arr,
            c="k",
            lw=0.5,
            label=f"${attr[0].upper()}_{{\\mathrm{{TT-{name.upper()}}}}}$",
        )
        kwargs = {"ls": "--", "c": _COLORS[c_idx], "lw": 1.0}
        varphi = f"\\varphi_{{{attr[0].upper()}}}"
        conv_so2 = f"\\ast S_{{\\mathrm{{TT-{name.upper()}}}}}"
        plot(
            dec.tau,
            rec_new,
            "-.",
            lw=1.0,
            c=_COLORS[1],
            label=f"${varphi}^{{\\mathrm{{INTERMEDIATE}}}}{conv_so2}$",
        )
        plot(
            dec.tau,
            rec_same,
            label=f"${varphi}^{{\\mathrm{{TT-{name.upper()}}}}}{conv_so2}$",
            **kwargs,
        )


def main() -> None:
    """Run the main script."""
    # check_waveform_responses(dec_cesm_2sep, dec_cesm_4sep)
    # CheckRecreatedWaveforms(dec_cesm_2sep, dec_cesm_4sep, scale_by_aod="log").run_loop()
    # CheckRecreatedWaveforms(
    #     dec_cesm_2sep, dec_cesm_4sep, scale_by_aod="log-inside"
    # ).run_loop()
    # CheckRecreatedWaveforms(
    #     dec_cesm_2sep, dec_cesm_4sep, scale_by_aod="root"
    # ).run_loop()
    CheckRecreatedWaveforms(dec_cesm_2sep, dec_cesm_4sep).run_loop()


if __name__ == "__main__":
    main()
