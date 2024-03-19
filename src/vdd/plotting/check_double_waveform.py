"""Check how well response functions are estimated from double waveforms.

We also check how well we are able to recreate the double waveform time series from the
response functions of other simulations.
"""

import matplotlib.pyplot as plt
import numpy as np
import volcano_base
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
    rf_f = plt.figure()
    rf_a = rf_f.gca()
    temp_f = plt.figure()
    temp_a = temp_f.gca()
    for dec in decs:
        # Check how well we are able to recreate the double waveform time series from the
        # response functions of other simulations.
        tau = dec.tau
        signal = dec.response_temp_so2
        kernel = dec.response_rf_so2
        name = "2sep" if "2sep" in dec.name else "4sep"
        rf_a.plot(tau, kernel, label=name)
        temp_a.plot(tau, signal, label=name)
    rf_a.set_xlabel("Time lag ($\\tau$) [yr]")
    rf_a.set_ylabel("Radiative forcing response to SO2")
    temp_a.set_xlabel("Time lag ($\\tau$) [yr]")
    temp_a.set_ylabel("Temperature response to SO2")
    rf_a.set_xlim((-2, 15))
    temp_a.set_xlim((-2, 15))
    rf_a.legend()
    temp_a.legend()
    rf_f.savefig(_SAVE_DIR / "responses_rf.png")
    temp_f.savefig(_SAVE_DIR / "responses_temp.png")
    plt.show()


def check_recreated_waveforms(*decs: vdd.load.DeconvolveCESM) -> None:
    """Check how well we are able to recreate the double waveform time series."""
    rf_f = plt.figure()
    rf_a = rf_f.gca()
    temp_f = plt.figure()
    temp_a = temp_f.gca()
    response_rf = dec_cesm_s.response_rf_so2
    response_temp = dec_cesm_s.response_temp_so2
    for i, dec in enumerate(decs):
        tau = dec.tau
        so2 = dec.so2
        temp = dec.temp
        rf = dec.rf
        diff_len = len(response_rf) - len(so2)
        r_rf = response_rf[diff_len // 2 : -diff_len // 2]
        r_rf = r_rf / np.max(r_rf) * np.max(dec.response_rf_so2)
        r_temp = response_temp[diff_len // 2 : -diff_len // 2]
        r_temp = r_temp / np.max(r_temp) * np.max(dec.response_temp_so2)
        recreated_rf = np.convolve(r_rf, so2, "same")
        recreated_temp = np.convolve(r_temp, so2, "same")
        name = "2sep" if "2sep" in dec.name else "4sep"
        pe = [
            patheffects.Stroke(linewidth=1, foreground=_COLORS[i]),
            patheffects.Normal(),
        ]
        rf_a.plot(tau, rf, c="k", lw=0.5, path_effects=pe, label=f"{name} old")
        temp_a.plot(tau, temp, c="k", lw=0.5, path_effects=pe, label=f"{name} old")
        rf_a.plot(tau, recreated_rf, c=_COLORS[i], label=f"{name} new")
        temp_a.plot(tau, recreated_temp, c=_COLORS[i], label=f"{name} new")
    rf_a.set_xlabel("Time [yr]")
    rf_a.set_ylabel("Radiative Forcing [W/m$^2$]")
    temp_a.set_xlabel("Time [yr]")
    temp_a.set_ylabel("Temperature [K]")
    rf_a.set_xlim((-2, 15))
    temp_a.set_xlim((-2, 15))
    rf_a.legend()
    temp_a.legend()
    rf_f.savefig(_SAVE_DIR / "recreated_waveforms_rf.png")
    temp_f.savefig(_SAVE_DIR / "recreated_waveforms_temp.png")
    plt.show()


def main() -> None:
    """Run the main script."""
    check_waveform_responses(dec_cesm_2sep, dec_cesm_4sep)
    check_recreated_waveforms(dec_cesm_2sep, dec_cesm_4sep)


if __name__ == "__main__":
    main()
