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
    aod_f.savefig(_SAVE_DIR / "responses_aod.png")
    rf_f.savefig(_SAVE_DIR / "responses_rf.png")
    temp_f.savefig(_SAVE_DIR / "responses_temp.png")
    plt.show()


def check_recreated_waveforms(*decs: vdd.load.DeconvolveCESM) -> None:
    """Check how well we are able to recreate the double waveform time series."""
    figs = {"aod": plt.figure(), "rf": plt.figure(), "temp": plt.figure()}
    axs = {k: v.gca() for k, v in figs.items()}
    response_rf = dec_cesm_s.response_rf_so2
    for i, dec in enumerate(decs):
        tau = dec.tau
        so2 = dec.so2
        diff_len = len(response_rf) - len(so2)
        name = "2sep" if "2sep" in dec.name else "4sep"
        pe = [
            patheffects.Stroke(linewidth=1, foreground=_COLORS[i]),
            patheffects.Normal(),
        ]
        for attr in ["aod", "rf", "temp"]:
            arr = getattr(dec, attr)
            dec_rec = getattr(dec, f"response_{attr}_so2")
            r_arr = getattr(dec_cesm_s, f"response_{attr}_so2")[
                diff_len // 2 : -diff_len // 2
            ]
            r_arr = r_arr / np.max(r_arr) * np.max(dec_rec)
            rec_same = np.convolve(dec_rec, so2, "same")
            rec_new = np.convolve(r_arr, so2, "same")
            plot = axs[attr].plot
            plot(tau, arr, c="k", lw=0.5, label=f"{name} original")
            kwargs = {"ls": "--", "c": "k", "lw": 0.5, "path_effects": pe}
            plot(tau, rec_same, label=f"{name} reconstruct self", **kwargs)  # type: ignore
            plot(tau, rec_new, c=_COLORS[i], label=f"{name} reconstruc other")
    [axs[k].set_xlabel("Time [yr]") for k in axs]
    [axs[k].set_xlim((-2, 25)) for k in axs]
    [axs[k].legend() for k in axs]
    axs["aod"].set_ylabel("Radiative Forcing [W/m$^2$]")
    axs["rf"].set_ylabel("Radiative Forcing [W/m$^2$]")
    axs["temp"].set_ylabel("Temperature [K]")
    figs["aod"].savefig(_SAVE_DIR / "recreated_waveforms_aod.png")
    figs["rf"].savefig(_SAVE_DIR / "recreated_waveforms_rf.png")
    figs["temp"].savefig(_SAVE_DIR / "recreated_waveforms_temp.png")
    plt.show()


def main() -> None:
    """Run the main script."""
    # check_waveform_responses(dec_cesm_2sep, dec_cesm_4sep)
    check_recreated_waveforms(dec_cesm_2sep, dec_cesm_4sep)


if __name__ == "__main__":
    main()