"""Check how well response functions are estimated from double waveforms.

We also check how well we are able to recreate the double waveform time series from the
response functions of other simulations.
"""

import cosmoplots
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
    files = [_SAVE_DIR / f"responses_{k}.png" for k in ["aod", "rf", "temp"]]
    aod_f.savefig(files[0])
    rf_f.savefig(files[1])
    temp_f.savefig(files[2])
    plt.show()


def check_recreated_waveforms(
    *decs: vdd.load.DeconvolveCESM, scale_by_aod: bool = False
) -> None:
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
            dec_resp = getattr(dec, f"response_{attr}_so2")
            resp_arr = getattr(dec_cesm_s, f"response_{attr}_so2")[
                diff_len // 2 : -diff_len // 2
            ]
            # Scale r_arr
            scale_arr = np.max(dec_resp) / np.max(resp_arr)
            resp_arr = resp_arr * scale_arr
            # Scale so2 for new array
            so2_new = so2.copy()
            idx = so2_new > 0
            extra_scale_arr = (dec.aod.max() - dec.aod[idx]) / dec.aod.max()
            so2_new[idx] = (
                so2_new[idx] * extra_scale_arr**0.5 if scale_by_aod else so2_new[idx]
            )
            rec_same = np.convolve(dec_resp, so2, "same")
            rec_new = np.convolve(resp_arr, so2_new, "same")
            plot = axs[attr].plot
            plot(tau, arr, c="k", lw=0.5, label=f"{name} original")
            kwargs = {"ls": "--", "c": "k", "lw": 0.5, "path_effects": pe}
            plot(tau, rec_same, label=f"{name} reconstruct self", **kwargs)  # type: ignore
            plot(tau, rec_new, c=_COLORS[i], label=f"{name} reconstruc other")
    [axs[k].set_xlabel("Time [yr]") for k in axs]
    [axs[k].set_xlim((-2, 25)) for k in axs]
    [axs[k].legend() for k in axs]
    axs["aod"].set_ylabel("Aerosol optical depth [1]")
    axs["rf"].set_ylabel("Radiative Forcing [W/m$^2$]")
    axs["temp"].set_ylabel("Temperature [K]")
    corrected = "-aod-corrected" if scale_by_aod else ""
    files = [
        _SAVE_DIR / f"recreated_waveforms_{k}{corrected}.png"
        for k in ["aod", "rf", "temp"]
    ]
    figs["aod"].savefig(files[0])
    figs["rf"].savefig(files[1])
    figs["temp"].savefig(files[2])
    cosmoplots.combine(*files).in_grid(1, 3).using(fontsize=50).save(
        _SAVE_DIR / f"responses_combined{corrected}.png"
    )
    for f in files:
        f.unlink()
    plt.show()


def main() -> None:
    """Run the main script."""
    # check_waveform_responses(dec_cesm_2sep, dec_cesm_4sep)
    check_recreated_waveforms(dec_cesm_2sep, dec_cesm_4sep, scale_by_aod=True)
    check_recreated_waveforms(dec_cesm_2sep, dec_cesm_4sep)


if __name__ == "__main__":
    main()
