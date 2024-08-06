"""Check how well response functions are estimated from double waveforms.

We also check how well we are able to recreate the double waveform time series from the
response functions of other simulations.
"""

from typing import Literal, Self

import matplotlib.pyplot as plt
import numpy as np
import volcano_base
import xarray as xr
from numpy.typing import NDArray

import vdd.load
import vdd.utils

plt.style.use(
    [
        "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle",
        "vdd.extra",
        "vdd.jgr",
    ],
)
_COLORS = plt.rcParams["axes.prop_cycle"].by_key()["color"]

_SAVE_DIR = volcano_base.config.SAVE_PATH / "waveform"
if not _SAVE_DIR.exists():
    _SAVE_DIR.mkdir(parents=False)

DataCESM = vdd.load.CESMData
DecCESM = vdd.load.DeconvolveCESM
# CESM2
padding = vdd.load.PaddingMethod.NOISE
dec_cesm_p4 = DecCESM(pad_before=padding, cesm=DataCESM(strength="tt-4sep"))
dec_cesm_p2 = DecCESM(pad_before=padding, cesm=DataCESM(strength="tt-2sep"))
dec_cesm_m4 = DecCESM(pad_before=padding, cesm=DataCESM(strength="medium-4sep"))
dec_cesm_m2 = DecCESM(pad_before=padding, cesm=DataCESM(strength="medium-2sep"))
dec_cesm_e = DecCESM(pad_before=padding, cesm=DataCESM(strength="size5000"))
dec_cesm_s = DecCESM(pad_before=padding, cesm=DataCESM(strength="strong"))
dec_cesm_p = DecCESM(pad_before=padding, cesm=DataCESM(strength="medium-plus"))
dec_cesm_m = DecCESM(pad_before=padding, cesm=DataCESM(strength="medium"))


def check_waveform_responses(*decs: vdd.load.DeconvolveCESM) -> None:
    """Check how well response functions are estimated from double waveforms."""
    aod_f = plt.figure()
    aod_a = aod_f.gca()
    rf_f = plt.figure()
    rf_a = rf_f.gca()
    temp_f = plt.figure()
    temp_a = temp_f.gca()
    for dec in decs:
        # Check how well we are able to recreate the double waveform time series from
        # the response functions of other simulations.
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


def curve_fit_aod(
    aod: NDArray[np.float64], alpha: float, beta: float
) -> NDArray[np.float64]:
    """Fit values using AOD as input."""
    return alpha * np.log(beta * aod + 1)


class CheckRecreatedWaveforms:
    """Check how well we are able to recreate the double waveform time series."""

    def __init__(
        self: Self,
        /,
        *decs: vdd.load.DeconvolveCESM,
        single_waveform: vdd.load.DeconvolveCESM | None = None,
        scale_by_aod: Literal["log", "log-inside", "root"] | bool = False,
        keys: dict[str, tuple[int, int]] | None = None,
    ) -> None:
        self.single_waveform = (
            dec_cesm_p if single_waveform is None else single_waveform
        )
        self.decs = decs
        self.scale_by_aod = scale_by_aod
        self.keys = (
            {"aod": (0, 1), "rf": (2, 3), "temp": (4, 5)} if keys is None else keys
        )
        self.figs, self.axs = vdd.utils.figure_multiple_rows_columns(
            len(self.keys),
            2,
            share_axes="x",
            columns_first=True,
        )

    def run_loop(self: Self) -> None:
        """Run the main loop."""
        max_len = 3
        for _, dec in enumerate(self.decs):
            # so2_new = self._get_so2_new(dec)
            if len(self.keys) == max_len:
                dec._data.initialise_data()  # noqa: SLF001
            so2 = dec.so2
            so2_new = so2.copy()
            for attr in list(self.keys.keys()):
                self._run_attr_loop(dec, attr, so2_new)
        [ax.set_xlabel("Time after first eruption [yr]") for ax in self.axs]
        if padding:
            [ax.set_xlim((-1, 21)) for ax in self.axs]
        [ax.legend(framealpha=0.5, numpoints=1) for ax in self.axs]
        if len(self.keys) == max_len:
            [
                self.axs[i].set_ylabel("Aerosol optical \ndepth [1]")
                for i in self.keys["aod"]
            ]
            [
                self.axs[i].set_ylabel("Radiative \nforcing [W/m$^2$]")
                for i in self.keys["rf"]
            ]
        [self.axs[i].set_ylabel("Temperature \nanomaly [K]") for i in self.keys["temp"]]
        corrected = f"-aod-{self.scale_by_aod}-corrected" if self.scale_by_aod else ""
        base = "small" if "medium" in self.decs[0].name else "int"
        self.figs.savefig(_SAVE_DIR / f"responses_combined_{base}{corrected}")
        plt.show()

    def _get_so2_new(self: Self, dec: vdd.load.DeconvolveCESM) -> xr.DataArray:
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

    def _prepare_plot_data(
        self: Self,
        dec: vdd.load.DeconvolveCESM,
        attr: str,
        so2_new: xr.DataArray,
    ) -> tuple[
        tuple[xr.DataArray, np.ndarray, np.ndarray, int],
        tuple[str, str, str],
        tuple[int, int],
    ]:
        base = "SMALL" if "medium" in dec.name else "INT"
        base_long = "SMALL" if "medium" in dec.name else "INTERMEDIATE"
        sep = "2sep" if "2sep" in dec.name else "4sep"
        idx = 0 if sep == "2sep" else 1
        c_idx = 4 if sep == "2sep" else 5
        so2 = dec.so2
        diff_len = len(self.single_waveform.response_temp_so2) - len(so2)
        arr: xr.DataArray = getattr(dec, attr)
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
        print(dec.name, attr, "*" * 50)
        print(dec.name)
        print(getattr(dec, f"response_{attr}_so2_err")[-1])
        print(self.single_waveform.name)
        print(getattr(self.single_waveform, f"response_{attr}_so2_err")[-1])
        # tmpfig = plt.figure()
        # tmpax = tmpfig.gca()
        # tmpax.semilogy(getattr(dec, f"response_{attr}_so2_err"))
        # tmpax.semilogy(getattr(self.single_waveform, f"response_{attr}_so2_err"))
        rec_same = np.convolve(dec_resp, so2, "same")
        rec_new = np.convolve(resp_arr, so2_new, "same")
        sign = 1 if attr == "aod" else -1
        return (arr, rec_same, rec_new, sign), (base, base_long, sep), (idx, c_idx)

    def _run_attr_loop(
        self: Self,
        dec: vdd.load.DeconvolveCESM,
        attr: str,
        so2_new: xr.DataArray,
    ) -> None:
        arrs, names, idx = self._prepare_plot_data(dec, attr, so2_new)
        *_, sign = arrs
        plot = self.axs[self.keys[attr][idx[0]]].plot
        plot(
            dec.temp.time,
            sign * arrs[0],
            c="k",
            lw=0.5,
            label=f"${attr[0].upper()}_{{\\text{{{names[0]}-{names[2].upper()}}}}}$",
        )
        kwargs = {"ls": "--", "c": _COLORS[idx[1]], "lw": 1.0}
        varphi = f"\\varphi_{{{attr[0].upper()}}}"
        conv_so2 = f"\\ast S_{{\\text{{{names[0]}-{names[2].upper()}}}}}"
        plot(
            dec.temp.time,
            sign * arrs[2],
            "-.",
            lw=1.0,
            c=_COLORS[1],
            label=f"${varphi}^{{\\text{{{names[1]}}}}}{conv_so2}$",
        )
        plot(
            dec.temp.time,
            sign * arrs[1],
            label=f"${varphi}^{{\\text{{{names[0]}-{names[2].upper()}}}}}{conv_so2}$",
            **kwargs,
        )


def main() -> None:
    """Run the main script."""
    CheckRecreatedWaveforms(
        dec_cesm_p2,
        dec_cesm_p4,
        single_waveform=dec_cesm_p,
    ).run_loop()
    CheckRecreatedWaveforms(
        dec_cesm_m2,
        dec_cesm_m4,
        single_waveform=dec_cesm_m,
    ).run_loop()


if __name__ == "__main__":
    main()
