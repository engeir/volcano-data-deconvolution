"""Analysis of how cutting off the response functions affect the result."""

# TODO: Using the ideas from `volcano-deconv`, implement a simple CutOff class that will
# repeatedly cut off estimated response functions, recreate, add noise (for example
# phase randomised input signal), and from this create an ensemble of response
# estimates.
#
# 1. Implement the CutOff class.
#    - [ ] Owns a Deconvolve object.
#    - [ ] Uses only one response function and corresponding signal (RF or temperature).

import pathlib

import cftime
import matplotlib.pyplot as plt
import numpy as np
import plastik
import volcano_base

import vdd.load
import vdd.utils

plt.style.use([
    "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle",
    "vdd.extra",
])
_SAVE_DIR = volcano_base.config.SAVE_PATH / "cut_off"

if not _SAVE_DIR.exists():
    _SAVE_DIR.mkdir(parents=False)

DataCESM = vdd.load.CESMData
DecCESM = vdd.load.DeconvolveCESM
# CESM2
dec_cesm_4sep = DecCESM(pad_before=True, cesm=DataCESM(strength="tt-4sep"))
dec_cesm_2sep = DecCESM(pad_before=True, cesm=DataCESM(strength="tt-2sep"))
dec_cesm_s = DecCESM(pad_before=True, cesm=DataCESM(strength="strong"))
dec_ob16_month = vdd.load.DeconvolveOB16(data="h0")
dec_ob16_month.name = "OB16 month"
all_decs = (dec_cesm_4sep, dec_cesm_2sep, dec_cesm_s, dec_ob16_month)


class PlotCutOff:
    """Plot the results of the CutOff class for any deconvolution object."""

    def __init__(self, *cut_offs: vdd.load.CutOff) -> None:
        self.cut_offs = cut_offs

    def call_cut_offs(self, method: str, *args, **kwargs) -> None:
        """Call a method on all CutOff objects."""
        match method, args, kwargs:
            case method, (arg1,), _:
                for co in self.cut_offs:
                    getattr(co, method)(arg1)

    def plot(self, remove_grid_parts: bool = True) -> None:
        """Plot the results of the CutOff class.

        Parameters
        ----------
        remove_grid_parts : bool
            If True, the individual images will be removed after combining them.

        Raises
        ------
        ValueError
            If no cuts have been made in the CutOff objects.
        """
        for co in self.cut_offs:
            if not co.cuts:
                raise ValueError(
                    "No cuts have been made. Run `call_cut_offs('cut_off', ...)`."
                )
            files = self._plot_single(co)
            f_l = list(files)
            f_l.sort(key=lambda x: x.name[-7:-4])
            f1: pathlib.Path = files[0]
            vdd.utils.combine(*f_l).in_grid(2, len(files) // 2).save(
                f1.parent / f"{vdd.utils.name_swap(f1.name[:-7])}combined.jpg"
            )
            if remove_grid_parts:
                for f in files:
                    f.unlink()

    @staticmethod
    def _plot_single(co: vdd.load.CutOff) -> tuple[pathlib.Path, ...]:
        """Plot the results of the CutOff class."""
        files: tuple[pathlib.Path, ...] = ()
        for k, v in co.cuts.items():
            resp_f = plt.figure()
            resp_a = resp_f.gca()
            temp_f = plt.figure()
            temp_a = temp_f.gca()
            resp_a.axvline(int(k) / 12, c="k", ls="--")
            resp_a.plot(co.dec.tau, co.response, label="Original")
            resp_a.plot(v.tau, v.response, label=f"Cut {int(k) // 12}")
            temp_a.plot(co.output.time, co.dec.temp_control, label="Control")
            temp_a.plot(co.output.time, co.output, label="Original")
            temp_a.plot(v.time, v.temp_rec, label=f"Cut {int(k) // 12}")
            percentiles = []
            for i, j in co.ensembles[k].items():
                if "response" in str(i):
                    percentiles.append(j)
                    # label = j.attrs["label"]
                    # resp_a.plot(j.tau, j, label=f"_{label}", c="grey", alpha=0.3)
            percs_np = np.asarray(percentiles)
            resp_a = plastik.percentiles(
                co.dec.tau, percs_np, plot_median=False, ax=resp_a
            )
            resp_a.set_xlabel("Time lag ($\\tau$) [yr]")
            resp_a.set_ylabel("Response [1]")
            resp_a.set_xlim((-2, 20))
            ymax = co.response.max()
            resp_a.set_ylim((ymax * (-0.05), ymax * 1.05))
            resp_a.legend(loc="upper right", framealpha=0.9)
            temp_a.set_xlabel("Time lag ($\\tau$) [yr]")
            temp_a.set_ylabel("Temperature anomaly [K]")
            temp_a.legend(loc="upper right", framealpha=0.9)
            match v.time.data[0]:
                case cftime._cftime.DatetimeNoLeap():
                    temp_a.set_xlim((-790 * 365, -650 * 365))
                case _:
                    pass
            num = "0" * (3 - len(k)) + k
            name = vdd.utils.name_swap(vdd.utils.clean_filename(co.dec.name))
            ts = vdd.utils.name_swap(
                vdd.utils.clean_filename("-".join(co.ts_specifier))
            )
            resp_name = _SAVE_DIR / f"{name}_resp_{ts}_{num}.jpg"
            resp_f.savefig(resp_name)
            temp_name = _SAVE_DIR / f"{name}_temp_{ts}_{num}.jpg"
            temp_f.savefig(temp_name)
            plt.close("all")
            files += (resp_name, temp_name)
        return files


def _use_cut_off() -> None:
    # OB16
    co_ob16_rf_so2 = vdd.load.CutOff(dec_ob16_month, ("rf", "so2"))
    co_ob16_temp_so2 = vdd.load.CutOff(dec_ob16_month, ("temp", "so2"))
    co_ob16_temp_rf = vdd.load.CutOff(dec_ob16_month, ("temp", "rf"))
    # CESM2 strong
    co_cesm_s_rf_so2 = vdd.load.CutOff(dec_cesm_s, ("rf", "so2"))
    co_cesm_s_temp_so2 = vdd.load.CutOff(dec_cesm_s, ("temp", "so2"))
    co_cesm_s_temp_rf = vdd.load.CutOff(dec_cesm_s, ("temp", "rf"))
    # CESM2 2sep
    co_2sep_rf_so2 = vdd.load.CutOff(dec_cesm_2sep, ("rf", "so2"))
    co_2sep_temp_so2 = vdd.load.CutOff(dec_cesm_2sep, ("temp", "so2"))
    co_2sep_temp_rf = vdd.load.CutOff(dec_cesm_2sep, ("temp", "rf"))
    # CESM2 4sep
    co_4sep_rf_so2 = vdd.load.CutOff(dec_cesm_4sep, ("rf", "so2"))
    co_4sep_temp_so2 = vdd.load.CutOff(dec_cesm_4sep, ("temp", "so2"))
    co_4sep_temp_rf = vdd.load.CutOff(dec_cesm_4sep, ("temp", "rf"))
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
    pco.call_cut_offs("cut_off", {12 * i for i in [6, 7, 8, 10]})
    pco.call_cut_offs("generate_ensembles", 50)
    pco.plot()


def _main() -> None:
    _use_cut_off()
    plt.show()


if __name__ == "__main__":
    _main()
