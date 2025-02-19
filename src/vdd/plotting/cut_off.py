"""Analysis of how cutting off the response functions affect the result.

Implementation of the basic CutOff class is done in `vdd.load`. This script implements
routines to plot the cut-off in a nice way.
"""

import datetime
import pathlib
from typing import Self

import cftime
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import plastik
import volcano_base

import vdd.load
import vdd.utils

plt.style.use(
    [
        "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle",
        "vdd.extra",
        "vdd.jgr",
    ],
)
_SAVE_DIR = volcano_base.config.SAVE_PATH / "cut_off"

if not _SAVE_DIR.exists():
    _SAVE_DIR.mkdir(parents=False)

DataCESM = vdd.load.CESMData
DecCESM = vdd.load.DeconvolveCESM
# CESM2
padding = vdd.load.PaddingMethod.NOISE
dec_cesm_m = DecCESM(pad_before=padding, cesm=DataCESM(strength="medium"))
dec_ob16_month = vdd.load.DeconvolveOB16(data="h0")
dec_ob16_month.name = "OB16 month"


class PlotCutOff:
    """Plot the results of the CutOff class for any deconvolution object."""

    def __init__(self: Self, *cut_offs: vdd.load.CutOff) -> None:
        self.cut_offs = cut_offs

    def call_cut_offs(self: Self, method: str, *args, **kwargs) -> None:  # noqa: ANN002,ANN003
        """Call a method on all CutOff objects."""
        match method, args, kwargs:
            case method, (arg1,), _:
                for co in self.cut_offs:
                    getattr(co, method)(arg1)

    def plot(self: Self, save_path: pathlib.Path | None = None) -> None:
        """Plot the results of the CutOff class.

        Parameters
        ----------
        save_path : pathlib.Path | None
            The directory to where the file should be saved.

        Raises
        ------
        ValueError
            If no cuts have been made in the CutOff objects.
        """
        save_path = save_path or _SAVE_DIR
        for co in self.cut_offs:
            fig, _ = vdd.utils.figure_multiple_rows_columns(
                len(co.cuts),
                2,
                share_axes="x",
            )
            if not co.cuts:
                msg = "No cuts have been made. Run `call_cut_offs('cut_off', ...)`."
                raise ValueError(
                    msg,
                )
            self._plot_single(fig, co)
            name = vdd.utils.name_swap(vdd.utils.clean_filename(co.dec.name))
            ts = vdd.utils.name_swap(
                vdd.utils.clean_filename("-".join(co.ts_specifier)),
            )
            fig.savefig(save_path / f"{name}_resp_{ts}_combined")

    @staticmethod
    def _plot_single(fig: mpl.figure.Figure, co: vdd.load.CutOff) -> mpl.figure.Figure:
        """Plot the results of the CutOff class."""
        colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        name = "OB16" if "OB16" in co.dec.name else "S26"
        for i, (k, v) in enumerate(co.cuts.items()):
            resp_a, temp_a = fig.axes[i * 2], fig.axes[i * 2 + 1]
            resp_a.axvline(int(k) / 12, c="k", ls="--")
            percentiles = [
                j for i_, j in co.ensembles[k].items() if "response" in str(i_)
            ]
            percs_np = np.asarray(percentiles)
            resp_a = plastik.percentiles(
                co.dec.tau,
                -1 * percs_np,
                plot_median=False,
                ax=resp_a,
                n=1,
                alpha=0.6,
                percentile_min=5,
                percentile_max=95,
            )
            resp_a.plot(
                co.dec.tau,
                -1 * co.response,
                c=colors[0],
                label=f"$\\varphi_T^{{\\text{{{name}}}}}$",
            )
            resp_a.plot(
                v.tau,
                -1 * v.response,
                c=colors[2],
                label=f"$\\varphi_{{T,{int(k) // 12}}}^{{\\text{{{name}}}}}$",
            )
            temp_a.plot(
                co.output.time,
                -1 * co.dec.temp_control,
                c=colors[1],
                label="$T_{\\text{CONTROL}}$",
            )
            temp_a.plot(
                co.output.time,
                -1 * co.output,
                c=colors[0],
                label=f"$T_{{\\text{{{name}}}}}$",
            )
            temp_a.plot(
                v.time,
                -1 * v.temp_rec,
                c=colors[2],
                label=f"$\\varphi_{{T,{int(k) // 12}}}^{{\\text{{{name}}}}}\\ast S_{{\\text{{{name}}}}}$",
            )
            resp_a.set_xlabel("Time lag [yr]")
            resp_a.set_ylabel("$\\varphi_T$")
            ymax = max(-1 * co.response, key=abs)
            match ymax:
                case val if val > 0:
                    resp_a.set_ylim((val * (-0.10), val * 1.10))
                case val if val < 0:
                    resp_a.set_ylim((val * (1.10), val * (-0.10)))
            resp_a.legend(loc="lower right", framealpha=0.9)
            xlab = (
                "Time [yr]"
                if "OB16" in co.dec.name
                else "Time after first eruption [yr]"
            )
            temp_a.set_xlabel(xlab)
            temp_a.set_ylabel("$T$ [K]")
            temp_a.legend(loc="lower right", framealpha=0.9)
            match v.time.data[0]:
                case cftime._cftime.DatetimeNoLeap():  # noqa: SLF001
                    temp_a.set_xlim(
                        (
                            vdd.utils.d2n(
                                datetime.datetime(1250, 1, 1, 0, 0, tzinfo=datetime.UTC)
                            ),
                            vdd.utils.d2n(
                                datetime.datetime(1310, 1, 1, 0, 0, tzinfo=datetime.UTC)
                            ),
                        ),
                    )
                    resp_a.set_xlim((-1, 20))
                case _:
                    temp_a.set_xlim((-1, 11))
                    resp_a.set_xlim((-1, 11))
        return fig


def _use_cut_off() -> None:
    # OB16
    # co_ob16_rf_so2 = vdd.load.CutOff(dec_ob16_month, ("rf", "so2"))
    co_ob16_temp_so2 = vdd.load.CutOff(dec_ob16_month, ("temp", "so2"))
    co_ob16_temp_so2.cut_off([12 * i for i in [5, 6, 7, 10]])
    # co_ob16_temp_rf = vdd.load.CutOff(dec_ob16_month, ("temp", "rf"))
    # CESM2 medium
    co_cesm_m_temp_so2 = vdd.load.CutOff(dec_cesm_m, ("temp", "so2"))
    co_cesm_m_temp_so2.cut_off([12 * i for i in [3, 4, 5, 7]])
    pco = PlotCutOff(
        # co_ob16_rf_so2,
        co_ob16_temp_so2,
        # co_ob16_temp_rf,
        co_cesm_m_temp_so2,
    )
    # pco.call_cut_offs("cut_off", [12 * i for i in [5, 6, 7, 10]])
    pco.call_cut_offs("generate_ensembles", 100)
    pco.plot()


def _main() -> None:
    _use_cut_off()
    plt.show()


if __name__ == "__main__":
    _main()
