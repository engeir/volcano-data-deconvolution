"""Plotting done for presenting as slides."""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pathlib

import contextlib
import datetime

import matplotlib.pyplot as plt
import numpy as np
import plastik
import volcano_base

import vdd.load
import vdd.utils
from vdd.load import PaddingMethod
from vdd.plotting.check_double_waveform import CheckRecreatedWaveforms
from vdd.plotting.cut_off import PlotCutOff
from vdd.plotting.deconv_ob16_cesm2 import PlotResponseFunctions
from vdd.plotting.reconstruction import PlotReconstruction

plt.style.use(
    [
        "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle",
        "vdd.extra",
        "vdd.jgr",
        {
            "savefig.format": "webp",
            "text.latex.preamble": r"\usepackage{amsmath,siunitx}",
        },
    ],
)

_SAVE_DIR = volcano_base.config.SAVE_PATH / "presentation"
if not _SAVE_DIR.exists():
    _SAVE_DIR.mkdir(parents=False)


def response_shape() -> None:
    """Run the main function."""
    data_cesm = vdd.load.CESMData
    dec_cesm = vdd.load.DeconvolveCESM
    padding = vdd.load.PaddingMethod.NOISE
    dec_cesm_int4 = dec_cesm(pad_before=padding, cesm=data_cesm(strength="tt-4sep"))
    dec_cesm_int2 = dec_cesm(pad_before=padding, cesm=data_cesm(strength="tt-2sep"))
    dec_cesm_med4 = dec_cesm(pad_before=padding, cesm=data_cesm(strength="medium-4sep"))
    dec_cesm_med2 = dec_cesm(pad_before=padding, cesm=data_cesm(strength="medium-2sep"))
    dec_cesm_e = dec_cesm(pad_before=padding, cesm=data_cesm(strength="size5000"))
    dec_cesm_s = dec_cesm(pad_before=padding, cesm=data_cesm(strength="strong"))
    dec_cesm_p = dec_cesm(pad_before=padding, cesm=data_cesm(strength="medium-plus"))
    dec_cesm_m = dec_cesm(pad_before=padding, cesm=data_cesm(strength="medium"))
    # Original
    dec_ob16 = vdd.load.DeconvolveOB16(data="h0", length=int(12 * 1000) + 1)
    dec_ob16.name = "OB16"
    all_decs = (
        dec_ob16,
        dec_cesm_m,
        dec_cesm_p,
        dec_cesm_s,
        dec_cesm_e,
        dec_cesm_med2,
        dec_cesm_med4,
        dec_cesm_int2,
        dec_cesm_int4,
    )

    plot = PlotResponseFunctions(*all_decs, norm=True)  # .run()
    fig, _ = vdd.utils.figure_multiple_rows_columns(1, 2)
    plot.norm_plot(fig.get_axes(), plot.plot_grayscale_highlight().get_axes())
    plot.plot_grayscale_highlight(fig, "temp-so2-gs")
    fig.savefig(_SAVE_DIR / "response_pulse.webp")
    plt.show()


def waveform() -> None:
    """Run the main script."""
    data_cesm = vdd.load.CESMData
    dec_cesm = vdd.load.DeconvolveCESM
    # CESM2
    padding = vdd.load.PaddingMethod.NOISE
    dec_cesm_p4 = dec_cesm(pad_before=padding, cesm=data_cesm(strength="tt-4sep"))
    dec_cesm_p2 = dec_cesm(pad_before=padding, cesm=data_cesm(strength="tt-2sep"))
    dec_cesm_m4 = dec_cesm(pad_before=padding, cesm=data_cesm(strength="medium-4sep"))
    dec_cesm_m2 = dec_cesm(pad_before=padding, cesm=data_cesm(strength="medium-2sep"))
    dec_cesm_p = dec_cesm(pad_before=padding, cesm=data_cesm(strength="medium-plus"))
    dec_cesm_m = dec_cesm(pad_before=padding, cesm=data_cesm(strength="medium"))

    name: pathlib.Path
    fig, axs, name = CheckRecreatedWaveforms(
        dec_cesm_p2, dec_cesm_p4, single_waveform=dec_cesm_p, keys={"temp": (0, 1)}
    ).run_loop(save_path=_SAVE_DIR, return_=True)
    fig.savefig(name.with_suffix(".webp"))
    fig, axs, name = CheckRecreatedWaveforms(
        dec_cesm_m2, dec_cesm_m4, single_waveform=dec_cesm_m, keys={"temp": (0, 1)}
    ).run_loop(save_path=_SAVE_DIR, return_=True)
    fig.savefig(name.with_suffix(".webp"))
    plt.show()


def _cut_off() -> None:
    data_cesm = vdd.load.CESMData
    dec_cesm = vdd.load.DeconvolveCESM
    # CESM2
    padding = vdd.load.PaddingMethod.NOISE
    dec_cesm_m = dec_cesm(pad_before=padding, cesm=data_cesm(strength="medium"))
    dec_ob16_month = vdd.load.DeconvolveOB16(data="h0")
    dec_ob16_month.name = "OB16 month"
    # OB16
    co_ob16_temp_so2 = vdd.load.CutOff(dec_ob16_month, ("temp", "so2"))
    co_ob16_temp_so2.cut_off(int(7 * 12))
    # CESM2 medium
    co_cesm_m_temp_so2 = vdd.load.CutOff(dec_cesm_m, ("temp", "so2"))
    co_cesm_m_temp_so2.cut_off(int(7 * 12))
    pco = PlotCutOff(
        co_ob16_temp_so2,
        co_cesm_m_temp_so2,
    )
    pco.call_cut_offs("generate_ensembles", 100)
    pco.plot(save_path=_SAVE_DIR)


def _ob16_timeseries() -> None:
    data_cesm = vdd.load.CESMData
    dec_cesm = vdd.load.DeconvolveCESM
    # CESM2
    padding = PaddingMethod.NOISE
    dec_cesm_m = dec_cesm(
        pad_before=padding, cesm=data_cesm(dims=["lat", "lon"], strength="medium")
    )
    # OB16
    dec_ob16_month = vdd.load.DeconvolveOB16(data="h0")
    dec_ob16_month.name = "OB16 month"
    # all_decs = (dec_cesm_m, dec_ob16_month)
    ob16 = dec_ob16_month.data
    rec_ob16_ = dec_ob16_month.dump_reconstructor()
    rec_small_ = dec_cesm_m.dump_reconstructor()
    zero_like = 1e-4
    valid_until = 7
    rec_ob16_.response_temp_so2[
        np.intersect1d(
            np.argwhere(rec_ob16_.response_temp_so2 < zero_like).flatten(),
            np.argwhere(rec_ob16_.tau > valid_until).flatten(),
        ).flatten()[0] :
    ] = 0
    with contextlib.suppress(Exception):
        # In case we never reach zero, for example when running with the large eruption
        # simulation output
        rec_small_.response_temp_so2[
            np.intersect1d(
                np.argwhere(rec_small_.response_temp_so2 < zero_like).flatten(),
                np.argwhere(rec_small_.tau > valid_until).flatten(),
            ).flatten()[0] :
        ] = 0
    rec_ob16 = PlotReconstruction(ob16, rec_ob16_)
    # rec_small = PlotReconstruction(ob16, rec_small_)

    figrec, axrec = plastik.figure_grid(1, 1, using={"share_axes": "x"})
    rec_ob16.plot_reconstruction_temp(axrec[0])
    # rec_small.plot_reconstruction_temp(axrec[1])
    xlim = (
        vdd.utils.d2n(datetime.datetime(850, 1, 1, 0, 0, tzinfo=datetime.UTC)),
        vdd.utils.d2n(datetime.datetime(1850, 1, 1, 0, 0, tzinfo=datetime.UTC)),
    )
    [ax.set_xlim(xlim) for ax in axrec]
    figrec.savefig(_SAVE_DIR / "compare-historical-size-temp-reconstructed")
    plt.show()


if __name__ == "__main__":
    # response_shape()
    # waveform()
    # cut_off()
    # ob16_timeseries()
    pass
