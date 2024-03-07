"""Plot the deconvolution comparison between OB16 and CESM2."""

import matplotlib.pyplot as plt
import volcano_base

import vdd.load
import vdd.utils
from vdd.deconvolve_methods import alternative_deconv

plt.rc("text.latex", preamble=r"\usepackage{amsmath}")
plt.style.use(
    "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle"
)

dec_cesm_s = vdd.load.DeconvolveCESM(pad_before=True)
dec_cesm_p = vdd.load.DeconvolveCESM(
    pad_before=True, cesm=vdd.load.CESMData("medium-plus")
)
dec_cesm_m = vdd.load.DeconvolveCESM(pad_before=True, cesm=vdd.load.CESMData("medium"))
dec_cesm_e = vdd.load.DeconvolveCESM(
    pad_before=True, cesm=vdd.load.CESMData("size5000")
)
# Original
dec_ob16 = vdd.load.DeconvolveOB16(data="h1")
# dec_ob16_month = vdd.load.DeconvolveOB16(data="h0")
# ob16_day = volcano_base.load.OttoBliesner(freq="h1", progress=True)
ob16_month = volcano_base.load.OttoBliesner(freq="h0", progress=True)
# dec_ob16 = vdd.load.DeconvolveOB16(data=ob16_day)
# dec_ob16.change_deconvolution_method(alternative_deconv)
dec_ob16_month = vdd.load.DeconvolveOB16(data=ob16_month)
dec_ob16_month.change_deconvolution_method(alternative_deconv)
rf_xlim = (-2, 10)
temp_xlim = (-2, 20)


def compare_ob16_with_cesm(norm=False) -> None:
    """Compare OB16 with CESM2."""
    ## Radiative forcing to SO2 response -----------------------------------------------
    rf = plt.figure()
    rf_a = rf.gca()
    # OB16 (CESM1)
    tau = volcano_base.manipulate.dt2float(dec_ob16.tau)
    rf_shape_ob16 = (
        vdd.utils.normalise(dec_ob16.response_rf_so2)
        if norm
        else dec_ob16.response_rf_so2
    )
    rf_a.plot(tau, rf_shape_ob16, label="OB16")
    # OB16 month (CESM1)
    tau = volcano_base.manipulate.dt2float(dec_ob16_month.tau)
    rf_shape_ob16 = (
        vdd.utils.normalise(dec_ob16_month.response_rf_so2)
        if norm
        else dec_ob16_month.response_rf_so2
    )
    rf_a.plot(tau, rf_shape_ob16, label="OB16 month")
    # CESM2
    rf_shape_cesm = (
        vdd.utils.normalise(dec_cesm_m.response_rf_so2)
        if norm
        else dec_cesm_m.response_rf_so2
    )
    rf_a.plot(dec_cesm_m.tau, rf_shape_cesm, label="CESM2 Small")
    rf_shape_cesm = (
        vdd.utils.normalise(dec_cesm_p.response_rf_so2)
        if norm
        else dec_cesm_p.response_rf_so2
    )
    rf_a.plot(dec_cesm_p.tau, rf_shape_cesm, label="CESM2 Intermediate")
    rf_shape_cesm = (
        vdd.utils.normalise(dec_cesm_s.response_rf_so2)
        if norm
        else dec_cesm_s.response_rf_so2
    )
    rf_a.plot(dec_cesm_s.tau, rf_shape_cesm, label="CESM2 Strong")
    rf_shape_cesm = (
        vdd.utils.normalise(dec_cesm_e.response_rf_so2)
        if norm
        else dec_cesm_e.response_rf_so2
    )
    rf_a.plot(dec_cesm_e.tau, rf_shape_cesm, label="CESM2 Extra Strong")
    rf_a.set_xlim(rf_xlim)
    rf_a.set_xlabel("Time lag ($\\tau$) [yr]")
    rf_a.set_ylabel(
        f"{"Normalised " if norm else ""}RF to {"\n" if norm else ""}SO2 response [1]"
    )
    rf_a.legend()
    rf.savefig(f"rf-so2{"-norm" if norm else ""}.png")
    ## Temperature to SO2 response -----------------------------------------------------
    temp = plt.figure()
    temp_a = temp.gca()
    # OB16 (CESM1)
    tau = volcano_base.manipulate.dt2float(dec_ob16.tau)
    temp_shape_ob16 = (
        vdd.utils.normalise(dec_ob16.response_temp_so2)
        if norm
        else dec_ob16.response_temp_so2
    )
    temp_a.plot(tau, temp_shape_ob16, label="OB16")
    # OB16 month (CESM1)
    tau = volcano_base.manipulate.dt2float(dec_ob16_month.tau)
    temp_shape_ob16 = (
        vdd.utils.normalise(dec_ob16_month.response_temp_so2)
        if norm
        else dec_ob16_month.response_temp_so2
    )
    temp_a.plot(tau, temp_shape_ob16, label="OB16 month")
    # CESM2
    temp_shape_cesm = (
        vdd.utils.normalise(dec_cesm_m.response_temp_so2)
        if norm
        else dec_cesm_m.response_temp_so2
    )
    temp_a.plot(dec_cesm_m.tau, temp_shape_cesm, label="CESM2 Small")
    temp_shape_cesm = (
        vdd.utils.normalise(dec_cesm_p.response_temp_so2)
        if norm
        else dec_cesm_p.response_temp_so2
    )
    temp_a.plot(dec_cesm_p.tau, temp_shape_cesm, label="CESM2 Intermediate")
    temp_shape_cesm = (
        vdd.utils.normalise(dec_cesm_s.response_temp_so2)
        if norm
        else dec_cesm_s.response_temp_so2
    )
    temp_a.plot(dec_cesm_s.tau, temp_shape_cesm, label="CESM2 Strong")
    temp_shape_cesm = (
        vdd.utils.normalise(dec_cesm_e.response_temp_so2)
        if norm
        else dec_cesm_e.response_temp_so2
    )
    temp_a.plot(dec_cesm_e.tau, temp_shape_cesm, label="CESM2 Extra Strong")
    temp_a.set_xlim(temp_xlim)
    temp_a.set_xlabel("Time lag ($\\tau$) [yr]")
    temp_a.set_ylabel(
        f"{"Normalised t" if norm else "T"}emperature to {"\n" if norm else ""}SO2 response [1]"
    )
    temp_a.legend()
    temp.savefig(f"temp-so2{"-norm" if norm else ""}.png")
    ## Temperature to RF response ------------------------------------------------------
    temp_rf = plt.figure()
    temp_rf_a = temp_rf.gca()
    # OB16 (CESM1)
    tau = volcano_base.manipulate.dt2float(dec_ob16.tau)
    temp_rf_shape_ob16 = (
        vdd.utils.normalise(dec_ob16.response_temp_rf)
        if norm
        else dec_ob16.response_temp_rf
    )
    temp_rf_a.plot(tau, temp_rf_shape_ob16, label="OB16")
    # OB16 month (CESM1)
    tau = volcano_base.manipulate.dt2float(dec_ob16_month.tau)
    temp_rf_shape_ob16 = (
        vdd.utils.normalise(dec_ob16_month.response_temp_rf)
        if norm
        else dec_ob16_month.response_temp_rf
    )
    temp_rf_a.plot(tau, temp_rf_shape_ob16, label="OB16 month")
    # CESM2
    temp_rf_shape_cesm = (
        vdd.utils.normalise(dec_cesm_m.response_temp_rf)
        if norm
        else dec_cesm_m.response_temp_rf
    )
    temp_rf_a.plot(dec_cesm_m.tau, temp_rf_shape_cesm, label="CESM2 Small")
    temp_rf_shape_cesm = (
        vdd.utils.normalise(dec_cesm_p.response_temp_rf)
        if norm
        else dec_cesm_p.response_temp_rf
    )
    temp_rf_a.plot(dec_cesm_p.tau, temp_rf_shape_cesm, label="CESM2 Intermediate")
    temp_rf_shape_cesm = (
        vdd.utils.normalise(dec_cesm_s.response_temp_rf)
        if norm
        else dec_cesm_s.response_temp_rf
    )
    temp_rf_a.plot(dec_cesm_s.tau, temp_rf_shape_cesm, label="CESM2 Strong")
    temp_rf_shape_cesm = (
        vdd.utils.normalise(dec_cesm_e.response_temp_rf)
        if norm
        else dec_cesm_e.response_temp_rf
    )
    temp_rf_a.plot(dec_cesm_e.tau, temp_rf_shape_cesm, label="CESM2 Extra Strong")
    temp_rf_a.set_xlim(rf_xlim)
    temp_rf_a.set_xlabel("Time lag ($\\tau$) [yr]")
    temp_rf_a.set_ylabel(
        f"{"Normalised t" if norm else "T"}emperature to {"\n" if norm else ""}RF response [1]"
    )
    temp_rf_a.legend()
    temp_rf.savefig(f"temp-rf{"-norm" if norm else ""}.png")
    plt.show()


if __name__ == "__main__":
    compare_ob16_with_cesm()
    compare_ob16_with_cesm(True)
