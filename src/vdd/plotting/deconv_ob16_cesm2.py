"""Plot the deconvolution comparison between OB16 and CESM2."""

import matplotlib.pyplot as plt
import volcano_base

import vdd.load
import vdd.utils

plt.rc("text.latex", preamble=r"\usepackage{amsmath}")
plt.style.use(
    "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle"
)

dec_cesm = vdd.load.DeconvolveCESM(pad_before=True)
dec_ob16 = vdd.load.DeconvolveOB16(freq="h1")
dec_ob16_month = vdd.load.DeconvolveOB16(freq="h0")
plt.figure()
dec_ob16_month.so2.plot()
dec_ob16_month.rf.plot()
dec_ob16_month.temp.plot()
plt.show()
xlim = (-5, 20)


def compare_ob16_with_cesm(norm=False) -> None:
    """Compare OB16 with CESM2."""
    ## Radiative forcing ---------------------------------------------------------------
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
        vdd.utils.normalise(dec_cesm.response_rf_so2)
        if norm
        else dec_cesm.response_rf_so2
    )
    rf_a.plot(dec_cesm.tau, rf_shape_cesm, label="CESM2")
    rf_a.set_xlim(xlim)
    rf_a.set_xlabel("Time lag ($\\tau$) [yr]")
    rf_a.set_ylabel(f"{"Normalised r" if norm else "R"}esponse [1]")
    rf_a.legend()
    rf.savefig(f"rf{"-norm" if norm else ""}.png")
    ## Temperature ---------------------------------------------------------------------
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
        vdd.utils.normalise(dec_cesm.response_temp_so2)
        if norm
        else dec_cesm.response_temp_so2
    )
    temp_a.plot(dec_cesm.tau, temp_shape_cesm, label="CESM2")
    temp_a.set_xlim(xlim)
    temp_a.set_xlabel("Time lag ($\\tau$) [yr]")
    temp_a.set_ylabel(f"{"Normalised r" if norm else "R"}esponse [1]")
    temp_a.legend()
    temp.savefig(f"temp{"-norm" if norm else ""}.png")
    plt.show()


# if __name__ == "__main__":
compare_ob16_with_cesm()
compare_ob16_with_cesm(True)
