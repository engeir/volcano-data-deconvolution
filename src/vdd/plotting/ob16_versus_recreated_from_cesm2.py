"""Compare original OB16 to the convolution from CESM2."""

import matplotlib.pyplot as plt
import numpy as np
import volcano_base

import vdd.load

plt.style.use(
    [
        "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle",
        "vdd.extra",
        "vdd.jgr",
    ],
)

ob16_month = volcano_base.load.OttoBliesner(freq="h0", progress=True)
dec_ob16_month = vdd.load.DeconvolveOB16(data=ob16_month, length=12000)
# dec_ob16_month = vdd.load.DeconvolveOB16(data=ob16_month)
# dec_ob16_month.change_deconvolution_method(alternative_deconv)
padding = vdd.load.PaddingMethod.ZEROS
cesm_e = vdd.load.DeconvolveCESM(
    pad_before=padding,
    cesm=vdd.load.CESMData(strength="size5000"),
)
cesm_s = vdd.load.DeconvolveCESM(
    pad_before=padding,
    cesm=vdd.load.CESMData(strength="strong"),
)
cesm_p = vdd.load.DeconvolveCESM(
    pad_before=padding,
    cesm=vdd.load.CESMData(strength="medium-plus"),
)
cesm_m = vdd.load.DeconvolveCESM(
    pad_before=padding,
    cesm=vdd.load.CESMData(strength="medium"),
)
ob16_time = dec_ob16_month.temp.time.data
ob16_so2 = dec_ob16_month.so2
ob16_rf = dec_ob16_month.rf
ob16_rt2_amplitude = dec_ob16_month.response_temp_so2.max()
ob16_rtr_amplitude = dec_ob16_month.response_temp_rf.max()
abso = plt.figure()
abso_a = abso.gca()
norm = plt.figure()
norm_a = norm.gca()
nor2 = plt.figure()
nor2_a = nor2.gca()
abso_a.plot(
    ob16_time,
    dec_ob16_month.temp.data,
    label="OB16",
)
norm_a.plot(
    ob16_time,
    vdd.utils.normalise(dec_ob16_month.temp.data),
    label="OB16",
)
nor2_a.plot(ob16_time, dec_ob16_month.temp.data, label="OB16")
for cesm, label in zip(
    [cesm_e, cesm_s, cesm_p, cesm_m],
    ["Extreme", "Strong", "Intermediate", "Small"],
    strict=True,
):
    # fmt: off
    conv_temp_so2 = np.convolve(ob16_so2.data, dec_ob16_month.response_temp_so2, "same")
    conv_temp_rf = np.convolve(ob16_rf.data, dec_ob16_month.response_temp_rf, "same")
    conv_norm_temp_so2 = np.convolve(ob16_so2.data, cesm.response_temp_so2 / cesm.response_temp_so2.max() * ob16_rt2_amplitude, "same")
    conv_norm_temp_rf = np.convolve(ob16_rf.data, cesm.response_temp_rf / cesm.response_temp_rf.max() * ob16_rtr_amplitude, "same")
    abso_a.plot(ob16_time, conv_temp_so2, label=f"Reconstructed from SO2 ({label})")
    abso_a.plot(ob16_time, conv_temp_rf, label=f"Reconstructed from RF ({label})")
    norm_a.plot(ob16_time, vdd.utils.normalise(conv_temp_so2), label=f"Reconstructed from SO2 ({label})")
    norm_a.plot(ob16_time, vdd.utils.normalise(conv_temp_rf), label=f"Reconstructed from RF ({label})")
    nor2_a.plot(ob16_time, conv_norm_temp_so2, label=f"Reconstructed from SO2 ({label})")
    nor2_a.plot(ob16_time, conv_norm_temp_rf, label=f"Reconstructed from RF ({label})")
    # fmt: on
abso_a.legend()
norm_a.legend()
nor2_a.legend()
plt.show()
