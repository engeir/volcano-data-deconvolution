"""Script showing how bad the results from the daily resolved deconvolution is."""

import fppanalysis
import matplotlib.pyplot as plt

import vdd.load

dec_ob16 = vdd.load.DeconvolveOB16(data="h1")
dec_ob16.name = "OB16"
rf = dec_ob16.rf
temp = dec_ob16.temp
# temp = temp[40:]
rrso2 = dec_ob16.response_rf_so2
rtso2 = dec_ob16.response_temp_so2
rtso2 = rtso2[30:]
rrso2[: len(rrso2) // 2] = 0
rtso2[: len(rtso2) // 2] = 0
srrso2 = fppanalysis.run_mean(rrso2, 50)
srtso2 = fppanalysis.run_mean(rtso2, 50)
dtr = fppanalysis.RL_gauss_deconvolve(srtso2, srrso2, 200)[0].flatten()
plt.figure()
plt.plot(rf)
plt.plot(temp)
plt.figure()
plt.plot(rrso2)
plt.plot(rtso2)
plt.figure()
plt.plot(srrso2)
plt.plot(srtso2)
plt.figure()
plt.plot(dtr)
plt.show()
