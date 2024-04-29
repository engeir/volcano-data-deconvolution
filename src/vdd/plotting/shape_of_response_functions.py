"""Simply plot the response function from CESM and compare with analytical response functions."""

import warnings

import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import gamma

import vdd.load

warnings.warn(
    "This script is deprecated and may be removed in the future.",
    category=DeprecationWarning,
    stacklevel=1,
)

plt.style.use([
    "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle",
    "vdd.extra",
])


def _plot():
    cesm = vdd.load.CESMData(strength="strong")
    dec = vdd.load.DeconvolveCESM(pad_before=True, cesm=cesm)
    rf = dec.response_rf_so2
    temp = dec.response_temp_so2
    rf /= rf.max()
    temp /= temp.max()
    time = dec.tau

    # Analytical response functions
    # We start by plotting a gamma function
    # Fit a gamma distribution to the data using curve_fit
    def gamma_pdf(x, sh, sc, pow):
        pdf = gamma.pdf(x, a=sh, scale=sc)
        pow_neg = np.zeros_like(x[x <= 0])
        pow_pos = (1 + x[x > 0]) ** (-pow)
        # pow_pos = np.exp(-x[x > 0] / pow) * np.exp(-x[x > 0] / pow2)
        out = pdf * np.concatenate([pow_neg, pow_pos])
        return out / out.max()

    rf_params, _ = curve_fit(gamma_pdf, time, rf)
    temp_params, _ = curve_fit(gamma_pdf, time, temp)
    print(rf_params, temp_params)
    plt.plot(time, gamma_pdf(time, *rf_params), label="Fitted RF Gamma")
    plt.plot(time, gamma_pdf(time, *temp_params), label="Fitted T Gamma")
    plt.plot(time, rf, label="RF")
    plt.plot(time, temp, label="T")
    plt.legend(framealpha=0.5)
    plt.figure()
    gap = 15
    x1 = np.linspace(0, 100, 1000)
    x2 = np.linspace(-gap * 1, 100 - gap * 1, 1000)
    x3 = np.linspace(-gap * 2, 100 - gap * 2, 1000)
    x4 = np.linspace(-gap * 3, 100 - gap * 3, 1000)
    temp_sum = (
        gamma_pdf(x1, *temp_params)
        + gamma_pdf(x2, *temp_params)
        + gamma_pdf(x3, *temp_params)
        + gamma_pdf(x4, *temp_params)
    )
    temp_single = gamma_pdf(x1, *temp_params)
    plt.plot(x1, temp_single, label="Fitted T Gamma")
    plt.plot(x1, temp_sum, label="Fitted T Gamma")
    print(temp_sum[-1], temp_single[-1])
    print(temp_sum.mean(), temp_single.mean())
    plt.show()


if __name__ == "__main__":
    _plot()
