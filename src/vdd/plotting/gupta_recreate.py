"""Implements a simple recreation of a plot from Gupta and Marshall (2018).

See https://doi.org/10.1175/JCLI-D-17-0703.s1 for the original paper.

Simply plot the response function from CESM and compare with analytical response
functions.
"""

import numpy as np
from matplotlib import pyplot as plt
from numpy.typing import NDArray
from scipy.optimize import curve_fit
from scipy.stats import gamma

import vdd.load

plt.style.use(
    [
        "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle",
        "vdd.extra",
    ],
)


def _plot() -> None:
    cesm = vdd.load.CESMData(strength="strong")
    dec = vdd.load.DeconvolveCESM(pad_before=vdd.load.PaddingMethod.ZEROS, cesm=cesm)
    rf = dec.response_rf_so2
    temp = dec.response_temp_so2
    rf /= rf.max()
    temp /= temp.max()
    time = dec.tau

    # Analytical response functions
    # We start by plotting a gamma function
    # Fit a gamma distribution to the data using curve_fit
    def gamma_pdf(
        x: NDArray[np.float64], sh: float, sc: float, pow_: float
    ) -> NDArray[np.float64]:
        pdf = gamma.pdf(x, a=sh, scale=sc)
        pow_neg = np.zeros_like(x[x <= 0])
        pow_pos = (1 + x[x > 0]) ** (-pow_)
        # pow_pos = np.exp(-x[x > 0] / pow) * np.exp(-x[x > 0] / pow2)
        out = pdf * np.concatenate([pow_neg, pow_pos])
        return out / out.max()

    rf_params, _ = curve_fit(gamma_pdf, time, rf)
    temp_params, _ = curve_fit(gamma_pdf, time, temp)
    plt.plot(time, gamma_pdf(time, *rf_params), label="Fitted RF Gamma")
    plt.plot(time, gamma_pdf(time, *temp_params), label="Fitted T Gamma")
    plt.plot(time, rf, label="RF")
    plt.plot(time, temp, label="T")
    plt.legend(framealpha=0.5)
    plt.figure()
    gap = 20
    x = [np.linspace(-gap * i, 1000 - gap * i, 1000) for i in range(40)]
    temp_sum = sum(gamma_pdf(xi, *temp_params) for xi in x)
    temp_single = gamma_pdf(x[0], *temp_params)
    plt.plot(x[0], temp_single, label="Fitted T Gamma")
    plt.plot(x[0], temp_sum, label="Fitted T Gamma")
    plt.show()


if __name__ == "__main__":
    _plot()
