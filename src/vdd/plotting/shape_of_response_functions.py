"""Simply plot the response function from CESM and compare with analytical response functions."""

from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import gamma

import vdd.load


def _plot():
    cesm = vdd.load.CESMData(strength="strong")
    dec = vdd.load.DeconvolveCESM(cesm=cesm)
    rf = dec.response_rf_so2
    temp = dec.response_temp_so2
    rf /= rf.max()
    temp /= temp.max()
    time = dec.tau

    # Analytical response functions
    # We start by plotting a gamma function
    # Fit a gamma distribution to the data using curve_fit
    def gamma_pdf(x, sh, sc):
        pdf = gamma.pdf(x, a=sh, scale=sc)
        return pdf / pdf.max()

    rf_params, _ = curve_fit(gamma_pdf, time, rf)
    temp_params, _ = curve_fit(gamma_pdf, time, temp)
    plt.plot(
        time,
        gamma_pdf(time, *rf_params),
        label="Fitted RF Gamma Distribution",
    )
    plt.plot(
        time,
        gamma_pdf(time, *temp_params),
        label="Fitted T Gamma Distribution",
    )
    plt.plot(time, rf, label="RF")
    plt.plot(time, temp, label="T")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    _plot()
