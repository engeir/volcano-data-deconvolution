"""Alternative deconvolution methods for the deconvolution of the forcing signal from the observed signal."""

import numpy as np
from scipy.optimize import minimize


def deconv_mean_2(signal, forcing):
    """Deconvolution with mean value and an unknown finite mean."""

    # This uses the last term in expression (24) with the mean value as an additional unknown.
    def costfun(x, signal, forcing):
        # Here, the last member of x is the unknown mean value.
        x_tmp = x[:-1]
        mu = x[-1]
        tmp = signal - np.convolve(forcing, x_tmp, "same") - mu * np.ones(x_tmp.size)
        return 0.5 * np.sum((tmp) ** 2)

    bounds = ((0.0, np.inf),) * signal.size
    bounds_mu = ((0.0, np.inf),)
    bounds += bounds_mu

    res = minimize(
        costfun, np.ones(signal.size + 1), args=(signal, forcing), bounds=bounds
    )
    return res
