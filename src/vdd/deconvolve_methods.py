"""Alternative deconvolution methods for a forcing signal and a observed signal."""

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import OptimizeResult, minimize


def deconv_mean_2(
    signal: NDArray[np.float64], forcing: NDArray[np.float64]
) -> OptimizeResult:
    """Deconvolution with mean value and an unknown finite mean."""

    # This uses the last term in expression (24) with the mean value as an additional
    # unknown.
    def costfun(
        x: NDArray[np.float64],
        signal: NDArray[np.float64],
        forcing: NDArray[np.float64],
    ) -> NDArray[np.float64]:
        # Here, the last member of x is the unknown mean value.
        x_tmp = x[:-1]
        mu = x[-1]
        tmp = signal - np.convolve(forcing, x_tmp, "same") - mu * np.ones(x_tmp.size)
        return 0.5 * np.sum((tmp) ** 2)

    bounds = ((0.0, np.inf),) * signal.size
    bounds_mu = ((0.0, np.inf),)
    bounds += bounds_mu

    return minimize(
        costfun,
        np.ones(signal.size + 1),
        args=(signal, forcing),
        bounds=bounds,
    )


def alternative_deconv(
    signal: np.ndarray,
    forcing: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Alternative deconvolution method.

    This method uses the last term in expression (24) with the mean value as an
    additional unknown.

    Parameters
    ----------
    signal : np.ndarray
        The observed signal.
    forcing : np.ndarray
        The forcing signal.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        A tuple containing the deconvolved forcing signal and the number of iterations
        performed by the optimizer.
    """
    res = deconv_mean_2(signal, forcing)
    return res.x[:-1], np.asarray(res.nit)
