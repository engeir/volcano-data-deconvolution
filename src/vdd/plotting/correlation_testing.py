"""Simple manipulation of time series to verify effect in correlation function."""

import fppanalysis
import matplotlib.pyplot as plt
import numpy as np

plt.style.use(
    [
        "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle",
        "vdd.extra",
        "vdd.jgr",
    ],
)


def shapes_different_timescales() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Create two shapes with different time scales."""
    t1 = np.linspace(0, 100, 1000)
    shape1 = np.linspace(
        -50,
        50,
        1000,
    )  # + np.random.default_rng().normal(0, 0.5, size=1000)
    shape2 = np.linspace(-50, 50, 1000)
    zero_time_idx = 50
    cutoff_time_idx = 51
    shape1[t1 > zero_time_idx] = np.exp(-shape1[t1 > zero_time_idx] / 1)
    shape1[t1 < zero_time_idx] = 0
    # shape2[t1 > zero_time_idx] = np.exp(-shape2[t1 > zero_time_idx] / 0.6)
    shape2[t1 > zero_time_idx] = np.exp(-shape2[t1 > zero_time_idx] / 1)
    shape2[t1 > cutoff_time_idx] = 0
    shape2[t1 < zero_time_idx] = 0
    return t1, shape1, shape2


def shapes_different_rise() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Create two shapes with different rise times."""
    t1 = np.linspace(0, 100, 1000)
    shape2 = np.linspace(
        -50,
        50,
        1000,
    )  # + np.random.default_rng().normal(0, 0.5, size=1000)
    shape1 = np.linspace(-50, 50, 1000)
    zero_time_idx = 50
    rise_time_idx = 49.1
    exp_region_start = np.max((zero_time_idx, rise_time_idx))
    zero_region_end = np.min((zero_time_idx, rise_time_idx))
    shape2[t1 > exp_region_start] = np.exp(-shape2[t1 > exp_region_start] / 1)
    shape2[t1 < zero_region_end] = 0
    shape2[(t1 < exp_region_start) & (t1 > zero_region_end)] = np.linspace(
        0,
        np.abs(shape2[int(exp_region_start * 10) - 0]),
        len(shape2[(t1 < exp_region_start) & (t1 > zero_region_end)]) + 1,
    )[:-1]
    shape1[t1 > zero_time_idx] = np.exp(-shape1[t1 > zero_time_idx] / 1)
    shape1[t1 < zero_time_idx] = 0
    return t1, shape1, np.roll(shape2, 3)


def correlation_test() -> (
    tuple[
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
    ]
):
    """Create a correlation function from time series using the shapes."""
    t1, shape1, shape2 = shapes_different_timescales()
    # t1, shape1, shape2 = shapes_different_rise()
    delta1 = np.zeros(1000)
    delta2 = np.zeros(1000)
    for i in np.random.default_rng().uniform(1, 900, size=30):
        delta1[int(i)] = np.random.default_rng().uniform(0, 0.5, size=1)[0]
        delta2[int(i)] = np.random.default_rng().uniform(0, 0.5, size=1)[0]
    f1 = -np.convolve(shape1, delta1, mode="same")
    f2 = -np.convolve(shape2, delta2, mode="same")
    t2, c = fppanalysis.corr_fun(f1, f2, 0.1)
    t3, auto = fppanalysis.corr_fun(f1, f1, 0.1)
    return t2, c, t3, auto, t1, f1, f2


if __name__ == "__main__":
    corrs = []
    auto = []
    t = None
    for _ in range(1000):
        t, c, t3, a, t1, f1, f2 = correlation_test()
        corrs.append(c)
        auto.append(a)
    c = np.mean(corrs, axis=0)
    a = np.mean(auto, axis=0)
    plt.plot(t1, f1, label="T")
    plt.plot(t1, f2, c="g", label="phi*f")
    plt.legend()
    plt.figure()
    plt.plot(t3, a, label="Auto")
    plt.plot(t, c, c="g", label="Deconvolved")  # type: ignore[arg-type]
    plt.xlim((-10, 10))
    plt.legend()
    plt.show()
