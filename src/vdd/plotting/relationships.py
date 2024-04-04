"""Plot the relationships between AOD, RF, T and more."""

import matplotlib.pyplot as plt
import numba as nb
import numpy as np
import sympy as sp
import xarray as xr

import vdd.load

plt.rc("text.latex", preamble=r"\usepackage{amsmath}")
plt.style.use(
    "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle"
)
DataCESM = vdd.load.CESMData
DecCESM = vdd.load.DeconvolveCESM
# CESM2
dec_cesm_4sep = DecCESM(pad_before=True, cesm=DataCESM(strength="double-overlap"))
dec_cesm_2sep = DecCESM(pad_before=True, cesm=DataCESM(strength="tt-2sep"))
dec_cesm_e = DecCESM(pad_before=True, cesm=DataCESM(strength="size5000"))
dec_cesm_s = DecCESM(pad_before=True, cesm=DataCESM(strength="strong"))
dec_cesm_p = DecCESM(pad_before=True, cesm=DataCESM(strength="medium-plus"))
dec_cesm_m = DecCESM(pad_before=True, cesm=DataCESM(strength="medium"))


class AnalyticSolution:
    """Attempt to find functions representing AOD from SO2."""

    # Assume A depend linearly on SO2 (at least for small/realistic values)
    # Assume RF depend non-linearly on AOD
    # Assume T superpose and is constructed from RF and a response function
    def ode_so2(self) -> None:
        """Sympy solution to SO2."""
        tau_s, t, k, n = sp.symbols("tau_s, t, k, n", real=True)
        s = sp.Function("s")
        t_k = sp.IndexedBase("t_k")
        s_k = sp.IndexedBase("s_k")
        dirac = sp.DiracDelta(t - t_k[k])
        rhs = -s(t) / tau_s + sp.Sum(s_k[k] * dirac, (k, 0, n))
        lhs = sp.diff(s(t), t)
        expr = sp.Eq(lhs, rhs)
        solution = sp.dsolve(expr, s(t))
        check = sp.checkodesol(expr, solution)
        print(expr)
        print(solution)
        # Eq((tau_s*Integral(exp(t/tau_s)*Sum(DiracDelta(t - t_k[_0])*s_k[_0], (_0, 0, n)), t) - Integral(s(t)*exp(t/tau_s), t))/tau_s, C3)
        print(check)

    def analytic_so2(self) -> None:
        """Second sympy solution to SO2."""
        tau_s, t, t_upper, k, n = sp.symbols("tau_s, t, t_upper, k, n", real=True)
        # s = sp.Function("s")
        t_k = sp.IndexedBase("t_k")
        s_k = sp.IndexedBase("s_k")
        dirac = sp.DiracDelta(t - t_k[k])
        exp = sp.exp(-(t_upper - t) / tau_s)
        sum = sp.Sum(s_k[k] * dirac, (k, 0, n))
        integral = sp.integrate(exp * sum, t)
        print(integral)


@nb.njit
def numerical_so2(time_axis: np.ndarray, delta_pulses: np.ndarray) -> np.ndarray:
    """Solution to SO2 exp decay time series."""
    out = np.zeros_like(time_axis)
    tau = 0.3
    scale = 3e-5
    for i, _ in enumerate(time_axis):
        t = time_axis[: i + 1]
        s_core = np.exp(-(t[-1] - t) / tau) * delta_pulses[: i + 1]
        hey = np.trapz(s_core, t)
        out[i] = hey
    return out * scale


@nb.njit
def numerical_aod(
    time_axis: np.ndarray, so2: np.ndarray, tau: float, scale: float
) -> np.ndarray:
    """Solution to AOD from SO2."""
    out = np.zeros_like(time_axis)
    tau_1 = tau
    # tau_2 = 1.8
    # scale = 30e3
    for i, _ in enumerate(time_axis):
        t = time_axis[: i + 1]
        a_core = (
            np.exp(-(t[-1] - t) / tau_1)
            # * np.exp(-(t[-1] - t) / tau_2)
            * scale
            * so2[: i + 1]
        )
        hey = np.trapz(a_core, t)
        out[i] = hey
    return out


@nb.njit
def numerical_rf(aod: np.ndarray) -> np.ndarray:
    """Solution to RF from AOD."""
    # We assume a non-linear relationship, and specifially a logaritmic one.
    scale = 24
    rf = scale * np.log(1 + aod)
    return rf


class PlotRelationship:
    """Plot the relationship between pairs of data arrays."""

    def __init__(self, *pairs: tuple[xr.DataArray, xr.DataArray, str]) -> None:
        self.pairs = pairs

    def plot(self) -> None:
        """Plot the relationships between the pairs of data arrays."""
        for i, (x, y, figname) in enumerate(self.pairs):
            f = plt.figure()
            a = f.gca()
            a.scatter(x, y)
            a.set_xlabel(str(x.name))
            a.set_ylabel(str(y.name))
            f.savefig(f"relationship_{figname}_{i}.png")
            plt.show()
            plt.close(f)


def _main() -> None:
    # Check how well response functions are estimated from double waveforms.
    decs = (dec_cesm_2sep, dec_cesm_e, dec_cesm_s, dec_cesm_p, dec_cesm_m)
    # Plot the relationships between pairs of data arrays.
    aod = dec_cesm_4sep.aod
    rf = dec_cesm_4sep.rf
    temp = dec_cesm_4sep.temp
    print(temp)
    for dec in decs:
        aod = xr.concat((aod, dec.aod), dim="time")
        rf = xr.concat((rf, dec.rf), dim="time")
        temp = xr.concat((temp, dec.temp), dim="time")
    print(temp)
    # PlotRelationship((aod, rf), (aod, temp), (rf, temp)).plot()


if __name__ == "__main__":
    # a = AnalyticSolution()
    # a.analytic_so2()
    dec = dec_cesm_2sep
    time_axis = dec.aod.time.data
    delta_pulses = dec.so2.data
    so2_true = dec.tmso2.data
    so2_fake = numerical_so2(time_axis, delta_pulses)
    aod_fake = numerical_aod(time_axis, so2_fake, 0.5, 18e3)
    aod_fake_2 = numerical_aod(time_axis, aod_fake, 0.4, 6)
    # aod_fake = numerical_aod(time_axis, so2=dec.tmso2.data)
    rf_fake = numerical_rf(aod_fake)
    plt.figure()
    plt.plot(time_axis, so2_fake)
    plt.plot(time_axis, so2_true)
    plt.figure()
    plt.plot(time_axis, aod_fake, label="fake")
    plt.plot(time_axis, aod_fake_2, label="fake2")
    plt.plot(time_axis, dec.aod, label="true")
    plt.legend()
    plt.figure()
    plt.plot(time_axis, rf_fake)
    plt.plot(time_axis, dec.rf)
    plt.show()
