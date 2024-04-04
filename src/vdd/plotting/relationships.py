"""Plot the relationships between AOD, RF, T and more."""

from collections.abc import Callable

import matplotlib.pyplot as plt
import numba as nb
import numpy as np
import sympy as sp
import volcano_base
import xarray as xr

import vdd.load

_SAVE_DIR = volcano_base.config.SAVE_PATH / "relationships"
if not _SAVE_DIR.exists():
    _SAVE_DIR.mkdir(parents=False)

plt.rc("text.latex", preamble=r"\usepackage{amsmath}")
plt.style.use(
    [
        "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle",
        "vdd.extra",
    ]
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
# OB16
dec_ob16_month = vdd.load.DeconvolveOB16(data="h0")
dec_ob16_month.name = "OB16 month"


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


class NumericalSolver:
    """Numerical solution to the SO2, AOD and RF relationships."""

    def __init__(self, time_axis: np.ndarray, delta_pulses: np.ndarray) -> None:
        self.time_axis = time_axis
        self.delta_pulses = delta_pulses


@nb.njit
def numerical_so2(
    time_axis: np.ndarray,
    delta_pulses: np.ndarray,
    tau: float = 0.3,
    scale: float = 3e-5,
) -> np.ndarray:
    """Solution to SO2 exp decay time series."""
    out = np.zeros_like(time_axis)
    for i, _ in enumerate(time_axis):
        t = time_axis[: i + 1]
        s_core = np.exp(-(t[-1] - t) / tau) * delta_pulses[: i + 1]
        out[i] = np.trapz(s_core, t)
    return out * scale


def numerical_so2_fit(base_array: np.ndarray) -> Callable:
    """Wrap so that I can do curve fitting."""

    def _wrapper(time_axis: np.ndarray, tau: float, scale: float) -> np.ndarray:
        return numerical_so2(time_axis, base_array, tau, scale)

    return _wrapper


@nb.njit
def numerical_aod(
    time_axis: np.ndarray, so2: np.ndarray, tau: float = 0.5, scale: float = 18e3
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
        out[i] = np.trapz(a_core, t)
    return out


def numerical_aod_fit(base_array: np.ndarray) -> Callable:
    """Wrap so that I can do curve fitting."""

    def _wrapper(time_axis: np.ndarray, tau: float, scale: float) -> np.ndarray:
        return numerical_aod(time_axis, base_array, tau, scale)

    return _wrapper


@nb.njit
def numerical_rf(
    aod: np.ndarray, scale_a: float = 24, scale_b: float = 1
) -> np.ndarray:
    """Solution to RF from AOD."""
    # We assume a non-linear relationship, and specifially a logaritmic one.
    rf = scale_a * np.log(1 + scale_b * aod)
    return rf


class PlotRelationship:
    """Plot the relationship between pairs of data arrays."""

    def __init__(
        self,
        *pairs: tuple[
            xr.DataArray | list[xr.DataArray],
            xr.DataArray | list[xr.DataArray],
            str | list[str],
        ],
    ) -> None:
        self.pairs = pairs

    def plot(self) -> None:
        """Plot the relationships between the pairs of data arrays."""
        for i, (x, y, labelname) in enumerate(self.pairs):
            match x, y, labelname:
                case xr.DataArray(), xr.DataArray(), str(lname):
                    self._plot_single(x, y, i, lname)
                case list(x_), list(y_), list(lname):
                    self._plot_multiple(x_, y_, i, lname)
                case _:
                    raise ValueError("Invalid input.")

    @staticmethod
    def _plot_single(x: xr.DataArray, y: xr.DataArray, i: int, labelname: str) -> None:
        f = plt.figure()
        a = f.gca()
        a.scatter(x, y, label=labelname)
        a.set_xlabel(str(x.name))
        a.set_ylabel(str(y.name))
        a.legend()
        f.savefig(_SAVE_DIR / f"relationship_{i}.png")
        plt.show()
        plt.close(f)

    @staticmethod
    def _plot_multiple(
        x: list[xr.DataArray], y: list[xr.DataArray], i: int, labelname: list[str]
    ) -> None:
        f = plt.figure()
        a = f.gca()
        for x_, y_, lab_ in zip(x, y, labelname, strict=True):
            a.scatter(x_, y_, label=lab_)
        a.set_xlabel(str(x[0].name))
        a.set_ylabel(str(y[0].name))
        a.legend()
        f.savefig(_SAVE_DIR / f"relationship_{i}.png")
        plt.show()
        plt.close(f)


def _scatterplot_comparison() -> None:
    # Check how well response functions are estimated from double waveforms.
    decs = (dec_cesm_2sep, dec_cesm_e, dec_cesm_s, dec_cesm_p, dec_cesm_m)
    # Plot the relationships between pairs of data arrays.
    # aod = dec_cesm_4sep.aod
    # rf = dec_cesm_4sep.rf
    # temp = dec_cesm_4sep.temp
    # for dec in decs[1:]:
    #     aod = xr.concat((aod, dec.aod), dim="time")
    #     rf = xr.concat((rf, dec.rf), dim="time")
    #     temp = xr.concat((temp, dec.temp), dim="time")
    # PlotRelationship(
    #     (aod, rf, "aod-rf"), (aod, temp, "aod-temp"), (rf, temp, "rf-temp")
    # ).plot()
    aod = []
    rf = []
    temp = []
    labels = []
    for dec in decs:
        aod.append(dec.aod)
        rf.append(dec.rf)
        temp.append(dec.temp)
        labels.append(dec.name)
    print(decs[0].rf.attrs)
    PlotRelationship((aod, rf, labels), (aod, temp, labels), (rf, temp, labels)).plot()


def _analytic() -> None:
    # a = AnalyticSolution()
    # a.analytic_so2()
    dec: vdd.load.DeconvolveCESM | vdd.load.DeconvolveOB16 = dec_ob16_month
    time_axis_ = dec.rf.dropna("time").time.data
    match dec:
        case vdd.load.DeconvolveOB16():
            time_axis = np.asarray(volcano_base.manipulate.dt2float(time_axis_))
            so2_ob16 = dec.data.so2.dropna("time")
            so2_ob16 = so2_ob16.assign_coords(
                {"time": volcano_base.manipulate.dt2float(so2_ob16.time.data)}
            )
            so2_ob16 = so2_ob16[int(349 * 12) + 2 :]
            so2_ob16 = so2_ob16[: len(time_axis)]
            so2_ob16 = so2_ob16.assign_coords({"time": time_axis})
            so2_ob16, *_ = xr.align(
                so2_ob16, dec.so2.assign_coords({"time": time_axis})
            )
            so2_true = np.roll(so2_ob16.data, 0) * 1e-6
            rf_true = dec.rf.dropna("time").data
        case vdd.load.DeconvolveCESM():
            time_axis = time_axis_
            so2_true = dec.tmso2.data
            aod_true = dec.aod.data
            rf_true = np.roll(dec.rf.data, -1)
    delta_pulses = dec.so2.dropna("time").data
    # print(aod_true)
    so2_fake = numerical_so2(time_axis, delta_pulses, 3.75111622e-01, 2.71076963e-05)
    # params_so2, _ = curve_fit(numerical_so2_fit(delta_pulses), time_axis, so2_true)
    # so2_fake_fit = numerical_so2(time_axis, delta_pulses, *params_so2)
    # print(params_so2)  # [3.75111622e-01 2.71076963e-05]
    aod_fake = numerical_aod(time_axis, so2_fake)
    aod_fake_2 = numerical_aod(time_axis, aod_fake, 0.40676405, 4.07634108)
    # params_aod, _ = curve_fit(numerical_aod_fit(aod_fake), time_axis, dec.aod.data)
    # aod_fake_fit = numerical_aod(time_axis, aod_fake, *params_aod)
    # print(params_aod)  # [0.40676405 4.07634108]
    # aod_fake = numerical_aod(time_axis, so2=dec.tmso2.data)
    rf_fake = numerical_rf(aod_fake, 18.50045044 * 0.8, 2.53849499)
    # params_rf, _ = curve_fit(numerical_rf, aod_fake, dec.rf.data)
    # rf_fake_fit = numerical_rf(aod_fake, *params_rf)
    # print(params_rf)  # [18.50045044  2.53849499]
    plt.figure()
    plt.plot(time_axis, so2_fake, label="fake")
    # plt.plot(time_axis, so2_fake_fit, label="fake fit")
    plt.plot(time_axis, so2_true, label="true")
    plt.legend()
    plt.figure()
    plt.plot(time_axis, aod_fake, label="fake")
    # plt.plot(time_axis, aod_fake_fit, label="fake fit")
    plt.plot(time_axis, aod_fake_2, label="fake2")
    # plt.plot(time_axis, aod_true, label="true")
    plt.legend()
    plt.figure()
    plt.plot(time_axis, rf_true, label="true")
    plt.plot(time_axis, rf_fake, label="fake")
    # plt.plot(time_axis, rf_fake_fit, label="fake fit")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    _analytic()
    # _scatterplot_comparison()
