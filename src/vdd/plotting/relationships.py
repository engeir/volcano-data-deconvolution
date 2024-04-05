"""Plot the relationships between AOD, RF, T and more."""

from collections.abc import Callable

import cosmoplots
import matplotlib.pyplot as plt
import numba as nb
import numpy as np
import sympy as sp
import volcano_base
import xarray as xr
from scipy.optimize import curve_fit

import vdd.load

_SAVE_DIR = volcano_base.config.SAVE_PATH / "relationships"
if not _SAVE_DIR.exists():
    _SAVE_DIR.mkdir(parents=False)

plt.rc("text.latex", preamble=r"\usepackage{amsmath}")
plt.style.use(
    [
        "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle",
        "vdd.extra",
        {"legend.fontsize": 6},
    ]
)
DataCESM = vdd.load.CESMData
DecCESM = vdd.load.DeconvolveCESM
# CESM2
dec_4sep = DecCESM(pad_before=True, cesm=DataCESM(strength="double-overlap"))
dec_2sep = DecCESM(pad_before=True, cesm=DataCESM(strength="tt-2sep"))
dec_e = DecCESM(pad_before=True, cesm=DataCESM(strength="size5000"))
dec_s = DecCESM(pad_before=True, cesm=DataCESM(strength="strong"))
dec_p = DecCESM(pad_before=True, cesm=DataCESM(strength="medium-plus"))
dec_m = DecCESM(pad_before=True, cesm=DataCESM(strength="medium"))
# OB16
dec_ob16 = vdd.load.DeconvolveOB16(data="h0")
dec_ob16.name = "OB16 month"


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

    def __init__(self, dec: vdd.load.Deconvolve) -> None:
        self.reset_all_switch = False
        self.params_so2 = (3.75111622e-01, 2.71076963e-05)
        # self.params_aod = (0.40676405, 4.07634108)
        self.params_aod = (0.5, 18e3)
        self.params_aod_aod = (0.4093926, 280.0960872, 0.40976931, 280.09609081)
        self.params_aod_rf = (0.5318085, 102.49535642, 29.52914528, 124.08601737)
        self.params_rf = (18.50045044, 2.53849499)
        self.delta_pulses = dec.so2.dropna("time").data
        match dec:
            case vdd.load.DeconvolveOB16():
                self.type_ = vdd.utils.clean_filename(f"ob16 {dec.name}")
                self._setup_ob16(dec)
            case vdd.load.DeconvolveCESM():
                self.type_ = vdd.utils.clean_filename(f"cesm {dec.name}")
                self._setup_cesm(dec)
            case _:
                raise ValueError("Invalid input.")

    def _setup_ob16(self, dec: vdd.load.DeconvolveOB16) -> None:
        time_axis_ = dec.rf.dropna("time").time.data
        self.time_axis = np.asarray(volcano_base.manipulate.dt2float(time_axis_))
        so2_ob16 = dec.data.so2.dropna("time")
        so2_ob16 = so2_ob16.assign_coords(
            {"time": volcano_base.manipulate.dt2float(so2_ob16.time.data)}
        )
        so2_ob16 = so2_ob16[int(349 * 12) + 2 :]
        so2_ob16 = so2_ob16[: len(self.time_axis)]
        so2_ob16 = so2_ob16.assign_coords({"time": self.time_axis})
        so2_ob16, *_ = xr.align(
            so2_ob16, dec.so2.assign_coords({"time": self.time_axis})
        )
        self.so2_true = np.roll(so2_ob16.data, 0) * 1e-6
        self.rf_true = dec.rf.dropna("time").data

    def _setup_cesm(self, dec: vdd.load.DeconvolveCESM) -> None:
        self.time_axis = dec.rf.dropna("time").time.data
        self.so2_true = dec.tmso2.data
        self.aod_true = dec.aod.data
        self.rf_true = np.roll(dec.rf.data, -1)

    @property
    def so2_fake(self) -> np.ndarray:
        """Fake SO2 data."""
        return self.numerical_so2(self.time_axis, self.delta_pulses, *self.params_so2)

    @property
    def aod_fake(self) -> np.ndarray:
        """Fake AOD data."""
        return self.numerical_aod(self.time_axis, self.so2_fake, *self.params_aod)

    @property
    def aod_aod_fake(self) -> np.ndarray:
        """Fake AOD data fitted from SO2 via AOD."""
        return self.numerical_aod_aod(
            self.time_axis, self.so2_fake, *self.params_aod_aod
        )

    @property
    def aod_rf_fake(self) -> np.ndarray:
        """Fake RF data fitted from SO2 via AOD."""
        return self.numerical_aod_rf(self.time_axis, self.so2_fake, *self.params_aod_rf)

    @property
    def rf_fake(self) -> np.ndarray:
        """Fake RF data."""
        return self.numerical_rf(self.aod_fake, *self.params_rf)

    def reset_params_so2(self) -> None:
        """Reset the SO2 parameters."""
        dp = self.delta_pulses
        ta = self.time_axis
        st = self.so2_true
        try:
            params_so2, _ = curve_fit(self.numerical_so2_fit(dp), ta, st)
        except RuntimeError as e:
            print(e)
            print("Using previous parameters")
            params_so2 = self.params_so2
        self.params_so2 = params_so2

    def reset_params_aod(self) -> None:
        """Reset the AOD parameters.

        Notes
        -----
        Keep in mind that when re-setting/updating the parameters, the fitting is done
        based on the SO2 estimate. Thus, re-setting/updating the SO2 estimate would
        invalidate this estimate.
        """
        ta = self.time_axis
        sf = self.so2_fake
        at = self.aod_true
        try:
            params_aod, _ = curve_fit(self.numerical_aod_fit(sf), ta, at)
        except RuntimeError as e:
            print(e)
            print("Using previous parameters")
            params_aod = self.params_aod
        self.params_aod = params_aod

    def reset_params_aod_aod(self) -> None:
        """Reset the AOD-AOD parameters.

        Notes
        -----
        Keep in mind that when re-setting/updating the parameters, the fitting is done
        based on the SO2 estimate. Thus, re-setting/updating of the SO2 estimate would
        invalidate this estimate.
        """
        ta = self.time_axis
        sf = self.so2_fake
        at = self.aod_true
        try:
            params_aod_aod, _ = curve_fit(self.numerical_aod_aod_fit(sf), ta, at)
        except RuntimeError as e:
            print(e)
            print("Using previous parameters")
            params_aod_aod = self.params_aod_aod
        self.params_aod_aod = params_aod_aod

    def reset_params_aod_rf(self) -> None:
        """Reset the AOD-RF parameters.

        Notes
        -----
        Keep in mind that when re-setting/updating the parameters, the fitting is done
        based on the SO2 estimate. Thus, re-setting/updating of the SO2 estimate would
        invalidate this estimate.
        """
        ta = self.time_axis
        sf = self.so2_fake
        rt = self.rf_true
        try:
            params_aod_rf, _ = curve_fit(self.numerical_aod_rf_fit(sf), ta, rt)
        except RuntimeError as e:
            print(e)
            print("Using previous parameters")
            params_aod_rf = self.params_aod_rf
        self.params_aod_rf = params_aod_rf

    def reset_params_rf(self) -> None:
        """Reset the RF parameters.

        Notes
        -----
        Keep in mind that when re-setting/updating the parameters, the fitting is done
        based on the AOD estimate. Thus, re-setting/updating of the AOD estimate would
        invalidate this estimate.
        """
        af = self.aod_fake
        rt = self.rf_true
        try:
            params_rf, _ = curve_fit(self.numerical_rf, af, rt)
        except RuntimeError as e:
            print(e)
            print("Using previous parameters")
            params_rf = self.params_rf
        self.params_rf = params_rf

    def reset_all(self) -> None:
        """Reset all parameters."""
        self.reset_all_switch = True
        print("")
        print(f"Resetting all parameters for {self.type_}")
        print("Resetting SO2 parameters")
        self.reset_params_so2()
        if self.type_.name.startswith("cesm"):
            print("Resetting AOD parameters")
            self.reset_params_aod()
            print("Resetting AOD-AOD parameters")
            self.reset_params_aod_aod()
        print("Resetting AOD-RF parameters")
        self.reset_params_aod_rf()
        print("Resetting RF parameters")
        self.reset_params_rf()

    @staticmethod
    @nb.njit
    def numerical_so2(
        time_axis: np.ndarray, delta_pulses: np.ndarray, tau: float, scale: float
    ) -> np.ndarray:
        """Solution to SO2 exp decay time series."""
        out = np.zeros_like(time_axis)
        for i, _ in enumerate(time_axis):
            t = time_axis[: i + 1]
            s_core = np.exp(-(t[-1] - t) / tau) * delta_pulses[: i + 1]
            out[i] = np.trapz(s_core, t)
        return out * scale

    def numerical_so2_fit(self, base_array: np.ndarray) -> Callable:
        """Wrap so that I can do curve fitting."""

        def _wrapped(time_axis: np.ndarray, tau: float, scale: float) -> np.ndarray:
            return self.numerical_so2(time_axis, base_array, tau, scale)

        return _wrapped

    @staticmethod
    @nb.njit
    def numerical_aod(
        time_axis: np.ndarray, so2: np.ndarray, tau: float, scale: float
    ) -> np.ndarray:
        """Solution to AOD from SO2."""
        out = np.zeros_like(time_axis)
        for i, _ in enumerate(time_axis):
            t = time_axis[: i + 1]
            a_core = np.exp(-(t[-1] - t) / tau) * scale * so2[: i + 1]
            out[i] = np.trapz(a_core, t)
        return out

    def numerical_aod_fit(self, base_array: np.ndarray) -> Callable:
        """Wrap so that I can do curve fitting."""

        def _wrapped(time_axis: np.ndarray, tau: float, scale: float) -> np.ndarray:
            return self.numerical_aod(time_axis, base_array, tau, scale)

        return _wrapped

    def numerical_aod_aod(
        self, time_axis: np.ndarray, so2: np.ndarray, *params: float
    ) -> np.ndarray:
        """Compute the AOD from an intermediate AOD of SO2, and SO2, simultaneously."""
        param_len = 4
        assert len(params) == param_len
        aod = self.numerical_aod(time_axis, so2, params[0], params[1])
        return self.numerical_aod(time_axis, aod, params[2], params[3])

    def numerical_aod_aod_fit(self, base_array: np.ndarray) -> Callable:
        """Wrap so that I can do curve fitting."""

        def _wrapped(
            time_axis: np.ndarray,
            tau1: float,
            scale1: float,
            tau2: float,
            scale2: float,
        ) -> np.ndarray:
            return self.numerical_aod_aod(
                time_axis, base_array, tau1, scale1, tau2, scale2
            )

        return _wrapped

    def numerical_aod_rf(
        self, time_axis: np.ndarray, so2: np.ndarray, *params: float
    ) -> np.ndarray:
        """Compute the RF from AOD and SO2 simultaneously."""
        param_len = 4
        assert len(params) == param_len
        aod = self.numerical_aod(time_axis, so2, params[0], params[1])
        return self.numerical_rf(aod, params[2], params[3])

    def numerical_aod_rf_fit(self, base_array: np.ndarray) -> Callable:
        """Wrap so that I can do curve fitting."""

        def _wrapped(
            time_axis: np.ndarray,
            tau: float,
            scale_aod: float,
            scale_rf1: float,
            scale_rf2: float,
        ) -> np.ndarray:
            return self.numerical_aod_rf(
                time_axis, base_array, tau, scale_aod, scale_rf1, scale_rf2
            )

        return _wrapped

    @staticmethod
    @nb.njit
    def numerical_rf(aod: np.ndarray, scale_a: float, scale_b: float) -> np.ndarray:
        """Solution to RF from AOD."""
        rf = scale_a * np.log(1 + scale_b * aod)
        return rf

    def _plot_so2(self, msg: str) -> None:
        plt.figure()
        plt.plot(self.time_axis, self.so2_true, label="Simulation output")
        plt.plot(self.time_axis, self.so2_fake, label="Numerical solution")
        plt.legend()
        plt.xlabel("Time [yr]")
        plt.ylabel("SO$_2$ [kg/m$^2$]")
        plt.savefig(_SAVE_DIR / f"numerical_so2_{self.type_}{msg}")

    def _plot_aod(self, msg: str) -> None:
        plt.figure()
        plt.plot(self.time_axis, self.aod_true, label="Simulation output")
        plt.plot(self.time_axis, self.aod_fake, label="Numerical soln (A(S))")
        plt.plot(self.time_axis, self.aod_aod_fake, label="Numerical soln (A(A(S)))")
        plt.legend()
        plt.xlabel("Time [yr]")
        plt.ylabel("Aerosol optical depth [1]")
        plt.savefig(_SAVE_DIR / f"numerical_aod_{self.type_}{msg}")

    def _plot_rf(self, msg: str) -> None:
        plt.figure()
        plt.plot(self.time_axis, self.rf_true, label="Simulation output")
        plt.plot(self.time_axis, self.rf_fake, label="Numerical soln (R(A))")
        plt.plot(self.time_axis, self.aod_rf_fake, label="Numerical soln (R(A(S)))")
        plt.legend()
        plt.xlabel("Time [yr]")
        plt.ylabel("Radiative forcing [W/m$^2$]")
        plt.savefig(_SAVE_DIR / f"numerical_rf_{self.type_}{msg}")

    def plot_available(self) -> None:
        """Plot the available data."""
        if self.reset_all_switch:
            msg = "_optimised"
        else:
            msg = ""
        self._plot_so2(msg)
        if self.type_.name.startswith("cesm"):
            self._plot_aod(msg)
        self._plot_rf(msg)


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
                case xr.DataArray(x), xr.DataArray(y), str(lname):
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
        cum_x = np.empty(0)
        cum_y = np.empty(0)
        for x_, y_, lab_ in zip(x, y, labelname, strict=True):
            if str(y_[0].name) == "RESTOM":
                noise_floor = 0.02
                idx = x_.data > noise_floor
                x_data = x_.data[idx]
                y_data = y_.data[idx]
                cum_x = np.append(cum_x, x_data)
                cum_y = np.append(cum_y, y_data)
                lin_regress = np.polyfit(cum_x, cum_y, 1)
                print(lin_regress[0])
                a.plot(
                    [cum_x.min(), cum_x.max()],
                    lin_regress[0] * np.asarray([cum_x.min(), cum_x.max()])
                    + lin_regress[1],
                    label="_he",
                )
            a.scatter(x_, y_, label=lab_)
        a.set_xlabel(str(x[0].name))
        a.set_ylabel(str(y[0].name))
        a.legend()
        f.savefig(_SAVE_DIR / f"relationship_{i}.png")
        plt.show()
        plt.close(f)


def _scatterplot_comparison() -> None:
    # Check how well response functions are estimated from double waveforms.
    decs = (dec_m, dec_p, dec_4sep, dec_2sep, dec_s, dec_e)
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
    PlotRelationship((aod, rf, labels), (aod, temp, labels), (rf, temp, labels)).plot()


def _analytic() -> None:
    a = AnalyticSolution()
    a.analytic_so2()


def _numerical_solver() -> None:
    decs = (dec_4sep, dec_2sep, dec_e, dec_s, dec_p, dec_m, dec_ob16)
    for dec in decs:
        ns = NumericalSolver(dec)
        ns.plot_available()
        ns.reset_all()
        ns.plot_available()
        for plot in ("so2", "aod", "rf"):
            f1 = _SAVE_DIR / f"numerical_{plot}_{ns.type_}.png"
            f2 = _SAVE_DIR / f"numerical_{plot}_{ns.type_}_optimised.png"
            try:
                cosmoplots.combine(f1, f2).in_grid(2, 1).using(fontsize=50).save(
                    _SAVE_DIR / f"numerical_{plot}_{ns.type_}_combined.png"
                )
            except FileNotFoundError:
                pass
            f1.unlink(missing_ok=True)
            f2.unlink(missing_ok=True)
        plt.show()
        plt.close("all")


if __name__ == "__main__":
    _numerical_solver()
    _scatterplot_comparison()
