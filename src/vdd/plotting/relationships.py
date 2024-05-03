"""Plot the relationships between AOD, RF, T and more."""

import json
import warnings
from collections.abc import Callable
from typing import Literal

import cosmoplots
import matplotlib.pyplot as plt
import numba as nb
import numpy as np
import rich.prompt
import sympy as sp
import volcano_base
import xarray as xr
from scipy.optimize import curve_fit

import vdd.load
from vdd.utils import name_swap as nsw

# TODO:
# - Create a plot where all simulations uses the parameters of a given simulation. It
#   would be nice if this was OB16, and in the case of RF we can do that, but since we
#   do not have any AOD data from OB16, we cannot recreate all plots from OB16.
# - Another alternative to using all parameters from the single simulation is to allow
#   the scaling to be optimised, but keep the time units as is.

warnings.warn(
    "This script is deprecated and may be removed in the future.",
    category=DeprecationWarning,
    stacklevel=1,
)

_SAVE_DIR = volcano_base.config.SAVE_PATH / "relationships"
if not _SAVE_DIR.exists():
    _SAVE_DIR.mkdir(parents=False)

plt.style.use([
    "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle",
    "vdd.extra",
    {
        "legend.fontsize": 6,
        "text.latex.preamble": r"\usepackage{amsmath,siunitx}",
    },
])
T_Params = Literal["SO2", "AOD", "AOD-AOD", "AOD-RF", "RF-as-AOD", "RF"]
DataCESM = vdd.load.CESMData
DecCESM = vdd.load.DeconvolveCESM
# CESM2
dec_4sep = DecCESM(pad_before=True, cesm=DataCESM(strength="tt-4sep"))
dec_2sep = DecCESM(pad_before=True, cesm=DataCESM(strength="tt-2sep"))
dec_e = DecCESM(pad_before=True, cesm=DataCESM(strength="size5000"))
dec_s = DecCESM(pad_before=True, cesm=DataCESM(strength="strong"))
dec_p = DecCESM(pad_before=True, cesm=DataCESM(strength="medium-plus"))
dec_m = DecCESM(pad_before=True, cesm=DataCESM(strength="medium"))
# OB16
dec_ob16 = vdd.load.DeconvolveOB16(data="h0", length=12000)

# decs = (dec_4sep, dec_2sep, dec_e, dec_s, dec_p, dec_m, dec_ob16)
decs = (dec_4sep, dec_2sep, dec_e, dec_s, dec_p, dec_m)
# decs = (dec_ob16,)


def s2n(num: float, decimal: int = 2) -> str:
    """Convert a number to scientific notation."""
    return f"\\num{{{num:.{decimal}e}}}"


class AnalyticSolution:
    """Attempt to find functions representing AOD from SO2."""

    # Assume A depend linearly on SO2 (at least for small/realistic values)
    # Assume RF depend non-linearly on AOD
    # Assume T superpose and is constructed from RF and a response function
    @staticmethod
    def ode_so2() -> None:
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

    @staticmethod
    def analytic_so2() -> None:
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

    _MAXFEV = 10000

    def __init__(self, dec: vdd.load.Deconvolve) -> None:
        self.reset_all_switch = ""
        self.estimate_so2 = True
        self.params_so2 = {"tau_s": 0.37511174, "scale_s": 13.82492275}
        self.params_aod = {"tau_a": 0.99455712, "scale_a": 0.02882828}
        self.params_aod_aod = {
            "tau_a1": 0.40939242,
            "tau_a2": 0.39193863,
            "scale_a": 0.40976906,
        }
        self.params_aod_rf = {
            "tau_a": 0.53180791,
            "scale_r": 0.15641879,
            "scale_a": 29.52907214,
        }
        self.params_rf_as_aod = {
            "tau_a": 1.40939242,
            "tau_r": 0.39193863,
            "scale_r": 0.40976906,
        }
        self.params_rf = {"scale_r": 18.50045044, "scale_a": 2.53849499}
        self.delta_pulses = dec.so2.dropna("time").data
        match dec:
            case vdd.load.DeconvolveOB16():
                self.type_ = nsw(vdd.utils.clean_filename(f"ob16 {dec.name}"))
                self._setup_ob16(dec)
            case vdd.load.DeconvolveCESM():
                self.type_ = nsw(vdd.utils.clean_filename(f"cesm {dec.name}"))
                self._setup_cesm(dec)
            case _:
                raise ValueError("Invalid input.")

    def _setup_ob16(self, dec: vdd.load.DeconvolveOB16) -> None:
        time_axis_ = dec.rf.dropna("time").time.data
        self.time_axis = np.asarray(volcano_base.manipulate.dt2float(time_axis_))
        so2_ob16 = dec.data.so2.dropna("time")
        so2_ob16 = so2_ob16.assign_coords({
            "time": volcano_base.manipulate.dt2float(so2_ob16.time.data)
        })
        so2_ob16 = so2_ob16[int(349 * 12) + 2 :]
        so2_ob16 = so2_ob16[: len(self.time_axis)]
        so2_ob16 = so2_ob16.assign_coords({"time": self.time_axis})
        so2_ob16, *_ = xr.align(
            so2_ob16, dec.so2.assign_coords({"time": self.time_axis})
        )
        self.so2_true = np.roll(so2_ob16.data, 0)
        self.rf_true = dec.rf.dropna("time").data

    def _setup_cesm(self, dec: vdd.load.DeconvolveCESM) -> None:
        self.time_axis = dec.rf.dropna("time").time.data
        self.so2_true = dec.tmso2.data * 510e3
        self.aod_true = dec.aod.data
        self.rf_true = np.roll(dec.rf.data, -1)

    def use_true_so2(self, use_true_so2: bool) -> None:
        """Decide whether to use the true or numerically estimated SO2 data."""
        self.estimate_so2 = not use_true_so2

    @property
    def so2_fake(self) -> np.ndarray:
        """Fake SO2 data."""
        if self.estimate_so2:
            return self._numerical_so2(
                self.time_axis, self.delta_pulses, *self.params_so2.values()
            )
        return self.so2_true

    @property
    def aod_fake(self) -> np.ndarray:
        """Fake AOD data."""
        return self._numerical_aod(
            self.time_axis, self.so2_fake, *self.params_aod.values()
        )

    @property
    def aod_aod_fake(self) -> np.ndarray:
        """Fake AOD data fitted from SO2 via AOD."""
        return self._numerical_aod_aod(
            self.time_axis, self.so2_fake, *self.params_aod_aod.values()
        )

    @property
    def aod_rf_fake(self) -> np.ndarray:
        """Fake RF data fitted from SO2 via AOD."""
        return self._numerical_aod_rf(
            self.time_axis, self.so2_fake, *self.params_aod_rf.values()
        )

    @property
    def rf_as_aod_fake(self) -> np.ndarray:
        """Fake AOD data fitted from SO2 via AOD."""
        return self._numerical_aod_aod(
            self.time_axis, self.so2_fake, *self.params_rf_as_aod.values()
        )

    @property
    def rf_fake(self) -> np.ndarray:
        """Fake RF data."""
        return self._numerical_rf(self.aod_fake, *self.params_rf.values())

    def reset_params_so2(self) -> None:
        """Reset the SO2 parameters."""
        dp = self.delta_pulses
        ta = self.time_axis
        st = self.so2_true
        try:
            params_so2, _ = curve_fit(
                self._numerical_so2_fit(dp), ta, st, maxfev=self._MAXFEV
            )
        except RuntimeError as e:
            print(e)
            print("Using previous parameters")
            params_so2 = self.params_so2.values()
        self.params_so2 = {
            key: value
            for key, value in zip(self.params_so2.keys(), params_so2, strict=False)
        }

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
            params_aod, _ = curve_fit(
                self._numerical_aod_fit(sf), ta, at, maxfev=self._MAXFEV
            )
        except RuntimeError as e:
            print(e)
            print("Using previous parameters")
            params_aod = self.params_aod.values()
        self.params_aod = {
            key: value
            for key, value in zip(self.params_aod.keys(), params_aod, strict=False)
        }

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
            params_aod_aod, _ = curve_fit(
                self._numerical_aod_aod_fit(sf), ta, at, maxfev=self._MAXFEV
            )
        except RuntimeError as e:
            print(e)
            print("Using previous parameters")
            params_aod_aod = self.params_aod_aod.values()
        self.params_aod_aod = {
            key: value
            for key, value in zip(
                self.params_aod_aod.keys(), params_aod_aod, strict=False
            )
        }

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
            params_aod_rf, _ = curve_fit(
                self._numerical_aod_rf_fit(sf), ta, rt, maxfev=self._MAXFEV
            )
        except RuntimeError as e:
            print(e)
            print("Using previous parameters")
            params_aod_rf = self.params_aod_rf.values()
        self.params_aod_rf = {
            key: value
            for key, value in zip(
                self.params_aod_rf.keys(), params_aod_rf, strict=False
            )
        }

    def reset_params_rf_as_aod(self) -> None:
        """Reset the R(A(S)) parameters, where R(A) uses the A(S) functional.

        This is equivalent to resetting the AOD-AOD parameters, but optimised for the RF
        data.

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
            params_rf_as_aod, _ = curve_fit(
                self._numerical_aod_aod_fit(sf), ta, rt, maxfev=self._MAXFEV
            )
        except RuntimeError as e:
            print(e)
            print("Using previous parameters")
            params_rf_as_aod = self.params_rf_as_aod.values()
        self.params_rf_as_aod = {
            key: value
            for key, value in zip(
                self.params_rf_as_aod.keys(), params_rf_as_aod, strict=False
            )
        }

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
            params_rf, _ = curve_fit(self._numerical_rf, af, rt, maxfev=self._MAXFEV)
        except RuntimeError as e:
            print(e)
            print("Using previous parameters")
            params_rf = self.params_rf.values()
        self.params_rf = {
            key: value
            for key, value in zip(self.params_rf.keys(), params_rf, strict=False)
        }

    def reset_all(
        self,
        custom_params: dict[T_Params, dict[str, float]] | None = None,
        from_: str = "self",
        from_json: bool = True,
        save_json: bool = False,
    ) -> None:
        """Reset all parameters."""
        self.reset_all_switch = from_
        if custom_params is not None:
            self._reset_all_custom(custom_params)
            return
        true_so2 = not self.estimate_so2
        filename = (
            f"numerical_params_{self.type_}_{"real" if true_so2 else "fake"}-so2.json"
        )
        if from_json and not save_json:
            with open(_SAVE_DIR / filename, encoding="locale") as f:
                params = json.load(f)
            self._reset_all_custom(params)
            return
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
        print("Resetting RF as AOD parameters")
        self.reset_params_rf_as_aod()
        print("Resetting RF parameters")
        self.reset_params_rf()
        if save_json:
            params = self._get_params()
            with open(_SAVE_DIR / filename, "w", encoding="locale") as f:
                json.dump(params, f, indent=2)

    def _reset_all_custom(
        self, custom_params: dict[T_Params, dict[str, float]]
    ) -> None:
        """Reset all parameters with custom parameters."""
        for param, values in custom_params.items():
            match param, values:
                case "SO2", {"tau_s": _, "scale_s": _}:
                    self.params_so2 = values
                case "AOD", {"tau_a": _, "scale_a": _}:
                    self.params_aod = values
                case "AOD-AOD", {"tau_a1": _, "tau_a2": _, "scale_a": _}:
                    self.params_aod_aod = values
                case "AOD-RF", {"tau_a": _, "scale_r": _, "scale_a": _}:
                    self.params_aod_rf = values
                case "RF-as-AOD", {"tau_a": _, "tau_r": _, "scale_r": _}:
                    self.params_rf_as_aod = values
                case "RF", {"scale_r": _, "scale_a": _}:
                    self.params_rf = values
                case _:
                    raise ValueError("Invalid input.")

    def _get_params(self) -> dict[T_Params, dict[str, float]]:
        d: dict[T_Params, dict[str, float]] = {
            "SO2": self.params_so2,
            "AOD": self.params_aod,
            "AOD-AOD": self.params_aod_aod,
            "AOD-RF": self.params_aod_rf,
            "RF-as-AOD": self.params_rf_as_aod,
            "RF": self.params_rf,
        }
        return d

    def print_params(self) -> dict[T_Params, dict[str, float]]:
        """Print the parameters."""
        print(f"{self.type_} parameters")
        d = self._get_params()
        print(d)
        return d

    @staticmethod
    @nb.njit
    def _numerical_so2(
        time_axis: np.ndarray, delta_pulses: np.ndarray, tau: float, scale: float
    ) -> np.ndarray:
        """Solution to SO2 exp decay time series."""
        out = np.zeros_like(time_axis)
        for i, _ in enumerate(time_axis):
            t = time_axis[: i + 1]
            s_core = np.exp(-(t[-1] - t) / tau) * delta_pulses[: i + 1]
            out[i] = np.trapz(s_core, t)
        return out * scale

    def _numerical_so2_fit(self, base_array: np.ndarray) -> Callable:
        """Wrap so that I can do curve fitting."""

        def _wrapped(time_axis: np.ndarray, tau: float, scale: float) -> np.ndarray:
            return self._numerical_so2(time_axis, base_array, tau, scale)

        return _wrapped

    @staticmethod
    @nb.njit
    def _numerical_aod(
        time_axis: np.ndarray,
        so2: np.ndarray,
        tau: float,
        scale: float,
    ) -> np.ndarray:
        """Solution to AOD from SO2."""
        out = np.zeros_like(time_axis)
        for i, _ in enumerate(time_axis):
            t = time_axis[: i + 1]
            a_core = np.exp(-(t[-1] - t) / tau) * so2[: i + 1]
            out[i] = np.trapz(a_core, t)
        return out * scale

    def _numerical_aod_fit(self, base_array: np.ndarray) -> Callable:
        """Wrap so that I can do curve fitting."""

        def _wrapped(time_axis: np.ndarray, tau: float, scale: float) -> np.ndarray:
            return self._numerical_aod(time_axis, base_array, tau, scale)

        return _wrapped

    def _numerical_aod_aod(
        self, time_axis: np.ndarray, so2: np.ndarray, *params: float
    ) -> np.ndarray:
        """Compute the AOD from an intermediate AOD of SO2, and SO2, simultaneously.

        Parameters
        ----------
        time_axis : np.ndarray
            Time axis.
        so2 : np.ndarray
            SO2 data.
        *params : float
            - tau_aod_1: float, time constant for initial AOD
            - tau_aod_2: float, time constant for second AOD
            - scale: float, scale of the SO2 inside the initial AOD

        Returns
        -------
        np.ndarray
            AOD data as a function of the true AOD as a function SO2.
        """
        param_len = 3
        assert len(params) == param_len
        aod = self._numerical_aod(time_axis, so2, params[0], 1.0)
        return self._numerical_aod(time_axis, aod, params[1], params[2])

    def _numerical_aod_aod_fit(self, base_array: np.ndarray) -> Callable:
        """Wrap so that I can do curve fitting."""

        def _wrapped(
            time_axis: np.ndarray,
            tau1: float,
            tau2: float,
            scale: float,
        ) -> np.ndarray:
            return self._numerical_aod_aod(time_axis, base_array, tau1, tau2, scale)

        return _wrapped

    def _numerical_aod_rf(
        self, time_axis: np.ndarray, so2: np.ndarray, *params: float
    ) -> np.ndarray:
        """Compute the RF from AOD and SO2 simultaneously.

        Parameters
        ----------
        time_axis : np.ndarray
            Time axis.
        so2 : np.ndarray
            SO2 data.
        *params : float
            - tau_aod: float, time constant for AOD
            - scale_rf: float, scale of the RF
            - scale_aod: float, scale of the AOD

        Returns
        -------
        np.ndarray
            RF data as a function of AOD as a function SO2.
        """
        param_len = 3
        assert len(params) == param_len
        aod = self._numerical_aod(time_axis, so2, params[0], 1.0)
        return self._numerical_rf(aod, params[1], params[2])

    def _numerical_aod_rf_fit(self, base_array: np.ndarray) -> Callable:
        """Wrap so that I can do curve fitting."""

        def _wrapped(
            time_axis: np.ndarray,
            tau: float,
            scale_aod: float,
            scale_rf: float,
        ) -> np.ndarray:
            return self._numerical_aod_rf(
                time_axis, base_array, tau, scale_aod, scale_rf
            )

        return _wrapped

    @staticmethod
    @nb.njit
    def _numerical_rf(aod: np.ndarray, scale_rf: float, scale_aod: float) -> np.ndarray:
        """Solution to RF from AOD."""
        rf = scale_rf * np.log(1 + scale_aod * aod)
        return rf

    def _plot_so2(self, msg: str) -> None:
        plt.figure()
        plt.plot(self.time_axis, self.so2_true, label="Simulation output")
        plt.plot(self.time_axis, self.so2_fake, label="Numerical solution")
        _so2_str = (
            f"$\\tau_S$: {s2n(self.params_so2["tau_s"])}"
            # f", $\\tau_{{A2}}$: {s2n(self.params_so2["scale_s"])}"
        )
        t = plt.gca().transAxes
        plt.text(0.98, 0.6, _so2_str, transform=t, size=6, ha="right")
        plt.legend()
        plt.xlabel("Time [yr]")
        plt.ylabel("SO$_2$ [kg/m$^2$]")
        self._set_xlim()
        plt.savefig(_SAVE_DIR / "single" / f"numerical_so2_{self.type_}{msg}")

    def _plot_aod(self, msg: str) -> None:
        plt.figure()
        plt.plot(self.time_axis, self.aod_true, label="Simulation output")
        plt.plot(self.time_axis, self.aod_fake, label="Numerical soln (A(S))")
        plt.plot(self.time_axis, self.aod_aod_fake, label="Numerical soln (A(A(S)))")
        _aod_str = (
            f"$\\tau_A$: {s2n(self.params_aod["tau_a"])}"
            # f", $C$: {s2n(self.params_aod["scale_a"])}"
        )
        _aod_aod_str = (
            f"$\\tau_{{A1}}$: {s2n(self.params_aod_aod["tau_a1"])}"
            f", $\\tau_{{A2}}$: {s2n(self.params_aod_aod["tau_a2"])}"
            # f", $C$: {s2n(self.params_aod_aod[2])}"
        )
        t = plt.gca().transAxes
        plt.text(0.98, 0.6, _aod_str, transform=t, size=6, ha="right")
        plt.text(0.98, 0.5, _aod_aod_str, transform=t, size=6, ha="right")
        plt.legend()
        plt.xlabel("Time [yr]")
        plt.ylabel("Aerosol optical depth [1]")
        self._set_xlim()
        plt.savefig(_SAVE_DIR / "single" / f"numerical_aod_{self.type_}{msg}")

    def _plot_rf(self, msg: str) -> None:
        plt.figure()
        plt.plot(self.time_axis, self.rf_true, label="Simulation output")
        plt.plot(self.time_axis, self.rf_fake, label="Numerical soln (R(A))")
        plt.plot(self.time_axis, self.aod_rf_fake, label="Numerical soln (R(A(S)))")
        plt.plot(self.time_axis, self.rf_as_aod_fake, label="Numerical soln (A(A(S)))")
        _aod_str_rf_str = (
            f"$\\tau_A$: {s2n(self.params_aod["tau_a"])}"
            # f", $C_A$: {s2n(self.params_aod["scale_a"])}"
            # f", $C_R$: {s2n(self.params_rf[0])}"
        )
        _aod_rf_str = (
            f"$\\tau_A$: {s2n(self.params_aod_rf["tau_a"])}"
            # f", $C_A$: {s2n(self.params_aod_rf[1])}"
            # f", $C_R$: {s2n(self.params_aod_rf[2])}"
        )
        _aod_aod_str = (
            f"$\\tau_{{A1}}$: {s2n(self.params_rf_as_aod["tau_a"])}"
            f", $\\tau_{{A2}}$: {s2n(self.params_rf_as_aod["tau_r"])}"
            # f", $C$: {s2n(self.params_rf_as_aod[2])}"
        )
        t = plt.gca().transAxes
        plt.text(0.98, 0.6, _aod_str_rf_str, transform=t, size=6, ha="right")
        plt.text(0.98, 0.5, _aod_rf_str, transform=t, size=6, ha="right")
        plt.text(0.98, 0.4, _aod_aod_str, transform=t, size=6, ha="right")
        plt.legend()
        plt.xlabel("Times [yr]")
        plt.ylabel("Radiative forcing [W/m$^2$]")
        self._set_xlim()
        plt.savefig(_SAVE_DIR / "single" / f"numerical_rf_{self.type_}{msg}")

    def _set_xlim(self) -> None:
        if self.type_.name.startswith("cesm"):
            plt.xlim((-2, 15))
        else:
            plt.xlim((1210, 1350))

    def plot_available(self) -> None:
        """Plot the available data."""
        if self.reset_all_switch:
            msg = f"_optimised_from_{self.reset_all_switch}"
        else:
            msg = ""
        if self.estimate_so2:
            msg += "_fake-so2"
        else:
            msg += "_true-so2"
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
        labels.append(nsw(dec.name))
    PlotRelationship((aod, rf, labels), (aod, temp, labels), (rf, temp, labels)).plot()


def _analytic() -> None:
    a = AnalyticSolution()
    a.analytic_so2()


class PlotNumerical:
    """Plot the numerical solutions for different parameters."""

    def __init__(
        self, ns: NumericalSolver, so2_timeseries: Literal["real", "fake"] = "fake"
    ) -> None:
        self.ns = ns
        self.from_json = True
        self.so2_ts: Literal["real", "fake"] = so2_timeseries

    def plot_with_params(
        self,
        params: dict[T_Params, dict[str, float]],
        from_: str,
    ) -> None:
        """Plot the numerical parameters.

        Parameters
        ----------
        params : dict[T_Params, dict[str, float]]
            Whether to plot using parameters saved in a JSON (default) file, or with the
            provided dictionary.
        from_ : str
            The source of the parameters.
        """
        filename = f"numerical_params_{self.ns.type_}_fake-so2.json"
        with open(_SAVE_DIR / filename, encoding="locale") as f:
            params_ = json.load(f)
        for key, value in params_.items():
            value.update(params[key])
        self.ns.use_true_so2(self.so2_ts == "real")
        self.ns.reset_all(params_, from_)
        self.ns.plot_available()
        plt.show()


def _create_single_plot_from_param_files() -> None:
    so2: Literal["fake", "real"] = "real"
    for base_str in ("cesm-cesm2-strong", "ob16-ob16-month"):
        with open(
            _SAVE_DIR / f"numerical_params_{base_str}_{so2}-so2.json",
            encoding="locale",
        ) as f:
            base_params = json.load(f)
        base_params_parts = base_params.copy()
        for key, value in base_params.items():
            old = value.copy()
            for k, _ in value.items():
                if "scale" in k:
                    del old[k]
            base_params_parts[key] = old
        empty: dict[T_Params, dict] = {key: {} for key in base_params.keys()}
        name = "cesm2-strong" if "cesm" in base_str else "ob16"
        for dec in decs:
            ns = NumericalSolver(dec)
            pns = PlotNumerical(ns, so2)
            for param in (
                (empty, "self"),
                (base_params, name),
                (base_params_parts, f"{name}-parts"),
            ):
                pns.plot_with_params(param[0], param[1])


def _reset_parameter_files() -> None:
    if not rich.prompt.Confirm.ask(
        "Are you sure you want to re-set the JSON parameter files?"
        " Rememember that this task takes several hours to complete to the OB16 dataset."
    ):
        return
    for dec in decs:
        ns = NumericalSolver(dec)
        for true_so2 in (True, False):
            ns.use_true_so2(true_so2)
            ns.reset_all(save_json=True)
            ns.plot_available()
            plt.show()


def _reconstruct_ob16_analysis() -> None:
    ns = NumericalSolver(dec_ob16)
    for so2 in ("fake", "real"):
        ns.use_true_so2(so2 == "real")
        ns.reset_all()
        compare = vdd.load.TSComparison(
            # dec_ob16.rf, ns.rf_true, dec_ob16.data.aligned_arrays["so2-rf"]
            # dec_ob16.rf, ns.rf_as_aod_fake, dec_ob16.data.aligned_arrays["so2-rf"]
            dec_ob16.rf,
            ns.aod_rf_fake,
            dec_ob16.data.aligned_arrays["so2-rf"][: len(ns.time_axis)],
        )
        compare.plot_reconstructions()
        compare.peak_difference_analysis()
        compare.correlation()
        compare.spectrum()
        plt.show()
        plt.close("all")


def _plot_combiner() -> None:
    for dec in decs:
        ns = NumericalSolver(dec)
        files = {}
        for plot in ("so2", "aod", "rf"):
            files[plot] = [
                _SAVE_DIR
                / "single"
                / f"numerical_{plot}_{ns.type_}_from_self_fake-so2.png",
                _SAVE_DIR
                / "single"
                / f"numerical_{plot}_{ns.type_}_optimised_from_self_fake-so2.png",
                _SAVE_DIR
                / "single"
                / f"numerical_{plot}_{ns.type_}_optimised_from_self_true-so2.png",
                _SAVE_DIR
                / "single"
                / f"numerical_{plot}_{ns.type_}_from_self_true-so2.png",
            ]
            try:
                cosmoplots.combine(files[plot][0], files[plot][1]).in_grid(2, 1).using(
                    fontsize=50
                ).save(
                    _SAVE_DIR / "combined" / f"numerical_{plot}_{ns.type_}_combined.png"
                )
            except FileNotFoundError:
                if not (f1 := files[plot][0]).exists():
                    print(f"Cannot find file {f1}")
                if not (f2 := files[plot][1]).exists():
                    print(f"Cannot find file {f2}")
            try:
                vdd.utils.combine(files[plot][1], files[plot][2]).in_grid(2, 1).save(
                    _SAVE_DIR
                    / "combined"
                    / f"numerical_{plot}_{ns.type_}_combined_so2.png"
                )
            except FileNotFoundError:
                if not (f1 := files[plot][1]).exists():
                    print(f"Cannot find file {f1}")
                if not (f2 := files[plot][2]).exists():
                    print(f"Cannot find file {f2}")
        try:
            cosmoplots.combine(files["aod"][1], files["rf"][1]).in_grid(2, 1).using(
                fontsize=50
            ).save(
                _SAVE_DIR / "combined" / f"numerical_aod_rf_{ns.type_}_combined_so2.png"
            )
        except FileNotFoundError:
            if not (f1 := files["aod"][1]).exists():
                print(f"Cannot find file {f1}")
            if not (f2 := files["rf"][1]).exists():
                print(f"Cannot find file {f2}")
        # response = dec.response_temp_rf
        # plt.plot(ns.time_axis, np.convolve(ns.rf_fake, response, mode="same"))
        # plt.plot(ns.time_axis, dec.temp)
        # plt.show()


if __name__ == "__main__":
    # _reset_parameter_files()
    _create_single_plot_from_param_files()
    # _reconstruct_ob16_analysis()
    # _plot_combiner()
    # _scatterplot_comparison()
