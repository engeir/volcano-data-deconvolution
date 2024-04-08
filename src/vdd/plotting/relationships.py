"""Plot the relationships between AOD, RF, T and more."""

from collections.abc import Callable
from typing import Literal

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
_MAXFEV = 0

plt.rc("text.latex", preamble=r"\usepackage{amsmath}")
plt.style.use(
    [
        "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle",
        "vdd.extra",
        {"legend.fontsize": 6},
    ]
)
Params_T = Literal["SO2", "AOD", "AOD-AOD", "AOD-RF", "RF"]
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

# cesm-cesm2-double-overlap parameters, fake so2
param_dict_4sep_fake: dict[Params_T, tuple[float, ...]] = {
    "SO2": (0.20648891265424113, 14.33051338439839),
    "AOD": (1.3701140219220893, 0.05344061819620313),
    "AOD-AOD": (
        0.28111196173184444,
        0.59365187011508,
        0.8667623279669701,
        0.4679529190185629,
    ),
    "AOD-RF": (
        0.7728077201434438,
        0.31837600012104006,
        20.783872003959562,
        0.3193658567153462,
    ),
    "RF": (18.50045044, 2.53849499),
}
# cesm-cesm2-double-overlap parameters, true so2
param_dict_4sep_true: dict[Params_T, tuple[float, ...]] = {
    "SO2": (0.20648891265424113, 14.33051338439839),
    "AOD": (1.4818306343193342, 0.04986290216778253),
    "AOD-AOD": (
        0.967108417490067,
        0.5469632323924328,
        0.24464680182046572,
        0.5238136820461587,
    ),
    "AOD-RF": (
        0.6523297356943707,
        0.5293992625811754,
        13.920698910261406,
        0.5532246284318106,
    ),
    "RF": (18.50045044, 2.53849499),
}

# cesm-cesm2-tt-2sep parameters, fake so2
param_dict_2sep_fake: dict[Params_T, tuple[float, ...]] = {
    "SO2": (0.2167868651072701, 13.53308031600478),
    "AOD": (1.203125151744134, 0.052691850392018943),
    "AOD-AOD": (
        0.49212675714195747,
        0.4716115580822787,
        0.5021080703992524,
        0.5119517599440576,
    ),
    "AOD-RF": (
        0.5685344215329661,
        0.5275774262454338,
        13.829249312705794,
        0.5283797893316774,
    ),
    "RF": (200.300474548285, 0.05311253870584331),
}
# cesm-cesm2-tt-2sep parameters, true so2
param_dict_2sep_true: dict[Params_T, tuple[float, ...]] = {
    "SO2": (0.2167868651072701, 13.53308031600478),
    "AOD": (1.3381979458180546, 0.048768716856846264),
    "AOD-AOD": (
        0.29439053899882894,
        0.5290502104203028,
        0.853603593442637,
        0.4640462596879809,
    ),
    "AOD-RF": (
        0.6289972280314285,
        -0.4227960794134276,
        16.01088137252591,
        -0.42435068061118325,
    ),
    "RF": (200.300474548285, 0.05311253870584331),
}

# cesm-cesm2-size5000 parameters, fake so2
param_dict_e_fake: dict[Params_T, tuple[float, ...]] = {
    "SO2": (0.49929255956304486, 13.666770936658683),
    "AOD": (0.8374064644849375, 0.019674030248830014),
    "AOD-AOD": (
        0.3903728505112651,
        0.3494130112162795,
        0.31730245658123996,
        0.3448631214287202,
    ),
    "AOD-RF": (
        0.3790555690876411,
        -0.13255768061180853,
        28.145013986772078,
        -0.1310272783432714,
    ),
    "RF": (18.50045044, 2.53849499),
}
# cesm-cesm2-size5000 parameters, true so2
param_dict_e_true: dict[Params_T, tuple[float, ...]] = {
    "SO2": (0.49929255956304486, 13.666770936658683),
    "AOD": (0.8907361662029977, 0.019143680887591163),
    "AOD-AOD": (
        0.5536372838129326,
        0.3611276179620605,
        0.2173158438760093,
        0.3560047476936677,
    ),
    "AOD-RF": (
        0.43208961426821774,
        0.13856045895467978,
        25.64485372331202,
        0.14032258394251887,
    ),
    "RF": (495.38696139905596, 0.008917135191374624),
}

# cesm-cesm2-strong parameters, fake so2
param_dict_s_fake: dict[Params_T, tuple[float, ...]] = {
    "SO2": (0.3751117448661117, 13.824922745411444),
    "AOD": (0.9945571186419672, 0.028828278613217097),
    "AOD-AOD": (
        0.40939242347216864,
        0.3919386260669227,
        0.409769055994987,
        0.3924878432544051,
    ),
    "AOD-RF": (
        0.531807911246878,
        0.15641879300382192,
        29.529072144938482,
        0.1594303090956324,
    ),
    "RF": (18.50045044, 2.53849499),
}
# cesm-cesm2-strong parameters, true so2
param_dict_s_true: dict[Params_T, tuple[float, ...]] = {
    "SO2": (0.3751117448661117, 13.824922745411444),
    "AOD": (1.0535013138978253, 0.02837383155494114),
    "AOD-AOD": (
        0.4264535902986881,
        0.3825137942433982,
        0.42710806772921195,
        0.38242429166124947,
    ),
    "AOD-RF": (
        0.5404768768976129,
        0.19056778896953694,
        25.37958322184627,
        0.18740378955177064,
    ),
    "RF": (18.50045044, 2.53849499),
}

# cesm-cesm2-medium-plus parameters, fake so2
param_dict_p_fake: dict[Params_T, tuple[float, ...]] = {
    "SO2": (0.24053102992514483, 14.154608712216556),
    "AOD": (1.3189023376025042, 0.047756232254152625),
    "AOD-AOD": (
        0.3384112260966441,
        0.5087015664959406,
        0.7291531985922045,
        0.4469191714817084,
    ),
    "AOD-RF": (
        0.7899896297431088,
        -0.2760069396404453,
        22.070795441096834,
        -0.27312241845184915,
    ),
    "RF": (18.50045044, 2.53849499),
}
# cesm-cesm2-medium-plus parameters, true so2
param_dict_p_true: dict[Params_T, tuple[float, ...]] = {
    "SO2": (0.24053102992514483, 14.154608712216556),
    "AOD": (1.3457829280226425, 0.04842795556722085),
    "AOD-AOD": (
        0.38698679506121814,
        0.4303465386484054,
        0.6794926156127085,
        0.5103959724549507,
    ),
    "AOD-RF": (
        0.7633365441277841,
        0.3142637772230818,
        19.7639144463403,
        0.3162380155223765,
    ),
    "RF": (18.50045044, 2.53849499),
}

# cesm-cesm2-medium parameters, fake so2
param_dict_m_fake: dict[Params_T, tuple[float, ...]] = {
    "SO2": (0.15205927660923824, 13.7343385672187),
    "AOD": (1.493642882429071, 0.0911080785832632),
    "AOD-AOD": (
        0.19122776892038992,
        0.6632526511506003,
        1.1240849391562433,
        0.8746048676041859,
    ),
    "AOD-RF": (
        0.872206432128419,
        0.6620391968401917,
        7.717841579740749,
        0.6632777677979476,
    ),
    "RF": (18.50045044, 2.53849499),
}
# cesm-cesm2-medium parameters, true so2
param_dict_m_true: dict[Params_T, tuple[float, ...]] = {
    "SO2": (0.15205927660923824, 13.7343385672187),
    "AOD": (1.5419746250465045, 0.09407498843543649),
    "AOD-AOD": (
        0.21903393729028367,
        0.8362481216678737,
        1.1058241798998705,
        0.653546777856078,
    ),
    "AOD-RF": (
        0.8539273546067695,
        0.8254226441630498,
        5.956262267697095,
        0.8254328108119243,
    ),
    "RF": (18.50045044, 2.53849499),
}

# ob16-ob16-month parameters, fake so2
param_dict_ob16_fake: dict[Params_T, tuple[float, ...]] = {
    "SO2": (1.0192882358961923, 12.08452750238792),
    "AOD": (0.99455712, 0.02882828),
    "AOD-AOD": (0.40939242, 0.39193863, 0.40976906, 0.39248784),
    "AOD-RF": (
        0.1821773347875807,
        0.03523696921973355,
        521.4159485898647,
        0.03506365160959408,
    ),
    "RF": (18.50045044, 2.53849499),
}
# ob16-ob16-month parameters, true so2
param_dict_ob16_true: dict[Params_T, tuple[float, ...]] = {
    "SO2": (1.0192882358961923, 12.08452750238792),
    "AOD": (0.99455712, 0.02882828),
    "AOD-AOD": (0.40939242, 0.39193863, 0.40976906, 0.39248784),
    "AOD-RF": (
        0.11168491836166736,
        0.17293165085203288,
        37.79409160216243,
        0.17295461630506495,
    ),
    "RF": (1599.2742610451014, 0.0026568467327948895),
}

decs = (dec_4sep, dec_2sep, dec_e, dec_s, dec_p, dec_m, dec_ob16)
# decs = (dec_4sep, dec_2sep, dec_e, dec_s, dec_p, dec_m)
# decs = (dec_m,)
param_dicts_fake = (
    param_dict_4sep_fake,
    param_dict_2sep_fake,
    param_dict_e_fake,
    param_dict_s_fake,
    param_dict_p_fake,
    param_dict_m_fake,
    param_dict_ob16_fake,
)
param_dicts_true = (
    param_dict_4sep_true,
    param_dict_2sep_true,
    param_dict_e_true,
    param_dict_s_true,
    param_dict_p_true,
    param_dict_m_true,
    param_dict_ob16_true,
)


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
        self.estimate_so2 = True
        self.params_so2 = (0.37511174, 13.82492275)
        # self.params_aod = (0.40676405, 4.07634108)
        self.params_aod = (0.99455712, 0.02882828)
        self.params_aod_aod = (0.40939242, 0.39193863, 0.40976906, 0.39248784)
        self.params_aod_rf = (0.53180791, 0.15641879, 29.52907214, 0.15943031)
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
        self.so2_true = np.roll(so2_ob16.data, 0)  # * 1e-6
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
            return self.numerical_so2(
                self.time_axis, self.delta_pulses, *self.params_so2
            )
        return self.so2_true

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
            params_so2, _ = curve_fit(
                self.numerical_so2_fit(dp), ta, st, maxfev=_MAXFEV
            )
        except RuntimeError as e:
            print(e)
            print("Using previous parameters")
            params_so2 = self.params_so2
        self.params_so2 = tuple(params_so2)

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
                self.numerical_aod_fit(sf), ta, at, maxfev=_MAXFEV
            )
        except RuntimeError as e:
            print(e)
            print("Using previous parameters")
            params_aod = self.params_aod
        self.params_aod = tuple(params_aod)

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
                self.numerical_aod_aod_fit(sf), ta, at, maxfev=_MAXFEV
            )
        except RuntimeError as e:
            print(e)
            print("Using previous parameters")
            params_aod_aod = self.params_aod_aod
        self.params_aod_aod = tuple(params_aod_aod)

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
                self.numerical_aod_rf_fit(sf), ta, rt, maxfev=_MAXFEV
            )
        except RuntimeError as e:
            print(e)
            print("Using previous parameters")
            params_aod_rf = self.params_aod_rf
        self.params_aod_rf = tuple(params_aod_rf)

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
            params_rf, _ = curve_fit(self.numerical_rf, af, rt, maxfev=_MAXFEV)
        except RuntimeError as e:
            print(e)
            print("Using previous parameters")
            params_rf = self.params_rf
        self.params_rf = tuple(params_rf)

    def reset_all(
        self, custom_params: dict[Params_T, tuple[float, ...]] | None = None
    ) -> None:
        """Reset all parameters."""
        self.reset_all_switch = True
        if custom_params is not None:
            self._reset_all_custom(custom_params)
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
        print("Resetting RF parameters")
        self.reset_params_rf()

    def _reset_all_custom(
        self, custom_params: dict[Params_T, tuple[float, ...]]
    ) -> None:
        """Reset all parameters with custom parameters."""
        for param, values in custom_params.items():
            match param, values:
                case "SO2", (tau, scale):
                    self.params_so2 = (tau, scale)
                case "AOD", (tau, scale):
                    self.params_aod = (tau, scale)
                case "AOD-AOD", (tau1, scale1, tau2, scale2):
                    self.params_aod_aod = (tau1, scale1, tau2, scale2)
                case "AOD-RF", (tau, scale_aod, scale_rf1, scale_rf2):
                    self.params_aod_rf = (tau, scale_aod, scale_rf1, scale_rf2)
                case "RF", (scale_a, scale_b):
                    self.params_rf = (scale_a, scale_b)
                case _:
                    raise ValueError("Invalid input.")

    def print_params(self) -> None:
        """Print the parameters."""
        print(f"{self.type_} parameters")
        print(f"SO2: {self.params_so2}")
        print(f"AOD: {self.params_aod}")
        print(f"AOD-AOD: {self.params_aod_aod}")
        print(f"AOD-RF: {self.params_aod_rf}")
        print(f"RF: {self.params_rf}")

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
        self._set_xlim()
        plt.savefig(_SAVE_DIR / f"numerical_so2_{self.type_}{msg}")

    def _plot_aod(self, msg: str) -> None:
        plt.figure()
        plt.plot(self.time_axis, self.aod_true, label="Simulation output")
        plt.plot(self.time_axis, self.aod_fake, label="Numerical soln (A(S))")
        plt.plot(self.time_axis, self.aod_aod_fake, label="Numerical soln (A(A(S)))")
        plt.legend()
        plt.xlabel("Time [yr]")
        plt.ylabel("Aerosol optical depth [1]")
        self._set_xlim()
        plt.savefig(_SAVE_DIR / f"numerical_aod_{self.type_}{msg}")

    def _plot_rf(self, msg: str) -> None:
        plt.figure()
        plt.plot(self.time_axis, self.rf_true, label="Simulation output")
        plt.plot(self.time_axis, self.rf_fake, label="Numerical soln (R(A))")
        plt.plot(self.time_axis, self.aod_rf_fake, label="Numerical soln (R(A(S)))")
        plt.legend()
        plt.xlabel("Time [yr]")
        plt.ylabel("Radiative forcing [W/m$^2$]")
        self._set_xlim()
        plt.savefig(_SAVE_DIR / f"numerical_rf_{self.type_}{msg}")

    def _set_xlim(self) -> None:
        if self.type_.name.startswith("cesm"):
            plt.xlim((-2, 15))
        else:
            plt.xlim((1210, 1350))

    def plot_available(self) -> None:
        """Plot the available data."""
        if self.reset_all_switch:
            msg = "_optimised"
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
        labels.append(dec.name)
    PlotRelationship((aod, rf, labels), (aod, temp, labels), (rf, temp, labels)).plot()


def _analytic() -> None:
    a = AnalyticSolution()
    a.analytic_so2()


def _numerical_solver() -> None:
    for dec, p_fake, p_true in zip(
        decs, param_dicts_fake, param_dicts_true, strict=True
    ):
        ns = NumericalSolver(dec)
        ns.plot_available()
        ns.reset_all(p_fake)
        ns.plot_available()
        ns.use_true_so2(True)
        ns.reset_all(p_true)
        ns.plot_available()
        for plot in ("so2", "aod", "rf"):
            f1 = _SAVE_DIR / f"numerical_{plot}_{ns.type_}_fake-so2.png"
            f2 = _SAVE_DIR / f"numerical_{plot}_{ns.type_}_optimised_fake-so2.png"
            f3 = _SAVE_DIR / f"numerical_{plot}_{ns.type_}_optimised_true-so2.png"
            f4 = _SAVE_DIR / f"numerical_{plot}_{ns.type_}_true-so2.png"
            try:
                cosmoplots.combine(f1, f2).in_grid(2, 1).using(fontsize=50).save(
                    _SAVE_DIR / f"numerical_{plot}_{ns.type_}_combined.png"
                )
            except FileNotFoundError:
                pass
            try:
                cosmoplots.combine(f2, f3).in_grid(2, 1).using(fontsize=50).save(
                    _SAVE_DIR / f"numerical_{plot}_{ns.type_}_combined_so2.png"
                )
            except FileNotFoundError:
                pass
            f1.unlink(missing_ok=True)
            f2.unlink(missing_ok=True)
            f3.unlink(missing_ok=True)
            f4.unlink(missing_ok=True)
        # plt.show()
        plt.close("all")


if __name__ == "__main__":
    # for dec in decs:
    #     ns = NumericalSolver(dec)
    #     ns.plot_available()
    #     # ns.use_true_so2(True)
    #     # ns.reset_all()
    #     # ns.plot_available()
    #     plt.show()
    _numerical_solver()
    # _scatterplot_comparison()
