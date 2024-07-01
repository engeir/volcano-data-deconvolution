"""Module for deconvolving data from Otto-Bliesner et at. (2016)."""

import matplotlib.pyplot as plt
import nc_time_axis  # noqa: F401
import numpy as np
import volcano_base
import xarray as xr

import vdd.compare
import vdd.load
import vdd.utils

plt.rc("text.latex", preamble=r"\usepackage{amsmath}")
plt.style.use(
    "https://raw.githubusercontent.com/uit-cosmo/cosmoplots/main/cosmoplots/default.mplstyle",
)


def check_cesm_output() -> None:
    """Check the CESM2 output data."""
    cesm = vdd.load.CESMData()
    x = cesm.aod.time.data
    if any(x != cesm.rf.time.data) or any(x != cesm.temp.time.data):
        msg = "Times do not match."
        raise ValueError(msg)
    plt.figure()
    vdd.utils.normalise(cesm.aod).plot()  # type: ignore[call-arg]
    vdd.utils.normalise(cesm.rf).plot()  # type: ignore[call-arg]
    vdd.utils.normalise(cesm.temp).plot()  # type: ignore[call-arg]
    # plt.plot(x, vdd.utils.normalise(cesm.aod))
    # plt.plot(x, vdd.utils.normalise(cesm.rf))
    # plt.plot(x, vdd.utils.normalise(cesm.temp))
    plt.show()


def extend_temp_shape() -> xr.DataArray:
    """Extend the shape of the temperature array."""
    t = vdd.load.CESMData().temp
    t *= -1
    count = 1
    x_base = np.concatenate((np.array([0]), t.time.data))
    if len(x_base) != int(12 * 20):
        raise ValueError
    x = x_base + 20 * count
    decay_year = 200
    while x[-1] < decay_year:
        count += 1
        x = np.append(x, x_base + 20 * count)
    decay = 1.5 * np.exp(-x / 30)
    decay += np.random.default_rng().normal(0, 0.1, len(x))
    signal = np.concatenate((t, decay))
    return xr.DataArray(
        -1 * signal,
        coords={"time": np.concatenate((t.time.data, x))},
        dims=["time"],
    )


def deconvolve_rf_and_temp_with_so2() -> None:
    """Deconvolve OB16 time series."""
    dec_class = vdd.load.DeconvolveOB16()
    dec_class.plot_dec_rf_with_so2()
    dec_class.plot_dec_temp_with_so2()
    plt.show()


def compare_ob16_with_cesm() -> None:
    """Compare OB16 with CESM2."""
    comp_class = vdd.compare.CompareOB16WithCESM()
    rf = comp_class.plot_rf()
    temp = comp_class.plot_temp()
    rf.gca().set_xlim((-5, 20))
    rf.savefig("rf.png")
    temp.savefig("temp.png")
    temp.gca().set_xlim((-5, 20))
    plt.show()


def _look_at_ob_data() -> None:
    """Deconvolve data."""
    # volcano_base.load.get_ob16_outputs()
    ob16_day = volcano_base.load.OttoBliesner(freq="h1", progress=True)
    temperature_day = ob16_day.temperature_median
    temperature = vdd.utils.weighted_monthly_avg(temperature_day)
    rf_day = ob16_day.rf_median
    rf = vdd.utils.weighted_monthly_avg(rf_day)
    so2_day = ob16_day.so2_delta
    # d1 = 190
    # so2_day = so2_day.assign_coords(
    #     time=so2_day.time.data - datetime.timedelta(days=d1)
    # )
    so2 = vdd.utils.weighted_monthly_avg(so2_day)

    so2, rf_fr, temp = xr.align(so2, rf, temperature)

    cesm_month = vdd.load.CESMData()
    frc = cesm_month.rf
    tmp = cesm_month.temp
    tmp_ext = vdd.utils.pad_before_convolution(extend_temp_shape())
    signal = np.convolve(so2, frc, mode="same")
    signal2 = np.convolve(so2, tmp, mode="same")
    signal3 = np.convolve(so2, tmp_ext, mode="same")
    plt.figure()
    vdd.utils.normalise(temp).plot()  # type: ignore[call-arg]
    plt.plot(so2.time.data, vdd.utils.normalise(signal3))
    plt.plot(so2.time.data, vdd.utils.normalise(signal2))
    plt.legend(["OB16 T", "CESM2 T extra", "CESM2 T"])
    plt.figure()
    vdd.utils.normalise(rf_fr).plot()  # type: ignore[call-arg]
    plt.plot(so2.time.data, vdd.utils.normalise(signal))
    plt.show()


if __name__ == "__main__":
    # deconvolve_rf_and_temp_with_so2()
    compare_ob16_with_cesm()
