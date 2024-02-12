"""Module for deconvolving data from Otto-Bliesner et at. (2016)."""

import fppanalysis
import matplotlib.pyplot as plt
import volcano_base


def main() -> None:
    """Deconvolve data."""
    temperature = volcano_base.load.get_ob16_temperature()
    rf = volcano_base.load.get_ob16_rf()
    _, so2 = volcano_base.load.get_so2_ob16_full_timeseries()
    t_so2, err_t_so2 = fppanalysis.deconvolution_methods.RL_gauss_deconvolve(
        temperature, so2, 200
    )
    t_rf, err_t_rf = fppanalysis.deconvolution_methods.RL_gauss_deconvolve(
        temperature, rf, 200
    )
    print(t_so2)
    plt.plot(t_so2)
    plt.plot(t_rf)
    plt.show()


if __name__ == "__main__":
    main()
