"""Implementation of comparisons between OB16 and CESM2 data."""

from functools import cached_property

import matplotlib.pyplot as plt
import volcano_base

import vdd.load


class CompareOB16WithCESM:
    """Class for comparing data from Otto-Bliesner et at. (2016) with my CESM2 data."""

    @cached_property
    def dec_ob16(self) -> vdd.load.DeconvolveOB16:
        """Get the deconvolved OB16 data.

        Returns
        -------
        vdd.load.DeconvolveOB16
            An instance of the DeconvolveOB16 class.
        """
        return vdd.load.DeconvolveOB16()

    @cached_property
    def dec_cesm(self) -> vdd.load.DeconvolveCESM:
        """Get the deconvolved CESM2 data.

        Returns
        -------
        vdd.load.DeconvolveCESM
            An instance of the DeconvolveCESM class.
        """
        return vdd.load.DeconvolveCESM(pad_before=True)

    def plot_rf(self) -> plt.Figure:
        """Plot the RF data."""
        plt.figure()
        # OB16 (CESM1)
        tau = volcano_base.manipulate.dt2float(self.dec_ob16.tau)
        plt.plot(tau, self.dec_ob16.response_rf_so2)
        # CESM2
        plt.plot(self.dec_cesm.tau, self.dec_cesm.response_rf_so2)
        return plt.gcf()

    def plot_temp(self) -> plt.Figure:
        """Plot the temperature data."""
        plt.figure()
        # OB16 (CESM1)
        tau = volcano_base.manipulate.dt2float(self.dec_ob16.tau)
        plt.plot(tau, self.dec_ob16.response_temp_so2)
        # CESM2
        plt.plot(self.dec_cesm.tau, self.dec_cesm.response_temp_so2)
        return plt.gcf()
