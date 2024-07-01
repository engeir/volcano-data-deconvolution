"""Small script to obtain useful parameters for the OB16 dataset."""

import json
from functools import cached_property
from typing import TYPE_CHECKING, Literal, Self, override

import closedexpressions
import fppanalysis
import matplotlib.pyplot as plt
import numpy as np
import rich.json
import scipy.optimize
from numpy.typing import NDArray
from rich.console import Console

import vdd
import vdd.load

if TYPE_CHECKING:
    import datetime

    import xarray as xr

console = Console()


class OB16Parameters:
    """Class to obtain useful parameters for the OB16 dataset.

    Parameters
    ----------
    resolution : Literal["h0", "h1"]
    """

    def __init__(self: Self, resolution: Literal["h0", "h1"]) -> None:
        self.data: vdd.load.DeconvolveOB16 = vdd.load.DeconvolveOB16(data=resolution)

    @override
    def __repr__(self: Self) -> str:
        """Return the representation of the class."""
        info: dict[str, str] = {
            f"{val[7:]}": getattr(self, val)
            for val in dir(self)
            if val.startswith("_param")
        }
        return json.dumps(info, indent=2)

    def __rich__(self: Self) -> rich.json.JSON:
        """Return a Rich representation of the class."""
        return rich.json.JSON(str(self), indent=2, sort_keys=True)

    @property
    def _param_delta_t(self: Self) -> str:
        """Return the sampling time."""
        time: xr.DataArray = self.data.temp.time
        delta_t: datetime.timedelta = time.to_numpy()[1] - time.to_numpy()[0]
        return f"{delta_t.total_seconds() / 3600 / 24 / 365} [yr] ({delta_t})"

    @property
    def _param_epsilon(self: Self) -> str:
        """Return the noise-to-signal ratio."""
        return str(self._gamma_epsilon[1])

    @property
    def _param_gamma(self: Self) -> str:
        """Return the intermittency parameter."""
        return str(self._gamma_epsilon[0])

    @cached_property
    def _gamma_epsilon(self: Self) -> tuple[float, float]:
        with console.status(
            "[bold yellow]Calculating gamma and epsilon ...",
            spinner="point",
        ):
            temp = (
                self.data.temp.to_numpy() - self.data.temp.to_numpy().mean()
            ) / self.data.temp.to_numpy().std()
            out = self._find_gamma_epsilon(temp)

        pdf, _, x = fppanalysis.distribution(temp, 30)
        pdf2 = closedexpressions.noisy_shot_noise(x, *out)
        plt.plot(x, pdf)
        plt.plot(x, pdf2)
        plt.show()
        return out

    @staticmethod
    def _pdf_from_gamma_epsilon(
        bins: NDArray[np.float64],
        gamma: float,
        epsilon: float,
    ) -> NDArray[np.float64]:
        """Return the PDF for the bins as a function of gamma and epsilon."""
        return closedexpressions.noisy_shot_noise(bins, gamma, epsilon)

    def _find_gamma_epsilon(
        self: Self,
        temp: NDArray[np.float64],
    ) -> tuple[float, float]:
        temp = (temp - temp.mean()) / temp.std()
        pdf, _, x = fppanalysis.distribution(temp, 30)
        params, _ = scipy.optimize.curve_fit(
            self._pdf_from_gamma_epsilon,
            x,
            pdf,
            p0=[0.1, 0.1],
        )
        return params

    @property
    def _param_num_events(self: Self) -> str:
        """Return the total number of volcanic events in the dataset."""
        so2: NDArray[np.float64] = self.data.so2.to_numpy()
        zero_like = 1e-3
        so2_peaks = np.where(so2 > zero_like)[0]
        num_events = len(so2_peaks)
        return str(num_events)

    @cached_property
    def _param_tau_d(self: Self) -> str:
        """Return the duration time of the volcanic events."""
        with console.status(
            "[bold yellow]Calculating the duration time ...",
            spinner="point",
        ):
            peak_window = (self.data.tau > -2) & (self.data.tau < 20)  # noqa: PLR2004
            tau = self.data.tau[peak_window]
            response = self.data.response_temp_so2[peak_window]
            max_val = response.max()
            idx = response > max_val / np.e
        plt.plot(tau, response)
        plt.plot(tau[idx], response[idx])
        plt.axvline(tau[idx][-1], color="k", linestyle="--")
        plt.text(
            tau[idx][-1],
            response[idx][-1],
            f"{tau[idx][-1]:.2f}",
        )
        plt.show()
        return f"{tau[idx][-1]:.2f} [yr]"

    @property
    def _param_total_time(self: Self) -> str:
        """Return the total time of the dataset."""
        time: xr.DataArray = self.data.temp.time
        total_time = (time.data[-1] - time.data[0]).total_seconds() / 3600 / 24 / 365
        return f"{total_time} [yr] (start: {time.data[0]}, end: {time.data[-1]})"


if __name__ == "__main__":
    month = OB16Parameters(resolution="h0")
    day = OB16Parameters(resolution="h1")
    console.print(month, soft_wrap=True)
    console.print(day)
