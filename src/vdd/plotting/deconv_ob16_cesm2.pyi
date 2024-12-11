import matplotlib as mpl
import pathlib
import vdd.load
from _typeshed import Incomplete as Incomplete

DataCESM = vdd.load.CESMData
DecCESM = vdd.load.DeconvolveCESM
console: Incomplete
padding: vdd.load.T_Padding
dec_cesm_int4: Incomplete
dec_cesm_int2: Incomplete
dec_cesm_med4: Incomplete
dec_cesm_med2: Incomplete
dec_cesm_e: Incomplete
dec_cesm_s: Incomplete
dec_cesm_p: Incomplete
dec_cesm_m: Incomplete
dec_ob16: Incomplete
all_decs: Incomplete

class PlotResponseFunctions:
    decs: Incomplete
    norm: Incomplete
    def __init__(self, *decs: vdd.load.Deconvolve, norm: bool = False) -> None: ...
    def plot_rf_so2(self, fig: mpl.figure.Figure | None = None, save_as: str = 'rf-so2') -> mpl.figure.Figure: ...
    def plot_rf_so2_decay(self, fig: mpl.figure.Figure | None = None, save_as: str = 'rf-so2_decay') -> mpl.figure.Figure: ...
    def plot_temp_so2(self, fig: mpl.figure.Figure | None = None, save_as: str = 'temp-so2') -> mpl.figure.Figure: ...
    def plot_temp_rf(self, fig: mpl.figure.Figure | None = None, save_as: str = 'temp-rf') -> mpl.figure.Figure: ...
    @staticmethod
    def plot_grayscale_highlight(fig: mpl.figure.Figure | None = None, save_as: str = 'temp-so2-gs', save_path: pathlib.Path | None = None) -> mpl.figure.Figure: ...
    def grayscale_plot(self, ax: mpl.axes.Axes, res_name: str, i_: int) -> None: ...
    def run(self, save_as: list[str] | None = None) -> None: ...

def main() -> None: ...
