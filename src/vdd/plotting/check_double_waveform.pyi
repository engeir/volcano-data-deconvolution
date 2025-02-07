import matplotlib as mpl
import numpy as np
import pathlib
import vdd.load
from _typeshed import Incomplete as Incomplete
from numpy.typing import NDArray as NDArray
from typing import Literal, overload

DataCESM = vdd.load.CESMData
DecCESM = vdd.load.DeconvolveCESM
padding: Incomplete
dec_cesm_p4: Incomplete
dec_cesm_p2: Incomplete
dec_cesm_m4: Incomplete
dec_cesm_m2: Incomplete
dec_cesm_e: Incomplete
dec_cesm_s: Incomplete
dec_cesm_p: Incomplete
dec_cesm_m: Incomplete

def check_waveform_responses(*decs: vdd.load.DeconvolveCESM) -> None: ...
def curve_fit_aod(aod: NDArray[np.float64], alpha: float, beta: float) -> NDArray[np.float64]: ...

class CheckRecreatedWaveforms:
    single_waveform: Incomplete
    decs: Incomplete
    scale_by_aod: Incomplete
    keys: Incomplete
    def __init__(self, /, *decs: vdd.load.DeconvolveCESM, single_waveform: vdd.load.DeconvolveCESM | None = None, scale_by_aod: Literal['log', 'log-inside', 'root'] | bool = False, keys: dict[str, tuple[int, int]] | None = None) -> None: ...
    @overload
    def run_loop(self, save_path: pathlib.Path | None, *, return_: Literal[True]) -> tuple[mpl.figure.Figure, list[mpl.axes.Axes], pathlib.Path]: ...
    @overload
    def run_loop(self, save_path: pathlib.Path | None, *, return_: Literal[False]) -> None: ...

def main() -> None: ...
