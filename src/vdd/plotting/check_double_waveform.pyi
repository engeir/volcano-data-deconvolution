import vdd.load
from _typeshed import Incomplete as Incomplete
from typing import Literal

DataCESM = vdd.load.CESMData
DecCESM = vdd.load.DeconvolveCESM
padding: vdd.load.T_Padding
dec_cesm_p4: Incomplete
dec_cesm_p2: Incomplete
dec_cesm_m4: Incomplete
dec_cesm_m2: Incomplete
dec_cesm_e: Incomplete
dec_cesm_s: Incomplete
dec_cesm_p: Incomplete
dec_cesm_m: Incomplete

def check_waveform_responses(*decs: vdd.load.DeconvolveCESM) -> None: ...
def curve_fit_aod(aod, alpha, beta) -> None: ...

class CheckRecreatedWaveforms:
    single_waveform: Incomplete
    decs: Incomplete
    scale_by_aod: Incomplete
    keys: Incomplete
    def __init__(self, /, *decs: vdd.load.DeconvolveCESM, single_waveform: vdd.load.DeconvolveCESM | None = None, scale_by_aod: Literal['log', 'log-inside', 'root'] | bool = False, keys: dict[str, tuple[int, int]] | None = None) -> None: ...
    def run_loop(self) -> None: ...

def main() -> None: ...
