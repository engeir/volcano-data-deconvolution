import pathlib
import vdd.load
from _typeshed import Incomplete as Incomplete

DataCESM = vdd.load.CESMData
DecCESM = vdd.load.DeconvolveCESM
padding: Incomplete
dec_cesm_m: Incomplete
dec_ob16_month: Incomplete

class PlotCutOff:
    cut_offs: Incomplete
    def __init__(self, *cut_offs: vdd.load.CutOff) -> None: ...
    def call_cut_offs(self, method: str, *args, **kwargs) -> None: ...
    def plot(self, save_path: pathlib.Path | None = None) -> None: ...
