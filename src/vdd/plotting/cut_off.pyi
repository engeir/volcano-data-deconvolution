import vdd.load
from _typeshed import Incomplete as Incomplete

DataCESM = vdd.load.CESMData
DecCESM = vdd.load.DeconvolveCESM
padding: vdd.load.T_Padding
dec_cesm_m: Incomplete
dec_ob16_month: Incomplete

class PlotCutOff:
    cut_offs: Incomplete
    def __init__(self, *cut_offs: vdd.load.CutOff) -> None: ...
    def call_cut_offs(self, method: str, *args, **kwargs) -> None: ...
    def plot(self, remove_grid_parts: bool = True) -> None: ...
