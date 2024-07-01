import vdd.load
from _typeshed import Incomplete as Incomplete

DataCESM = vdd.load.CESMData
DecCESM = vdd.load.DeconvolveCESM
dec_cesm_4sep: Incomplete
dec_cesm_2sep: Incomplete
dec_cesm_m4: Incomplete
dec_cesm_m2: Incomplete
dec_cesm_e: Incomplete
dec_cesm_s: Incomplete
dec_cesm_p: Incomplete
dec_cesm_m: Incomplete
dec_ob16_month: Incomplete
all_decs: Incomplete

class PlotParametrisation:
    decs: Incomplete
    results: Incomplete
    def __init__(self, *decs: vdd.load.Deconvolve) -> None: ...
    def strategy(self) -> None: ...
    def plot_simulations(self) -> None: ...
    def plot_method_comparisons(self) -> None: ...

def main() -> None: ...
