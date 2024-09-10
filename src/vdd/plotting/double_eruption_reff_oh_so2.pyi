import xarray as xr
from _typeshed import Incomplete as Incomplete
from typing import NamedTuple
from volcano_base.load import FindFiles as FindFiles

COLORS: Incomplete

class PlotArgs(NamedTuple):
    ls: str
    c: str
    lab: str

class ReffPlot:
    FINDER: Incomplete
    sad_: Incomplete
    reff_: Incomplete
    temp_: Incomplete
    def print(self) -> None: ...
    def ens2median(self, arr: list[xr.DataArray]) -> xr.DataArray: ...
    def compute(self) -> xr.Dataset: ...
    def load(self) -> xr.Dataset: ...
    def plot(self) -> None: ...

class OHPlot:
    oh_c: FindFiles
    oh_m: FindFiles
    oh_p: FindFiles
    oh_s: FindFiles
    oh_e: FindFiles
    oh_m2: FindFiles
    oh_m4: FindFiles
    oh_p2: FindFiles
    oh_p4: FindFiles
    def ens2median(self, arr: list[xr.DataArray]) -> xr.DataArray: ...
    def print_available(self) -> None: ...
    def print(self) -> None: ...
    def compute(self) -> xr.Dataset: ...
    def load(self) -> xr.Dataset: ...
    def plot(self) -> None: ...

class SO2BurdenPlot:
    FINDER: Incomplete
    so2_c: Incomplete
    so2_m: Incomplete
    so2_m2: Incomplete
    so2_m4: Incomplete
    so2_p: Incomplete
    so2_s: Incomplete
    so2_e: Incomplete
    so2_p2: Incomplete
    so2_p4: Incomplete
    @staticmethod
    def ens2median(arr: list[xr.DataArray]) -> xr.DataArray: ...
    def print(self) -> None: ...
    def compute(self) -> xr.Dataset: ...
    def load(self) -> xr.Dataset: ...
    def plot(self) -> None: ...
