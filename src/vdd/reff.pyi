import xarray as xr
from _typeshed import Incomplete as Incomplete

ILEV: Incomplete
LEV: Incomplete

class Reff:
    reff: Incomplete
    temp: Incomplete
    sad: Incomplete
    h: Incomplete
    coeff: Incomplete
    def __init__(self, reff: xr.DataArray, temp: xr.DataArray, sad: xr.DataArray) -> None: ...
    def calculate_reff(self) -> xr.DataArray: ...
