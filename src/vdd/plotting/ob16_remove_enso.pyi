import xarray as xr

def combine_times(first: xr.DataArray, second: xr.DataArray) -> xr.DataArray: ...
def remove_quadratic_fit(da: xr.DataArray) -> xr.DataArray: ...
def remove_all_freqs(da: xr.DataArray) -> xr.DataArray: ...
def select_enso_area(da: xr.DataArray) -> xr.DataArray: ...
def regress_a_on_b(a: xr.DataArray, b: xr.DataArray) -> xr.DataArray: ...
def _compute_enso() -> xr.DataArray: ...
def main() -> None: ...
