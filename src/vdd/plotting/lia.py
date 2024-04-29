"""Analysis to see if the little ice age is reproduced by the model."""

# TODO:
# - [ ] Grab data from OB16 of the relevant period
# - [ ] Convolve forcing data with responses from CESM2 single eruption simulations
# - [ ] Compare with a simple gamma distribution with similar tail as the temperatures

import warnings

import volcano_base

import vdd.load

warnings.warn(
    "This script is deprecated and may be removed in the future.",
    category=DeprecationWarning,
    stacklevel=1,
)


class LittleIceAge:
    """Analysis of the little ice age."""

    def __init__(self, dec: vdd.load.DeconvolveCESM) -> None:
        self.ob16 = volcano_base.load.OttoBliesner(freq="h0", progress=True)
        self.dec = dec
