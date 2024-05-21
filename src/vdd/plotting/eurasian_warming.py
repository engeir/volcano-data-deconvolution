"""Script that checks is the Eruasian warming pattern is present in the data.

From Schneider et al. (2009), they refer to Shindell et al. (2003), who found a
difference between 2P and 3P volcanic events (n times Pinatubo). 2P showed a Eurasian
winter warming, 3P did not. Schneider et al. (2009) finds the same pattern, where
Samalas 1258 has no pattern, while composites of three large but comparatively small
eruptions do show such a winter warming pattern.
"""

import warnings

import cartopy.crs as ccrs
import cmcrameri.cm as cmc
import matplotlib.pyplot as plt
import numpy as np
import volcano_base

import vdd.load

warnings.warn(
    "This script is deprecated and may be removed in the future.",
    category=DeprecationWarning,
    stacklevel=1,
)


def _main() -> None:
    d_ = vdd.load.CESMData(strength="size5000", dims=["lat", "lon"])
    d_.temp.plot()
    plt.figure()
    d = vdd.load.CESMData(strength="medium", dims=[])
    # time_sel = np.array([11, 12, 13, 23, 24, 25]) + 12 * 5
    p1 = volcano_base.manipulate.mean_flatten(
        d.temp.isel(time=np.array([11, 12, 13]) - 6), dims=["time"]
    )
    p2 = volcano_base.manipulate.mean_flatten(
        d.temp.isel(time=np.array([23, 24, 25]) - 6), dims=["time"]
    )
    prat = p1 - p2
    p = prat.plot(
        transform=ccrs.PlateCarree(),
        vmin=-5,
        vmax=7,
        cmap=cmc.batlow,
        subplot_kws={"projection": ccrs.PlateCarree()},
    )
    # p = volcano_base.manipulate.mean_flatten(
    #     d.temp.isel(time=time_sel), dims=["time"]
    # ).plot(
    #     transform=ccrs.PlateCarree(),
    #     vmin=-20,
    #     vmax=15,
    #     cmap=cmc.batlow,
    #     subplot_kws={"projection": ccrs.PlateCarree()},
    # )
    p.axes.coastlines()
    plt.show()


if __name__ == "__main__":
    _main()
