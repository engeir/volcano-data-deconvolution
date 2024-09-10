"""Calculate the effective radius.

This calculation follows the derivation presented in Clyne et al. (2021) [1], Appendix
A.

[1]: Clyne, M., Lamarque, J-F., et al. (2021), Model physics and chemistry causing
     intermodel disagreement within the VolMIP-Tambora Interactive Stratospheric Aerosol
     ensemble. DOI: 10.5194/acp-21-3317-2021
"""

import numpy as np
import volcano_base
import xarray as xr

# These are extracted from the dimensions of the output files from CESM2, and represent
# the vertical level borders (ILEV) and the middle value at each level (LEV).
# fmt: off
ILEV = [
    4.50049997269275e-06, 7.42010008991656e-06, 1.22337002750328e-05, 2.01700007806949e-05,
    3.32545013748131e-05, 5.48274989853326e-05, 9.03979966437873e-05, 0.000149040005226198,
    0.000245720002567396, 0.000405125007318929, 0.000667940014409396, 0.00110126495656004,
    0.00181564996637462, 0.00299349994747899, 0.00496299981023185, 0.00815065141068771,
    0.0134769998112461, 0.0223189999815077, 0.0367965003533754, 0.060664999182336,
    0.0991565029835328, 0.157389993546531, 0.238849999732338, 0.34520000917837,
    0.475134991575032, 0.631804985459894, 0.829154974780977, 1.08274002559483,
    1.40684994403273, 1.81885005440563, 2.33979988843203, 2.99505004659295,
    3.81470005959272, 4.83444985002279, 6.09635002911091, 7.64934998005629,
    9.55010019242764, 11.8640000000596, 14.6655002608895, 18.0380009114742,
    22.0755003392696, 26.8824994564056, 32.5734987854958, 39.2730012536049,
    47.1144989132881, 56.2404990196228, 66.8004974722862, 80.7014182209969,
    94.9410423636436, 111.69321089983, 131.401270627975, 154.586806893349,
    181.863352656364, 213.952820748091, 251.704417169094, 296.117216348648,
    348.366588354111, 409.83521938324, 482.149928808212, 567.22442060709,
    652.332969009876, 730.445891618729, 796.363070607185, 845.353666692972,
    873.715866357088, 900.324631482363, 924.964462406933, 947.432334534824,
    967.538624536246, 985.112190246582, 1000,
]
LEV = [
    5.96030003130466e-06, 9.82690018247467e-06, 1.62018505278638e-05, 2.6712251077754e-05,
    4.40410001800728e-05, 7.261274781456e-05, 0.000119719000934992, 0.000197380003896797,
    0.000325422504943162, 0.000536532510864163, 0.00088460248548472, 0.00145845746146733,
    0.00240457495692681, 0.00397824987885542, 0.00655682561045978, 0.0108138256109669,
    0.0178979998963769, 0.0295577501674416, 0.0487307497678557, 0.0799107510829344,
    0.128273248265032, 0.198119996639434, 0.292025004455354, 0.410167500376701,
    0.553469988517463, 0.730479980120435, 0.955947500187904, 1.24479498481378,
    1.61284999921918, 2.07932497141883, 2.66742496751249, 3.40487505309284,
    4.32457495480776, 5.46539993956685, 6.8728500045836, 8.59972508624196,
    10.7070500962436, 13.2647501304746, 16.3517505861819, 20.0567506253719,
    24.4789998978376, 29.7279991209507, 35.9232500195503, 43.1937500834465,
    51.6774989664555, 61.5204982459545, 73.7509578466415, 87.8212302923203,
    103.317126631737, 121.547240763903, 142.994038760662, 168.225079774857,
    197.908086702228, 232.828618958592, 273.910816758871, 322.241902351379,
    379.100903868675, 445.992574095726, 524.687174707651, 609.778694808483,
    691.389430314302, 763.404481112957, 820.858368650079, 859.53476652503,
    887.020248919725, 912.644546944648, 936.198398470879, 957.485479535535,
    976.325407391414, 992.556095123291,
]
# fmt: on


class Reff:
    r"""Compute effective radius.

    Parameters
    ----------
    reff : xr.DataArray
        The 3D aerosol effective radius in cm (`REFF_AERO`).
    temp : xr.DataArray
        The 3D temperature in kelvin (K) (`T`).
    sad : xr.DataArray
        The 3D aerosol surface area density in cm2/cm3 (`SAD_AERO`).

    Examples
    --------
    Computing :math:`R_{\text{eff}}` using the format of output files provided in the accompanying NIRD
    archive can be done in conjunction with the `volcano_base` package.

    >>> import volcano_base

    >>> FINDER = (
    ...     volcano_base.load.FindFiles()
    ...     .find({"SAD_AERO", "REFF_AERO", "T"}, {"ens1", "ens3"}, {"tt-2sep", "tt-4sep"})
    ...     .sort("ensemble", "sim")
    ...     .copy
    ... )
    >>> sad_ = FINDER().keep("SAD_AERO")
    >>> reff_ = FINDER().keep("REFF_AERO")
    >>> temp_ = FINDER().keep("T")
    >>> sad_xr = sad_.load()
    >>> s1 = sad_xr[0]
    >>> reff_xr = reff_.load()
    >>> r1 = reff_xr[0]
    >>> temp_xr = temp_.load()
    >>> t1 = temp_xr[0]
    >>> reff = Reff(r1, t1, s1).calculate_reff()
    >>> # reff.to_netcdf("reff.nc")
    >>> reff.plot()
    >>> plt.show()
    """

    def __init__(
        self, reff: xr.DataArray, temp: xr.DataArray, sad: xr.DataArray
    ) -> None:
        # Convert from cm to µm
        self.reff = (self._flatten_if_3d(reff) * 10_000).assign_attrs(reff.attrs)
        self.temp = self._flatten_if_3d(temp)
        self.sad = self._flatten_if_3d(sad)
        self.h = self._calculate_h()
        self.coeff = self.sad * self.h

    def _calculate_h(self) -> xr.DataArray:
        r"""Calculate the vertical thickness.

        Returns
        -------
        h : xr.DataArray
            The vertical layer thickness at each vertical coordinate.

        Notes
        -----
        This implements equation A1 of Clyne et al. (2021).

        .. math::

            h = \frac{RT}{g}\ln{\frac{P_2}{P_1}}

        where :math:`R` is the dry air gas constant, :math:`g` is acceleration due to
        gravity at sea level, and :math:`T` is the average temperature of the layer
        between pressure levels :math:`P_2<P_1`.
        """
        dry_air_gas_const = 287  # J kg^-1 K^-1
        g = 9.81  # m s^-2
        for p1, p2, temp_2d in zip(ILEV[1:], ILEV[:-1], self.temp.T, strict=True):
            # p1 is closest to the surface (higher pressure) than p2. Both go top-down.
            h: xr.DataArray = dry_air_gas_const * temp_2d / g * np.log(p2 / p1)
        return h.compute()

    def _flatten_if_3d(self, arr: xr.DataArray) -> xr.DataArray:
        """If the input arrays include lat/lon coordinates, average them out."""
        max_dims = 4
        return (
            volcano_base.manipulate.mean_flatten(arr, dims=["lat", "lon"])
            if len(arr.dims) == max_dims
            else arr
        )

    def calculate_reff(self) -> xr.DataArray:
        r"""Calculate the aerosol effective radius.

        Returns
        -------
        reff : xr.DataArray
            The global stratospheric aerosol effective radius.

        Notes
        -----
        This implements equation A3 of Clyne et al. (2021).

        .. math::

            R_{\text{eff}} =
            \frac{\sum_{\tau=\tau_1} (\text{SAD}\cdot h\cdot R_{\text{eff}})_{\tau}}
            {\sum_{\tau=\tau_1} (\text{SAD}\cdot h)_{\tau}}

        where :math:`\{text{SAD}}` is the surface aerosol density, and :math:`h` is the
        vertical layer thickness in meters from A1/`self._calculate_h`.
        """
        num = self.coeff * self.reff
        reff = num.sum(dim="lev") / self.coeff.sum(dim="lev")
        return reff.assign_attrs(self.reff.attrs)
