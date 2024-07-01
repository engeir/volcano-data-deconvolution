"""Script that answers the question "How close are the eruptions in OB16?"."""

import matplotlib.pyplot as plt
import numpy as np
import volcano_base
import xarray as xr

ob16 = volcano_base.load.OttoBliesner(freq="h0", progress=True)
eruption_distances = np.asarray(
    volcano_base.manipulate.dt2float(
        ob16.so2_delta.where(ob16.so2_delta > 0, drop=True).time.data,
    )
    .diff()
    .dropna(),
)
arr = (
    xr.DataArray(
        volcano_base.manipulate.dt2float(
            ob16.so2_delta.where(ob16.so2_delta > 0, drop=True).time.data,
        ),
        dims=("time",),
    )
    .diff("time")
    .dropna("time")
)
print(arr.mean().data)
for call in ["min", "max"]:
    # Print the statistic as well as the index and value of the eruption that caused it
    print(
        f"{call}: a {getattr(arr, call)().data} year gap ended in the decimal year {getattr(arr, f"idx{call}")().data}",
    )

ob16.so2_delta.plot()
plt.show()

# The answer is:
# Minimum eruption gap: 1.915 decimal years that ended in the decimal year 1476.619
# Maximum eruption gap: 107.0 decimal years that ended in the decimal year 854.704
