"""Check which eruptions are close and very large or small."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import volcano_base

import vdd.load

_SAVE_DIR = volcano_base.config.SAVE_PATH / "ob16-so2-info"
if not _SAVE_DIR.exists():
    _SAVE_DIR.mkdir(parents=False)
PEAK_LIMIT = 0.1  # Tg -- in order to remove all the zeros
SMALL_LIMIT = 20  # Tg -- limit defining small eruptions
BIG_LIMIT = 80  # Tg -- limit defining big eruptions
CLOSE_LIMIT = 4  # yr -- limit defining close eruptions

_ = volcano_base.load.OttoBliesner(freq="h0", progress=True)
ob16 = vdd.load.DeconvolveOB16(data="h0", length=12001)
so2 = ob16.so2  # .plot()  # type: ignore[call-arg]
so2 = so2.assign_coords(time=volcano_base.manipulate.dt2float(so2.time.data))
so2 = so2.assign_coords(time=np.asarray([round(val) for val in so2.time.data]))
so2_peaks = so2.where(so2 > PEAK_LIMIT).dropna("time")
so2_peaks_series = so2_peaks.to_series()
peak_times = np.array(so2_peaks.time.data)
print(peak_times)
peaks_sort: np.ndarray = np.concatenate(([100], np.diff(peak_times)))
print(peaks_sort)
so2_peaks_close = so2_peaks.where(peaks_sort < CLOSE_LIMIT).dropna("time")
so2_peaks_close_series = so2_peaks_close.to_series()
so2_peaks_big = so2_peaks.where(so2_peaks > BIG_LIMIT).dropna("time")
so2_peaks_big_series = so2_peaks_big.to_series()
so2_peaks_small = so2_peaks.where(so2_peaks < SMALL_LIMIT).dropna("time")
so2_peaks_small_series = so2_peaks_small.to_series()
isolated_normal = []
i = 0
for val in so2_peaks:
    match val.time.data:
        case item if item in so2_peaks_big.time.data:
            isolated_normal.append(False)
        case item if item in so2_peaks_small.time.data:
            isolated_normal.append(False)
        case item if item in so2_peaks_close.time.data:
            isolated_normal.append(False)
        case _:
            isolated_normal.append(True)
peaks_normal = so2_peaks.where(isolated_normal).dropna("time")
peaks_normal_series = peaks_normal.to_series()
print(so2_peaks_close_series)
peak_series = pd.DataFrame(
    {
        # "peaks": so2_peaks_series,
        "Close": so2_peaks_close_series,
        "Big": so2_peaks_big_series,
        "Small": so2_peaks_small_series,
        "Normal": peaks_normal_series,
    },
)
print(
    f"There are {len(so2_peaks)} peaks, of which {len(so2_peaks_close)} are close,"
    f" {len(so2_peaks_big)} are big, and {len(so2_peaks_small)} are small. Overlapping"
    f" traits happen"
    f" {len(peaks_normal) - (len(so2_peaks) - len(so2_peaks_close) - len(so2_peaks_big) - len(so2_peaks_small))}"
    f" times. This leaves {len(peaks_normal)} normal peaks.",
)
# There are 71 peaks, of which 3 are close, 5 are big, and 47 are small. Overlapping
# traits happen 2 times. This leaves 18 normal peaks.

# Calculate the number of overlaps for each element. Adjust the values by dividing by
# the number of overlaps to make the stacked version have the correct value.
overlap_counts = (peak_series > 0).sum(axis=1)
adjusted_peak_series = peak_series.div(overlap_counts, axis=0).fillna(0)

plt.figure()
plt.semilogy()
plt.axhline(SMALL_LIMIT, c="gray", zorder=-10, lw=0.5)
plt.axhline(BIG_LIMIT, c="gray", zorder=-10, lw=0.5)
adjusted_peak_series.plot.bar(ax=plt.gca(), fontsize=3, stacked=True)
plt.xlabel("Time [yr]")
plt.ylabel(r"$\mathrm{SO_{2}} [\mathrm{Tg}]$")
plt.ylim((0.6, 1000))
plt.savefig(_SAVE_DIR / "so2-peak-info")
plt.show()
