# Manuscript overview

## Structure and planned talking points

- The deconvolution works well with SO2 delta pulses as forcing, both for daily and
  monthly resolved dataset
- Deconvolution of temperature to RF is worse due to the noise in the forcing signal
  (But this should be attempted with the alternative method)
- The development of the temperature response to SO2 delta forcing is consistent between
  the Otto-Bliesner et al. (2016) dataset and the CESM2 small volcanic eruption forcing
- The amplitude of the temperature response to SO2 delta forcing within OB16 is smaller,
  but similar to the intermediate CESM2 forcing (same as seen on previous paper)
- The radiative forcing response to SO2 delta forcing within OB16 is relatively stronger
  than that of temperature, with an amplitude between the smallest and second smallest
  CESM2 forcing

## To do

1. Compute the residual as $\mathrm{conv}(T, R_{\mathrm{RF~or~SO2}})-T$, that is, both
   reconstructed from T---RF and T---SO2 response
2. Look at the residual spectrum compared to that of a control simulation
   - Gaussian?
   - Any difference between the proper reconstruction (T---SO2) and the noisy one
     (T---RF)?
   - What happen at low frequency?
3. Look at the correlation between the residual and the reconstructed
   - Any extra peaks in the residual?
4. Look at the difference in peak values between reconstructions and the original
   - Is the peak differences symmetric or skewed? If skewed, how?

## Response functions

### RF to SO2 response

Let us first look at the response function for the radiative forcing, both the true
values and normalised for comparing the evolution.

![Radiative forcing to SO2 response functions](./rf-so2.png)

![Normalised radiative forcing to SO2 response functions](./rf-so2-norm.png)

As expected, the RF response functions are similar in shape across CESM2 simulations,
but differ in amplitude. Compared to the RF response function of Otto-Bliesner et al.
(2016), the peak is reached at similar times, but the Otto-Bliesner et al. (2016)
response has a sharper peak and thus an earlier decay. The amplitude of the RF response
of Otto-Bliesner et al. (2016) is right between the smallest and second smallest CESM2
RF response functions, as expected considering this is where the most prominent volcanic
eruptions in the Otto-Bliesner et al. (2016) dataset are located.

### Temperature to SO2 response

Now we can look at the temperature response functions. Again, we show both the true
responses and the normalised responses.

![Temperature SO2 response functions](./temp-so2.png)

![Normalised temperature to SO2 response functions](./temp-so2-norm.png)

For temperature, we find the Otto-Bliesner et al. (2016) response amplitude to lie
between the second smallest and second largest CESM2 response amplitudes. The shape is
however very similar to the smallest CESM2 response function, which in turn differs from
the three larger CESM2 response functions.

### Temperature to RF response

Finally, we look at the temperature response to the radiative forcing.

![Temperature response functions](./temp-rf.png)

![Normalised temperature response functions](./temp-rf-norm.png)

The temperature to RF response functions are much harder to compute using the
deconvolution algorithm. Since each of the two input time series to the algorithm are
noisy, the performance suffer from poor alignment in time as well as the general noise
level. The daily resolved Otto-Bliesner et al. (2016) dataset fails to produce a
meaningful result due to the time misalignment, while the monthly resolved dataset
produces a result that is noisy yet showing a reasonable shape. The CESM2 response
functions are all very noisy, but the general shape, at least of the three larger
response functions, is reasonable.

## Reconstruction

We next reconstruct the OB16 dataset to compare the residual to a control simulation as
well as the reconstructions themselves against the original.

### Time series

Let us first have a look at the temperature from both a control simulation, and the SO2
forced simulation, as well as the reconstructed temperature from the SO2 response and
the RF response.

![Reconstructed temperature time series](./temp-reconstructed.png)

We notice that the SO2 response reconstruction is almost without noise, which is to be
expected as it is the convolution between the response function, and a train of delta
pulses. The RF response reconstruction on the other hand is very noisy, and follows the
original temperature time series more closely as more of the variability is captured
from using the RF time series rather than the SO2 time series.

### Residuals

Next, we look at the residuals from the reconstructions, and specifically the
correlation between the residuals and the reconstructions.

![Residuals](./correlation-residual-reconstructed.png)

The residuals are here defined as the difference between the original temperature time
series and the reconstructed time series. As such, the correlation function between
reconstructed and residual from using the RF response function give strong negative
correlation at small time lags, with a weaker positive correlation at larger time lags.

The correlation function between the SO2 response reconstruction and the residual is
much more flat, but with spuriously strong correlations at all time lags.

### Spectrum

Below is a plot showing the power spectral density of the two residual time series, and
the control temperature time series.

![Spectrum](./spectrum-residual-control_temp.png)

Based on how the reconstructed temperature time series from the SO2 and RF response
functions are constructed, we expect the residual from the RF response reconstruction to
feature more power at high frequencies, as the reconstructed time series contain noise
from both the response function and the radiative forcing time series.

We would expect the control simulation temperature to be close to Gaussian noise, and
this is indeed what the spectrum shows.

### Peak differences

We want to investigate how well we are able to resolve the true peaks in the temperature
from the reconstructions. Ideally, the only differences in peak values should be due to
noise in the temperature time series, and not due to the reconstruction method.

We plot both the PDF, and the CDF of the peak difference time series.

![Peak difference PDF](./peak-difference-pdf.png)

![Peak differences CDF](./peak-difference-cdf.png)
