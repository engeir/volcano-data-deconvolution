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
