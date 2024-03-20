---
header-includes: |
   \usepackage{tcolorbox}
   \usepackage{tikz}
   \usetikzlibrary{calc}
   \tcbuselibrary{skins,xparse,breakable}
   \tcbset{%
       % enhanced jigsaw,  % Makes the breakpoint sharp
       colback=white,
       colframe=black,
       % title=Answer,  % Title of the block
   %   bookmark={Q\arabic{\tcbcounter}}
   }
   \newtcolorbox{myquote}{%
       breakable,
       pad at break*=1.5pc,
       colback=gray!15!white,
       colframe=gray!15!white,
       overlay first and middle={
         \coordinate (A1) at ($(interior.south east) + (-10pt,5pt)$);
         \coordinate (C1) at ($(interior.south east) + (-6pt,7.5pt)$);
         \draw[fill=gray] (A1) -- +(0,5pt) -- (C1) -- cycle;
       }
   }
   \renewenvironment{quote}{\begin{NoHyper}\begin{myquote}}{\end{myquote}\end{NoHyper}}
---

# Manuscript overview

## Structure and planned talking points

- The deconvolution works well with SO2 delta pulses as forcing, both for daily and
  monthly resolved dataset
- Deconvolution of temperature to RF is worse due to the noise in the forcing signal
  (But this may be attempted with the alternative method)
- The development of the temperature response to SO2 delta forcing is consistent between
  the Otto-Bliesner et al. (2016) dataset and the CESM2 small volcanic eruption forcing
- The amplitude of the temperature response to SO2 delta forcing within OB16 is smaller,
  but similar to the intermediate CESM2 forcing (same as seen on previous paper)
- The radiative forcing response to SO2 delta forcing within OB16 is relatively stronger
  than that of temperature, with an amplitude between the smallest and second smallest
  CESM2 forcing
- The temperature response to radiative forcing normalised to OB16 has the best
  correspondence with the 2-year separated double waveform response function, suggesting
  that (1) the temperature does not combine linearly with the radiative forcing when
  perturbed and (2) the deconvolution may not be appropriate to fully describe the
  response function.

## To do

1. ~~Compute the residual as $\mathrm{conv}(T, R_{\mathrm{RF~or~SO2}})-T$, that is, both
   reconstructed from T---RF and T---SO2 response~~
2. ~~Look at the residual spectrum compared to that of a control simulation~~
   - Gaussian? **Yes**
   - Any difference between the proper reconstruction (T---SO2) and the noisy one
     (T---RF)? **Slightly, see section below**
   - What happen at low frequency? **RF looses power**
3. ~~Look at the correlation between the residual and the reconstructed~~
   - Any extra peaks in the residual? **Stronger correlation in the RF case**
4. ~~Look at the difference in peak values between reconstructions and the original~~
   - Is the peak differences symmetric or skewed? If skewed, how? **Only about 22%
     confidence that the population mean is different from zero**

---

## Table of contents

<!-- vim-markdown-toc GFM -->

- [Definitions](#definitions)
- [Should we expect linear temperature dependence on radiative forcing?](#should-we-expect-linear-temperature-dependence-on-radiative-forcing)
- [Response functions](#response-functions)
  - [RF to SO2 response](#rf-to-so2-response)
  - [Temperature to SO2 response](#temperature-to-so2-response)
  - [Temperature to RF response](#temperature-to-rf-response)
- [Reconstruction](#reconstruction)
  - [Time series](#time-series)
  - [Residuals](#residuals)
  - [Spectrum](#spectrum)
  - [Peak differences](#peak-differences)
- [Parametrisation](#parametrisation)
  - [Experiment results](#experiment-results)
- [Double waveform](#double-waveform)
  - [Double waveform response functions](#double-waveform-response-functions)
  - [Reconstructed double waveforms](#reconstructed-double-waveforms)
- [Cut off response function](#cut-off-response-function)
  - [Inspecting the noise floor](#inspecting-the-noise-floor)

<!-- vim-markdown-toc -->

## Definitions

- **SO2**: Sulphur dioxide --- $S$
- **RF**: Radiative forcing --- $R$
- **T**: Temperature --- $T$
- Convolution is denoted as $\ast$
- Deconvolution is denoted as $\tilde\ast$
- The response function of $A$ from $B$ is denoted $\phi_{AB}$, where $A$ is the signal
  and $B$ is the kernel

> Thus, $\phi_{AB}=A\tilde\ast B$, and $A=B\ast\phi_{AB}$.

## Should we expect linear temperature dependence on radiative forcing?

A lot of these results suggest that even though there are non-linearities in the
conversion from SO2 to both radiative forcing and temperature, the temperature response
may still be linearly dependent on the radiative forcing. This is a strong result in our
previous paper investigating single waveform volcanic eruption simulations, but if this
is either trivial or at least to be expected, there is less to be gained from finding a
good estimate of the temperature to radiative forcing response function.

## Response functions

### RF to SO2 response

Let us first look at the response function for the radiative forcing, both the true
values and normalised for comparing the evolution.

![Radiative forcing to SO2 response functions](./generated_files/deconv_ob16_cesm2/rf-so2.png)

![Normalised radiative forcing to SO2 response functions](./generated_files/deconv_ob16_cesm2/rf-so2-norm.png)

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

![Temperature SO2 response functions](./generated_files/deconv_ob16_cesm2/temp-so2.png)

![Normalised temperature to SO2 response functions](./generated_files/deconv_ob16_cesm2/temp-so2-norm.png)

For temperature, we find the Otto-Bliesner et al. (2016) response amplitude to lie
between the second smallest and second largest CESM2 response amplitudes. The shape is
however very similar to the smallest CESM2 response function, which in turn differs from
the three larger CESM2 response functions.

### Temperature to RF response

Finally, we look at the temperature response to the radiative forcing.

![Temperature response functions](./generated_files/deconv_ob16_cesm2/temp-rf.png)

![Normalised temperature response functions](./generated_files/deconv_ob16_cesm2/temp-rf-norm.png)

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
well as the reconstructions themselves against the original. We use both the OB16
response functions themselves, and the CESM2 strong simulation response functions in a
pairwise comparison.

### Time series

Let us first have a look at the temperature from both a control simulation, and the SO2
forced simulation, as well as the reconstructed temperature from the SO2 response and
the RF response.

![Reconstructed temperature time series OB16](./generated_files/reconstruction/ob16-month-temp-reconstructed.png){width=49%}
![Reconstructed temperature time series CESM2](./generated_files/reconstruction/cesm2-strong-temp-reconstructed.png){width=49%}

We notice that the SO2 response reconstruction is almost without noise, which is to be
expected as it is the convolution between the response function, and a train of delta
pulses. The RF response reconstruction is very noisy, yet follows the original
temperature time series more closely as more of the variability is captured from using
the RF time series rather than the SO2 time series.

### Residuals

Next, we look at the residuals from the reconstructions, and specifically the
correlation between the residuals and the reconstructions.

![Residuals OB16](./generated_files/reconstruction/ob16-month-correlation-residual-reconstructed.png){width=49%}
![Residuals CESM2](./generated_files/reconstruction/cesm2-strong-correlation-residual-reconstructed.png){width=49%}

The residuals are here defined as the difference between the original temperature time
series and the reconstructed time series. As such, the correlation function between
reconstructed and residual from using the RF response function give strong negative
correlation at small time lags, with a weaker positive correlation at larger time lags.

The correlation function between the SO2 response reconstruction and the residual is
much more flat, but with spuriously strong correlations at all time lags.

### Spectrum

Below is a plot showing the power spectral density of the two residual time series, and
the control temperature time series.

![Spectrum OB16](./generated_files/reconstruction/ob16-month-spectrum-residual-control_temp.png){width=49%}
![Spectrum CESM2](./generated_files/reconstruction/cesm2-strong-spectrum-residual-control_temp.png){width=49%}

Based on how the reconstructed temperature time series from the SO2 and RF response
functions are constructed, we expect the residual from the RF response reconstruction to
feature more power at high frequencies, as the reconstructed time series contain noise
from both the response function and the radiative forcing time series. Likewise, the
RF-reconstructed temperature follows the original temperature time series more closely,
thus reducing the slow variability in the residual time series.

We would expect the control simulation temperature to be close to Gaussian noise, and
this is indeed what the spectrum shows.

### Peak differences

We want to investigate how well we are able to resolve the true peaks in the temperature
from the reconstructions. Ideally, the only differences in peak values should be due to
noise in the temperature time series, and not due to the reconstruction method.

We plot both the PDF, and the CDF of the peak difference time series.

![Peak difference PDF OB16](./generated_files/reconstruction/ob16-month-peak-difference-pdf.png){width=49%}
![Peak difference PDF CESM2](./generated_files/reconstruction/cesm2-strong-peak-difference-pdf.png){width=49%}

![Peak differences CDF OB16](./generated_files/reconstruction/ob16-month-peak-difference-cdf.png){width=49%}
![Peak differences CDF CESM2](./generated_files/reconstruction/cesm2-strong-peak-difference-cdf.png){width=49%}

## Parametrisation

We have good estimates of the RF to SO2 and T to SO2 response functions. Both from CESM2
simulations, but also from the Otto-Bliesner et al. (2016) dataset. The next step is to
obtain a good estimate of the T to RF response function. This is ideally as simple as
deconvolving the RF to SO2 response function from the T to SO2 response function, but as
both response functions have a width as well as being noisy, this is not a trivial task.

We first list up some useful relations in order to better understand the problem at
hand.

<!-- dprint-ignore-start -->
$$
\begin{aligned}
R&:=S\ast\phi_{RS}\\
T&:=S\ast\phi_{TS}\\
T&:=R\ast\phi_{TR}=S\ast\phi_{RS}\ast\phi_{TR}\\
\Rightarrow \phi_{TS}&=\phi_{RS}\ast\phi_{TR}\\
\Rightarrow \phi_{TR}&=\phi_{TS}\tilde\ast\phi_{RS}=(T\tilde\ast S)\tilde\ast (R\tilde\ast S)=T\tilde\ast R\quad[=T\tilde\ast(S\ast\phi_{RS})]\\
\end{aligned}
$$ {#eq:label}
<!-- dprint-ignore-end -->

- [x] $\phi_{RS}=R\tilde\ast S$
- [x] $\phi_{TS}=T\tilde\ast S$
- [ ] $\phi_{TR}=T\tilde\ast R = (T\tilde\ast S)\tilde\ast (R\tilde\ast S) =
      \phi_{TS}\tilde\ast\phi_{RS}$

We conclude that we may compute $\phi_{TR}$ in two ways, (1) either by **deconvolving
$T$ with $R$**, or (2) by __deconvolving $\phi_{TS}$ with $\phi_{RS}$__, which in turn
are the results of two further deconvolutions. A third option is to use the
deconvolution between $T$ and the convolution between $S$ and $\phi_{RS}$, but this is a
reconstruction of the original $R$, and perhaps not what we are looking for.

### Experiment results

Below, we calculate the response function representing the temperature response to RF in
three different ways in all figures. The responses are

1. Deconvolving $T$ with $R$ (`dec(T, RF)`)
2. Deconvolving $\phi_{TS}$ with $\phi_{RS}$ (`dec(dec(T, SO2), dec(R, SO2))`)
3. Deconvolving $\phi_{TS}$ with $\phi_{RS}$ while enforcing negative time lags to be
   zero (`dec(dec(T, SO2), dec(R, SO2)), corrected`)

![Otto-Bliesner daily](./generated_files/parametrisation/parametrisation_ob16.png)

![Otto-Bliesner monthly](./generated_files/parametrisation/parametrisation_ob16-month.png)

![CESM2 small](./generated_files/parametrisation/parametrisation_cesm2-medium.png)

![CESM2 intermediate](./generated_files/parametrisation/parametrisation_cesm2-medium-plus.png)

![CESM2 large](./generated_files/parametrisation/parametrisation_cesm2-strong.png)

![CESM2 extra large](./generated_files/parametrisation/parametrisation_cesm2-size5000.png)

![CESM2 double 2-sep](./generated_files/parametrisation/parametrisation_cesm2-tt-2sep.png)

![CESM2 double 4-sep](./generated_files/parametrisation/parametrisation_cesm2-double-overlap.png)

An average across the three estimation methods, with all datasets shown in the same
plot:

![Average](./generated_files/parametrisation/parametrisation_ensemble.png)

An average across a subset of all simulation cases for all three estimation methods:

![Method means](./generated_files/parametrisation/parametrisation_method.png)

From this, we find that the overall best estimates are form directly deconvolving $T$
with $R$. Evidently, it seems that noise is not the main issue, but rather the
misalignment between the signal and kernel that is used in the deconvolution. This is
backed up by the very poor results form using daily resolved data within the
Otto-Bliesner dataset, which just gets worse when smoothed out.

## Double waveform

We want to see how well the deconvolution is suited for estimating the response
functions in simulations where volcanic eruptions occur close in time.

### Double waveform response functions

We first compare the response functions estimated from the CESM2 double waveform
simulations, where a two year (Fig \ref{fig:responses_rf}) and four year (Fig
\ref{fig:responses_temp}) separation between the volcanic eruptions have been used.

![Response functions RF](./generated_files/waveform/responses_rf.png){#fig:responses_rf}

![Response functions T](./generated_files/waveform/responses_temp.png){#fig:responses_temp}

### Reconstructed double waveforms

Let us use the CESM2 strong eruption ensemble as a good estimate of the true response
functions, and then compare the climate model output of $R$ and $T$ with the
reconstructed ones from the CESM2 strong ensemble response functions.

![True (old) versus reconstructed (new) AOD](./generated_files/waveform/recreated_waveforms_aod.png)

![True (old) versus reconstructed (new) RF](./generated_files/waveform/recreated_waveforms_rf.png)

![True (old) versus reconstructed (new) T](./generated_files/waveform/recreated_waveforms_temp.png)

The response functions have here been rescaled to have amplitude equal to the amplitude
of the corresponding response functions generated from the simulation in question. The
plots give good indications that for the four year separated double waveform, the
reconstruction is indistinguishable from the true response functions at the current
noise level. However, the two year separated double waveform shows a clear non-linear
dependence between the injected SO2 and resulting $R$ and $T$ time series. A
significantly weaker response in both $R$ and $T$ is seen in the second eruption, most
notably in the $T$ plot.

## Cut off response function

To check how far into the response functions the noise is substantial, we cut the
response function at a certain time lag, before we generate an ensemble of temperature
output signals by first convolving the cut response function with the forcing, and then
adding phase shifted noise to the temperature signal. The precedure can be summarise in
the following steps:

1. Cut the response function at a given index.
2. Convolve the cut response function with the forcing signal to obtain a temperature
   estimate $T_{\mathrm{est}}$.
3. Add noise represented by a phase shifted temperature control signal to obtain an
   ensemble of temperature estimates $T_i$.
4. Deconvolve the ensemble of temperature estimates with the forcing signal to obtain an
   ensemble of response functions for the given cut off index.

### Inspecting the noise floor

We get estimates of the noise present in the response functions from doing the
bootstrapping method outlined above. The noise is then illustrated as percentile plots
that show the percentiles as layers of shaded filled in area plots.

Let us look at how the cut off ensemble looks like for the Otto-Bliesner dataset when
cutting the response function off at years 2, 4, 8 and 16.

![Cut off response function OB16](./generated_files/cut_off/resp_024.png)

![Cut off response function OB16](./generated_files/cut_off/resp_048.png)

![Cut off response function OB16](./generated_files/cut_off/resp_096.png)

![Cut off response function OB16](./generated_files/cut_off/resp_192.png)
