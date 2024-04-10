---
title: "Manuscript overview"
author:
- Eirik Rolland Enger
date: \today{}
geometry: margin=3cm
numbersections: true
toc: true
shiftheadinglevelby: -1
header-includes: |
  ```{=latex}
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
       arc=0pt,
       outer arc=0pt,
       breakable,
       boxrule=0.5pt,
       pad at break*=1.5pc,
       colback=gray!15!white,
       colframe=black,
       overlay first and middle={
         \coordinate (A1) at ($(interior.south east) + (-10pt,5pt)$);
         \coordinate (C1) at ($(interior.south east) + (-6pt,7.5pt)$);
         \draw[fill=gray] (A1) -- +(0,5pt) -- (C1) -- cycle;
       }
   }
   \renewenvironment{quote}{\begin{NoHyper}\begin{myquote}}{\end{myquote}\end{NoHyper}}
  ```
---

## Preface {-}

### Structure and planned talking points {-}

<!-- FIXME: update the talking points
-->

- The deconvolution works well with SO~2~ delta pulses as forcing, both for daily and
  monthly resolved dataset (section \ref{response-functions})
- Deconvolution of temperature to RF is worse due to the noise in the forcing signal,
  but most likely sufficient from the strongest CESM2 simulations as well as the OB16
  data (may attempt with the alternative method, but probably too computationally
  expensive, at least with OB16) (section \ref{parametrisation})
- The development of the temperature response to SO~2~ delta forcing is consistent
  between the OB16 dataset and the CESM2 small volcanic eruption forcing
- The amplitude of the temperature response to SO~2~ delta forcing within OB16 is
  smaller, but similar to the intermediate CESM2 forcing (same as seen on previous
  paper)
- The radiative forcing response to SO~2~ delta forcing within OB16 is relatively
  stronger than that of temperature, with an amplitude between the smallest and second
  smallest CESM2 forcing
- The temperature response to radiative forcing normalised to OB16 has the best
  correspondence with the 2-year separated double waveform response function, suggesting
  that (1) the temperature does not combine linearly with the radiative forcing when
  perturbed and (2) the deconvolution may not be appropriate to fully describe the
  response function (section \ref{reconstruction})
- Cut off analysis shows that the response of RF to SO~2~ is significant only up to
  about 4 years, while the temperature response to both SO~2~ and RF is significant up
  to at least 16 years (section \ref{cut-off-response-function})
- The data suggest that as long as AOD is not perturbed, the temperature response to
  volcanic eruptions combine linearly, and so does the radiative forcing (section
  \ref{double-waveform}, double waveform analysis, Fig. \ref{waveform-comparison})

### To do {-}

1. ~~Compute the residual as $T-\mathrm{conv}(T, \phi_{TS\mathrm{\,or\,}TR})$, that is,
   both reconstructed from T---RF and T---$\mathrm{SO_2}$ response~~
2. ~~Look at the residual spectrum compared to that of a control simulation~~
   - Gaussian? **Yes**
   - Any difference between the proper reconstruction (T---SO~2~) and the noisy one
     (T---RF)? **Slightly, see section below**
   - What happen at low frequency? **RF looses power**
3. ~~Look at the correlation between the residual and reconstructed~~
   - Any extra peaks in the residual? **Stronger correlation in the RF case**
4. ~~Look at the difference in peak values between reconstructions and the original~~
   - Is the peak differences symmetric or skewed? If skewed, how? **Only about 22%
     confidence that the population mean is different from zero**

---

## Intro

### Definitions

- **SO~2~**: Sulphur dioxide --- $S$
- **RF**: Radiative forcing --- $R$
- **T**: Temperature --- $T$
- Convolution is denoted as $\ast$
- Deconvolution is denoted as $\tilde\ast$
- The response function of $A$ from $B$ is denoted $\phi_{AB}$, where $A$ is the signal
  and $B$ is the kernel
- @ottobliesner2016 is denoted OB16

> Thus, $\phi_{AB}=A\tilde\ast B$, and $A=B\ast\phi_{AB}$.

### Should we expect linear temperature dependence on radiative forcing?

1. From simulations using ESMs, should we expect a linear dependence between temperature
   and radiative forcing? Yes, provided that the forcing does not last very long, and
   can be seen as a perturbation to the system.
2. In nature, should we expect a linear dependence between temperature and radiative
   forcing?
3. The linear dependence is only valid for peak values, is this a counter-argument that
   a linear dependence is not to be expected?

A lot of these results suggest that even though there are non-linearities in the
conversion from SO~2~ to both radiative forcing and temperature, the temperature
response may still be linearly dependent on the radiative forcing. This is a strong
result in our previous paper investigating single waveform volcanic eruption
simulations, but if this is either trivial or at least to be expected, there is less to
be gained from finding a good estimate of the temperature to radiative forcing response
function.

Simple energy balance models [(EBMs) often assume a linear dependence between
temperature and radiative forcing]{.underline}. It is common to discuss the climate
feedback parameter to be a sum of individual feedbacks (clouds, aerosols, CO~2~, and
more). However, this assumes that the climate sensitivity is a constant, while _from
most AOGCM simulations of constant $4\times\mathrm{CO_2}$ forcing, the climate
sensitivity is found to vary_ [@gregory2016].

However, these processes occur on much longer time scales than the volcanic eruptions we
simulate. Possible non-linear effects such as sea ice growth and subsequent stop due to
lack of sea to grow into, and changes in ocean currents, both to the poles but also
across the Pacific, are not fast enough to play a role over the three years until the
anomalies peak.

So, this topic revolves around whether the climate sensitivity is a constant for any
given forcing, and further whether the climate sensitivity is the same for all forcings.
The latter is very likely no, and the former seems to be incorrect for at least some
forcings. Thus, we should not automatically expect the temperature response to radiative
forcing to be linear, or dismiss it as trivial. But it also matters how fast the changes
are and for how long the changes are sustained, for which the single volcanic eruptions
are likely too short-lived to yield non-linear effects.

## Response functions

### RF to SO~2~ response

Let us first look at the response function for the radiative forcing, both the true
values and normalised for comparing the evolution.

![Absolute and normalised radiative forcing to SO~2~ response functions](./generated_files/deconv_ob16_cesm2/rf-so2.png)

As expected, the RF response functions are similar in shape across CESM2 simulations,
but differ in amplitude. Compared to the RF response function of OB16, the peak is
reached at similar times, but the OB16 response has a sharper peak and thus an earlier
decay. The amplitude of the RF response of OB16 is right between the smallest and second
smallest CESM2 RF response functions, as expected considering this is where the most
prominent volcanic eruptions in the OB16 dataset are located.

### Temperature to SO~2~ response

Now we can look at the temperature response functions. Again, we show both the true
responses and the normalised responses.

![Absolute and normalised temperature SO~2~ response functions](./generated_files/deconv_ob16_cesm2/temp-so2.png)

For temperature, we find the OB16 response amplitude to lie between the second smallest
and second largest CESM2 response amplitudes. We also note that the initial rise of the
responses across all intermediate and smaller eruption sizes is similar. However, they
all reach the peak at different times, and decay at different rates. The shape of the
OB16 response is similar to the smallest CESM2 response function during the decaying
phase (not shown), although there is significant noise in the smallest CESM2 response.
Both OB16 and the smallest CESM2 differs from the three larger CESM2 response functions.

### Temperature to RF response

Finally, we look at the temperature response to the radiative forcing.

![Absolute and normalised temperature response functions](./generated_files/deconv_ob16_cesm2/temp-rf.png)

The temperature to RF response functions are much harder to compute using the
deconvolution algorithm. Since each of the two input time series to the algorithm are
noisy, the performance suffer from poor alignment in time as well as the general noise
level. The daily resolved OB16 dataset fails to produce a meaningful result due to the
time misalignment, while the monthly resolved dataset produces a result that is noisy
yet showing a reasonable shape. The CESM2 response functions are all very noisy, but the
general shape, at least of the three larger response functions, is reasonable.

## Reconstruction

We next reconstruct the OB16 dataset to compare the residual to a control simulation as
well as the reconstructions themselves against the original.

> **NOTE**: In the images that follow, we use both the OB16 response functions
> themselves, and the CESM2 strong simulation response functions in a pairwise
> comparison.

### Time series

Let us first have a look at the temperature from both a control simulation, and the
SO~2~ forced simulation, as well as the reconstructed temperature from the SO~2~
response and the RF response.

![Reconstructed temperature time series OB16](./generated_files/reconstruction/ob16-month-temp-reconstructed.png){width=49%}
![Reconstructed temperature time series CESM2](./generated_files/reconstruction/cesm2-strong-temp-reconstructed.png){width=49%}

We notice that the SO~2~ response reconstruction is almost without noise, which is to be
expected as it is the convolution between the response function, and a train of delta
pulses. The RF response reconstruction is very noisy, yet follows the original
temperature time series more closely as more of the variability is captured from using
the RF time series rather than the SO~2~ time series.

### Residuals

Next, we look at the residuals from the reconstructions, and specifically the
correlation between the residuals and the reconstructions.

![Residuals OB16](./generated_files/reconstruction/ob16-month-correlation-residual-reconstructed.png){width=49%}
![Residuals CESM2](./generated_files/reconstruction/cesm2-strong-correlation-residual-reconstructed.png){width=49%}

The residuals are here defined as the difference between the original temperature time
series and the reconstructed time series. As such, the correlation function between
reconstructed and residual from using the RF response function give strong negative
correlation at small time lags, with a weaker positive correlation at larger time lags.

The correlation function between the SO~2~ response reconstruction and the residual is
much more flat, but with spuriously strong correlations at all time lags.

### Spectrum

Below is a plot showing the power spectral density of the two residual time series, and
the control temperature time series.

![Spectrum OB16](./generated_files/reconstruction/ob16-month-spectrum-residual-control_temp.png){width=49%}
![Spectrum CESM2](./generated_files/reconstruction/cesm2-strong-spectrum-residual-control_temp.png){width=49%}

Based on how the reconstructed temperature time series from the SO~2~ and RF response
functions are constructed, we expect the residual from the RF response reconstruction to
feature more power at high frequencies, as the reconstructed time series contain noise
from both the response function and the radiative forcing time series. Likewise, the
RF-reconstructed temperature follows the original temperature time series more closely,
thus reducing the slow variability in the residual time series.

We would expect the control simulation temperature to be close to Gaussian noise, and
this is indeed what the spectrum shows.

### Peak differences

We want to investigate how well we can resolve the true peaks in the temperature from
the reconstructions. Ideally, the only differences in peak values should be due to noise
in the temperature time series, and not due to the reconstruction method.

We plot both the PDF, and the CDF of the peak difference time series.

![Peak difference PDF OB16](./generated_files/reconstruction/ob16-month-peak-difference-pdf.png){width=49%}
![Peak difference PDF CESM2](./generated_files/reconstruction/cesm2-strong-peak-difference-pdf.png){width=49%}

![Peak differences CDF OB16](./generated_files/reconstruction/ob16-month-peak-difference-cdf.png){width=49%}
![Peak differences CDF CESM2](./generated_files/reconstruction/cesm2-strong-peak-difference-cdf.png){width=49%}

## Parametrisation

We have good estimates of the RF to SO~2~ and T to SO~2~ response functions. Both from
CESM2 simulations, but also from the OB16 dataset. The next step is to get a good
estimate of the T to RF response function. This is ideally as simple as deconvolving the
RF to SO~2~ response function from the T to SO~2~ response function, but as both
response functions have a width and are noisy, this is not a trivial task.

We first list up some useful relations to better understand the problem at hand.

<!-- dprint-ignore-start -->
$$
\begin{aligned}
R&:=S\ast\phi_{RS}\\
T&:=S\ast\phi_{TS}\\
T&:=R\ast\phi_{TR}=S\ast\phi_{RS}\ast\phi_{TR}\\
\Rightarrow \phi_{TS}&=\phi_{RS}\ast\phi_{TR}\\
\Rightarrow \phi_{TR}&=\phi_{TS}\tilde\ast\phi_{RS}
                      =(T\tilde\ast S)\tilde\ast (R\tilde\ast S)
                      =T\tilde\ast R\quad[=T\tilde\ast(S\ast\phi_{RS})]\\
\end{aligned}
$$ {#eq:label}
<!-- dprint-ignore-end -->

- [x] $\phi_{RS}=R\tilde\ast S$
- [x] $\phi_{TS}=T\tilde\ast S$
- [x] $\phi_{TR}=T\tilde\ast R = (T\tilde\ast S)\tilde\ast (R\tilde\ast S) =
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

![In the figures above, response functions are generated from (a) OB16, (b) CESM2
intermediate, (c) CESM2 large, (d) CESM2 extra large, (e) CESM2 double 2-year, and (f)
CESM2 double 4-year
simulations.](./generated_files/parametrisation/parametrisation_combined.png)

An average across the three estimation methods, with all datasets shown in the same
plot; an average across a subset of all simulation cases for all three estimation
methods:

![Average](./generated_files/parametrisation/parametrisation_ensemble.png){width=49%}
![Method means](./generated_files/parametrisation/parametrisation_method.png){width=49%}

From this, we find that the overall best estimates are form directly deconvolving $T$
with $R$. Evidently, it seems that noise is not the main issue, but rather the
misalignment between the signal and kernel that is used in the deconvolution. This is
backed up by the very poor results form using daily resolved data within the OB16
dataset, which just gets worse when smoothed out.

## Double waveform

We want to see how well the deconvolution is suited for estimating the response
functions in simulations where volcanic eruptions occur close in time.

### Double waveform response functions

We first compare the response functions estimated from the CESM2 double waveform
simulations, where a two year and four year separation between the volcanic eruptions
have been used.

![Response functions RF](./generated_files/waveform/responses_rf.png){#fig:responses_rf
width=49%}
![Response functions T](./generated_files/waveform/responses_temp.png){#fig:responses_temp
width=49%}

_The RF to SO~2~ and T to SO~2~ response functions for the double waveform simulations._

### Reconstructed double waveforms

Let us use the CESM2 strong eruption ensemble as a good estimate of the true response
functions, and then compare the climate model output of $R$ and $T$ with the
reconstructed ones from the CESM2 strong ensemble response functions.

![True versus reconstructed waveforms of (a) AOD, (b) RF, and (c) T](./generated_files/waveform/responses_combined.png){#waveform-comparison}

The response functions have here been rescaled to have amplitude equal to the amplitude
of the corresponding response functions generated from the simulation in question. The
plots give good indications that for the four year separated double waveform, the
reconstruction is indistinguishable from the true response functions at the current
noise level. However, the two year separated double waveform shows a clear non-linear
dependence between the injected SO~2~ and resulting $R$ and $T$ time series. A
significantly weaker response in both $R$ and $T$ is seen in the second eruption, most
notably in the plot of $T$.

### Reconstructed double waveforms, scaled

In our previous paper, we found a non-linear dependence on the peak amplitude of
temperature and radiative forcing to aerosol optical depth. Let us assume that the RF
and temperature combine depending on the level of AOD, and as such scale the forcing in
the convolution accordingly.

How? When the AOD is equal to zero, or when there is no perturbation and deviation from
equilibrium, the response is scaled by $1$. So if we scale by the difference between
unity and the ratio between the AOD at the time of the eruption and the maximum AOD, we
should get a better estimate of $R$ and $T$ (but of course not of AOD, which we assume
is roughly linearly dependent on injected SO~2~).

So, while we previously scaled the single eruption response function, say, from the
CESM2 strong simulation, to get estimates of $R$ and $T$ as (where $a$ is $R$ or $T$)

<!-- dprint-ignore-start -->
$$
a=\phi_{aS,\mathrm{strong}}
  \frac{\phi_{aS,\mathrm{double}}^{\max}}{\phi_{aS,\mathrm{strong}}^{\max}}
  *S,
$$
<!-- dprint-ignore-end -->

we now scale the SO~2~ value as well (which is equivalent to scaling the response
function a second time). This yields

<!-- dprint-ignore-start -->
$$
a^{\dagger}=\phi_{aS,\mathrm{strong}}
            \frac{\phi_{aS,\mathrm{double}}^{\max}}{\phi_{aS,\mathrm{strong}}^{\max}}
            *S\left(1-\frac{A}{A^{\max}}\right)^{1/2}.
$$
<!-- dprint-ignore-end -->

Since $S$ is a train of delta pulses, multiplying by the AOD is the same as multiplying
by the AOD values at the time of the eruptions. Fig. \ref{waveform-comparison-scaled}
shows the results of this scaling, which are far better than I expected, albeit with an
additional square root that needs to be addressed.

![True versus reconstructed and AOD square root corrected waveforms of (a) AOD, (b) RF,
and (c)
T](./generated_files/waveform/responses_combined-aod-root-corrected.png){#waveform-comparison-scaled}

Can we do better? Let us try to use a logarithmic scaling instead, and see if this gives
a better result or at least a better understanding of the scaling. We again use the "one
minus AOD ratio" as the starting point of our scaling factor, so we need a mapping that
preserves the endpoint values of $0$ and $1$. The logarithmic scaling is then given by
$\log(1+1-A/A^{\max})/\log(2)$. (A third option could be $1-\log(1+A/A^{\max})/\log(2)$,
but this yields a too small scaling.)

<!-- dprint-ignore-start -->
$$
a^{\dagger}=\phi_{aS,\mathrm{strong}}
            \frac{\phi_{aS,\mathrm{double}}^{\max}}{\phi_{aS,\mathrm{strong}}^{\max}}
            *S\frac{\log\left(1+1-\frac{A}{A^{\max}}\right)}{\log(2)}.
$$
<!-- dprint-ignore-end -->

![True versus reconstructed and AOD logarithmic corrected waveforms of (a) AOD, (b) RF,
and (c)
T](./generated_files/waveform/responses_combined-aod-log-corrected.png){#waveform-comparison-scaled}

### Unwrapping the scaling

There are a few things to note from this procedure. We not only scale the SO~2~
injections, we also scale the single waveform response function to have the same maximum
as the response function from the double waveform simulation. Then, after this initial
scaling, we compare the scale of the SO~2~ time series according to the AOD time series
to manipulate the amplitude of the response function depending on the state of the AOD
at the time of the eruption.

So, if we assume for a moment that the peak temperature response goes as the square root
of the injected SO~2~. Then the response function from eruption $\alpha$ would be scaled
as $\phi_{aS}^\alpha=\sqrt{S^\beta/S^\alpha}\phi_{aS}^\beta$ to represent the response
function from eruption $\beta$. We get the inverse expression in the ratio of SO~2~s
since the deconvolution mimic a division in the Fourier domain. Combining this with the
above scaling, we get

#### Scale the response function

$$ \phi_{aS}^\alpha=\sqrt{S^\beta/S^\alpha}\phi_{aS}^\beta $$

or more pragmatically

$$ \phi_{aS}^\alpha=(\phi_{aS}^{\alpha,\max}/\phi_{aS}^{\beta,\max})\phi_{aS}^\beta. $$

#### Scale the time series

<!-- dprint-ignore-start -->
$$
a^{\alpha,\dagger}=
  \phi_{aS}^\beta
  \frac{\phi_{aS}^{\alpha,\max}}{\phi_{aS}^{\beta,\max}}
  *[S^\alpha]_i\frac{\log\left(1+1-\frac{[A^\alpha]_i}{A^{\alpha,\max}}\right)}{\log(2)}
  \approx
  \phi_{aS}^\beta \sqrt{\frac{S^\beta}{S^\alpha}}
  *[S^\alpha]_i
  \frac{\log\left(1+1-\frac{[A^\alpha]_i}{A^{\alpha,\max}}\right)}{\log(2)}.
$$
<!-- dprint-ignore-end -->

### Numerical solution

We here assume that we can get the SO~2~ content from a simple exponential decay model
as

<!-- dprint-ignore-start -->
$$
S(t)=C_1\int_{0}^{t}\exp\left(-\frac{t-t'}{\tau_S}\right)\sum_k S_k\delta(t'-t_k)dt',
$$
<!-- dprint-ignore-end -->

and further that the AOD is given by a second exponential decay model as

<!-- dprint-ignore-start -->
$$
A(t)=\int_{0}^{t}C_1\exp\left(-\frac{t-t'}{\tau_A}\right)S(t)dt'.
$$
<!-- dprint-ignore-end -->

We assume a logarithmic relation between the RF and the AOD, yielding

<!-- dprint-ignore-start -->
$$
R(A(t))=C_1\log(1+C_2A(t)),
$$
<!-- dprint-ignore-end -->

where $C_i$ are constants.

### But wait, there is more

We actually get pretty good results _if we optimize the RF signal based on SO~2~ input
via the AOD estimation_, which is not really what we wanted! Ideally, we would get a
good AOD estimate from SO~2~ input, and then a good RF estimate from using the AOD
estimate as input. Similarly, we get a good AOD estimate when optimising the AOD based
on SO~2~ as input, via an intermediate AOD estimate.

Thus, what we are left with are that the best estimate RF and AOD are given as $R$ and
$A^\dagger$ as

$$ R(t)=R(A(S(t))), $$

and

$$ A^\dagger(t)=A^\dagger(A(S(t))). $$

![SO2 fit for CESM small](./generated_files/relationships/numerical_so2_cesm-cesm2-medium_combined.png)

![AOD fit for CESM small](./generated_files/relationships/numerical_aod_cesm-cesm2-medium_combined.png)

![RF fit for CESM small](./generated_files/relationships/numerical_rf_cesm-cesm2-medium_combined.png)

![SO2 fit for CESM intermediate](./generated_files/relationships/numerical_so2_cesm-cesm2-medium-plus_combined.png)

![AOD fit for CESM intermediate](./generated_files/relationships/numerical_aod_cesm-cesm2-medium-plus_combined.png)

![RF fit for CESM intermediate](./generated_files/relationships/numerical_rf_cesm-cesm2-medium-plus_combined.png)

## Cut off response function

To check how far into the response functions the noise is substantial, we cut the
response function at a certain time lag, before we generate an ensemble of temperature
output signals by first convolving the cut response function with the forcing, and then
adding phase shifted noise to the temperature signal. The procedure can be summarise in
the following steps:

1. Cut the response function at a given index.
2. Convolve the cut response function with the forcing signal to get a temperature
   estimate $T_{\mathrm{est}}$.
3. Add noise represented by a phase shifted temperature control signal to create an
   ensemble of temperature estimates $T_i$.
4. Deconvolve the ensemble of temperature estimates with the forcing signal to create an
   ensemble of response functions for the given cut off index.

### Inspecting the noise floor

We get estimates of the noise present in the response functions from doing the
bootstrapping method outlined above. The noise is then illustrated as percentile plots
that show the percentiles as layers of shaded filled in area plots.

Let us look at how the cut off ensemble looks like for the OB16 dataset when cutting the
response function off at years 2, 4, 8 and 16 (specified in months in the plot labels).

#### Noise of OB16

We find that when estimating the temperature time series back from the radiative forcing
and the temperature to radiative forcing response function, most of the features are
resolved even when cutting off at two years. Keeping the first four years also captures
more irregularities during the decay. However, there might be considerable contributions
of the decay from the largest eruption in the dataset, which also features a large bump
in during the decay, between roughly 2 and 4 years.

For example, the eruption in about the year $1340$ is well captured in all scenarios,
while the largest as well as the two eruption following closely in time deviates more.

![Cut off response function OB16, $T$ against $R$](./generated_files/cut_off/ob16-month_resp_temp-rf_combined.png)

#### Noise of CESM2 strong

The time series from the CESM2 simulations are not long enough to fully capture the
temperature decay. It is evident that the temperature response to RF is not fully
decayed after $20$ years, but the shape seems rather well captured up to the point where
the time series ends. Noise come mostly from having a kernel in the deconvolution that
is non-delta like.

![Cut off response function CESM2 strong, $T$ against $R$](./generated_files/cut_off/cesm2-strong_resp_temp-rf_combined.png)

#### Noise of CESM2 2-sep

Similarly to the strong CESM2 simulation, the response function is not fully decayed
after the $12$ years available in the time series. The noise is more substantial in the
double waveform simulations due to a smaller ensemble. We do however get very different
response function shapes, where the $2$ year separated double waveform response has a
much faster decay.

![Cut off response function CESM2 2-sep, $T$ against $R$](./generated_files/cut_off/cesm2-tt-2sep_resp_temp-rf_combined.png)

#### Noise of CESM2 4-sep

As the $4$ year separated double waveform response function contain only a single
ensemble member, the noise is substantial. However, we can recreate both the RF and
temperature time series using the raw SO~2~ time series, and the response functions from
CESM2 strong of the RF and temperature to SO~2~. Therefore, we expect the shape of the
response of temperature to RF to be similar to what obtained from CESM2 strong.

![Cut off response function CESM2 4-sep, $T$ against $R$](./generated_files/cut_off/cesm2-double-overlap_resp_temp-rf_combined.png)

### Cut-off reconstructions

Let us finally try to use the cut off response functions to recreate the temperature
time series, and do the same reconstruction analysis as we started off with. The cuts
are the same as before, at years 2, 4, 8, and 16.

![Reconstructed from cut-offs at (a) 2, (b) 4, (c) 8, and (d) 16 years; temperature](./generated_files/reconstruction/cut-off-temp-combined.png)

![Reconstructed from cut-offs at (a) 2, (b) 4, (c) 8, and (d) 16 years; power spectra](./generated_files/reconstruction/cut-off-power-combined.png)

![Reconstructed from cut-offs at (a) 2, (b) 4, (c) 8, and (d) 16 years; PDF](./generated_files/reconstruction/cut-off-pdf-combined.png)

![Reconstructed from cut-offs at (a) 2, (b) 4, (c) 8, and (d) 16 years; CDF](./generated_files/reconstruction/cut-off-cdf-combined.png)

![Reconstructed from cut-offs at (a) 2, (b) 4, (c) 8, and (d) 16 years; correlation](./generated_files/reconstruction/cut-off-correlation-combined.png)

\clearpage{}

## References
