---
title: 2D TOCSY with DIPSI-2 Isotropic Mixing
order: 3
summary: Simulating a phase-sensitive 2D TOCSY spectrum of a three-spin linear chain using an explicit DIPSI-2 RF mixing sequence in Signal2D.
---

Total Correlation Spectroscopy (TOCSY) transfers longitudinal or transverse magnetization across an entire coupled spin network via isotropic mixing. Unlike COSY, which produces cross-peaks only between directly coupled spins, TOCSY generates both direct and relayed cross-peaks across multi-step coupling pathways.

## Physical Model

Consider a linear three-spin-1/2 (AMX) chain with chemical-shift offsets $$\nu_1 = -120\ \mathrm{Hz}$$, $$\nu_2 = 0\ \mathrm{Hz}$$, and $$\nu_3 = +160\ \mathrm{Hz}$$, coupled sequentially by $$J_{12} = 10\ \mathrm{Hz}$$ and $$J_{23} = 8\ \mathrm{Hz}$$, with **no direct coupling** between spins 1 and 3 ($$J_{13} = 0$$):

$$
\hat{H}_0 = 2\pi \left(
  \nu_1 \hat{I}_{1z} + \nu_2 \hat{I}_{2z} + \nu_3 \hat{I}_{3z}
  + J_{12}\,\hat{\mathbf{I}}_1 \cdot \hat{\mathbf{I}}_2
  + J_{23}\,\hat{\mathbf{I}}_2 \cdot \hat{\mathbf{I}}_3
\right).
$$

During the mixing period, an explicit **DIPSI-2** radiofrequency sequence (Shaka, Lee, and Pines, *J. Magn. Reson.* **77**, 274&ndash;293, 1988) is applied at nutation frequency $$\nu_{\mathrm{RF}} = 2.5\ \mathrm{kHz}$$. Its basic composite element $$R$$ comprises nine contiguous pulses of alternating $$+x / -x$$ phase:

$$
R = 320^\circ_x \; 410^\circ_{-x} \; 290^\circ_x \; 285^\circ_{-x} \; 30^\circ_x \; 245^\circ_{-x} \; 375^\circ_x \; 265^\circ_{-x} \; 370^\circ_x,
$$

assembled into the four-step supercycle $$R\,\bar{R}\,\bar{R}\,R$$ (duration $$\approx 11.51\ \mathrm{ms}$$ per cycle at $$2.5\ \mathrm{kHz}$$). Repeating the supercycle five times yields a mixing time of $$\tau_m \approx 57.6\ \mathrm{ms}$$, suppressing chemical-shift differences while preserving the scalar-coupling Hamiltonian so that magnetization relays from spin 1 through spin 2 to spin 3.

## Mathematica Implementation

{% raw %}
```mathematica
Needs["SpinDynamica`"];

(* 1. Define a linear three-spin-1/2 chain (1-2-3) with J13 = 0 *)
SetSpinSystem[{1, 2, 3}];

nu1 = -120.; nu2 = 0.; nu3 = 160.;
J12 = 10.; J23 = 8.;

H0 = 2 Pi (
  nu1 opI[1, "z"] + nu2 opI[2, "z"] + nu3 opI[3, "z"] +
  J12 opI[1] . opI[2] + J23 opI[2] . opI[3]
);

(* 2. Construct the explicit DIPSI-2 RF mixing supercycle R Rbar Rbar R *)
nuRF = 2500.;
wRF = 2 Pi nuRF;
dipsi2Angles = {320, 410, 290, 285, 30, 245, 375, 265, 370} * Degree;

dipsi2Block[sign_] := MapIndexed[
  {sign * (-1)^(First[#2] - 1) * wRF * opI["x"], #1 / wRF} &,
  dipsi2Angles
];

dipsi2Cycle = Join[
  dipsi2Block[1], dipsi2Block[-1], dipsi2Block[-1], dipsi2Block[1]
];
nCycles = 5;  (* Total mixing time ~ 57.6 ms *)

(* 3. Simulate phase-sensitive 2D States-TOCSY with z-filtered DIPSI-2 mixing *)
sw = 2 Pi * 500.;
sig2D = Signal2D[
  {{{sw, 64}}, {{sw, 128}}},
  BackgroundGenerator -> H0,
  Preparation -> Pulse[{Pi/2, "x"}],
  Mixing -> {Pulse[{Pi/2, Pi}], Repeat[dipsi2Cycle, nCycles], Pulse[{Pi/2, "x"}]},
  Observable -> opI["-"],
  QuadratureMethod -> "States",
  LineBroadening -> {8., 8.}
];

(* 4. 2D hypercomplex Fourier transform with pure-absorption phasing *)
spec2D = FT[sig2D,
  ZeroFill -> {2, 2},
  Apodization -> {{"SquaredSineBell", Pi/2}, {"SquaredSineBell", Pi/2}},
  PhaseCorrection -> {Pi/2, 0}
];

(* 5. Plot the 2D contour spectrum with 1D skyline projections *)
PlotSpectrum2D[spec2D,
  ReverseAxes -> True,
  BaseContourLevel -> 0.05,
  ShowProjections -> True]
```
{% endraw %}

## Interpreting the Spectrum

- **Diagonal peaks** appear at $$(\nu_1, \nu_1)$$, $$(\nu_2, \nu_2)$$, and $$(\nu_3, \nu_3)$$.
- **Direct cross-peaks** connect directly coupled partners at $$(\nu_1, \nu_2)$$, $$(\nu_2, \nu_1)$$, $$(\nu_2, \nu_3)$$, and $$(\nu_3, \nu_2)$$.
- **Relayed cross-peaks** appear cleanly in-phase at $$(\nu_1, \nu_3)$$ and $$(\nu_3, \nu_1)$$ despite $$J_{13} = 0$$, demonstrating multi-step isotropic mixing along the $$1 \leftrightarrow 2 \leftrightarrow 3$$ chain.

## Further Reading

- [Signal2D_Tutorial_01_COSY.nb]({{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/Examples/Signal2D/Signal2D_Tutorial_01_COSY.nb) — Introduction to `Signal2D`, `Spectrum2D`, phasing, and `PlotSpectrum2D`.
- [Signal2D_Tutorial_03_DQF_COSY.nb]({{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/Examples/Signal2D/Signal2D_Tutorial_03_DQF_COSY.nb) — Phase cycling and coherence-pathway selection in 2D NMR.
