---
title: Single-Spin Trajectory under RF Nutation
order: 1
summary: Simulating spin-1/2 magnetization evolution under a constant radiofrequency field using Trajectory.
---

This example computes the time evolution of longitudinal and transverse magnetization for an isolated spin-1/2 nucleus subjected to a resonant radiofrequency (RF) field along the $$x$$-axis.

## Physical Model

Consider a single spin-1/2 initialized in longitudinal magnetization $$\hat{\rho}(0) = \hat{I}_z$$ and evolving under a time-independent RF Hamiltonian along $$x$$:

$$
\hat{H} = \omega_1 \hat{I}_x, \qquad \omega_1 = 2\pi\ \mathrm{rad\,s^{-1}}.
$$

In Hilbert space, the density operator evolves according to the Liouville&ndash;von Neumann equation:

$$
\frac{\mathrm{d}\hat{\rho}(t)}{\mathrm{d}t} = -i \left[\hat{H}, \hat{\rho}(t)\right],
\qquad
\hat{\rho}(t) = \exp\!\left(-i \hat{H} t\right) \hat{\rho}(0) \exp\!\left(i \hat{H} t\right).
$$

For an observable operator $$\hat{Q}$$, SpinDynamica's `Trajectory` returns the normalized expectation value formed via the Liouville bracket $$(\hat{Q} \mid \hat{\rho}(t)) / (\hat{Q} \mid \hat{Q})$$. Under a $$1\ \mathrm{Hz}$$ nutation field ($$\omega_1 = 2\pi\ \mathrm{rad\,s^{-1}}$$), the normalized observables follow:

$$
\langle \hat{I}_z \rangle(t) = \cos(\omega_1 t),
\qquad
\langle \hat{I}_y \rangle(t) = -\sin(\omega_1 t).
$$

## Mathematica Implementation

```mathematica
Needs["SpinDynamica`"];

(* Define an isolated spin-1/2 system *)
SetSpinSystem[1];

(* Compute TrajectoryFunction objects for Iz and Iy over 1 s under a 1 Hz RF field *)
{tz, ty} = Trajectory[
  opI["z"] -> {opI["z"], opI["y"]},
  {2 Pi opI["x"], 1}
];

(* Plot both expectation-value trajectories over one nutation period *)
Plot[{tz[t], ty[t]}, {t, 0, 1},
  PlotLegends -> {"<Iz>(t)", "<Iy>(t)"},
  FrameLabel -> {"t / s", "Normalized expectation value"}]
```

## Further Reading

- [Part 2: Simulations (`SDv4.0.0-part2-Simulations.nb`)]({{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/documentation/SDv4.0.0-part2-Simulations.nb) covers `Trajectory`, `TransformationAmplitude`, and `Signal1D` in detail.
