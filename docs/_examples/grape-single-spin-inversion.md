---
title: Single-Spin Inversion via GRAPE Optimal Control
order: 2
summary: Designing an amplitude-bounded state-to-state inversion pulse with analytical GRAPE gradients and independent verification.
---

SpinDynamica includes an optimal-control module (`OptimalControl` / `GRAPE`) for designing shaped pulses using analytical gradients and quasi-Newton or second-order optimizers. This example reproduces the workflow from the Optimal Control tutorial series for inverting an off-resonance spin-1/2 state.

## Problem Formulation

We consider a single spin-1/2 subject to an off-resonance drift Hamiltonian $$\hat{H}_0$$ and two orthogonal RF control channels along $$x$$ and $$y$$ over duration $$T = 10\ \mathrm{ms}$$ divided into $$N = 30$$ piecewise-constant time slices:

$$
\hat{H}(t) = \hat{H}_0 + u_x(t)\hat{I}_x + u_y(t)\hat{I}_y,
\qquad
\frac{\hat{H}_0}{2\pi} = 50\ \mathrm{Hz} \cdot \hat{I}_z.
$$

Each quadrature control amplitude is bounded by $$\lvert u_x(t)\rvert, \lvert u_y(t)\rvert
\le 120\ \mathrm{Hz}$$. The optimization drives the Zeeman basis state $$\lvert\alpha\rangle$$ to $$\lvert\beta\rangle$$, maximizing the pure-state fidelity:

$$
F = \left|\langle \beta \mid \hat{U}(T) \mid \alpha \rangle\right|^2.
$$

## Mathematica Implementation

```mathematica
Needs["SpinDynamica`"];

(* 1. Define a single spin-1/2 system and retrieve Zeeman basis kets {|alpha>, |beta>} *)
SetSpinSystem[1];
b1 = ZeemanBasis[];
ks = BasisKets[b1];

(* 2. Define +50 Hz off-resonance drift and x/y RF control operators in Hz *)
H0 = 50. opI[1, "z"];
Hc = {opI[1, "x"], opI[1, "y"]};

(* 3. Specify the state-to-state target over 10 ms divided into 30 slices *)
target = StateTarget[ks[[1]] -> ks[[2]]];
prob = OptimalControlProblem[H0, Hc, target, {0.01, 30},
  FrequencyUnits -> "Hz", ControlUnits -> "Hz"];

(* 4. Verify analytical GRAPE gradients against central finite differences *)
uInit = ConstantArray[{25., -15.}, 30];
gradChk = VerifyOptimalControlGradient[prob, uInit];

(* 5. Optimize with bounded L-BFGS up to target fidelity 0.9999 *)
res = GRAPEOptimize[prob,
  AmplitudeBound -> 120.,
  TargetFidelity -> 0.9999,
  MaxIterations -> 100,
  RandomSeed -> 101];

(* 6. Independently verify against SpinDynamica's NPropagationOperator *)
sim = VerifyOptimalControlSimulation[res];
wf = OptimalControlWaveform[res];

(* 7. Visualize waveform, convergence, and Bloch-vector trajectory *)
PlotOptimalControlWaveform[res]
PlotOptimalControlWaveform[res, Method -> "AmplitudePhase"]
PlotOptimalControlConvergence[res, Method -> "Infidelity"]
PlotOptimalControlTrajectory[res, ks[[1]], {opI[1, "x"], opI[1, "y"], opI[1, "z"]}]
```

## Notes

- `VerifyOptimalControlGradient` compares exact analytical GRAPE gradients against central finite differences before optimization.
- `VerifyOptimalControlSimulation` exports the optimized pulse via `OptimalControlEventList` and re-propagates it independently with `NPropagationOperator`.
- See the full interactive notebook [01_GRAPE_and_SingleSpin_Rotations.nb]({{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/Examples/OptimalControl/Notebooks/01_GRAPE_and_SingleSpin_Rotations.nb) on GitHub for additional explanations and variations.
