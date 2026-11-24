---
title: Getting Started and Conventions
order: 1
summary: Loading SpinDynamica, setting up spin systems, inspecting options, and chronology conventions.
---

This guide summarizes the core conventions used throughout SpinDynamica, based on
Part 1 (*Introduction to SpinDynamica*) and Part 2 (*Simulations*) of the Mathematica documentation notebooks.

## Loading SpinDynamica

Once the `SpinDynamica` directory has been added to `$Path` (see [Installation]({{ '/installation/' | relative_url }})), load the package with:

```mathematica
Needs["SpinDynamica`"]
```

You can confirm the loaded package version at any time using:

```mathematica
$SpinDynamicaVersion
```

## Defining a Spin System

Before building spin operators or running simulations, define the active spin system with `SetSpinSystem`. Shortcuts are provided for systems of spins-1/2, while explicit pairs `{label, spinQuantumNumber}` support arbitrary spin quantum numbers:

{% raw %}
```mathematica
(* Single spin-1/2 with label 1 *)
SetSpinSystem[1]

(* Two spins-1/2 with labels 1 and 2 *)
SetSpinSystem[{1, 2}]

(* Spin-1/2 coupled to a spin-1 nucleus *)
SetSpinSystem[{{1, 1/2}, {2, 1}}]
```
{% endraw %}

Evaluate `SpinSystem[]` to inspect the currently active spin system.

## Spin Operators

Cartesian and shift angular momentum operators are constructed with `opI`:

- `opI[1, "x"]`, `opI[1, "y"]`, `opI[1, "z"]` construct Cartesian angular momentum operators for spin `1`.
- `opI[1, "+"]` and `opI[1, "-"]` construct raising and lowering shift operators.
- `opI["z"]` (omitting the spin label) sums over all spins in the current spin system.
- `opI[1] . opI[2]` denotes the scalar operator product between spins `1` and `2`.

## Inspecting Routines and Options

SpinDynamica exposes clickable symbol listings and standard Wolfram Language help queries:

```mathematica
(* Clickable overview of all top-level SpinDynamica symbols *)
SpinDynamicaSymbols[]

(* Usage message and default options for a specific routine *)
?Trajectory
Options[Trajectory]
```

Optional arguments to SpinDynamica functions are passed as standard Mathematica replacement rules such as `BackgroundGenerator -> ...` or `EnsembleAverage -> ...`.

## Chronology Convention

When a pulse sequence or event list `{e1, e2, ...}` is evaluated, the temporal ordering is controlled by the `Chronology` option and the global setting `$Chronology`:

- `"LeftToRight"` (the default): `e1` is applied first, followed by `e2`, matching standard pulse-sequence diagrams.
- `"RightToLeft"`: `e2` acts before `e1`, matching the order in which matrix operators act from the right.

## Next Steps

- Browse the [Mathematica Documentation Notebooks]({{ '/documentation/' | relative_url }}) for full module-by-module coverage.
- Review the [Worked Examples]({{ '/examples/' | relative_url }}) for executable simulations and optimal-control workflows.
