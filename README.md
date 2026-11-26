# SpinDynamica

**Symbolic and numerical spin dynamics in Wolfram Mathematica.**

SpinDynamica is an open-source Wolfram Mathematica package suite for symbolic and numerical calculations in magnetic resonance, combining Hilbert- and Liouville-space formalisms with Mathematica's symbolic algebra, numerical solvers, and graphics.

- **Website:** https://christianbengs.github.io/SpinDynamica-pub/
- **Releases:** https://github.com/ChristianBengs/SpinDynamica-pub/releases
- **Issues:** https://github.com/ChristianBengs/SpinDynamica-pub/issues

## Highlights in v4.0.0

- **Bruker DataLoader:** Import, inspect, plot, and convert processed 1D Bruker NMR spectra (`ImportBrukerSpectrum`, `PlotBrukerSpectrum`, `ToSignal`).
- **Optimal Control:** Analytical-gradient GRAPE pulse optimization (`OptimalControlProblem`, `GRAPEOptimize`), robust ensemble design, constrained/second-order optimizers, and time-optimal control.
- **Signal2D & Performance:** High-performance 2D NMR simulation (`Signal2D`) with Hilbert/Liouville blocked double-diagonalization, hypercomplex/TPPI quadrature, periodic MAS `COMPUTE2D`, compiled kernels, and 2D processing/plotting (`FT`, `PlotSpectrum2D`).
- **Wolfram Paclet Distribution:** Packaged as a `.paclet` archive for straightforward installation via `PacletInstall`.

## Installation

Requires **Wolfram Mathematica 8.0 or later** (tested under Mathematica 14.1–15.0).

### Option 1: Wolfram Paclet (Recommended)

Download `SpinDynamica-4.0.0.paclet` from [GitHub Releases](https://github.com/ChristianBengs/SpinDynamica-pub/releases) and evaluate in Mathematica:

```mathematica
PacletInstall["/path/to/SpinDynamica-4.0.0.paclet"]
Needs["SpinDynamica`"]
```

### Option 2: Manual `$Path` Setup

Clone or unpack the repository and add the inner `SpinDynamica/SpinDynamica` directory (containing `SpinDynamica.m`) to `$Path`:

```mathematica
PrependTo[$Path, "/path/to/SpinDynamica/SpinDynamica"];
Needs["SpinDynamica`"]
```

Verify the installation with:

```mathematica
$SpinDynamicaVersion
SpinDynamicaSymbols[]
```

See the [Installation Guide](https://christianbengs.github.io/SpinDynamica-pub/installation/) for persistent `init.m` setup and additional details.

## Repository Structure

- `SpinDynamica/SpinDynamica/` — Core Mathematica package source (`SpinDynamica.m`, `PacletInfo.wl`, and submodules).
- `SpinDynamica/documentation/` — Core Mathematica documentation notebooks (`SDv4.0.0-part1-Intro.nb` through `SDv4.0.0-part12-Groups+Algebras.nb`).
- `SpinDynamica/Examples/` — Interactive tutorial notebooks covering `Signal2D` (COSY, HSQC, DQF-COSY), `OptimalControl` (GRAPE, ensembles, time-optimal control), and `Bruker_DataLoader`.
- `docs/` — Source for the GitHub Pages website.

## How to Cite

If you use SpinDynamica in published research, please cite:

> Christian Bengs and Malcolm H. Levitt, *"SpinDynamica: Symbolic and numerical magnetic resonance in a Mathematica environment"*, *Magnetic Resonance in Chemistry*, **56**(6), 374–414 (2018). DOI: [10.1002/mrc.4642](https://doi.org/10.1002/mrc.4642)

```bibtex
@article{Bengs2018SpinDynamica,
  author  = {Bengs, Christian and Levitt, Malcolm H.},
  title   = {SpinDynamica: Symbolic and numerical magnetic resonance in a Mathematica environment},
  journal = {Magnetic Resonance in Chemistry},
  volume  = {56},
  number  = {6},
  pages   = {374--414},
  year    = {2018},
  doi     = {10.1002/mrc.4642}
}
```

## Authors & Acknowledgements

- **Authors:** Malcolm H. Levitt and Christian Bengs
- **Acknowledgements:** Recent changes and additions were initiated by Mohamed Sabba and Christian Bengs and developed with the assistance of Astra and Antigravity.

## License

SpinDynamica is open-source software. A formal `LICENSE` file will be added to the repository.
