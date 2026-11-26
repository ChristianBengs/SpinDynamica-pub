---
layout: default
title: Home
section: home
permalink: /
---
<section class="hero">
  <h1 class="hero-title">{{ site.title }}</h1>
  <p class="hero-tagline">{{ site.tagline }}</p>
  <p class="hero-summary">
    SpinDynamica is an open-source set of Wolfram Mathematica packages for symbolic and numerical
    spin dynamics calculations in magnetic resonance. It combines exact quantum-mechanical formalisms
    in Hilbert and Liouville space with Mathematica's symbolic algebra, numerical solvers, and graphics.
  </p>

  <div class="button-group">
    {% include download-button.html %}
    <a class="btn btn-secondary" href="{{ site.repository_url }}">GitHub Repository</a>
    <a class="btn btn-secondary" href="{{ '/installation/' | relative_url }}">Installation Guide</a>
  </div>

  {% assign rel = site.data.release %}
  {% if rel.published and rel.version != "" %}
    <p class="release-meta">
      Latest stable release: <strong>{{ rel.version }}</strong>
      {% if rel.release_date != "" %}({{ rel.release_date }}){% endif %}
      &middot; Requires Mathematica {{ rel.min_mathematica }}+ (tested under {{ rel.tested_mathematica }})
      &middot; <a href="{{ '/download/' | relative_url }}">Release details</a>
    </p>
  {% else %}
    <p class="release-meta">
      First versioned GitHub Release pending &middot; Requires Mathematica {{ rel.min_mathematica }}+
      (tested under {{ rel.tested_mathematica }}) &middot;
      <a href="{{ '/download/' | relative_url }}">Download options</a>
    </p>
  {% endif %}
</section>

<section>
  <h2>Introduction</h2>
  <p>
    SpinDynamica provides a unified computational environment for magnetic resonance theory and simulation.
    Users can define spin systems, construct spin Hamiltonians and Liouvillians symbolically or numerically,
    propagate quantum states and density operators through pulse sequences, simulate one- and two-dimensional
    spectra, and design optimal-control radiofrequency waveforms within a single Mathematica notebook.
  </p>
</section>

<section>
  <h2>Capabilities</h2>
  <div class="capability-grid">
    <div class="capability-card">
      <h3>Hilbert &amp; Liouville Space</h3>
      <p>
        Flexible spin-system definitions, Zeeman and custom Hilbert-space bases, spin angular momentum
        and shift operators, operator bases, and Liouville-space superoperators.
      </p>
    </div>
    <div class="capability-card">
      <h3>Interactions &amp; Relaxation</h3>
      <p>
        Symbolic and numerical construction of spin interactions, Hamiltonians, and relaxation
        superoperators in the Redfield regime.
      </p>
    </div>
    <div class="capability-card">
      <h3>Propagation &amp; Pulse Sequences</h3>
      <p>
        Analytical and numerical propagation using matrix exponentials, diagonalization, or ODE integration
        over arbitrary event lists and pulse sequences, alongside perturbation approximations.
      </p>
    </div>
    <div class="capability-card">
      <h3>1D &amp; 2D NMR Simulations</h3>
      <p>
        Expectation-value trajectories, transformation amplitudes, 1D and 2D signal simulation,
        polarization dynamics, ensemble averaging, and orientational sampling schemes.
      </p>
    </div>
    <div class="capability-card">
      <h3>Signals &amp; Bruker Import</h3>
      <p>
        Time-domain signal processing, Fourier transformation, lineshape functions, and direct import
        and plotting of processed 1D Bruker spectra.
      </p>
    </div>
    <div class="capability-card">
      <h3>Optimal Control &amp; Optimization</h3>
      <p>
        GRAPE optimal control with analytical gradients, time-optimal and robust ensemble pulse design,
        and compiled transformation-amplitude optimization across pulse-sequence parameters.
      </p>
    </div>
    <div class="capability-card">
      <h3>Rotations, Tensors &amp; Symmetry</h3>
      <p>
        Euler angles, Wigner D-matrices, quaternions, Clebsch&ndash;Gordan coefficients, irreducible
        spherical tensors, Lie theory, commutants, and spin permutation groups.
      </p>
    </div>
    <div class="capability-card">
      <h3>Visualization</h3>
      <p>
        Integrated plotting of spin trajectories, optimal-control waveforms and convergence histories,
        and 3D spherical tensor ovaloids using Mathematica graphics.
      </p>
    </div>
  </div>
</section>

<section id="whats-new">
  <h2>What&rsquo;s New in v4.0.0</h2>
  <p>
    Version 4.0.0 introduces major additions for experimental data import, pulse-sequence optimal control,
    and high-performance two-dimensional NMR simulations:
  </p>
  {% include whats-new.html %}
  <p>
    Explore executable walkthroughs of these features in the <a href="{{ '/examples/' | relative_url }}">Examples</a> section.
  </p>
</section>

<section>
  <h2>Getting Started</h2>
  <p>
    To get started, obtain the package from the <a href="{{ '/download/' | relative_url }}">Download</a> page,
    configure Mathematica's <code>$Path</code> following the <a href="{{ '/installation/' | relative_url }}">Installation</a> instructions,
    and explore the <a href="{{ '/documentation/' | relative_url }}">Documentation</a> and <a href="{{ '/examples/' | relative_url }}">Examples</a>.
  </p>
</section>

<section id="how-to-cite">
  <h2>How to Cite SpinDynamica</h2>
  <p>
    If you use SpinDynamica in published research, please cite the following article:
  </p>
  {% include citation.html %}
</section>
