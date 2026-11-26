---
layout: page
title: Examples
summary: Worked SpinDynamica web examples and interactive Mathematica tutorial notebooks.
section: examples
permalink: /examples/
---

## Worked Web Examples

{% assign examples = site.examples | sort: 'order' %}
{% if examples.size > 0 %}
<ul class="card-list">
  {% for ex in examples %}
  <li class="card-list-item">
    <a class="card-list-title" href="{{ ex.url | relative_url }}">{{ ex.title }}</a>
    {% if ex.summary %}
    <p class="card-list-summary">{{ ex.summary }}</p>
    {% endif %}
  </li>
  {% endfor %}
</ul>
{% endif %}

<div class="placeholder-note">
  <p>
    <span class="placeholder-label">More examples planned:</span>
    Additional worked web examples covering coupled spin systems, relaxation, and pulse-sequence optimization will be added here.
  </p>
</div>

## Interactive Mathematica Tutorial Notebooks

The repository contains interactive Mathematica tutorial notebooks under <code>SpinDynamica/Examples/</code>, organized by area:

### 2D NMR Simulations (`SpinDynamica/Examples/Signal2D/`)

<ul>
  <li><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/Examples/Signal2D/Signal2D_Tutorial_01_COSY.nb"><code>Signal2D_Tutorial_01_COSY.nb</code></a> &mdash; 2D COSY simulation, hypercomplex Fourier transform, and contour plotting</li>
  <li><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/Examples/Signal2D/Signal2D_Tutorial_02_HSQC.nb"><code>Signal2D_Tutorial_02_HSQC.nb</code></a> &mdash; Heteronuclear 2D HSQC simulation and indirect quadrature encoding</li>
  <li><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/Examples/Signal2D/Signal2D_Tutorial_03_DQF_COSY.nb"><code>Signal2D_Tutorial_03_DQF_COSY.nb</code></a> &mdash; Double-quantum-filtered COSY (DQF-COSY) with phase cycling and coherence-pathway selection</li>
</ul>

### Optimal Control (`SpinDynamica/Examples/OptimalControl/`)

<ul>
  <li><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/Examples/OptimalControl/OptimalControl_Tutorial_01_SingleSpin_GRAPE.nb"><code>OptimalControl_Tutorial_01_SingleSpin_GRAPE.nb</code></a> &mdash; Single-spin GRAPE optimization, analytical gradients, and verification</li>
  <li><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/Examples/OptimalControl/OptimalControl_Tutorial_02_TwoSpin_Ensembles_and_TimeOptimal.nb"><code>OptimalControl_Tutorial_02_TwoSpin_Ensembles_and_TimeOptimal.nb</code></a> &mdash; Two-spin transfer, robust ensemble optimization, and time-optimal control</li>
</ul>

### Bruker DataLoader (`SpinDynamica/Examples/Bruker_DataLoader/`)

<ul>
  <li><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/Examples/Bruker_DataLoader/BrukerDataLoader_Tutorial_01_1DImport.nb"><code>BrukerDataLoader_Tutorial_01_1DImport.nb</code></a> &mdash; Importing processed 1D Bruker spectra, inspecting metadata, plotting, and converting to SpinDynamica signals</li>
</ul>
