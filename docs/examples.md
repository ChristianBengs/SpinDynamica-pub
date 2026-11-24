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
    Additional worked web examples covering coupled spin systems, relaxation, 2D NMR, and general pulse-sequence optimization will be added here.
  </p>
</div>

## Interactive Mathematica Tutorial Notebooks

The repository contains tutorial and walkthrough notebooks under <code>SpinDynamica/Examples/</code>, organized by area:

### Optimal Control (`SpinDynamica/Examples/OptimalControl/`)

<ul>
  <li><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/Examples/OptimalControl/OptimalControl_Walkthrough.nb"><code>OptimalControl_Walkthrough.nb</code></a> &mdash; End-to-end optimal-control walkthrough</li>
  <li><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/Examples/OptimalControl/Notebooks/00_OptimalControl_Tutorials_Index.nb"><code>00_OptimalControl_Tutorials_Index.nb</code></a> &mdash; Tutorial series index</li>
  <li><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/Examples/OptimalControl/Notebooks/01_GRAPE_and_SingleSpin_Rotations.nb"><code>01_GRAPE_and_SingleSpin_Rotations.nb</code></a> &mdash; GRAPE and single-spin rotations</li>
  <li><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/Examples/OptimalControl/Notebooks/02_Analytical_Gradients_and_Propagation.nb"><code>02_Analytical_Gradients_and_Propagation.nb</code></a> &mdash; Analytical gradients and propagation</li>
  <li><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/Examples/OptimalControl/Notebooks/03_Polarisation_and_Coherence_Transfer.nb"><code>03_Polarisation_and_Coherence_Transfer.nb</code></a> &mdash; Polarisation and coherence transfer</li>
  <li><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/Examples/OptimalControl/Notebooks/04_Cooperative_and_Advanced_Targets.nb"><code>04_Cooperative_and_Advanced_Targets.nb</code></a> &mdash; Cooperative and advanced targets</li>
  <li><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/Examples/OptimalControl/Notebooks/05_Robust_Ensemble_Optimisation.nb"><code>05_Robust_Ensemble_Optimisation.nb</code></a> &mdash; Robust ensemble optimisation</li>
  <li><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/Examples/OptimalControl/Notebooks/06_Optimisers_Constraints_and_Regularisation.nb"><code>06_Optimisers_Constraints_and_Regularisation.nb</code></a> &mdash; Optimisers, constraints, and regularisation</li>
  <li><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/Examples/OptimalControl/Notebooks/07_Waveform_Parameterisations.nb"><code>07_Waveform_Parameterisations.nb</code></a> &mdash; Waveform parameterisations</li>
  <li><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/Examples/OptimalControl/Notebooks/08_TimeDependent_and_MAS_OptimalControl.nb"><code>08_TimeDependent_and_MAS_OptimalControl.nb</code></a> &mdash; Time-dependent and MAS optimal control</li>
  <li><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/Examples/OptimalControl/Notebooks/09_TimeOptimal_Control.nb"><code>09_TimeOptimal_Control.nb</code></a> &mdash; Time-optimal control</li>
</ul>

### General Optimization (`SpinDynamica/Examples/Optimization/`)

<ul>
  <li><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/Examples/Optimization/Notebooks/00_GeneralOptimization_Tutorials_Index.nb"><code>00_GeneralOptimization_Tutorials_Index.nb</code></a> &mdash; General optimization tutorials index</li>
  <li><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/Examples/Optimization/Notebooks/01_Transformation_Amplitude_Optimization_Basics.nb"><code>01_Transformation_Amplitude_Optimization_Basics.nb</code></a> &mdash; Transformation-amplitude optimization basics</li>
  <li><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/Examples/Optimization/Notebooks/02_Bounded_MultiScale_and_SecondOrder_Optimization.nb"><code>02_Bounded_MultiScale_and_SecondOrder_Optimization.nb</code></a> &mdash; Bounded, multi-scale, and second-order optimization</li>
  <li><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/Examples/Optimization/Notebooks/03_Ensembles_TimeDependent_and_Liouville_Optimization.nb"><code>03_Ensembles_TimeDependent_and_Liouville_Optimization.nb</code></a> &mdash; Ensembles, time-dependent, and Liouville optimization</li>
</ul>

### Signals &amp; Bruker Data (`SpinDynamica/Examples/Signals/`)

<ul>
  <li><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/Examples/Signals/Notebooks/01_Bruker_1D_Data_Import_and_Signal_Workflow.nb"><code>01_Bruker_1D_Data_Import_and_Signal_Workflow.nb</code></a> &mdash; Bruker 1D data import and signal workflow</li>
  <li><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/Examples/Signals/Notebooks/02_2D_NMR_Simulation_Processing_and_MAS.nb"><code>02_2D_NMR_Simulation_Processing_and_MAS.nb</code></a> &mdash; 2D NMR simulation, processing, and MAS</li>
</ul>
