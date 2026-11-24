---
layout: page
title: Documentation
summary: Web guides, Mathematica documentation notebooks, and interactive help commands.
section: documentation
permalink: /documentation/
---

## Web Guides

{% assign docs = site.documentation | sort: 'order' %}
{% if docs.size > 0 %}
<ul class="card-list">
  {% for doc in docs %}
  <li class="card-list-item">
    <a class="card-list-title" href="{{ doc.url | relative_url }}">{{ doc.title }}</a>
    {% if doc.summary %}
    <p class="card-list-summary">{{ doc.summary }}</p>
    {% endif %}
  </li>
  {% endfor %}
</ul>
{% endif %}

<div class="placeholder-note">
  <p>
    <span class="placeholder-label">Documentation roadmap:</span>
    Additional web-native reference pages covering individual SpinDynamica modules and symbols will be added here in future updates.
  </p>
</div>

## Mathematica Documentation Notebooks

The repository includes eleven Mathematica documentation notebooks under <code>SpinDynamica/documentation/</code> introducing core concepts, conventions, and top-level routines:

<div class="table-wrap">
  <table>
    <thead>
      <tr>
        <th>Part</th>
        <th>Topic</th>
        <th>Notebook</th>
      </tr>
    </thead>
    <tbody>
      <tr>
        <td>Part 1</td>
        <td>Introduction, setup, and conventions</td>
        <td><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/documentation/SD.3.10-part1-Intro.nb"><code>SD.3.10-part1-Intro.nb</code></a></td>
      </tr>
      <tr>
        <td>Part 2</td>
        <td>Numerical and symbolic simulations</td>
        <td><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/documentation/SD.3.10-part2-Simulations.nb"><code>SD.3.10-part2-Simulations.nb</code></a></td>
      </tr>
      <tr>
        <td>Part 3</td>
        <td>Spin operators and Hilbert space</td>
        <td><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/documentation/SD.3.10-part3-Operators.nb"><code>SD.3.10-part3-Operators.nb</code></a></td>
      </tr>
      <tr>
        <td>Part 4</td>
        <td>Superoperators and Liouville space</td>
        <td><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/documentation/SD.3.10-part4-Superoperators.nb"><code>SD.3.10-part4-Superoperators.nb</code></a></td>
      </tr>
      <tr>
        <td>Part 5</td>
        <td>Spin interactions and Hamiltonians</td>
        <td><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/documentation/SD.3.10-part5-Interactions.nb"><code>SD.3.10-part5-Interactions.nb</code></a></td>
      </tr>
      <tr>
        <td>Part 6</td>
        <td>Relaxation superoperators (Redfield regime)</td>
        <td><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/documentation/SD.3.10-part6-Relaxation.nb"><code>SD.3.10-part6-Relaxation.nb</code></a></td>
      </tr>
      <tr>
        <td>Part 7</td>
        <td>Signals and Fourier transformation</td>
        <td><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/documentation/SD.3.10-part7-Signal+FT.nb"><code>SD.3.10-part7-Signal+FT.nb</code></a></td>
      </tr>
      <tr>
        <td>Part 8</td>
        <td>3D geometry and graphics</td>
        <td><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/documentation/SD.3.10-part8-3D.nb"><code>SD.3.10-part8-3D.nb</code></a></td>
      </tr>
      <tr>
        <td>Part 9</td>
        <td>Euler angles and Wigner matrices</td>
        <td><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/documentation/SD.3.10-part9-Euler+Wigner.nb"><code>SD.3.10-part9-Euler+Wigner.nb</code></a></td>
      </tr>
      <tr>
        <td>Part 10</td>
        <td>Polarization and thermal equilibrium</td>
        <td><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/documentation/SD.3.10-part10-Polarization.nb"><code>SD.3.10-part10-Polarization.nb</code></a></td>
      </tr>
      <tr>
        <td>Part 12</td>
        <td>Groups, algebras, and spin permutations</td>
        <td><a href="{{ site.repository_url }}/blob/{{ site.repository_branch }}/SpinDynamica/documentation/SD.3.10-part12-Groups+Algebras.nb"><code>SD.3.10-part12-Groups+Algebras.nb</code></a></td>
      </tr>
    </tbody>
  </table>
</div>

## Interactive Help in Mathematica

Once SpinDynamica is loaded, built-in Wolfram Language help commands display module listings, usage messages, and default options directly in your notebook:

```mathematica
SpinDynamicaSymbols[]
?Trajectory
Options[Trajectory]
```
