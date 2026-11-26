---
layout: page
title: About
summary: Authors, scope, citation, and licensing information for SpinDynamica.
section: about
permalink: /about/
---

## Overview

SpinDynamica is a Wolfram Mathematica framework for symbolic and numerical calculations in magnetic resonance spin dynamics. It is designed to provide a flexible, transparent platform in which spin Hamiltonians, Liouvillians, pulse sequences, trajectories, and optimal-control tasks can be expressed close to the underlying quantum-mechanical notation.

As noted in the introductory documentation (`SDv4.0.0-part1-Intro.nb`), exact Hilbert- and Liouville-space calculations are general for any spin system, though memory requirements and runtime grow rapidly beyond roughly six coupled spins. Relaxation is supported in the Redfield regime, while slower fluctuating processes such as chemical exchange can be modeled through ensemble averaging of time-domain trajectories.

## Authors

- **Malcolm H. Levitt**
- **Christian Bengs**

## Acknowledgements

Recent changes and additions in SpinDynamica were initiated by **Mohamed Sabba** and **Christian Bengs** and developed with the assistance of **Astra** and **Antigravity**.

<div class="placeholder-note">
  <p>
    <span class="placeholder-label">Affiliation &amp; Contact (Placeholder):</span>
    Institutional affiliation and direct contact details will be added here.
  </p>
</div>

For bug reports, feature requests, or questions about the public repository, please open an issue on the
<a href="{{ site.repository_url }}/issues">GitHub issue tracker</a>.

<h2 id="how-to-cite">How to Cite</h2>

If you use SpinDynamica in published work, please cite:

{% include citation.html %}

## License

SpinDynamica is open-source software, as stated in the introductory documentation (`SDv4.0.0-part1-Intro.nb`).

<div class="placeholder-note">
  <p>
    <span class="placeholder-label">License Details (Placeholder):</span>
    No formal <code>LICENSE</code> file is currently present in the public repository. The exact open-source license text and SPDX identifier will be documented here once added.
  </p>
</div>
