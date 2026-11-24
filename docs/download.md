---
layout: page
title: Download
summary: Stable versioned GitHub Releases and access to the public source repository.
section: download
permalink: /download/
---
{% assign rel = site.data.release %}

## Latest Stable Release

{% if rel.published and rel.tag != "" %}
<div class="info-box">
  <p>
    <strong>Version:</strong> {{ rel.version }}
    {% if rel.release_date != "" %}&middot; <strong>Released:</strong> {{ rel.release_date }}{% endif %}<br>
    <strong>Git Tag:</strong> <code>{{ rel.tag }}</code>
  </p>
  {% if rel.notes_summary != "" %}
  <p>{{ rel.notes_summary }}</p>
  {% endif %}
  <p>
    {% include download-button.html %}
    <a class="btn btn-secondary" href="{{ site.repository_url }}/releases/tag/{{ rel.tag }}">Release Notes on GitHub</a>
  </p>
</div>
{% else %}
<div class="placeholder-note">
  <p>
    <span class="placeholder-label">Release pending:</span>
    The first versioned stable GitHub Release has not yet been published. Once tagged and published
    on GitHub, direct archive download links and release metadata will appear here automatically.
  </p>
  <p>
    {% include download-button.html %}
  </p>
</div>
{% endif %}

## Previous Releases

<p>
  All versioned distributions are archived on the
  <a href="{{ site.repository_url }}/releases">GitHub Releases page</a>.
</p>

{% if rel.previous_releases.size > 0 %}
<div class="table-wrap">
  <table>
    <thead>
      <tr>
        <th>Version</th>
        <th>Release Date</th>
        <th>Notes</th>
        <th>Archive</th>
      </tr>
    </thead>
    <tbody>
      {% for prev in rel.previous_releases %}
      <tr>
        <td><code>{{ prev.version }}</code></td>
        <td>{{ prev.release_date }}</td>
        <td>{{ prev.notes }}</td>
        <td><a href="{{ site.repository_url }}/releases/tag/{{ prev.tag }}">View release</a></td>
      </tr>
      {% endfor %}
    </tbody>
  </table>
</div>
{% else %}
<p>
  No previous versioned releases are listed yet. Historical releases will be tabulated here once published.
</p>
{% endif %}

## Source Code Repository

<p>
  The public source code is hosted on GitHub at
  <a href="{{ site.repository_url }}">{{ site.repository_url }}</a>.
  For reproducible calculations, we recommend using a versioned GitHub Release rather than the active development branch.
</p>

## System Requirements

<ul>
  <li><strong>Wolfram Mathematica:</strong> version {{ rel.min_mathematica }} or higher (tested under Mathematica {{ rel.tested_mathematica }}).</li>
  <li><strong>Platform:</strong> cross-platform (Windows, macOS, and Linux, wherever Mathematica is supported).</li>
</ul>

<p>
  After downloading and unpacking the archive, follow the <a href="{{ '/installation/' | relative_url }}">Installation</a> guide to configure Mathematica and load the package.
</p>
