---
layout: page
title: Installation
summary: How to install SpinDynamica, configure Mathematica's $Path, and verify the installation.
section: installation
permalink: /installation/
---
{% assign rel = site.data.release %}

## Requirements

- **Wolfram Mathematica:** version {{ rel.min_mathematica }} or later (documented and tested under Mathematica {{ rel.tested_mathematica }}).
- **SpinDynamica distribution archive** from the [Download]({{ '/download/' | relative_url }}) page.

## 1. Unpack the Package

Download the release archive and extract it to a permanent location on your computer. Keep the `SpinDynamica` directory (containing `SpinDynamica.m` and its module subdirectories) intact and unmodified.

## 2. Inform Mathematica Where SpinDynamica Is Located

Mathematica searches the directories listed in `$Path` when loading packages. You can add your `SpinDynamica` directory either per notebook or persistently for all sessions.

### Option A: Per-Notebook Setup (Recommended for Reproducibility)

Place one of the following lines at the top of your Mathematica notebook, replacing the path with the directory that contains `SpinDynamica.m`:

```mathematica
(* macOS / Linux *)
PrependTo[$Path, "/path/to/SpinDynamica"];

(* Windows *)
PrependTo[$Path, "C:\\path\\to\\SpinDynamica"];
```

### Option B: Persistent Setup via `init.m`

To configure `$Path` automatically whenever a Mathematica kernel starts, add the same `PrependTo[$Path, ...]` statement to `init.m` in Mathematica's `Autoload` or user `Kernel` directory (as described in Part 1 of the [documentation notebooks]({{ '/documentation/' | relative_url }})). Restart the kernel and evaluate `$Path` to verify that the directory appears in the list.

## 3. Load SpinDynamica

Execute the following command in a Mathematica notebook:

```mathematica
Needs["SpinDynamica`"]
```

When loaded, SpinDynamica displays a startup banner showing the active version, Mathematica version, and operating system.

<div class="info-box">
  <p>
    <strong>Note on built-in symbols:</strong> SpinDynamica extends several built-in Mathematica symbols
    (including <code>Dot</code> and <code>Times</code>) to support spin-operator algebra. In the same kernel
    session, this can slow down unrelated Mathematica calculations that rely heavily on those functions.
  </p>
</div>

## 4. Verify the Installation

Evaluate `$SpinDynamicaVersion` and `SpinDynamicaSymbols[]` to confirm that the package loaded properly and inspect the active routines:

```mathematica
$SpinDynamicaVersion
SpinDynamicaSymbols[]
```

- `$SpinDynamicaVersion` returns the version string of the loaded SpinDynamica package.
- `SpinDynamicaSymbols[]` displays a clickable overview of top-level SpinDynamica symbols grouped by module.
