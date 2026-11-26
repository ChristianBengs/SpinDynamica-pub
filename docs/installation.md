---
layout: page
title: Installation
summary: How to install SpinDynamica via Wolfram Paclet or manual $Path setup, and verify the installation.
section: installation
permalink: /installation/
---
{% assign rel = site.data.release %}

## Requirements

- **Wolfram Mathematica:** version {{ rel.min_mathematica }} or later (tested under Mathematica {{ rel.tested_mathematica }}).
- **SpinDynamica distribution** (`.paclet` file or source archive) from the [Download]({{ '/download/' | relative_url }}) page.

## Method 1: Install as a Wolfram Paclet (Recommended)

Download the `.paclet` file (for example, `SpinDynamica-{{ rel.version }}.paclet`) from the [Download]({{ '/download/' | relative_url }}) page and install it from a Mathematica notebook:

```mathematica
PacletInstall["/path/to/SpinDynamica-{{ rel.version }}.paclet"]
```

On Windows, use escaped backslashes or forward slashes in the file path (for example, `"C:/path/to/SpinDynamica-{{ rel.version }}.paclet"`).

Once installed as a paclet, SpinDynamica can be loaded directly in any new Mathematica session without modifying `$Path`:

```mathematica
Needs["SpinDynamica`"]
```

To remove or replace an installed paclet later, evaluate `PacletUninstall["SpinDynamica"]`.

## Method 2: Manual Installation via `$Path`

If you prefer not to install a paclet (or are working directly from an unpacked source archive), extract the archive to a permanent location while keeping the inner `SpinDynamica` directory (containing `SpinDynamica.m` and its module subdirectories) intact.

### Per-Notebook Setup

Place one of the following lines at the top of your Mathematica notebook, replacing the path with the directory that contains `SpinDynamica.m`:

```mathematica
(* macOS / Linux *)
PrependTo[$Path, "/path/to/SpinDynamica"];

(* Windows *)
PrependTo[$Path, "C:\\path\\to\\SpinDynamica"];

Needs["SpinDynamica`"]
```

### Persistent Setup via `init.m`

To configure `$Path` automatically whenever a Mathematica kernel starts, add the same `PrependTo[$Path, ...]` statement to `init.m` in Mathematica's `Autoload` or user `Kernel` directory (as described in Part 1 of the [documentation notebooks]({{ '/documentation/' | relative_url }})). Restart the kernel and evaluate `$Path` to verify that the directory appears in the list.

## Notes on Loading

When loaded, SpinDynamica displays a startup banner showing the active version, Mathematica version, and operating system.

<div class="info-box">
  <p>
    <strong>Note on built-in symbols:</strong> SpinDynamica extends several built-in Mathematica symbols
    (including <code>Dot</code> and <code>Times</code>) to support spin-operator algebra. In the same kernel
    session, this can slow down unrelated Mathematica calculations that rely heavily on those functions.
  </p>
</div>

## Verify the Installation

Evaluate `$SpinDynamicaVersion` and `SpinDynamicaSymbols[]` to confirm that the package loaded properly and inspect the active routines:

```mathematica
$SpinDynamicaVersion
SpinDynamicaSymbols[]
```

- `$SpinDynamicaVersion` returns the version string of the loaded SpinDynamica package.
- `SpinDynamicaSymbols[]` displays a clickable overview of top-level SpinDynamica symbols grouped by module.
