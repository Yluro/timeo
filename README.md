# Timeo

Timeo is an [Olex2](https://www.olexsys.org/olex2/)$^1$ plugin that integrates several symmetry and shape analysis tools into a single graphical interface within Olex2.

_Timeo was called Timaeus before version 1.2.0. If you have Timaeus installed, see [Upgrading from Timaeus](#upgrading-from-timaeus)._
### Features

All analysis methods in Timeo will read the atomic coordinates directly from the `OlexRefinementModel`; there is no need for manual writing of any i/o files needed for external programs.

**Olex2 GUI panel:**

_Almost_ all functions of Timeo are accessible through the GUI panel.

**[Cosmochlore](https://github.com/Yluro/cosmochlore) integration:**

Cosmochlore was developed in parallel to Timeo to organically integrate into it. All options of cosmochlore can be accessed through the Continuous Shape Measures, Continuous Symmetry Operation Measures and Octahedral Distortion Parameters tabs.

**SHAPE 2.1 wrapper:**

SHAPE 2.1 i/o files are supported by Timeo can run and parse SHAPE 2.1 files

**Octahedral distortion parameters:**

The plugin comes with its own simpler implementation of the OctaDist algorithm (separate from cosmochlore's one) based on topological arguments of the octahedral geometry.

**$\eta^n$ pi-bonded ligands:**

Timeo can detect pi-bonded ligands and the option is given to average the coordinates of all donor atoms of the ligand to its centroid 
as per the Cirera _et al_$^2$ paper. 

**Disorder handling:**

Automatic disorder handling. If more than two parts are found in the selected structure, calculations will be run on each part separately.

By default, all measurements run on part 0 + part N. So if the selected atom is in part 0 and has some connected atoms in part 1 and 2, the measures will run on parts 0 + 1 and 0 + 2. 

Part codes can be specified: `0` will include only atoms in part 0 in the measurement. `01` will include atoms part 0 and 1. `12 34` will separate the atoms in part 1 and 2 into one measurement and the atoms in part 3 and 4 into another.


## Requirements
- Olex2 1.5.
- [Cosmochlore](https://github.com/Yluro/cosmochlore) 1.2.0 or newer. Cosmochlore is my own implementation of the measures algorithms developed specially to work with Timeo.
- _(Optional)_ SHAPE 2.1
  - Download SHAPE 2.1 from the [Electronic Structure Group's webpage](https://www.ee.ub.edu/downloads/)

The plugin was developed/tested using a Windows 10/11 machine. The plugin should be system agnostic but please report any bugs found in any other operating systems. 

_**Note:** It is known that SHAPE 2.1 gives trouble in Mac machines with operating systems newer than 2022 and some Linux systems. Unfortunately, I can't do anything about that since the ESG hasn't published a precompiled SHAPE version since 2013. In the future I might introduce support for cosymlib by ESG which is the updated version of their shape and symmetry measures program._

## Installation
1. [Download the files](https://github.com/Yluro/timeo/releases) from the lastest GitHub release. 
2. Go to the base directory of your Olex2 installation. You can open the directory by typing `shell BaseDir()` in the Olex2 console.
3. Create a file called `plugins.xld` and write the following contents inside `plugins.xld`:
```xml
<Plugin
 <Timeo>
>
``` 
4. Go to the `<BaseDir>\util\pyUtil\PluginLib\` folder and create a folder called `plugin-Timeo`.
5. Extract the downloaded zip and drop the contents of the `timeo-X.X.X` folder, into `plugin-Timeo` so that `Timeo.py` sits directly at `<BaseDir>\util\pyUtil\PluginLib\plugin-Timeo\Timeo.py`.
6. On restarting Olex2, a Timeo panel should appear under the Tools tab.

**Troubleshooting:** 
 - _If the panel does not appear, the Olex2 console prints `Failed to load plugin 'X'`. `No module named 'X'` means the name in `plugins.xld` does not match `<BaseDir>\util\pyUtil\PluginLib\plugin-X\X.py` letter for letter (including case), or `Timeo.py` is nested one folder too deep._
 - _Make sure `plugins.xld` is not saved as `plugins.xld.txt`._
 - _If `plugins.xld` already exists, do not overwrite it. Add the `<Timeo>` line between `<Plugin` and `>`._
 - _`ModuleNotFoundError: No module named 'timeo'` occurs if the name of the plugin is written in lowercase._

### Upgrading from Timaeus
1. Delete the `<BaseDir>\util\pyUtil\PluginLib\plugin-Timaeus` folder and install Timeo as described above.
2. In `plugins.xld`, replace the `<Timaeus>` line with `<Timeo>`.
3. Console commands are now `spy.Timeo.*` instead of `spy.Timaeus.*`.

_**Note:** I push development changes constantly to the master branch. You could git clone this repository to automatically keep the plugin updated. Usually, if changes are pushed it means that the plugin is in a usable state. But it does not guarantee that things won't break._

## Cosmochlore
[cosmochlore](https://github.com/Yluro/cosmochlore) is a separate Rust program that calculates shape and symmetry measures. Timeo wraps its three subcommands; each writes its output next to the structure and prints the results table to the console. Every option below has a phil parameter, editable in `Options` > `Settings` or from the panel itself. Cosmochlore was intentionally developed to function as a background process to Timeo. Most options of the program can be accessed through the GUI.

| Function               | Description                                                                                                                              |
|------------------------|------------------------------------------------------------------------------------------------------------------------------------------|
| `spy.Timeo.autoCSHM()` | Continuous Shape Measures against the built-in and user-defined reference polyhedra.                                                     |
| `spy.Timeo.autoCSOM()` | Continuous Symmetry Operation Measures against a list of point groups.                                                                   |
| `spy.Timeo.autoODIS()` | Octahedral distortion parameters, independent of the OctaDist reimplementation below. Select the central atom of a 6-coordinate complex. |

### Usage

1. Open a structure in Olex2.
2. Select the atoms you want to include in the measurement.
 - If one atom is selected, the neighbouring atoms will be taken into account to form a centered shape (e.g. a coordination structure).
 - If multiple atoms are selected, they will be interpreted as a non-centered shape (e.g. a borane cluster).
4. Run any of the cosmochlore functions from the Olex2 console or from its own tab in the Tools/Timeo panel.
5. Results are printed to the console and saved in `<FilePath>/cosmochlore/`.



### CShM runs table
The bottom of the Continuous Shape Measures tab lists every CShM run saved for the structure: one row per measured centre and disorder component, found as `_ideal.xyz` or `_cshm_table.csv` files in `<FilePath>/cosmochlore/`. Runs that saved neither are not listed. Per run:
- `Shape`: the run's shapes, best (lowest CShM) first. `CShM` shows the value for the one picked, which is also the shape drawn.
- `Delete`: deletes the run's `_ideal.xyz` and `_cshm_table.csv` after asking, plus its input `.xyz` unless a CSoM or ODis output of the same centre still uses it.

### Ideal polyhedra in the model
Tick `Show Ideal Polyhedra` in the CShM options to draw the ideal polyhedra from CShM over the structure. Every run in the table with an `_ideal.xyz` is drawn at once and redrawn after each new run. While the box is ticked the drawing settings appear, and the runs table gets a `Draw` column to leave runs out:
- `Opacity %`: face transparency. Set it to 0 to draw only edges and vertices.
- `Edges` and `Vertex spheres`: draw the polyhedron edges and mark the ideal vertex positions.
- `Colour by`: `shape` gives every shape its own colour, so centres are easy to compare; `centre` gives every run its own colour; `element` uses the colour Olex2 draws the central atom in.

Opacity, edges, vertex spheres and colour are saved between sessions; the runs table's shape and Draw choices last for the session.

| Function                                | Description                                                                         |
|-----------------------------------------|-------------------------------------------------------------------------------------|
| `spy.Timeo.draw_ideal()`                | Draws the chosen shape of every ticked run with the current settings.               |
| `spy.Timeo.clear_ideal()`               | Hides every drawn ideal polyhedron.                                                 |

The overlay uses Olex2's built-in `AddObject` user objects, which disappear when another structure is loaded.

### User-defined shapes
Reference shapes beyond the built-in 90 can be added as `.yaml` files in the `user_shapes/` folder of the plugin (`Open user defined shapes folder` in the Continuous Shape Measures tab). A mismatched vertex count is reported as an error by `cosmochlore` itself.

## AutoSHAPE
[SHAPE 2.1](https://www.ee.ub.edu/continuous-shape-and-symmetry-measures/)$^2$ is a software published by ESG used to calculate Continuous Shape Measures (CShM's). `autoSHAPE` is a collection of personal Python scripts I developed and used to run and parse SHAPE i/o files. Timeo contains an implementation of `autoSHAPE` to:
- Generate the necessary `.dat` input files for SHAPE automatically.
- Run SHAPE and parse the resulting `.tab` output.
- Output a summary table from the `.out` and `.tab` files.
SM's autoSHAPE does not overwrite previous runs as it stores each run in a dedicated folder: `<FilePath>\autoSHAPE\<FileName>_<part>_<atoms>\<run>`.


### Usage
 1. `autoSHAPE` uses the same selection criteria as `cosmochlore`. Select a single atom for a centered measurement or multiple for a non-centered one.
 2. Run `spy.Timeo.autoSHAPE()` from the Olex2 command line or from the Tools/Timeo/Other programs panel.

## Octahedral Distortion Parameters.
Timeo includes a reimplementation of the [OctaDist](https://octadist.github.io/)$^3$ algorithm. Unlike the original implementation, this version identifies opposite faces and vertices of an octahedron using topological criteria. It relies on constructing a convex hull that is topologically equivalent to an octahedron. As a result, the algorithm may fail when more than three vertices are coplanar, causing the convex hull to degenerate into a different polyhedral shape.

Several features of the original OctaDist program have been omitted to simplify integration with Olex2 and to remove redundant functionality.

The Octahedral Distortion Parameters module computes all distortion parameters reported by OctaDist, together with an additional $\tau$ and $\mu$ parameters that quantify the deviation of trans angles from the ideal 180 $^\circ$ and the deviation of the metal centre to the centroid of the octahedron, respectively. The calculated values are generally consistent with those produced by OctaDist, although discrepancies may occur for highly distorted octahedra or for structures approaching an ideal trigonal prism.


- Bond length distortion:
```math
\zeta = \sum_{i=1}^6 |d_i - d_{mean}|
```
where $d_i$ are the M-X bond distances and $d_{mean}$ is the mean M-X bond distance.
- Octahedral tilting parameter:
```math
\Delta = \frac{1}{6}\sum_{i=1}^6 \left(\frac{d_i - d_{mean}}{d_{mean}}\right)^2
```
where $d_i$ are the M-X bond distances and $d_{mean}$ is the mean M-X bond distance.
- Cis angle distortion:
```math
\Sigma = \sum_{i=1}^{12} |90^\circ - \phi_i|
```
where $\phi_i$ are the cis angles.
- Octahedral twisting distortion:
```math
\Theta = \sum_{i=1}^{24}  |60^\circ - \theta_i|
```
where $\theta_i$ are twisting angles between vectors of two opposite faces.

#### Extra parameters:
- Trans angle distortion:
```math
\tau = \sum_{i=1}^3 |180^\circ - \psi_i|
```
where $\psi_i$ are the trans angles.
- Deviation of the metal from the centroid:
```math
\vec{c} = \frac{1}{6}\sum^6_{i=1} \vec{r}_i
```
```math
\mu = |\vec{r}_M - \vec{c}|
```
where $\vec{r}_M$ is the position of the metal and $\vec{r}_i$ are the position of the donor atoms.


### Usage
1. Open a structure in Olex2.
2. Select the central atom of a 6-coordinate complex.
3. Run `spy.Timeo.autoOCTADIST()` from the Olex2 console or from the Tools/Timeo/Other programs panel.
4. Results are printed in the console. A graph will be saved in `<FilePath>/OH_distortion` showing the extracted octahedron. 


## Known limitations/upcoming features.

| Status | Features                                                                                                                 |
|--------|--------------------------------------------------------------------------------------------------------------------------|
| ✔️     | Suppport for disordered structures.                                                                                      |
| ✔️     | Centorid search for pi-bonded ligands. *Doesn't work with ligands outside ASU.                                           |
| ✔️     | Custom reference shapes. *Via cosmochlore, see `user_shapes/`.                                                           | 
| ✔️     | Smarter program logic (automatic coordination site detection, multiple selections, etc.).                                | 
| ✔️     | Support for other measurement programs (cosmochlore: CShM, CSoM, ODis).                                                  |
| ✔️     | Reimplementation of octahedral distortion parameters (Zeta, Sigma, Theta) — relevant for spin-crossover (SCO) complexes. |
| ✔️     | Non centered shapes.                                                                                                     |
| ✔️     | HTML UI.                                                                                                                 |

## License
Copyright (C) 2026 José Serrano Guarinos.

Timeo is free software, released under the [GNU General Public License v3.0](LICENSE.MD) or (at your option) any later version.

## Citations

1. _J. Appl. Cryst._ (2009). 42, 339–341. DOI: https://doi.org/10.1107/S0021889808042726
2. _Organometallics_ (2005), 24, **7**, 1556–1562. DOI: https://doi.org/10.1021/om049150z
3. _Dalton Trans._ (2021) 50, **3**, 1086–1096. DOI: https://doi.org/10.1039/d0dt03988h
