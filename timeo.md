#Timeo Plugin
Shape and symmetry analysis: a SHAPE 2.1 wrapper, an octahedral distortion calculator, and a cosmochlore integration (CShM, CSoM, ODis). URL[https://github.com/Yluro/cosmochlore]

#Other programs
Tools that do not use cosmochlore: the SHAPE 2.1 wrapper and the OctaDist reimplementation.

#SHAPE 2.1
SHAPE 2.1 wrapper: Continuous Shape Measures against the built-in reference polyhedra.

#SHAPE Status
Whether **shape.exe**/**shape_2.1.bat** (SHAPE 2.1) was found on **PATH**.

#SHAPE 1
##SHAPE analysis
Runs SHAPE 2.1 on the current selection. One atom selected = centered shape (its neighbours are added); several atoms = non-centered shape. Results print to the console and save under **autoSHAPE/**.

#OctaDist
This plugin's own reimplementation of the OctaDist distortion parameters.

#OctaDist 1
##Octahedral distortion
Runs this plugin's own OctaDist reimplementation. Select the central atom of a 6-coordinate complex; its neighbours are added automatically. Saves a graph under **Oh_distortion/**.

#Options
General plugin-wide options.

#Options 1
##Merge pi-bonded ligands
When checked, pi-bonded ligand fragments are merged into one centroid before measuring (Cirera *et al.*, *Organometallics* 2005, 24, 1556).
##Exclude hydrogens
Exclude any selected or connected hydrogen/deuterium atoms from the measure calculation.

#Options 2
##Disorder parts
Part codes for the measurements. Code **01 02** will analyse parts 0 and 1 together in a single measurement and parts 0 and 2 in another one. 

#Options 3
##SHAPE executable
Path to the SHAPE 2.1 executable (**shape.exe**/**shape_2.1.bat**). Saved between sessions. Leave empty to search **PATH**.

#Options 4
##cosmochlore executable
Path to the cosmochlore executable (1.2.0+). Saved between sessions. Leave empty to search **PATH**.

#Options 5
##Reload ALL
Reloads the plugin's Python code.

##Settings
Opens this plugin's phil parameter editor.

##Open Folder
Opens the plugin's install folder in the file browser.

#Timeo Debug
Low-level diagnostics for the selection/refinement-model machinery.

#Debug 1
##Sel
Prints the current selection string.
##Neighbours
Prints the bonded neighbours of the selected atoms.
##Selection Class
Builds and splits an **AtomSelection**, printing each step.

#Debug 2
##XYZ-SEL
Prints the coordinates of a selected atom.
##ORM
Writes the refinement model's atom list to **orm.txt**.
##SHAPE
Prints whether **shape** was found on **PATH**.

#Cosmochlore Status
Whether a **cosmochlore** executable (1.2.0+) was found. The Run button below needs this to be green. URL[https://github.com/Yluro/cosmochlore]

#Continuous Shape Measures
CShM: compares the selection to built-in and user-defined reference polyhedra. Same selection rule as SHAPE analysis.

#CShM 1
##Run CShM
Runs the measure and prints the results table.

#CShM 2
Options for this CShM run.

#CShM 3
##Results CSV
Also saves a **_cshm_table.csv** file.
##Ideal Polyhedra XYZ
Also saves a **_ideal.xyz** file with the aligned reference shapes, needed to draw them from the runs table below.

#CShM 4
Custom reference shapes beyond the built-in 90.

#CShM 5
One checkbox per **.yaml** file in **user_shapes/**. Only checked ones are used; a mismatched vertex count is reported by cosmochlore itself.

#CShM 6
##Open user defined shapes folder
Opens the **user_shapes** folder, creating it if needed.

#CShM 7
Every CShM run saved for this structure (in **cosmochlore/**), one per centre and disorder component: each **_ideal.xyz** or **_cshm_table.csv** found. Runs that saved neither are not listed.

#CShM 8
**Draw** draws the run's ideal polyhedron over the structure; the appearance settings appear below the table while at least one run is ticked. A run with only a **.csv** has no **_ideal.xyz** to draw: ticking it prints an error asking to rerun CShM with **Ideal polyhedra .xyz** checked. **Shape** browses the run's shapes, best (lowest CShM) first, and **CShM** shows the value for the one picked, which is also the shape drawn. Both choices last this session only.
##Delete
Deletes the run's **_ideal.xyz** and **_cshm_table.csv**, after asking. The input **.xyz** goes too unless a CSoM or ODis output of the same centre still uses it.

#Ideal polyhedra 1
Appearance of the drawn polyhedra, saved between sessions.
##Opacity
Face opacity in percent. **0** draws only edges and vertices.
##Edges
Draws the polyhedron edges.
##Vertices
Marks the ideal vertex positions with small spheres.
##Colour
**shape**: one colour per shape, to compare centres. **centre**: one colour per run, to tell centres and disorder components apart. **element**: the colour Olex2 draws the central atom in.
##Redraw
Draws the polyhedra again, e.g. after **_ideal.xyz** files changed outside Timeo.

#Continuous Symmetry Operation Measures
CSoM: how well the selection matches given point groups. Slower than CShM/ODis, more so with more point groups
or atoms.

#CSoM 1
##Run CSoM
Runs the measure for every listed point group.

#CSoM 2
Options for this CSoM run.

#CSoM 3
##Point groups
Space-separated, case-sensitive Schoenflies symbols (e.g. **Oh D4h D3d**).
At least one is required.

#CSoM 4
##Centering mode
**auto** centers on the first atom (centered selection) or centroid (non-centered); 

**first**/**centroid** always use one or the other; 

**manual** uses the vector on the right.
##Manual vector (x y z)
Used only when the centering mode is set to manual, e.g. **5.1 3.3 -4.9**. In orthogonalised coordinates in Angstroms. 

#CSoM 5
##Summary CSV
Also saves a **_csom_table.csv** summary file.
##Per-operation details CSV
Also saves a per-operation **_details.csv** file, per point group.
##Operated .xyz/.mol2
Also saves the structure with each symmetry operation applied.

#CSoM 6
##Ignore atom labels
Match atoms by position only, ignoring element/label.

#CSoM 7
Axis-search settings (cosmochlore's **--seeds**, **--iterations** and **--tolerance**). It is recommended not to change them.
##Seeds
Number of candidate symmetry-axis directions, spread over a Fibonacci hemisphere, that seed the search. Each one is refined at every in-plane orientation the point group tells apart. Default **20**.
##Iterations
Maximum number of Nelder-Mead iterations spent refining each seed. Default **200**.
##Tolerance
Convergence tolerance of the Nelder-Mead refinement, as a power of ten: the box holds **N** for a tolerance of **1e-N**. Default **6** (1e-6).

#CSoM 8
Every CSoM run saved for this structure (in **cosmochlore/**), one per centre and disorder component: each **_csom_table.csv**, **_details.csv**, **_operated.xyz** or **_merged.mol2** found.

#CSoM 9
**Point group** browses the run's point groups, best (lowest CSoM) first, and **CSoM** shows the value for the one picked. The choice lasts this session only. Only the **_csom_table.csv** holds the CSoM, so a run saved without **Summary .csv** lists its point groups and shows **no .csv**.
##Delete
Deletes all of the run's CSoM files, after asking. The input **.xyz** goes too unless a CShM or ODis output of the same centre still uses it.

#Octahedral Distortion Parameters
ODis: a second, independent octahedral distortion implementation, via cosmochlore. Select the central atom of a 6-coordinate complex.

#ODis 1
##Run ODis
Runs the analysis and prints the results table.

#ODis 2
Options for this odis run.

#ODis 3
##Results CSV
Also saves a **_odis_table.csv** file.
##Full analysis
Also computes CShM (vs. octahedron/trigonal prism) and CSoM (vs. common distortions).
