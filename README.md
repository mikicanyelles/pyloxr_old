# pyLOXr README

pyLOXr is a set of Python scripts based on [MDAnalysis](https://www.mdanalysis.org), which are useful for analysing Molecular Dynamics trajectories calculated using any of the most widespread MD simulation programs ([list of supported coordinates files](https://www.mdanalysis.org/docs/documentation_pages/coordinates/init.html#supported-coordinate-formats)) and to prepare the required files for carrying out QM/MM calculations using [ChemShell](https://www.chemshell.org/). It is presented as an option-based program (a very simple command-line GUI), but it is also possible to use it as a set of modules to integrate it in any other program.

## Modules

pyLOXr is divided into 6 different modules, some of them depending on other ones.

1. **Trajectory summariser**: generates a short *txt* file where the basic characteristics of the trajectroy are specified.
2. **Plots and histograms of distances**: generates the plots and histograms of the distances between an atom of the protein (commonly an O) and 1, 2 or 3 C(s) and the corresponding H(s).
3. **RMSD**: (*requires `cpptraj` (embedded in [AmberTools](http://www.ambermd.org/AmberTools.php)) to work*) generates a plot of the RMSD of the backbone of the protein or of the substrate.
4. **Frame selector**: selects the frames of the dynamics which satisfy a given *cut-off* distance (specified by user) between two specified atoms.
5. **QM/MM models**: from the frames saved by the *"Frame selector"*, it generates a cropped model where only a drop of solvent molecules are kept. As a results, the module gives a *pdb* file which can be used on ChemShell.
6. ***set act* creator**: from the *QM/MM models* module, it generates a list of the atoms which can be moved during QM/MM calculations on ChemShell.

## Dependencies

pyLOXr depends on the following Python packages:

- Python >= 3.6
- [MDAnalysis](https://www.mdanalysis.org)
- [Matplotlib](https://matplotlib.org/)
- [Numpy](http://www.numpy.org/)
- [pandas](https://pandas.pydata.org/)