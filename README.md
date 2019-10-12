# pyLOXr README

pyLOXr is a set of Python scripts based on [MDAnalysis](https://www.mdanalysis.org) Python package, which are useful for analysing Molecular Dynamics trajectories calculated using AMBER and to prepare the required files for carrying out QM/MM calculations using [ChemShell](https://www.chemshell.org/). It is presented as an option-based program (as-simple-as-possible command-line UI), but it is also possible to use it as a set of modules in order to integrate it in any other program.

## Modules

pyLOXr is divided into 6 different modules, some of them depending on other ones.

1. **Trajectory summariser**: generates a short *txt* file where the basic characteristics of the trajectroy are specified.
2. **Plots and histograms of distances**: generates the plots and histograms of the distances between groups of two atoms.
3. **RMSD**: generates a plot of the RMSD of the backbone of the protein or of the substrate.
4. **Frame selector**: selects the frames of the trajectory which satisfy a given criteria(on). This(These) criterion can be either a *cut-off* distance (specified by user) between two specified atoms or the difference of distance between two different bonds.
5. **QM/MM models**: it generates a cropped model where only a pseudospherical drop of solvent molecules are kept. As a results, the module gives a *pdb* file which can be used on ChemShell. Moreover, it lets crop also the topology and parameters and coordinates files using the ParmEd package.

## Dependencies

pyLOXr depends on the following Python packages:

- Python >= 3.6
- [MDAnalysis](https://www.mdanalysis.org)
- [Matplotlib](https://matplotlib.org/)
- [Numpy](http://www.numpy.org/)
- ParmEd ([Documentation](http://parmed.github.io/ParmEd/html/index.html), [Repository](https://github.com/ParmEd/ParmEd))