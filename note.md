# Arguments:

### Mandatory

- **Topology and parameters, coordinates**: -t, --topology
- **Single-frame structure**: -f, --frame

### Optional

- **LaTeX**: -l, --latex; _-> Enables LaTeX typography for plots_. The default is _False_.
- **Parallel calculations**: -p, --parallel; _-> Enables parallel calculation of distances_. The default is _False_.
- **Timer**: -tm, --timer; _-> Enables the time counter for each module_. The default is _False_.
- **csv**: -c, --csv; _-> Enables the creation of csv charts of the plotted data._ The default is _False_.

# Coses a afegir

- [X] AFEGIR QUE O NECESSITA -t O NECESSITA -f
- [X] Comprovar si fa falta posar primer la topologia i després la trajectoria o és indiferent. -> s'han de posar en ordre: primer la top i despres les crds.
- ~~[ ] Afegir barra de progrés a RMSD~~ -> és una funció de MDAnalysis la que calcula el RMSD, no hi ha cap tipus de cicle.
- [X] Corregir títols gràfics dos RMSD
