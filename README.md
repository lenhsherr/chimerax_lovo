# ChimeraX pyLovo

Protein sub structure alignment

## Installation
First install ChimeraX [here](https://www.rbvi.ucsf.edu/chimerax/)

Then download the alphaspace and store it in your preferred location.

Open ChimeraX and in the command line, type:

```devel install pylovo_installation_folder```

To reinstall after any change in code, simply repeat the command.

## Usage


### For protein with different topology

`lovo match_atoms to ref_atoms n_cluster 3 plot true`

For example:

`lovo #2 to #1` would align all atoms of model 2 to model 1
`lovo #2@Ca to #1@Ca` would align just Ca carbon atoms of model 2 to model 1


### For trajectory
When two or more structures are from a trajectory, meaning they have exactly same topology, you can align them with this
command:

`lovo traj match_atoms to ref_atom n_cluster 3 plot true`

Use *to* option to specify which snapshot to align to, if not provided it defaults to first snapshots

Use n_clusters to specify how many clusters to produce, default is 3, with core, intermediate and flexible regions

You can also produce plots with `plot true` command


This program is based on the MDlovo paper:

 [L. Martínez, Automatic identification of mobile and rigid substructures
 in molecular dynamics simulations and fractional structural fluctuation
 analysis. PLoS One 10(3): e0119264, 2015.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0119264)

