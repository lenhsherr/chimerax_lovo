# ChimeraX pyLovo

Protein sub structure alignment

## Installation
First install ChimeraX [here](https://www.rbvi.ucsf.edu/chimerax/)

Then download the alphaspace and store it in your prefered location.

Open ChimeraX and in the command line, type:

```devel install pylovo_installation_folder```

To reinstall after any alteration, simply repeat the command.

## Usage

### Snapshot to snapshot alignment

1. Identify stable substructure
    Use lovoScan command to generate plot of RMSD vs percentage
    of total structure to be aligned with.
    `lovo #1@Ca to_atom #2@Ca scan_phi 1`

2. Pick a percentage according to your system.
3. Perform alignment of those two structure with the previous defined
percentage

    `lovo #1@Ca to_atom #2@Ca phi 0.8`


This program is based on the MDlovo paper:

 [L. Martínez, Automatic identification of mobile and rigid substructures
 in molecular dynamics simulations and fractional structural fluctuation
 analysis. PLoS One 10(3): e0119264, 2015.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0119264)


Automatically selecting cut off value for system

