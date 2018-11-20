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
2. Pick a percentage according to your system.
3. Perform alignment of those two structure with the previous defined
percentage


`lovoScan #1&protein #2&protein`

`lovoAlign #1&protein #2&protein 0.2`

### Things you should know

By default, lovo only use Protein Ca atom for alignment.
