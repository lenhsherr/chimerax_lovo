# ==========================================================================
# Functions and descriptions for registering using ChimeraX bundle API
# ==========================================================================

import matplotlib.pyplot as plt
from chimerax.atomic import AtomsArg, Atoms
from chimerax.core.commands import CmdDesc, BoolArg, EnumOf  # Command description
from chimerax.core.commands import FloatArg
from chimerax.std_commands.align import *

from .lovo import *


def lovoScan(session, ref: Atoms, moving):
    """

    """

    log = session.logger


    ref_coords = ref.coords
    moving_coords = moving.coords

    print(ref_coords[0])
    print(moving_coords[0])
    print()

    phis = []
    partial_rmsds = []

    for phi in np.arange(0.1, 1.0, 0.05):
        partial_rmsd, align_indices = getPartialRMSD(moving=moving_coords, ref=ref_coords, phi=phi)
        phis.append(phi)
        partial_rmsds.append(partial_rmsd)
        # print("{}".format(partial_rmsd))

    plt.plot(phis, partial_rmsds, linestyle='-', marker='o')
    plt.show()


lovoScan_desc = CmdDesc(
    required=[("ref", AtomsArg), ("moving", AtomsArg)])




def lovoAlign(session, atoms, to_atoms=None, move=None, each=None,
          match_chain_ids=False, match_numbering=False, match_atom_names=False,
          sequence=None, report_matrix=False, phi = None, scan_phi = None):
    """Move atoms to minimize RMSD with to_atoms.
    Returns matched atoms and matched to_atoms, matched atom rmsd, paired atom rmsd, and transform.
    The matched atoms can be fewer than the paired atoms if cutoff distance is specified.
    If "each" is not None then nothing is returned.

    If 'move' is 'structures', superimpose the models by changing the model positions.
    If it is 'atoms', 'residues', 'chains' or 'structure atoms', then atoms promoted extended
    to this level are moved.  If move is False move nothing or True move structures.
    If move is an Atoms collection then move only the specified atoms.

    If 'each' is "structure" then each structure of atoms is separately
    aligned to the to_atoms.  If 'each' is "chain" then each chain is
    aligned separately.  If 'each' is "coordset" then each coordinate set
    of the first set of atoms (which must belong to a single structure)
    is aligned. Default is that all atoms are aligned as one group.

    If 'match_chain_ids' is true then only atoms with matching chain identifiers are paired.
    Unpaired atoms or to_atoms are not used for alignment.

    If 'match_numbering' is true then only atoms with matching residue numbers are paired.
    It is assumed the atoms are in residue number order.  Unpaired atoms or to_atoms
    are not used.

    If 'match_atom_names' is true then only atoms with matching names are paired.
    Unpaired atoms or to_atoms are not used for alignment.

    If 'sequence' names a reference sequence then the previously calculated alignment
    of atoms and to_atoms to that reference sequence is used to pair the atoms.

    If 'report_matrix' is True, report the transformation matrix to
    the Reply Log.

    If 'scan_phi' is True, output the rmsd vs phi graph.

    """

    import numpy as np
    from chimerax.core.geometry import align_points

    if move is None:
        move = 'structures'

    log = session.logger
    if sequence is None:
        patoms, pto_atoms = paired_atoms(atoms, to_atoms, match_chain_ids,
                                         match_numbering, match_atom_names)
        da, dra = len(atoms) - len(patoms), len(to_atoms) - len(pto_atoms)
        if da > 0 or dra > 0:
            log.info('Pairing dropped %d atoms and %d reference atoms' % (da, dra))
    else:
        patoms, pto_atoms = sequence_alignment_pairing(atoms, to_atoms, sequence)

    npa, npta = len(patoms), len(pto_atoms)
    if npa != npta:
        raise UserError('Must align equal numbers of atoms, got %d and %d' % (npa, npta))
    elif npa == 0:
        raise UserError('No atoms paired for alignment')

    xyz_to = pto_atoms.scene_coords

    xyz_from = patoms.scene_coords




    if scan_phi:
        phis = []
        partial_rmsds = []


        for phi in np.arange(0.1, 1.0, 0.05):
            partial_rmsd, align_indices = getPartialRMSD(moving=xyz_from, ref=xyz_to, phi=phi)
            phis.append(phi)
            partial_rmsds.append(partial_rmsd)

        plt.plot(phis, partial_rmsds, linestyle='-', marker='o')
        plt.xlabel("Fraction of atoms to use")
        plt.ylabel("RMSD of used atoms")
        plt.show()
        return
    else:
        if phi is None:
            phi = 1.0
        partial_rmsd, stable_indices = (getPartialRMSD(xyz_from, xyz_to, phi))
        tf, rmsd = align_points(xyz_from[stable_indices], xyz_to[stable_indices])

        full_rmsd = rmsd
        matched_patoms, matched_pto_atoms = patoms, pto_atoms
        msg = 'RMSD between %d atom pairs is %.3f angstroms' % (len(stable_indices), rmsd)

        if report_matrix:
            log.info(matrix_text(tf, atoms, to_atoms))

        log.status(msg)
        log.info(msg)

        move_atoms(atoms, to_atoms, tf, move)

        return matched_patoms, matched_pto_atoms, rmsd, full_rmsd, tf





lovoAlign_desc = CmdDesc(required = [('atoms', AtomsArg)],
               keyword = [('to_atoms', AtomsArg),
                          ('move', EnumOf(('atoms', 'residues', 'chains', 'structures',
                                           'structure atoms', 'nothing'))),
                          ('each', EnumOf(('chain', 'structure', 'coordset'))),
                          ('match_chain_ids', BoolArg),
                          ('match_numbering', BoolArg),
                          ('match_atom_names', BoolArg),
                          ('report_matrix', BoolArg),
                          ('phi', FloatArg),
                          ('scan_phi', BoolArg)],
               required_arguments = ['to_atoms'],
               synopsis = 'Align one set of atoms to another')