# ==========================================================================
# Functions and descriptions for registering using ChimeraX bundle API
# ==========================================================================

from chimerax.atomic import AtomsArg
from chimerax.core.colors import Color
from chimerax.core.commands import CmdDesc, BoolArg, EnumOf, IntArg  # Command description
from chimerax.core.commands import FloatArg
from chimerax.core.geometry import align_points
from chimerax.core.objects import Objects
from chimerax.std_commands.align import *
from chimerax.std_commands.color import color
from chimerax.std_commands.select import select

from .lovo import *


def lovo(session, atoms, to=None, n_clusters = 3,plot = False):
    from chimerax.match_maker.match import cmd_match
    from .lovo import normalize_msd

    for ret in cmd_match(session, match_atoms=atoms, to=to, cutoff_distance=None, compute_ss=False):
        ref_atoms, match_atoms, paired_RMSD, overall_RMSD, transformation_matrix = ret

        ref_atoms = ref_atoms.residues.atoms
        match_atoms = match_atoms.residues.atoms


        match_atoms_traj = np.array([match_atoms.scene_coords])

        # Cluster based on indices
        section_indices, msd = alovo_traj(match_atoms_traj, ref=ref_atoms.scene_coords, n_clusters=n_clusters)

        b_factor = rmsf(match_atoms_traj)

        section_colors = [[0, 240, 0, 1],  # Green RGBA
                          [240, 130, 0, 1],  # Yellow
                          [240, 0, 0, 1],  # Red
                          ]

        if plot:
            # Histogram Plot
            plt.figure(0)
            bins = np.linspace(0, max(msd), int(max(msd)) * 10)
            for i, indices in enumerate(section_indices):
                plt.hist(msd[indices], bins=bins, alpha=0.5, label=str(i), color=[np.array(section_colors[i]) / 240])
            plt.xlabel('RMSD')
            plt.ylabel('Count')
            plt.show()

            # Scatter Plot
            plt.figure(1)
            for i, indices in enumerate(section_indices):
                plt.scatter(msd[indices], normalize_msd(msd)[indices], color=[np.array(section_colors[i]) / 240],
                            alpha=1)
            plt.xlabel('RMSD')
            plt.ylabel('nRMSD')
            plt.show()

            # # rmsf vs msd
            # plt.figure(2)
            # for i, indices in enumerate(section_indices):
            #     plt.scatter(b_factor[indices], msd[indices], color=[np.array(section_colors[i]) / 240], alpha=1)
            # plt.xlabel('B-factor')
            # plt.ylabel('RMSD')
            # plt.show()

        for i in range(len(section_indices)):
            color(session, Objects(match_atoms[section_indices[i]]), Color(rgba=section_colors[i]))
            color(session, Objects(ref_atoms[section_indices[i]]), Color(rgba=section_colors[i]))

        match_atoms.scene_coords = match_atoms.scene_coords - (centroid(match_atoms.scene_coords) - centroid(ref_atoms.scene_coords))
        tf, rmsd = align_points(match_atoms.scene_coords[section_indices[0]], ref_atoms.scene_coords[section_indices[0]])
        move_atoms(match_atoms, ref_atoms, tf, move='structure')

    return


lovo_desc = CmdDesc(required=[('atoms', AtomsArg)],
                    keyword=[('to', AtomsArg),
                             ('move', EnumOf(('atoms', 'residues', 'chains', 'structures',
                                              'structure atoms', 'nothing'))),
                             ('n_cluster', IntArg),
                             ('plot', BoolArg),

                             ],
                    required_arguments=['to'],
                    synopsis='Align one set of atoms to another')


def lovo_traj(session, match_atoms, to=None, n_clusters=3, plot=False):
    match_structures = match_atoms.structures.unique()

    if to is None:
        ref_atoms = match_structures[0].atoms.intersect(match_atoms)
        match_atoms_list = [structure.atoms.intersect(match_atoms) for structure in match_structures[1:]]
    else:
        ref_atoms = to
        match_atoms_list = [structure.atoms.intersect(match_atoms) for structure in match_structures]

    ref_atoms = ref_atoms[np.argsort(ref_atoms.serial_numbers)]
    match_atoms_list = [match_atoms[np.argsort(match_atoms.serial_numbers)] for match_atoms in match_atoms_list]

    match_atoms_list[0] = match_atoms_list[0][np.argsort(match_atoms_list[0].serial_numbers)]

    match_atoms_traj = np.array([atoms.scene_coords for atoms in match_atoms_list])

    # Cluster based on indices
    section_indices, msd = alovo_traj(match_atoms_traj, ref=ref_atoms.scene_coords, n_clusters=n_clusters)

    b_factor = rmsf(match_atoms_traj)

    section_colors = [[0, 240, 0, 1],  # Green RGBA
                      [240, 130, 0, 1],  # Yellow
                      [240, 0, 0, 1],  # Red
                      ]

    if plot:
        # Histogram Plot
        plt.figure(0)
        bins = np.linspace(0, max(msd), int(max(msd)) * 10)
        for i, indices in enumerate(section_indices):
            plt.hist(msd[indices], bins=bins, alpha=0.5, label=str(i), color=[np.array(section_colors[i]) / 240])
        plt.xlabel('RMSD')
        plt.ylabel('Count')
        plt.show()

        # Scatter Plot
        plt.figure(1)
        for i, indices in enumerate(section_indices):
            plt.scatter(msd[indices], normalize_msd(msd)[indices], color=[np.array(section_colors[i]) / 240], alpha=1)
        plt.xlabel('RMSD')
        plt.ylabel('nRMSD')
        plt.show()

        #rmsf vs msd
        # plt.figure(2)
        # for i, indices in enumerate(section_indices):
        #     plt.scatter(b_factor[indices], msd[indices], color=[np.array(section_colors[i]) / 240], alpha=1)
        # plt.xlabel('B-factor')
        # plt.ylabel('RMSD')
        # plt.show()

    for i in range(len(section_indices)):
        for atoms in match_atoms_list:
            color(session, Objects(atoms[section_indices[i]]), Color(rgba=section_colors[i]))
        color(session, Objects(ref_atoms[section_indices[i]]), Color(rgba=section_colors[i]))

    for atoms in match_atoms_list:
        atoms.scene_coords = atoms.scene_coords - (centroid(atoms.scene_coords) - centroid(ref_atoms.scene_coords))
        tf, rmsd = align_points(atoms.scene_coords[section_indices[0]], ref_atoms.scene_coords[section_indices[0]])
        move_atoms(atoms, ref_atoms, tf, move='structure')

    return


lovo_traj_desc = CmdDesc(required=[('match_atoms', AtomsArg)],
                         keyword=[('to', AtomsArg),
                                  ('n_clusters', IntArg),
                                  ('plot', BoolArg),
                                  ],
                         synopsis='Align one set of atoms to another')