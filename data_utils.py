import io
import re
import math
import pathlib
import zipfile

import numpy as np


###############################################################################
# Basic geometry manipulation utilities
###############################################################################
def split_to_atoms(geom: str):
    atoms = []
    coords = []
    for line in geom.strip().split('\n'):
        atom = line.split()[0]
        coord = [float(coord) for coord in line.split()[1:]]
        atoms.append(atom)
        coords.append(coord)

    return atoms, coords


###############################################################################
# ORCA data parsing utilities
###############################################################################
def matrices_from_orca(out_file: str, dummy: bool = False):
    if dummy:
        dm = np.array([[1., 1., 1., 1.],
                       [1., 1., 1., 1.],
                       [1., 1., 1., 1.],
                       [1., 1., 1., 1.]])
        ovrlp = np.array([[0.1, 0.1, 0.1, 0.1],
                          [0.1, 0.1, 0.1, 0.1],
                          [0.1, 0.1, 0.1, 0.1],
                          [0.1, 0.1, 0.1, 0.1]])
        atoms = [[0.09], [0.09]]
        coords = [[0., 0., 0.], [0., 0., 1.]]
        nbas_tot = 4
        natoms = 2
        return dm, ovrlp, atoms, coords, nbas_tot, natoms

    supported_orca_versions = [
        '2.9.1',
        '5.0.0',
        '5.0.1',
        '5.0.2',
        '5.0.3',
        '5.0.4',
        '6.0.0',
        '6.0.1']
    charg = {
        'H': 1.0,
        'B': 5.0,
        'C': 6.0,
        'N': 7.0,
        'O': 8.0,
        'F': 9.0,
        'Si': 14.0,
        'P': 15.0,
        'S': 16.0,
        'Cl': 17.0}
    nbas_tot = 0
    natoms = 0
    basis_pattern = \
        (r'Group +1 +Type +[A-Z][a-z]? +: '
         r'\d+s\d+p\d*[d]?\d*[f]? contracted to '
         r'(\d+)s(\d+)p(\d*)[d]?(\d*)[f]? pattern')
    shells = [1, 3, 5, 7]
    geom = ''
    if out_file.endswith('.zip'):
        z = zipfile.ZipFile(out_file)
        f = io.TextIOWrapper(z.open(f'{pathlib.Path(out_file).stem}.out'),
                             encoding='utf-8')
    else:
        f = open(out_file)
    for line in f:
        if 'Program Version' in line:
            orca_version = line.split()[2]
            if orca_version not in supported_orca_versions:
                raise NotImplementedError(f'Orca version {orca_version}'
                                          f' outputs are not supported')
        if line == 'CARTESIAN COORDINATES (ANGSTROEM)\n':
            f.__next__()
            line = f.__next__()
            while len(line) > 1:
                geom += line
                line = f.__next__()
        if line == 'BASIS SET INFORMATION\n':
            for i in range(4):
                line = f.__next__()
            m = re.search(basis_pattern, line)
            for i in range(len(shells)):
                shells[i] = \
                    int(m.groups()[i]) * shells[i] if m.groups()[i] else 0
        if 'Basis Dimension' in line and nbas_tot == 0:
            nbas_tot = int(line.split()[-1])
            dm = np.zeros((nbas_tot, nbas_tot))
            ovrlp = np.zeros((nbas_tot, nbas_tot))
            natoms = nbas_tot // sum(shells)
        if line == 'DENSITY\n':
            f.__next__()
            for block in range(math.ceil(nbas_tot / 6)):
                f.__next__()
                for row in range(nbas_tot):
                    line = f.__next__()
                    cols = len(line.split()) - 1
                    column = 0
                    while column < cols:
                        dm[row, column + 6 * block] = \
                            float(line.split()[column + 1])
                        column += 1
        if line == 'OVERLAP MATRIX\n':
            f.__next__()
            for block in range(math.ceil(nbas_tot / 6)):
                f.__next__()
                for row in range(nbas_tot):
                    line = f.__next__()
                    cols = len(line.split()) - 1
                    column = 0
                    while column < cols:
                        ovrlp[row, column + 6 * block] = \
                            float(line.split()[column + 1])
                        column += 1
    f.close()
    if out_file.endswith('.zip'):
        z.close()
    atoms, coords = split_to_atoms(geom)
    atoms = [[charg[atom] / 100.0] for atom in atoms]

    return dm, ovrlp, atoms, coords, nbas_tot, natoms


def blocks_from_orca(out_file: str, overlap_thresh: float,
                     dummy: bool = False):
    dm, ovrlp, atoms, coords, nbas_tot, natoms = \
        matrices_from_orca(out_file, dummy=dummy)
    nbas = nbas_tot // natoms
    diagonal_densities = []
    off_diagonal_densities = []
    off_diagonal_overlaps = []
    adjacency_atom2link_sources = []
    adjacency_atom2link_targets = []
    adjacency_link2atom_sources = []
    adjacency_link2atom_targets = []
    nlinks = 0
    for atom in range(natoms):
        diagonal_density = \
            np.hstack([dm[atom * nbas + bas: atom * nbas + bas + 1,
                       atom * nbas + bas: (atom + 1) * nbas]
                       for bas in range(nbas)])
        diagonal_densities.append(diagonal_density.flatten().tolist())
        for other_atom in range(natoms):
            if atom == other_atom:
                continue
            off_diagonal_overlap = \
                ovrlp[atom * nbas: (atom + 1) * nbas,
                      other_atom * nbas: (other_atom + 1) * nbas]
            if np.abs(off_diagonal_overlap).max() < overlap_thresh:
                continue
            off_diagonal_overlaps.append(
                off_diagonal_overlap.flatten().tolist())
            adjacency_atom2link_sources.append(atom)
            adjacency_atom2link_targets.append(nlinks)
            adjacency_link2atom_sources.append(nlinks)
            adjacency_link2atom_targets.append(other_atom)
            nlinks += 1
            off_diagonal_density = \
                dm[atom * nbas: (atom + 1) * nbas,
                   other_atom * nbas: (other_atom + 1) * nbas]
            off_diagonal_densities.append(
                off_diagonal_density.flatten().tolist())

    return (diagonal_densities,
            off_diagonal_densities,
            off_diagonal_overlaps,
            adjacency_atom2link_sources,
            adjacency_atom2link_targets,
            adjacency_link2atom_sources,
            adjacency_link2atom_targets,
            atoms, natoms, nlinks, nbas)


def serv_arg_name():
    idx = 0
    while True:
        yield f'flat_args_0_{idx}' if idx else 'flat_args_0'
        idx += 1


def served_input_from_orca(out_file: str,
                           overlap_thresh: float,):
    diagonal_densities,               \
        off_diagonal_densities,       \
        off_diagonal_overlaps,        \
        adjacency_atom2link_sources,  \
        adjacency_atom2link_targets,  \
        adjacency_link2atom_sources,  \
        adjacency_link2atom_targets,  \
        atoms, natoms, nlinks, nbas = \
        blocks_from_orca(out_file, overlap_thresh)

    output_tensors = {}
    idx = serv_arg_name()

    for i in range(6):
        output_tensors[next(idx)] = \
            np.array([0.], dtype=np.float32).reshape((1, 1))
    output_tensors[next(idx)] = np.array([1], dtype=np.int32)
    output_tensors[next(idx)] = \
        np.array(adjacency_atom2link_sources, dtype=np.int32)
    output_tensors[next(idx)] = \
        np.array(adjacency_atom2link_targets, dtype=np.int32)
    output_tensors[next(idx)] = \
        np.array(off_diagonal_overlaps, dtype=np.float32)
    output_tensors[next(idx)] = \
        np.array([len(off_diagonal_overlaps)], dtype=np.int32)
    output_tensors[next(idx)] = \
        np.array(adjacency_link2atom_sources, dtype=np.int32)
    output_tensors[next(idx)] = \
        np.array(adjacency_link2atom_targets, dtype=np.int32)
    output_tensors[next(idx)] = \
        np.array(off_diagonal_overlaps, dtype=np.float32)
    output_tensors[next(idx)] = \
        np.array([len(off_diagonal_overlaps)], dtype=np.int32)
    output_tensors[next(idx)] = \
        np.array(diagonal_densities, dtype=np.float32)
    output_tensors[next(idx)] = \
        np.array(atoms, dtype=np.float32)
    output_tensors[next(idx)] = \
        np.array([len(diagonal_densities)], dtype=np.int32)
    output_tensors[next(idx)] = \
        np.array(off_diagonal_densities, dtype=np.float32)
    output_tensors[next(idx)] = \
        np.array([len(off_diagonal_densities)], dtype=np.int32)

    return output_tensors
