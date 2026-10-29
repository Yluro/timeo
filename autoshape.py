# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 José Serrano Guarinos

"""Generating SHAPE 2.1 input, running it, and parsing what it writes back."""
import os
import shutil
import subprocess
from typing import List

import numpy as np

from constants import REF_SHAPE_DICT
from selection import MolecularStructure

# Olex2 has no console, so Windows would pop up a console window for every SHAPE run.
_NO_WINDOW = getattr(subprocess, 'CREATE_NO_WINDOW', 0)


class ShapeCalculation:
    """Builds the .dat input file SHAPE 2.1 expects for one structure."""

    def __init__(self, structure: MolecularStructure, file_name: str, centered: bool,
                 keywords: List[str]):

        self.coords = np.array(structure.coords, dtype=np.float64)
        self.labels = structure.labels

        self.is_centered = centered
        self.file_name = file_name
        self.dat_keywords = '\n'.join(keywords) + '\n'

        if centered:
            self.centre = 1
            self.ligands = len(self.coords) - 1
            self.dat_title = f'$ {self.file_name}_{self.labels[0]}\n'
            self.subtitle = f'{self.labels[0]}\n'.upper()  # First label contains the metal atom.
        else:
            self.centre = 0
            self.ligands = len(self.coords)
            self.dat_title = f"$ {self.file_name}_{self.labels[0]}-{self.labels[-1]}\n"
            self.subtitle = f'{self.labels[0]}-{self.labels[-1]}\n'.upper()

        self.subfolder_name = self.dat_title.strip().split(' ')[-1]
        self.positions = f'{self.ligands} {self.centre}\n'

        # The SHAPE 2.1 documentation gives these strings of numbers.
        if self.ligands not in REF_SHAPE_DICT:
            raise ValueError(f'No reference geometries defined for {self.ligands} vertices. '
                             f'Check the structure for extra bonds or atoms.')
        self.geometries = f'{REF_SHAPE_DICT[self.ligands]}\n'

        # Builds the Label x y z table.
        self.table = '\n'.join(
            f'{label} {" ".join(f"{x:g}" for x in xyz)}'
            for label, xyz in zip(self.labels, self.coords)
        )  # Joins all rows of the table

        self.dat_file_contents = (self.dat_title + self.dat_keywords + self.positions
                                  + self.geometries + self.subtitle + self.table)

    def write_tab(self, file_path):
        """Writes the .dat file into a fresh numbered run folder.

        Returns the folder it wrote to, or None if all 99 run slots are taken.
        """
        for i in range(99):
            file_name = f'{self.subfolder_name}_{i}.dat'
            file_dir = os.path.join(file_path, 'autoSHAPE', self.subfolder_name, str(i))
            dat_file_path = os.path.join(file_dir, file_name)

            if os.path.exists(dat_file_path):
                continue

            os.makedirs(file_dir, exist_ok=True)

            with open(dat_file_path, 'w') as f:
                print(f'Writing {file_name} at {dat_file_path}...')
                f.write(self.dat_file_contents)

            return file_dir

        print(f'Could not write {self.subfolder_name}: all 99 run folders already exist.')
        return None


class ShapeError(Exception):
    """Raised when no usable SHAPE executable can be found; the message is fit for the console."""


def locate_shape(configured_path=None) -> str:
    """Returns the SHAPE 2.1 executable: `configured_path` if given, else the one on PATH."""
    if not configured_path:
        exe = shutil.which('shape')
        if exe is None:
            raise ShapeError('Unable to find shape.exe in the system path.')
        return exe

    if not os.path.isfile(configured_path):
        raise ShapeError(f'timeo.shape.exe_path is set to "{configured_path}", '
                         f'but that file does not exist.')
    if 'shape' not in os.path.basename(configured_path).lower():
        raise ShapeError(f'timeo.shape.exe_path is set to "{configured_path}", '
                         f'which does not look like the SHAPE executable.')
    return configured_path


def can_find_shape_msg(configured_path=None, silent=True) -> bool:
    try:
        exe = locate_shape(configured_path)
    except ShapeError as e:
        print(e)
        return False

    if not silent:
        print(f"SHAPE executable found at: {exe}")
    return True


def run_shape(folder, exe='shape') -> List[str]:
    """
    Runs a SHAPE instance on all dat files in a specified folder.
    :param folder:
    :param exe: SHAPE executable to run.
    :return files: Name of the output files without file extension.
    """
    if not folder:
        return []

    print("Running SHAPE...")
    dat_files = [f for f in os.listdir(folder) if f.endswith('.dat')]
    for file in dat_files:
        process = subprocess.Popen([exe], stdin=subprocess.PIPE,
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                   text=True, cwd=folder, creationflags=_NO_WINDOW)
        out, err = process.communicate(input=f'{file}\n')  # Send the file name and an Enter key
        print(out)

    return [f[:-4] for f in dat_files]


def parse_shape_tab(tab_path):
    """Reads a SHAPE .tab file. Returns (atom_label, shape_names, shape_labels, values)."""
    with open(tab_path, 'r') as f:
        lines = f.readlines()

    # Find the polyhedra table (Label  N  Symm  Name).
    # The slice skips the header of the tab file, which would otherwise confuse the search.
    shape_labels = []
    shape_names = []
    for tab_line in lines[5:]:
        if len(tab_line.split()) >= 4 and tab_line.split()[1].isdigit():
            shape_labels.append(tab_line.split()[0])
            shape_names.append(' '.join(tab_line.split()[3:]))

    # Second pass through the tab file to find the data row with the CShMs.
    data_row = None
    for tab_line in lines[5:]:
        if ',' in tab_line and any(c.isdigit() for c in tab_line):
            data_row = tab_line

    if data_row is None:
        print(f"Could not parse data row in {tab_path}")
        return None

    parts = [c.strip() for c in data_row.split(',')]
    atom_label = parts[0]
    values = [float(v) for v in parts[1:]]

    return atom_label, shape_names, shape_labels, values


def print_shape_table(tab_path):
    result = parse_shape_tab(tab_path)
    if result is None:
        return False

    atom_label, shape_names, shape_labels, values = result
    max_name_length = max(len(name) for name in shape_names) + 2
    min_val = min(values)

    print('\n' + '=' * 65)
    print(f'SHAPE 2.1 results for {atom_label} in {os.path.basename(tab_path)}:')
    print('-' * 65)
    print(f"{'Polyhedron':<{max_name_length}}{'Symbol':<10}{'CShM':>10}")
    print('-' * 65)
    for name, label, val in zip(shape_names, shape_labels, values):
        marker = ' <---- best fit' if val == min_val else ''
        print(f'{name:<{max_name_length}}{label:<10}{val:>10.3f}{marker:>10}')
    print('=' * 65)

    if min_val > 10:
        print(f'WARNING: extremely distorted geometries found for {atom_label}. '
              f'Make sure there is no unnatural bonds in the model.')
    return True
