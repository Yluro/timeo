# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 José Serrano Guarinos

"""Thin wrapper around the cosmochlore CLI (cshm / csom / odis).

cosmochlore (https://github.com/Yluro/cosmochlore) is a separate Rust project
to calculate some shape and symmetry measures.

This module finds the executable, version checks it,
turning a MolecularStructure into the .xyz file it expects,
building argv for each of its three subcommands, and running it.

It knows nothing about Olex2's phil/param system (Timeo.py) will
read whatever settings they need and pass plain values to cosmochlore.py.
"""
import os
import re
import subprocess
import shutil

from selection import MolecularStructure

MIN_VERSION = (1, 2, 0)

# Point groups accepted by csom's -p/--pg flag, in Schoenflies notation.
# Names are case-sensitive and must match one of these exactly.
POINT_GROUPS = {
    'C2', 'C2h', 'C2v', 'C3', 'C3h', 'C3v', 'C4', 'C4h', 'C4v', 'C5', 'C5h', 'C5v',
    'C6', 'C6h', 'C6v', 'C7', 'C8', 'Ci', 'Cs', 'D2', 'D2d', 'D2h', 'D3', 'D3d',
    'D3h', 'D4', 'D4d', 'D4h', 'D5', 'D5d', 'D5h', 'D6', 'D6h', 'D7h', 'D8h', 'E',
    'I', 'Ih', 'O', 'Oh', 'S10', 'S4', 'S6', 'S8', 'T', 'Td', 'Th',
}

CENTERING_MODES = {'auto', 'first', 'centroid', 'manual'}

# Olex2 has no console, so Windows would pop up a console window for every cosmochlore run.
_NO_WINDOW = getattr(subprocess, 'CREATE_NO_WINDOW', 0)

# (exe_path, mtime) -> version, so redrawing the GUI does not rerun `cosmochlore --version`.
_version_cache = {}


class CosmochloreError(Exception):
    """Raised for anything that stops a cosmochlore run: not found, too old,
    bad arguments, or a non-zero exit. `str(error)` is fit to print as-is."""


def get_version(exe_path):
    """Returns the (major, minor, patch) version of `exe_path`, or None if it
    could not be run or its output could not be parsed."""
    try:
        key = (exe_path, os.path.getmtime(exe_path))
    except OSError:
        return None
    if key in _version_cache:
        return _version_cache[key]

    try:
        result = subprocess.run([exe_path, '--version'], capture_output=True,
                                text=True, timeout=10, creationflags=_NO_WINDOW)
    except (OSError, subprocess.SubprocessError):
        return None

    for line in reversed(result.stdout.splitlines()):
        match = re.search(r'(\d+)\.(\d+)\.(\d+)\s*$', line.strip())
        if match:
            _version_cache[key] = tuple(int(g) for g in match.groups())
            return _version_cache[key]
    return None


_DOWNLOAD_HINT = 'Get the latest version of cosmochlore from: https://github.com/Yluro/cosmochlore'


def _format_version(version):
    return '.'.join(str(v) for v in version)


def locate_cosmochlore(configured_path=None):
    """Returns the cosmochlore executable (MIN_VERSION or newer), or raises CosmochloreError."""
    if configured_path:
        exe = configured_path
        if not os.path.isfile(exe):
            raise CosmochloreError(f'timeo.cosmochlore.exe_path is set to "{exe}", '
                                   f'but that file does not exist.')
    else:
        exe = shutil.which('cosmochlore')
        if exe is None:
            raise CosmochloreError(f'cosmochlore executable not found on PATH. {_DOWNLOAD_HINT}')

    version = get_version(exe)
    if version is None:
        raise CosmochloreError(f'Could not determine the version of the cosmochlore '
                               f'executable at "{exe}".')

    if version < MIN_VERSION:
        raise CosmochloreError(
            f'cosmochlore at "{exe}" is version {_format_version(version)}, but '
            f'{_format_version(MIN_VERSION)} or newer is required. {_DOWNLOAD_HINT}')

    return exe


def can_find_cosmochlore_msg(configured_path=None, silent=True) -> bool:
    """Console helper mirroring autoshape.can_find_shape_msg()."""
    try:
        exe = locate_cosmochlore(configured_path)
    except CosmochloreError as e:
        print(e)
        return False

    if not silent:
        print(f'cosmochlore {_format_version(get_version(exe))} found at: {exe}')
    return True


def write_xyz(structure: MolecularStructure, path, comment=''):
    """Writes `structure` to a .xyz file, atoms in the order MolecularStructure
    holds them. A centered structure's centre (index 0).
    """
    with open(path, 'w') as f:
        f.write(f'{len(structure)}\n')
        f.write(f'{comment}\n')
        for label, (x, y, z) in zip(structure.labels, structure.coords):
            f.write(f'{label} {x:.6f} {y:.6f} {z:.6f}\n')


def run_cosmochlore(exe, args, cwd=None, timeout=300) -> str:
    """
    Runs `exe` with `args`. Returns captured stdout on success. Raises CosmochloreError carrying
    cosmochlore's own stderr message when there is one.
    """
    try:
        result = subprocess.run([exe] + args, cwd=cwd, capture_output=True,
                                text=True, timeout=timeout, creationflags=_NO_WINDOW)
    except subprocess.TimeoutExpired:
        raise CosmochloreError(f'cosmochlore timed out after {timeout}s (args: {args}).')
    except OSError as e:
        raise CosmochloreError(f'Could not run cosmochlore at "{exe}": {e}')

    if result.returncode != 0:
        message = result.stderr.strip() or (
            f'cosmochlore exited with code {result.returncode} and no error message.')
        raise CosmochloreError(message)

    return result.stdout


def trim_banner(output: str) -> str:
    """Drops cosmochlore's ASCII banner/preamble from captured stdout"""
    lines = output.splitlines()
    for i, line in enumerate(lines):
        if line.startswith('Input file:'):
            return '\n'.join(lines[i:])
    return output


def _centering_args(centered: bool, center_index: int = 1):
    """The -c/-n argv fragment shared by cshm and csom."""
    return ['-c', str(center_index)] if centered else ['-n']


def build_cshm_args(xyz_path, centered, shapes=None, user_shapes=None,
                    table=False, ideal=False):
    """Builds argv for `cosmochlore cshm`. `shapes` is a list of built-in
    reference-shape indices; `user_shapes` a list of YAML file paths."""
    args = ['cshm', xyz_path] + _centering_args(centered)

    if shapes:
        args += ['-s'] + [str(i) for i in shapes]
    if user_shapes:
        args += ['-r'] + list(user_shapes)
    if table:
        args.append('-t')
    if ideal:
        args.append('-i')

    return args


def build_csom_args(xyz_path, centered, point_groups, mode='auto', vector=None,
                    full=False, table=False, operated=False, seeds=None,
                    iterations=None, tolerance=None, ignore_labels=False):
    """Builds argv for `cosmochlore csom`. Validates `point_groups` and `mode`
    itself and raises CosmochloreError on anything invalid."""
    if not point_groups:
        raise CosmochloreError('csom requires at least one point group.')

    unknown = sorted(pg for pg in point_groups if pg not in POINT_GROUPS)
    if unknown:
        raise CosmochloreError(f'Unknown point group(s): {", ".join(unknown)}. Point group '
                               f'names are case-sensitive (e.g. "Oh", not "oh" or "OH").')

    if mode not in CENTERING_MODES:
        raise CosmochloreError(f'Unknown centering mode "{mode}". Expected one of: '
                               f'{", ".join(sorted(CENTERING_MODES))}.')
    if mode == 'manual' and not vector:
        raise CosmochloreError('centering mode "manual" requires a 3-value vector.')

    # cosmochlore accepts 0 for any of these and quietly returns an unrefined
    # (wrong) deviation instead of erroring.
    for name, value in (('seeds', seeds), ('iterations', iterations), ('tolerance', tolerance)):
        if value is not None and not value > 0:
            raise CosmochloreError(f'csom {name} must be greater than 0, got {value!r}.')

    args = ['csom', xyz_path] + _centering_args(centered)
    args += ['-p'] + list(point_groups)
    args += ['-m', mode]

    if vector:
        args += ['-u'] + [str(v) for v in vector]
    if full:
        args.append('-f')
    if table:
        args.append('-t')
    if operated:
        args.append('-o')
    if seeds is not None:
        args += ['-s', str(seeds)]
    if iterations is not None:
        args += ['-i', str(iterations)]
    if tolerance is not None:
        args += ['-T', f'{tolerance:g}']
    if ignore_labels:
        args.append('-g')

    return args


def build_odis_args(xyz_path, full=False, table=False):
    """Builds argv for `cosmochlore odis`. odis has no --nc: it always expects
    a centered, 7-atom structure with the centre at position 1."""
    args = ['odis', xyz_path, '-c', '1']
    if full:
        args.append('-f')
    if table:
        args.append('-t')
    return args


def list_user_shapes(shapes_dir):
    """Returns the .yaml reference-shape files found in `shapes_dir`.

    Each entry is (display_name, full_path): display_name is the filename
    without its .yaml extension, full_path is what cshm's -r/--ref expects.
    Returns [] if the directory doesn't exist.
    """
    if not os.path.isdir(shapes_dir):
        return []
    names = sorted(f for f in os.listdir(shapes_dir) if f.lower().endswith('.yaml'))
    return [(os.path.splitext(f)[0], os.path.join(shapes_dir, f)) for f in names]