# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 José Serrano Guarinos

"""Draws cosmochlore's ideal polyhedra (cshm --ideal) over the model in Olex2's 3D view.

Every CShM run with an _ideal.xyz or a _cshm_table.csv in <FilePath>/cosmochlore/
for the current structure is a "run", listed in the CShM runs table; all ticked
runs with an _ideal.xyz are drawn at once, styled by the timeo.overlay phil params.

Built on Olex2's own AddObject macro (olex/xpprocs.cpp, TMainForm::macAddObject):
  AddObject poly <name> 0 "x y z,x y z,...;i j k,i j k,..."   Cartesian vertices, 0-based triangles
  AddObject sphere <name> 0 <radius> <x> <y> <z>               Cartesian unless -c is given
  SetMaterial <name>.Object <material>                         colour/transparency
  kill <name>                                                  hides the object again
"""
import csv
import json
import os
import re
from itertools import combinations

import numpy as np
import olex
import olx
from olexFunctions import OlexFunctions
from scipy.spatial import ConvexHull

from helper_functions import as_bool, get_id_from_label, get_orm_atoms

OV = OlexFunctions()

# Below this ratio of singular values the vertices are treated as coplanar / collinear.
_FLAT_TOL = 1e-3
# Faces whose unit normals agree this closely are one face split by the triangulation.
_COPLANAR_TOL = 1e-3
# Gap between the two copies of a flat shape's faces, one facing each side.
_SIDE_GAP = 0.005

# Face colours (r, g, b), handed out per shape or per run in drawing order.
_PALETTE = [(70, 130, 230), (230, 120, 40), (60, 180, 90), (200, 60, 160),
            (220, 200, 40), (60, 190, 200), (150, 90, 220), (200, 70, 60)]

MARKER_RADIUS = 0.12
EDGE_RADIUS = 0.035
EDGE_SIDES = 8

_IDEAL = '_ideal.xyz'
_RUN_FILE = re.compile(r'^(?P<label>.+)_(?P<part>\d+)(?P<kind>_ideal\.xyz|_cshm_table\.csv)$')
# Maps each run's file stem to whether it was centered, written by record_run().
_REGISTRY = 'ideal_runs.json'

# Names of the objects drawn so far, so clear_ideal() can hide them again.
_drawn = []
# Olex2 reuses the geometry of an existing collection with the same name, so every
# draw needs fresh object names even after the old ones are hidden.
_draw_count = 0
# Face objects currently drawn and their colours, so opacity can change in place.
_faces = {}
# Runs table state, this session only: stems of runs with Draw ticked, and the shape
# picked for each run when it isn't the best match.
_ticked = set()
_run_shapes = {}


def read_ideal_xyz(path):
    """Parses a multi-frame _ideal.xyz into dicts with shape, cshm, labels and coords."""
    with open(path) as f:
        lines = [line.rstrip() for line in f]

    frames = []
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i])
        comment = lines[i + 1]
        atoms = [line.split() for line in lines[i + 2:i + 2 + n]]
        i += 2 + n

        match = re.search(r'CShM\s*=\s*([-+0-9.eE]+)', comment)
        frames.append({
            'shape': comment.split()[0] if comment.split() else f'shape{len(frames)}',
            'cshm': float(match.group(1)) if match else float('inf'),
            'labels': [a[0] for a in atoms],
            'coords': np.array([[float(c) for c in a[1:4]] for a in atoms]),
        })
    return frames


def triangulate(points):
    """Returns (vertices, triangles) covering `points`' convex hull, or None if collinear."""
    points = np.asarray(points, dtype=np.float64)
    if len(points) < 3:
        return None

    centroid = points.mean(axis=0)
    _, s, vt = np.linalg.svd(points - centroid)
    if s[1] < _FLAT_TOL * s[0]:
        return None

    if s[2] >= _FLAT_TOL * s[0]:
        return points, ConvexHull(points).simplices

    # Planar (SP-4, TP-3, HP-6...): fan the in-plane hull outline around its centroid.
    outline = ConvexHull((points - centroid) @ vt[:2].T).vertices
    vertices = np.vstack([points, centroid])
    hub = len(points)
    triangles = [(hub, outline[k], outline[(k + 1) % len(outline)]) for k in range(len(outline))]
    return vertices, np.array(triangles)


def face_meshes(points):
    """(vertices, triangles) of each face object for `points`' hull, flat shapes from both sides."""
    points = np.asarray(points, dtype=np.float64)
    mesh = triangulate(points)
    if mesh is None:
        return []
    vertices, triangles = mesh
    _, s, vt = np.linalg.svd(points - points.mean(axis=0))
    if s[2] >= _FLAT_TOL * s[0]:
        return [mesh]

    # Olex2 culls the back of transparent faces, and AddObject poly turns every triangle
    # to face away from the mean of the object's vertices, a coin-flip for a flat shape.
    # So a flat shape gets a copy a hair to either side of its plane, around a mean in it.
    centroid = vertices.mean(axis=0)
    flat = centroid + (vertices - centroid) @ vt[:2].T @ vt[:2]
    offset = _SIDE_GAP / 2 * vt[2]
    return [(np.vstack([flat + offset, flat - offset]),
             np.vstack([triangles, triangles + len(flat)]))]


def polyhedron_edges(points):
    """The (a, b) point pairs forming the polyhedron's edges, without triangulation diagonals."""
    points = np.asarray(points, dtype=np.float64)
    mesh = triangulate(points)
    if mesh is None:
        # Collinear: join the points in order along the line.
        if len(points) < 2:
            return []
        axis = np.linalg.svd(points - points.mean(axis=0))[2][0]
        order = np.argsort(points @ axis)
        return [(points[i], points[j]) for i, j in zip(order, order[1:])]

    vertices, triangles = mesh
    faces_of = {}
    for t, tri in enumerate(triangles):
        for i, j in combinations(sorted(tri), 2):
            faces_of.setdefault((i, j), []).append(t)

    def normal(t):
        a, b, c = vertices[triangles[t]]
        n = np.cross(b - a, c - a)
        return n / np.linalg.norm(n)

    edges = []
    for (i, j), faces in faces_of.items():
        # Kept: outline edges (one face) and folds between two non-coplanar faces.
        if len(faces) == 2 and abs(np.dot(normal(faces[0]), normal(faces[1]))) > 1 - _COPLANAR_TOL:
            continue
        edges.append((vertices[i], vertices[j]))
    return edges


def tube(a, b, radius=EDGE_RADIUS, sides=EDGE_SIDES):
    """(vertices, triangles) of an open-ended cylinder from `a` to `b`."""
    axis = (b - a) / np.linalg.norm(b - a)
    helper = np.array([1.0, 0.0, 0.0]) if abs(axis[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    p1 = np.cross(axis, helper)
    p1 /= np.linalg.norm(p1)
    p2 = np.cross(axis, p1)

    angles = 2 * np.pi * np.arange(sides) / sides
    ring = radius * (np.outer(np.cos(angles), p1) + np.outer(np.sin(angles), p2))
    vertices = np.vstack([a + ring, b + ring])
    triangles = []
    for k in range(sides):
        k1 = (k + 1) % sides
        triangles += [(k, k1, sides + k), (k1, sides + k1, sides + k)]
    return vertices, np.array(triangles)


def _rgba(r, g, b, a):
    """Packs a colour the way Olex2's TGlOption reads a uint32 (R in the low byte)."""
    return r | (g << 8) | (b << 16) | (a << 24)


def _material(rgb, alpha):
    """Front and back material string: ambient, diffuse, specular, shininess, transparent."""
    r, g, b = rgb
    a = int(round(alpha * 255))
    ambient = _rgba(r // 2, g // 2, b // 2, a)
    diffuse = _rgba(r, g, b, a)
    specular = _rgba(128, 128, 128, 255)
    # AmbientF|AmbientB|DiffuseF|DiffuseB|SpecularF|ShininessF, plus Transparent if alpha < 1.
    flags = 0x1 | 0x2 | 0x4 | 0x8 | 0x10 | 0x40 | (0x400 if a < 255 else 0)
    return f'{flags};{ambient};{ambient};{diffuse};{diffuse};{specular};32'


def _poly_command(name, vertices, triangles):
    verts = ','.join(f'{x:.4f} {y:.4f} {z:.4f}' for x, y, z in vertices)
    tris = ','.join(f'{i} {j} {k}' for i, j, k in triangles)
    return f'AddObject poly {name} 0 "{verts};{tris}"'


def options():
    """The timeo.overlay drawing options, normalised from whatever the GUI stored."""
    return {
        'opacity': min(100, max(0, int(OV.GetParam('timeo.overlay.opacity', 35) or 0))),
        'edges': as_bool(OV.GetParam('timeo.overlay.edges', True)),
        'markers': as_bool(OV.GetParam('timeo.overlay.markers', True)),
        'colour_by': str(OV.GetParam('timeo.overlay.colour_by', 'shape') or 'shape').strip().lower(),
    }


def shape_commands(frame, centered, rgb, base, opacity=35, edges=True, markers=True):
    """Returns (macro commands, object names, face object names) drawing one ideal shape."""
    commands, names = [], []
    donors = frame['coords'][1:] if centered else frame['coords']

    def add(name, vertices, triangles, alpha, colour=rgb):
        commands.append(_poly_command(name, vertices, triangles))
        commands.append(f'SetMaterial {name}.Object {_material(colour, alpha)}')
        names.append(name)

    faces = []
    if opacity > 0:
        for k, mesh in enumerate(face_meshes(donors)):
            add(f'{base}f{k}', *mesh, opacity / 100)
            faces.append(f'{base}f{k}')

    if edges:
        darker = tuple(int(c * 0.6) for c in rgb)
        for k, (a, b) in enumerate(polyhedron_edges(donors)):
            add(f'{base}e{k}', *tube(a, b), 1.0, darker)

    if markers:
        name = f'{base}v'
        for x, y, z in donors:
            commands.append(f'AddObject sphere {name} 0 {MARKER_RADIUS} {x:.4f} {y:.4f} {z:.4f}')
        commands.append(f'SetMaterial {name}.Object {_material(rgb, 1.0)}')
        names.append(name)
    return commands, names, faces


def _pick_frame(frames, choice):
    """The frame for shape `choice`, or the lowest-CShM one for 'best' or a shape not in `frames`."""
    # Olex2's combo can hand the value back with its case changed (e.g. 'Best').
    wanted = str(choice or 'best').strip().lower()
    match = next((f for f in frames if f['shape'].lower() == wanted), None)
    return match or min(frames, key=lambda f: f['cshm'])


def _element_rgb(label, orm):
    """The colour Olex2 draws the element of atom `label` in, or None if it can't be looked up."""
    try:
        element = olx.xf.au.GetAtomType(get_id_from_label(label, orm))
        colour = int(olx.GetMaterial(f'{element}.Sphere').split(';')[1])
    except (RuntimeError, TypeError, ValueError, IndexError):
        return None
    return colour & 255, (colour >> 8) & 255, (colour >> 16) & 255


def _workdir():
    return os.path.join(olx.FilePath(), 'cosmochlore')


def _read_registry():
    try:
        with open(os.path.join(_workdir(), _REGISTRY)) as f:
            return json.load(f)
    except (OSError, ValueError):
        return {}


def _write_registry(registry):
    try:
        with open(os.path.join(_workdir(), _REGISTRY), 'w') as f:
            json.dump(registry, f, indent=1)
    except OSError as e:
        print(f'Could not save {_REGISTRY}: {e}')


def record_run(paths, centered):
    """Remembers whether a CShM run's _ideal.xyz files are centered."""
    registry = _read_registry()
    for path in paths:
        registry[os.path.basename(path)[:-len(_IDEAL)]] = bool(centered)
    _write_registry(registry)


def find_runs():
    """Every CShM run on the current structure, as dicts sorted by stem.

    A run is an _ideal.xyz, a _cshm_table.csv or both; 'ideal' and 'table' hold their
    paths, or None for the one that wasn't written.
    """
    folder = _workdir()
    if not os.path.isdir(folder):
        return []

    prefix = f'{olx.FileName()}_'
    registry = _read_registry()
    runs = {}
    for name in os.listdir(folder):
        match = _RUN_FILE.match(name[len(prefix):]) if name.startswith(prefix) else None
        if match is None:
            continue
        stem = name[:-len(match.group('kind'))]
        run = runs.setdefault(stem, {
            'stem': stem,
            'ideal': None,
            'table': None,
            'label': match.group('label'),
            'part': int(match.group('part')),
            # Runs from before the registry existed: Timeo's usual single-atom selection.
            'centered': registry.get(stem, True),
        })
        run['ideal' if match.group('kind') == _IDEAL else 'table'] = os.path.join(folder, name)
    runs = sorted(runs.values(), key=lambda r: r['stem'])

    # 'Fe1', or 'Fe1 (2)' when that centre was split into several disorder components.
    for run in runs:
        split = sum(r['label'] == run['label'] for r in runs) > 1
        run['title'] = f'{run["label"]} ({run["part"] + 1})' if split else run['label']
    return runs


def read_results(run):
    """The run's shapes as dicts with shape and cshm, best first: from its _ideal.xyz, else its .csv."""
    try:
        if run['ideal']:
            results = read_ideal_xyz(run['ideal'])
        else:
            with open(run['table'], newline='') as f:
                results = [{'shape': row['Symbol'], 'cshm': float(row['CShM'])}
                           for row in csv.DictReader(f)]
    except (OSError, ValueError, IndexError, KeyError, TypeError):
        return []
    return sorted(results, key=lambda r: r['cshm'])


def draw_runs(runs, verbose=False):
    """Replaces the overlay with each run's chosen shape; `verbose` prints a one-line summary."""
    global _draw_count
    clear_ideal()

    opts = options()
    orm = get_orm_atoms() if opts['colour_by'] == 'element' else None
    colours = {}
    summary = []
    for i, run in enumerate(runs):
        try:
            frames = read_ideal_xyz(run['ideal'])
        except (OSError, ValueError, IndexError) as e:
            print(f'Could not read {os.path.basename(run["ideal"])}: {e}')
            continue
        if not frames:
            continue
        frame = _pick_frame(frames, _run_shapes.get(run['stem']))

        rgb = None
        if opts['colour_by'] == 'element':
            key = run['label']
            rgb = _element_rgb(run['label'], orm)
        else:
            key = run['stem'] if opts['colour_by'] == 'centre' else frame['shape']
        rgb = rgb or colours.setdefault(key, _PALETTE[len(colours) % len(_PALETTE)])

        base = f'timeo{_draw_count}_{i}_{re.sub(r"[^A-Za-z0-9]", "", frame["shape"])}'
        commands, names, faces = shape_commands(frame, run['centered'], rgb, base, opts['opacity'],
                                                opts['edges'], opts['markers'])
        for command in commands:
            olex.m(command)
        _drawn.extend(names)
        for face in faces:
            _faces[face] = rgb
        summary.append(f'{run["title"]} {frame["shape"]} ({frame["cshm"]:.3f})')

    _draw_count += 1
    OV.Refresh()
    if verbose and summary:
        print('Ideal polyhedra: ' + ', '.join(summary))
    return bool(summary)


def draw_ideal(verbose=True):
    """Draws the chosen ideal shape of every ticked run, or clears the overlay if none is."""
    runs = find_runs()
    ticked_runs = [run for run in runs if ticked(run)]
    if not ticked_runs:
        if verbose and any(run['ideal'] for run in runs):
            print('No CShM run is ticked to draw: tick Draw in the CShM runs table.')
        elif verbose:
            print('No ideal polyhedra for this structure: run CShM with '
                  '"Ideal polyhedra .xyz" checked.')
        clear_ideal()
        return False
    return draw_runs(ticked_runs, verbose)


def ticked(run):
    """Whether the run's Draw box is ticked; never for a run without an _ideal.xyz."""
    return bool(run['ideal']) and run['stem'] in _ticked


def any_ticked():
    """Whether any run on this structure is ticked to draw."""
    return any(ticked(run) for run in find_runs())


def refresh():
    """Quietly redraws the ticked runs, or clears the overlay if none is ticked."""
    draw_ideal(verbose=False)


def option_changed(param):
    """Applies a changed timeo.overlay option: opacity in place, anything else by redrawing."""
    opacity = options()['opacity']
    if param.endswith('.opacity') and _faces and opacity > 0:
        for name, rgb in _faces.items():
            olex.m(f'SetMaterial {name}.Object {_material(rgb, opacity / 100)}')
        OV.Refresh()
    else:
        refresh()


def set_ideal_run_visible(stem, state):
    """Ticks or unticks one run's Draw box, then redraws. False for a run without an _ideal.xyz."""
    run = next((r for r in find_runs() if r['stem'] == stem), None)
    if run is None or not run['ideal']:
        if as_bool(state):
            print(f'Error: no ideal polyhedra to draw for {run["title"] if run else stem}. '
                  f'Rerun cosmochlore (Run CShM) on it with "Ideal polyhedra .xyz" checked.')
        return False

    if as_bool(state):
        _ticked.add(stem)
    else:
        _ticked.discard(stem)
    refresh()
    return True


def set_ideal_run_shape(stem, shape):
    """Picks the shape shown and drawn for one run, then redraws."""
    _run_shapes[stem] = shape
    refresh()


def delete_run(stem):
    """Deletes one CShM run's files after asking, and redraws the overlay without it."""
    run = next((r for r in find_runs() if r['stem'] == stem), None)
    if run is None:
        return False

    folder = _workdir()
    paths = [p for p in (run['ideal'], run['table']) if p]
    # The input .xyz too, unless a CSoM or ODis output of the same centre is still next to it.
    xyz = os.path.join(folder, f'{stem}.xyz')
    others = [n for n in os.listdir(folder)
              if n.startswith(f'{stem}_') and os.path.join(folder, n) not in paths]
    if os.path.isfile(xyz) and not others:
        paths.append(xyz)

    files = '\n'.join(os.path.basename(p) for p in paths)
    if olx.Alert('Please confirm', f'Delete the CShM run of {run["title"]}?\n\n{files}', 'YCQ') != 'Y':
        return False

    for path in paths:
        try:
            os.remove(path)
        except OSError as e:
            print(f'Could not delete {os.path.basename(path)}: {e}')
    registry = _read_registry()
    if registry.pop(stem, None) is not None:
        _write_registry(registry)
    _ticked.discard(stem)
    _run_shapes.pop(stem, None)
    refresh()
    return True


def cshm_runs_html():
    """The CShM runs table: per run a checkbox to draw it, its centre, a shape picker, that
    shape's CShM and a Delete button."""
    runs = find_runs()
    if not runs:
        return ('<i>No CShM results saved for this structure yet: run CShM with '
                '"Results .csv" or "Ideal polyhedra .xyz" checked.</i>')

    bg = olx.GetVar('HtmlTableBgColour', '#ffffff')
    combo_height = olx.GetVar('HtmlComboHeight', '20')
    button_height = olx.GetVar('HtmlInputHeight', '20')
    columns = [('Draw', 10), ('Centre', 20), ('Shape', 35), ('CShM', 15), ('', 20)]
    rows = [''.join(f'<td align="center" width="{width}%"><b>{title}</b></td>'
                    for title, width in columns)]
    for run in runs:
        stem = run['stem']
        key = re.sub(r'[^A-Za-z0-9_-]', '_', stem)
        # Unticked without an _ideal.xyz: ticking it then only prints how to get one.
        checked = 'true' if ticked(run) else 'false'
        cells = [
            f'<input type="checkbox" name="CShMDraw_{key}" checked="{checked}" '
            f'onclick="spy.Timeo.set_ideal_run_visible(\'{stem}\', html.GetState(\'~name~\'))">',
            f'<b>{run["title"]}</b>',
        ]

        results = read_results(run)
        if results:
            picked = _pick_frame(results, _run_shapes.get(stem))
            # Single-quoted inside, as in Olex2's own combos.
            cells.append(
                f'<input type="combo" name="CShMShape_{key}" width="100%" height="{combo_height}" '
                f'items="\'{";".join(r["shape"] for r in results)}\'" value="{picked["shape"]}" '
                f'readonly="true" '
                f'onchange="spy.Timeo.set_ideal_run_shape(\'{stem}\', html.GetValue(\'~name~\'))">')
            cells.append(f'{picked["cshm"]:.3f}')
        else:
            cells += ['<i>unreadable</i>', '']

        cells.append(
            f'<input type="button" name="CShMDelete_{key}" value="Delete" width="100%" '
            f'height="{button_height}" bgcolor="GetVar(linkButton.bgcolor)" '
            f'fgcolor="GetVar(linkButton.fgcolor)" flat="GetVar(linkButton.flat)" '
            f'custom="GetVar(custom_button)" onclick="spy.Timeo.delete_cshm_run(\'{stem}\')">')
        rows.append(''.join(f'<td align="center">{cell}</td>' for cell in cells))

    return ('<table bgcolor="#ababab" border="0" cellspacing="1" cellpadding="2" width="100%">'
            + ''.join(f'<tr bgcolor="{bg}">{row}</tr>' for row in rows) + '</table>')


def clear_ideal():
    """Hides every ideal polyhedron drawn so far."""
    for name in _drawn:
        olex.m(f'kill {name}')
    _drawn.clear()
    _faces.clear()
    OV.Refresh()
