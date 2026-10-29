# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 José Serrano Guarinos

"""The CSoM runs table: every csom output in <FilePath>/cosmochlore/ for the current structure.

A run is one measured centre and disorder component: its _csom_table.csv and the per point
group _details.csv, _operated.xyz and _merged.mol2 files cosmochlore wrote next to its .xyz.
"""
import csv
import os
import re

import olx

_RUN_FILE = re.compile(r'^(?P<label>.+)_(?P<part>\d+)(?P<kind>_csom_table\.csv'
                       r'|_(?P<pg>[A-Za-z0-9]+)_(?:details\.csv|operated\.xyz|merged\.mol2))$')

# Point group picked for each run when it isn't the best match, this session only.
_run_point_groups = {}


def _workdir():
    return os.path.join(olx.FilePath(), 'cosmochlore')


def find_runs():
    """Every CSoM run on the current structure, as dicts sorted by stem.

    'table' holds the run's _csom_table.csv path or None, 'files' all its output files and
    'point_groups' the point groups its per point group files are named after.
    """
    folder = _workdir()
    if not os.path.isdir(folder):
        return []

    prefix = f'{olx.FileName()}_'
    runs = {}
    for name in os.listdir(folder):
        match = _RUN_FILE.match(name[len(prefix):]) if name.startswith(prefix) else None
        if match is None:
            continue
        stem = name[:-len(match.group('kind'))]
        run = runs.setdefault(stem, {
            'stem': stem,
            'table': None,
            'files': [],
            'point_groups': [],
            'label': match.group('label'),
            'part': int(match.group('part')),
        })
        path = os.path.join(folder, name)
        run['files'].append(path)
        if match.group('pg') is None:
            run['table'] = path
        elif match.group('pg') not in run['point_groups']:
            run['point_groups'].append(match.group('pg'))
    runs = sorted(runs.values(), key=lambda r: r['stem'])

    # 'Fe1', or 'Fe1 (2)' when that centre was split into several disorder components.
    for run in runs:
        split = sum(r['label'] == run['label'] for r in runs) > 1
        run['title'] = f'{run["label"]} ({run["part"] + 1})' if split else run['label']
    return runs


def read_results(run):
    """The run's point groups as dicts with point_group and csom, best first, from its .csv."""
    if not run['table']:
        return []
    try:
        with open(run['table'], newline='') as f:
            results = [{'point_group': row['PointGroup'], 'csom': float(row['Dev'])}
                       for row in csv.DictReader(f)]
    except (OSError, ValueError, KeyError, TypeError):
        return []
    return sorted(results, key=lambda r: r['csom'])


def _pick(results, choice):
    """The result for point group `choice`, or the lowest-CSoM one for 'best' or one not in `results`."""
    # Olex2's combo can hand the value back with its case changed (e.g. 'Best').
    wanted = str(choice or 'best').strip().lower()
    match = next((r for r in results if r['point_group'].lower() == wanted), None)
    return match or results[0]


def set_point_group(stem, point_group):
    """Picks the point group shown for one run."""
    _run_point_groups[stem] = point_group


def delete_run(stem):
    """Deletes one CSoM run's files after asking."""
    run = next((r for r in find_runs() if r['stem'] == stem), None)
    if run is None:
        return False

    folder = _workdir()
    paths = sorted(run['files'])
    # The input .xyz too, unless a CShM or ODis output of the same centre is still next to it.
    xyz = os.path.join(folder, f'{stem}.xyz')
    others = [n for n in os.listdir(folder)
              if n.startswith(f'{stem}_') and os.path.join(folder, n) not in paths]
    if os.path.isfile(xyz) and not others:
        paths.append(xyz)

    files = '\n'.join(os.path.basename(p) for p in paths)
    if olx.Alert('Please confirm', f'Delete the CSoM run of {run["title"]}?\n\n{files}', 'YCQ') != 'Y':
        return False

    for path in paths:
        try:
            os.remove(path)
        except OSError as e:
            print(f'Could not delete {os.path.basename(path)}: {e}')
    _run_point_groups.pop(stem, None)
    return True


def csom_runs_html():
    """The CSoM runs table: per run its centre, a point group picker, that point group's CSoM
    and a Delete button."""
    runs = find_runs()
    if not runs:
        return ('<i>No CSoM results saved for this structure yet: run CSoM with '
                '"Summary .csv" checked.</i>')

    bg = olx.GetVar('HtmlTableBgColour', '#ffffff')
    combo_height = olx.GetVar('HtmlComboHeight', '20')
    button_height = olx.GetVar('HtmlInputHeight', '20')
    columns = [('Centre', 25), ('Point group', 40), ('CSoM', 15), ('', 20)]
    rows = [''.join(f'<td align="center" width="{width}%"><b>{title}</b></td>'
                    for title, width in columns)]
    for run in runs:
        stem = run['stem']
        key = re.sub(r'[^A-Za-z0-9_-]', '_', stem)
        cells = [f'<b>{run["title"]}</b>']

        results = read_results(run)
        if results:
            picked = _pick(results, _run_point_groups.get(stem))
            # Single-quoted inside, as in Olex2's own combos.
            cells.append(
                f'<input type="combo" name="CSoMPointGroup_{key}" width="100%" '
                f'height="{combo_height}" '
                f'items="\'{";".join(r["point_group"] for r in results)}\'" '
                f'value="{picked["point_group"]}" readonly="true" '
                f'onchange="spy.Timeo.set_csom_run_point_group(\'{stem}\', html.GetValue(\'~name~\'))">')
            cells.append(f'{picked["csom"]:.3f}')
        elif run['table']:
            cells += ['<i>unreadable</i>', '']
        else:
            # Only per point group files: their names give the point groups, not the CSoM.
            cells += [' '.join(run['point_groups']), '<i>no .csv</i>']

        cells.append(
            f'<input type="button" name="CSoMDelete_{key}" value="Delete" width="100%" '
            f'height="{button_height}" bgcolor="GetVar(linkButton.bgcolor)" '
            f'fgcolor="GetVar(linkButton.fgcolor)" flat="GetVar(linkButton.flat)" '
            f'custom="GetVar(custom_button)" onclick="spy.Timeo.delete_csom_run(\'{stem}\')">')
        rows.append(''.join(f'<td align="center">{cell}</td>' for cell in cells))

    return ('<table bgcolor="#ababab" border="0" cellspacing="1" cellpadding="2" width="100%">'
            + ''.join(f'<tr bgcolor="{bg}">{row}</tr>' for row in rows) + '</table>')
