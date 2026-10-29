# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 José Serrano Guarinos

__version__ = '1.3.1'

import os
import re

import olex
import olx
from olexFunctions import OlexFunctions

OV = OlexFunctions()

from reload_all import reload_all

from autoshape import (ShapeCalculation, ShapeError, can_find_shape_msg, locate_shape,
                       print_shape_table, run_shape)
from constants import octadist_citation, shape21_citation
import cosmochlore
import csom_runs
import ideal_overlay
from helper_functions import (get_neighbours, get_neighbours_on_sel, get_selected_atoms,
                              get_xyz_sel, print_console_bs, print_orm, test_selection_class, as_bool)
from octahedral_distortion import CalcDistortion
from selection import AtomSelection, split_by_parts
import update_check

debug = bool(OV.GetParam("olex2.debug", False))

instance_path = OV.DataDir()

try:
    from_outside = False
    p_path = os.path.dirname(os.path.abspath(__file__))
except:
    from_outside = True
    p_path = os.path.dirname(os.path.abspath("__file__"))

l = open(os.sep.join([p_path, 'def.txt'])).readlines()
d = {}
for line in l:
    line = line.strip()
    if not line or line.startswith("#"):
        continue
    d[line.split("=")[0].strip()] = line.split("=")[1].strip()

p_name = d['p_name']
p_htm = d['p_htm']
p_img = eval(d['p_img'])
p_scope = d['p_scope']

OV.SetVar('Timeo_plugin_path', p_path)

from PluginTools import PluginTools as PT


def _exclude_h():
    return as_bool(OV.GetParam('timeo.exclude_h', False))


def _split_selection(selection, centered=True):
    """Splits `selection` into structures using the timeo.part_codes param."""
    part_codes = OV.GetParam('timeo.part_codes', '') or ''
    return split_by_parts(selection, part_codes, centered)


def _prepare_structures(sel_string, merge=None):
    """Builds one MolecularStructure per disorder component from `sel_string`.

    A single selected atom is grown into a centered coordination polyhedron
    (its neighbours are added, and pi-bonded ligands are merged into a
    centroid if `merge`, or failing to the merge_ligands phil param.

    Multiple selected atoms are taken as-is as a non-centered shape.

    Returns (structures, centered), or (None, None) if the selection is empty.
    """
    if sel_string == '':
        print('Invalid atom selection: no atoms selected.')
        return None, None

    selection = AtomSelection(sel_string)

    if len(selection) > 1:
        # Multiple atoms selected: treat them as the vertices of a non-centered shape.
        selection.remove_duplicates()
        centered = False
    else:
        # A single atom selected: grow it into a coordination polyhedron.
        selection.add_neighbours(_exclude_h())
        if merge is None:
            merge = as_bool(OV.GetParam('timeo.merge_ligands', False))
        if merge:
            selection.merge_ligands()
        centered = True

    # Returns one structure per disorder component, or a single structure when the strucutre is not disordered
    return _split_selection(selection, centered), centered


def _cosmochlore_exe_path():
    return OV.GetParam('timeo.cosmochlore.exe_path', '') or None


def check_cosmochlore(silent=True):
    return cosmochlore.can_find_cosmochlore_msg(_cosmochlore_exe_path(), silent)


def save_params():
    """Saves the changed timeo params to <DataDir>/timeo.phil."""
    user_phil_file = os.path.join(OV.DataDir(), f'{p_scope}.phil')
    olx.phil_handler.save_param_file(file_name=user_phil_file, scope_name=p_scope, diff_only=True)


def SetOption(param, value):
    """Sets a timeo param and saves it so it survives a restart."""
    OV.SetParam(param, value)
    save_params()


def _set_exe_param(param, path):
    path = str(path).strip().strip('"\'').replace('\\', '/')
    SetOption(param, path)
    if path and not os.path.exists(path):
        print(f'Warning: "{path}" does not exist. The path was saved anyway.')


def SetCosmochloreExe(path=''):
    """Sets and saves timeo.cosmochlore.exe_path; empty falls back to PATH."""
    _set_exe_param('timeo.cosmochlore.exe_path', path)


def SetShapeExe(path=''):
    """Sets and saves timeo.shape.exe_path; empty falls back to PATH."""
    _set_exe_param('timeo.shape.exe_path', path)


def _shape_exe_path():
    return OV.GetParam('timeo.shape.exe_path', '') or None


def check_shape(silent=True):
    return can_find_shape_msg(_shape_exe_path(), silent)


def _cosmochlore_workdir():
    """The scratch folder cosmochlore .xyz inputs and outputs are written to."""
    path = os.path.join(olx.FilePath(), 'cosmochlore')
    os.makedirs(path, exist_ok=True)
    return path


def _user_shapes_dir():
    """Where autoCSHM()'s -r/--ref candidates are read from: a `user_shapes`
    folder living next to the plugin's own code, alongside SHAPE's `custom_shapes`."""
    return os.path.join(p_path, 'user_shapes')


def _safe_key(name: str) -> str:
    """Sanitises a shape's display name into a safe HTML control name / phil
    param suffix - the display name itself (a .yaml filename the user chose)
    isn't guaranteed to be free of spaces, quotes, or other characters that
    would break the generated markup or a dotted phil path."""
    return re.sub(r'[^A-Za-z0-9_-]', '_', name)


def _user_shape_param(name: str) -> str:
    return f'timeo.cosmochlore.cshm.user_shape.{_safe_key(name)}'


def user_shapes_checkboxes_html():
    """Returns one checkbox per .yaml file found in the user_shapes folder, for
    picking which ones autoCSHM() passes to cshm's -r/--ref. Selection state is
    stored per file (see _user_shape_param) rather than as one phil-declared
    list, since the set of files - and so the set of possible param names - is
    only known at runtime; this means it is not phil-persisted across Olex2
    restarts, only for the current session.

    Mirrors the table/tr/td structure gui/snippets/input-checkbox-td expands
    to: the checkbox's own `label` attribute is not what renders visible text
    in Olex2's control - the label is a separate <td><b>...</b></td> cell next
    to it.
    """
    shapes = cosmochlore.list_user_shapes(_user_shapes_dir())
    if not shapes:
        return '<i>No .yaml files found.</i>'

    rows = []
    for name, _ in shapes:
        param = _user_shape_param(name)
        checked = 'true' if as_bool(OV.GetParam(param, False)) else 'false'
        label = name.replace('"', "'")
        rows.append(
            f'<tr><td><input type="checkbox" name="UserShape_{_safe_key(name)}" '
            f'checked="{checked}" '
            f'onclick="spy.SetParam(\'{param}\', html.GetState(\'~name~\'))">'
            f'</td><td align="left"><b>{label}</b></td></tr>'
        )
    return '<table cellpadding="0" cellspacing="0">' + ''.join(rows) + '</table>'


def _selected_user_shapes():
    """Full paths of the user_shapes .yaml files currently checked in the GUI."""
    return [path for name, path in cosmochlore.list_user_shapes(_user_shapes_dir())
           if as_bool(OV.GetParam(_user_shape_param(name), False))]


def open_user_shapes_folder():
    """Opens the user_shapes folder in the OS file browser. Creates it first
    if it doesn't exist yet, so the button always has somewhere to open."""
    path = _user_shapes_dir()
    os.makedirs(path, exist_ok=True)
    olx.Shell(path)


def reload_plugin():
    """Re-imports every plugin module from wherever the plugin is installed,
    so the Reload ALL button doesn't depend on a hardcoded install path."""
    reload_all()


def open_plugin_folder():
    """Opens the plugin's own folder in the OS file browser."""
    olx.Shell(p_path)


# MAIN LOGIC FUNCTIONS.
def autoSHAPE():
    print('\n' + '-' * 50)
    print('Simple continuous Shape Analysis Using autoSHAPE')
    try:
        shape_exe = locate_shape(_shape_exe_path())
    except ShapeError as e:
        print(e)
        return False

    structures, centered = _prepare_structures(olex.f('sel()'))
    if structures is None:
        return False

    for i, structure in enumerate(structures):
        try:
            shape_measurement = ShapeCalculation(structure, f'{olx.FileName()}_{i}',
                                                 centered, ['%fullout'])
        except ValueError as e:
            print(f'Skipping part {i}: {e}')
            continue

        folder = shape_measurement.write_tab(olx.FilePath())
        if folder is None:
            continue

        for f in run_shape(folder, shape_exe):
            print_shape_table(os.path.join(folder, f'{f}.tab'))

    print(shape21_citation)
    return True


def autoOCTADIST():
    # Get the selected atoms.
    sel_string = olex.f('sel()')
    if sel_string == '':
        print('Invalid atom selection: no atoms selected.')
        return False

    selection = AtomSelection(sel_string)

    # Exit if the selection is anything other than the single central atom.
    if len(selection) != 1:
        print(f'Invalid atom selection: expected 1 atom, found {len(selection)}.')
        return False

    # Add coordinated atoms to the current selection.
    selection.add_neighbours(_exclude_h())

    for structure in _split_selection(selection):
        try:
            calculation = CalcDistortion(structure)
        except ValueError as e:
            print(f'Invalid polyhedra for {structure.labels[0]}: {e}')
            continue

        calculation.print_results(olx.FileName())
        calculation.draw_octahedron()

    # Citation
    print('\nThis calculations were made using a reimplementation of the OctaDist '
          'algorithm by David J. Harding et al.')
    print(octadist_citation)
    return True


_COSMOCHLORE_CREDIT = 'Computed using cosmochlore (https://github.com/Yluro/cosmochlore).'


def autoCSHM(shapes=None, user_shapes=None, table=None, ideal=None):
    """Continuous Shape Measures via cosmochlore, on the current selection.

    `shapes`: built-in reference-shape indices to restrict to (None = all
    applicable for the detected vertex count). `user_shapes`: paths to
    user-defined shape .yaml files; when not given, falls back to whichever
    files are checked in the user_shapes_checkboxes_html() GUI list. Passed to
    cshm's -r/--ref exactly as given - no vertex-count pre-filtering, cosmochlore
    reports a mismatch itself (aborting that part's run) rather than this
    having any fallback logic of its own to hide the error.
    `table`/`ideal`: also write the corresponding cosmochlore output files
    next to the .xyz; default to the matching phil params when not given.
    """
    print('\n' + '-' * 50)
    print('Continuous Shape Measures using cosmochlore')
    try:
        exe = cosmochlore.locate_cosmochlore(_cosmochlore_exe_path())
    except cosmochlore.CosmochloreError as e:
        print(e)
        return False

    if user_shapes is None:
        user_shapes = _selected_user_shapes() or None
    if table is None:
        table = as_bool(OV.GetParam('timeo.cosmochlore.cshm.table', False))
    if ideal is None:
        ideal = as_bool(OV.GetParam('timeo.cosmochlore.cshm.ideal', False))

    structures, centered = _prepare_structures(olex.f('sel()'))
    if structures is None:
        return False

    workdir = _cosmochlore_workdir()
    ran_any = False
    ideal_paths = []
    for i, structure in enumerate(structures):
        xyz_path = os.path.join(workdir, f'{olx.FileName()}_{structure.labels[0]}_{i}.xyz')
        cosmochlore.write_xyz(structure, xyz_path, comment=f'{olx.FileName()} part {i}')

        try:
            args = cosmochlore.build_cshm_args(xyz_path, centered, shapes, user_shapes, table,
                                               ideal)
            print(cosmochlore.trim_banner(cosmochlore.run_cosmochlore(exe, args, cwd=workdir)))
            ran_any = True
        except cosmochlore.CosmochloreError as e:
            print(f'cshm failed for part {i}: {e}')
            continue

        # cosmochlore names it after the input: <stem>_ideal.xyz, next to the .xyz.
        ideal_path = os.path.splitext(xyz_path)[0] + '_ideal.xyz'
        if ideal and os.path.isfile(ideal_path):
            ideal_paths.append(ideal_path)

    if ideal_paths:
        ideal_overlay.record_run(ideal_paths, centered)
        # Redraws ticked runs from their new files.
        ideal_overlay.refresh()
    if ran_any:
        # Refreshes the CShM runs table.
        olex.m('html.Update()')
        print(_COSMOCHLORE_CREDIT)
    return ran_any


def SetOverlayOption(param, value):
    """Sets and saves a timeo.overlay drawing option and applies it to the drawn overlay."""
    SetOption(param, value)
    ideal_overlay.option_changed(param)


# Registered instead of the ideal_overlay functions themselves: reload_all() reloads
# this module first, so those would keep running the code from before a Reload ALL.
def draw_ideal():
    """Draws the chosen ideal shape of every ticked CShM run."""
    return ideal_overlay.draw_ideal()


def clear_ideal():
    """Hides every drawn ideal polyhedron."""
    ideal_overlay.clear_ideal()


def set_ideal_run_visible(stem, state):
    """Runs table Draw checkbox: draws or hides one CShM run."""
    ideal_overlay.set_ideal_run_visible(stem, state)
    # Shows or hides the appearance settings, and unticks a run that has nothing to draw.
    olex.m('html.Update()')


def set_ideal_run_shape(stem, shape):
    """Runs table combo: the shape shown and drawn for one CShM run."""
    ideal_overlay.set_ideal_run_shape(stem, shape)
    # Shows the new shape's CShM.
    olex.m('html.Update()')


def delete_cshm_run(stem):
    """Runs table Delete button: deletes one CShM run's files after asking."""
    if ideal_overlay.delete_run(stem):
        olex.m('html.Update()')


def cshm_runs_html():
    return ideal_overlay.cshm_runs_html()


def ideal_settings_shown():
    """Whether any CShM run is ticked to draw, for the appearance settings' #includeif."""
    return ideal_overlay.any_ticked()


def autoCSOM(point_groups=None, mode=None, vector=None, full=None, table=None,
            operated=None, seeds=None, iterations=None, tolerance=None,
            ignore_labels=None):
    """Continuous Symmetry Operation Measures via cosmochlore, on the current selection.

    point_groups: a space-separated string or list of Schoenflies, falls back to the
    timeo.cosmochlore.csom.point_groups phil param when not given.

    mode: centering mode (auto/first/centroid/manual), defaults to the matching phil param.

    `vector`: required 3-value centering vector when mode is 'manual'.

    `full`/`table`/`operated`/`ignore_labels`: default to the matching phil params.

    `seeds`/`iterations`/`tolerance`: the axis-search knobs (cosmochlore's
    -s/--seeds, -i/--iterations and -T/--tolerance), defaulting to the matching phil
    params. `tolerance` is the value itself (e.g. 1e-7); the phil param behind
    it - and the GUI's spin box - hold the exponent N of a 1e-N tolerance.
    """
    print('\n' + '-' * 50)
    print('Continuous Symmetry Operation Measures using cosmochlore')
    try:
        exe = cosmochlore.locate_cosmochlore(_cosmochlore_exe_path())
    except cosmochlore.CosmochloreError as e:
        print(e)
        return False

    if point_groups is None:
        point_groups = OV.GetParam('timeo.cosmochlore.csom.point_groups', '')
    if isinstance(point_groups, str):
        point_groups = point_groups.split()

    if not point_groups:
        print('No point groups given. Pass e.g. '
              "spy.Timeo.autoCSOM('Oh D4h D3d'), or set "
              "timeo.cosmochlore.csom.point_groups.")
        return False

    if mode is None:
        mode = OV.GetParam('timeo.cosmochlore.csom.mode', 'auto')
    # Olex2's combo control capitalises the value it hands back (e.g. 'auto'
    # -> 'Auto') regardless of the case used in the combo's own item list, so
    # normalise here rather than trust whatever case arrives from the GUI or a
    # console caller.
    mode = str(mode).strip().lower()

    if vector is None:
        vector_str = OV.GetParam('timeo.cosmochlore.csom.vector', '')
        if vector_str.strip():
            try:
                vector = [float(v) for v in vector_str.split()]
            except ValueError:
                print(f'Could not parse the manual centering vector "{vector_str}": '
                      f'expected three numbers separated by spaces, e.g. "0.0 0.0 0.0".')
                return False
            if len(vector) != 3:
                print(f'The manual centering vector needs exactly 3 numbers, '
                      f'found {len(vector)} in "{vector_str}".')
                return False

    if mode == 'manual' and not vector:
        print('Centering mode is "manual" but no centering vector was given. Set it in '
              'the Cosmochlore section (x y z, space-separated), or pass '
              "spy.Timeo.autoCSOM(vector=[x, y, z]).")
        return False

    if full is None:
        full = as_bool(OV.GetParam('timeo.cosmochlore.csom.full', False))
    if table is None:
        table = as_bool(OV.GetParam('timeo.cosmochlore.csom.table', False))
    if operated is None:
        operated = as_bool(OV.GetParam('timeo.cosmochlore.csom.operated', False))
    if ignore_labels is None:
        ignore_labels = as_bool(OV.GetParam('timeo.cosmochlore.csom.ignore_labels', False))

    if seeds is None:
        seeds = OV.GetParam('timeo.cosmochlore.csom.seeds', 20)
    if iterations is None:
        iterations = OV.GetParam('timeo.cosmochlore.csom.iterations', 200)
    if tolerance is None:
        tolerance = 10.0 ** -OV.GetParam('timeo.cosmochlore.csom.tolerance_exponent', 6)

    structures, centered = _prepare_structures(olex.f('sel()'))
    if structures is None:
        return False

    print('Running cosmochlore, this may take a few seconds...')

    workdir = _cosmochlore_workdir()
    ran_any = False
    for i, structure in enumerate(structures):
        xyz_path = os.path.join(workdir, f'{olx.FileName()}_{structure.labels[0]}_{i}.xyz')
        cosmochlore.write_xyz(structure, xyz_path, comment=f'{olx.FileName()} part {i}')

        try:
            args = cosmochlore.build_csom_args(xyz_path, centered, point_groups, mode, vector, full,
                                               table, operated, seeds, iterations, tolerance,
                                               ignore_labels)
            print(cosmochlore.trim_banner(cosmochlore.run_cosmochlore(exe, args, cwd=workdir)))
            ran_any = True
        except cosmochlore.CosmochloreError as e:
            print(f'csom failed for part {i}: {e}')

    if ran_any:
        # Refreshes the CSoM runs table.
        olex.m('html.Update()')
        print(_COSMOCHLORE_CREDIT)
    return ran_any


# Registered instead of the csom_runs functions themselves, for the same reason as above.
def set_csom_run_point_group(stem, point_group):
    """CSoM runs table combo: the point group shown for one run."""
    csom_runs.set_point_group(stem, point_group)
    # Shows the new point group's CSoM.
    olex.m('html.Update()')


def delete_csom_run(stem):
    """CSoM runs table Delete button: deletes one CSoM run's files after asking."""
    if csom_runs.delete_run(stem):
        olex.m('html.Update()')


def csom_runs_html():
    return csom_runs.csom_runs_html()


def autoODIS(full=None, table=None):
    """Octahedral distortion analysis via cosmochlore, on the current selection.

    Requires a single selected atom with exactly six neighbours (a 7-atom
    polyhedron once its neighbours are added). `full`/`table` default to the
    matching phil params when not given.
    """
    print('\n' + '-' * 50)
    print('Octahedral distortion analysis using cosmochlore')
    try:
        exe = cosmochlore.locate_cosmochlore(_cosmochlore_exe_path())
    except cosmochlore.CosmochloreError as e:
        print(e)
        return False

    if full is None:
        full = as_bool(OV.GetParam('timeo.cosmochlore.odis.full', False))
    if table is None:
        table = as_bool(OV.GetParam('timeo.cosmochlore.odis.table', False))

    sel_string = olex.f('sel()')
    if sel_string == '':
        print('Invalid atom selection: no atoms selected.')
        return False

    selection = AtomSelection(sel_string)
    if len(selection) != 1:
        print(f'Invalid atom selection: expected 1 atom, found {len(selection)}.')
        return False

    selection.add_neighbours(_exclude_h())

    workdir = _cosmochlore_workdir()
    ran_any = False
    for i, structure in enumerate(_split_selection(selection)):
        if len(structure) != 7:
            print(f'Skipping part {i}: expected 7 atoms (centre + six donors), '
                  f'found {len(structure)}.')
            continue

        xyz_path = os.path.join(workdir, f'{olx.FileName()}_{structure.labels[0]}_{i}.xyz')
        cosmochlore.write_xyz(structure, xyz_path, comment=f'{olx.FileName()} part {i}')

        try:
            args = cosmochlore.build_odis_args(xyz_path, full, table)
            print(cosmochlore.trim_banner(cosmochlore.run_cosmochlore(exe, args, cwd=workdir)))
            ran_any = True
        except cosmochlore.CosmochloreError as e:
            print(f'odis failed for part {i}: {e}')

    if ran_any:
        print(_COSMOCHLORE_CREDIT)
    return ran_any


def shape_status_html():
    try:
        exe = locate_shape(_shape_exe_path())
    except ShapeError as e:
        return f"<font color='{OV.GetParam('gui.grey')}'>{e}</font>"

    return f"<font color='{OV.GetParam('gui.green')}'>SHAPE executable found at: {exe}</font>"


def cosmochlore_status_html():
    try:
        exe = cosmochlore.locate_cosmochlore(_cosmochlore_exe_path())
    except cosmochlore.CosmochloreError as e:
        return f"<font color='{OV.GetParam('gui.grey')}'>{e}</font>"

    version = '.'.join(str(v) for v in cosmochlore.get_version(exe))
    text = f'cosmochlore {version} found at: {exe}'
    return f"<font color='{OV.GetParam('gui.green')}'>{text}</font>"


def _update_cache_file():
    return os.path.join(OV.DataDir(), 'timeo_update.json')


def _update_check_enabled():
    return as_bool(OV.GetParam('timeo.check_for_updates', True))


def _refresh_update_cache(verbose):
    # The check runs in a thread, so Schedule hands the result back to Olex2's main thread,
    # as threads_imp.py does for the news image.
    update_check.refresh_in_background(
        _update_cache_file(),
        on_done=lambda: olx.Schedule(1, f'spy.Timeo.announce_update({verbose})'))


def _print_update_notice(verbose=False):
    cache = update_check.read_cache(_update_cache_file())
    newer = update_check.newer_release(__version__, cache)
    if newer:
        tag, url = newer
        print(f'Timeo {tag} is available (you have {__version__}). Get it from: {url}')
    elif verbose and cache.get('ok'):
        print(f'Timeo {__version__} is up to date.')
    elif verbose:
        print('Could not reach GitHub to check for Timeo updates.')


def start_update_check():
    """Checks GitHub for a newer Timeo release, at most once a day and once per session."""
    # Timeo.py runs twice on load (reload_all re-imports it), so a module global would not do.
    if OV.GetVar('timeo_update_check_started', '') == 'true':
        return
    OV.SetVar('timeo_update_check_started', 'true')
    if not _update_check_enabled() or not OV.canNetwork(show_msg=False):
        return
    if update_check.is_due(update_check.read_cache(_update_cache_file())):
        _refresh_update_cache(verbose=False)
    else:
        _print_update_notice()


def announce_update(verbose=False):
    """Prints the result of a finished update check and redraws the panel."""
    _print_update_notice(as_bool(verbose))
    olex.m('html.Update()')


def CheckForUpdates():
    """Asks GitHub for the latest Timeo release now, ignoring the daily cache."""
    if not OV.canNetwork():
        return
    print('Checking GitHub for Timeo updates...')
    _refresh_update_cache(verbose=True)


def update_status_html():
    grey = OV.GetParam('gui.grey')
    if not _update_check_enabled():
        return f"<font color='{grey}'>Timeo {__version__}. The update check is off.</font>"

    cache = update_check.read_cache(_update_cache_file())
    newer = update_check.newer_release(__version__, cache)
    if newer:
        tag, url = newer
        return (f"<font color='{OV.GetParam('gui.orange')}'>Timeo {tag} is available "
                f"(you have {__version__}): <a href='shell {url}'>download</a></font>")
    if 'tag' not in cache:
        return f"<font color='{grey}'>Timeo {__version__}. Not checked for updates yet.</font>"
    return f"<font color='{OV.GetParam('gui.green')}'>Timeo {__version__} is up to date.</font>"


class Timeo(PT):
    def __init__(self):
        super(Timeo, self).__init__()
        self.p_name = p_name
        self.p_path = p_path
        self.p_scope = p_scope
        self.p_htm = p_htm
        self.p_img = p_img
        self.deal_with_phil(operation='read')
        self.print_version_date()
        if not from_outside:
            self.setup_gui()

        # Main entry points.
        OV.registerFunction(autoSHAPE, True, "Timeo")
        OV.registerFunction(autoOCTADIST, True, "Timeo")
        OV.registerFunction(check_shape, True, "Timeo")
        OV.registerFunction(SetShapeExe, True, "Timeo")
        OV.registerFunction(SetOption, True, "Timeo")
        OV.registerFunction(shape_status_html, False, 'Timeo')

        # cosmochlore entry points.
        OV.registerFunction(autoCSHM, True, "Timeo")
        OV.registerFunction(autoCSOM, True, "Timeo")
        OV.registerFunction(autoODIS, True, "Timeo")
        OV.registerFunction(check_cosmochlore, True, "Timeo")
        OV.registerFunction(cosmochlore_status_html, False, "Timeo")
        OV.registerFunction(SetCosmochloreExe, True, "Timeo")
        OV.registerFunction(save_params, True, "Timeo")
        OV.registerFunction(user_shapes_checkboxes_html, False, "Timeo")
        OV.registerFunction(open_user_shapes_folder, True, "Timeo")

        # Ideal polyhedra overlay in the 3D view.
        OV.registerFunction(draw_ideal, True, "Timeo")
        OV.registerFunction(clear_ideal, True, "Timeo")
        OV.registerFunction(set_ideal_run_visible, True, "Timeo")
        OV.registerFunction(set_ideal_run_shape, True, "Timeo")
        OV.registerFunction(delete_cshm_run, True, "Timeo")
        OV.registerFunction(cshm_runs_html, False, "Timeo")
        OV.registerFunction(ideal_settings_shown, False, "Timeo")
        OV.registerFunction(SetOverlayOption, True, "Timeo")

        # CSoM runs table.
        OV.registerFunction(set_csom_run_point_group, True, "Timeo")
        OV.registerFunction(delete_csom_run, True, "Timeo")
        OV.registerFunction(csom_runs_html, False, "Timeo")

        # Options panel.
        OV.registerFunction(reload_plugin, True, "Timeo")
        OV.registerFunction(open_plugin_folder, True, "Timeo")

        # Update check.
        OV.registerFunction(CheckForUpdates, True, "Timeo")
        OV.registerFunction(announce_update, True, "Timeo")
        OV.registerFunction(update_status_html, False, "Timeo")

        # Debug panel helpers.
        OV.registerFunction(get_selected_atoms, True, "Timeo")
        OV.registerFunction(get_neighbours, True, "Timeo")
        OV.registerFunction(get_neighbours_on_sel, True, "Timeo")
        OV.registerFunction(get_xyz_sel, True, "Timeo")
        OV.registerFunction(print_console_bs, False, 'Timeo')
        OV.registerFunction(print_orm, False, 'Timeo')
        OV.registerFunction(test_selection_class, False, 'Timeo')

        # After registering, so announce_update exists when the check finishes.
        start_update_check()
    # END Generated =======================================


Timeo_instance = Timeo()
print("Loading Timeo modules.")
reload_all()
print(f"Timeo {__version__} by JSG loaded.")
