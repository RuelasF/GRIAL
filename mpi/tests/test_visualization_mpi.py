#!/usr/bin/env python3
"""C04 actual writer, independent reconstruction, VTK consumer and segment controls."""
import argparse
from collections import Counter
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
from parallel_vtk import FIELDS, collect_frames, read_frame, read_piece
from test_observer_mpi import cases, compare, receipt


def compare_fields(a, b):
    maximum = 0.0
    for field in FIELDS:
        x, y = a['fields'][field], b['fields'][field]
        assert len(x) == len(y)
        # Legacy and local loops can select different scalar/vector exp kernels
        # for log-radius visualization. The native solution is checked separately.
        if field not in ['Velocity_X', 'Shift_beta_r']:
            assert x.tobytes() == y.tobytes(), field
        error = max(abs(u-v)/max(1.0, abs(u), abs(v)) for u, v in zip(x, y))
        assert error < 1e-11, (field, error)
        maximum = max(maximum, error)
    return maximum


def geometry(directory):
    cells = {}
    for path in sorted(directory.glob('piece_*.vtu')):
        p = read_piece(path)
        a, start = p['arrays'], 0
        for gid, end in zip(a['GlobalCellId'], a['offsets']):
            cells[gid] = sorted(tuple(a['Points'][3*v:3*v+3]) for v in a['connectivity'][start:end])
            start = end
    return cells


def vtk_check(directory, required):
    try:
        import vtk
    except ImportError:
        assert not required, 'VTK consumer requested but not installed'
        return dict(checked=False)
    from vtk.util.numpy_support import vtk_to_numpy
    reader = vtk.vtkXMLPUnstructuredGridReader()
    errors = []
    reader.AddObserver('ErrorEvent', lambda *args: errors.append(str(args)))
    reader.SetFileName(str(directory / 'frame.pvtu'))
    reader.Update()
    mesh = reader.GetOutput()
    expected = read_frame(directory, reconstruct=True)
    assert not errors and mesh.GetNumberOfCells() == expected['cells']
    identifiers = vtk_to_numpy(mesh.GetCellData().GetArray('GlobalCellId'))
    for field in FIELDS:
        actual = vtk_to_numpy(mesh.GetCellData().GetArray(field))
        assert all(actual[i] == expected['fields'][field][gid] for i, gid in enumerate(identifiers)), field
    size = vtk.vtkCellSizeFilter()
    size.SetInputData(mesh)
    size.Update()
    dimension = 1 + (expected['shape'][1] > 1) + (expected['shape'][2] > 1)
    measure = size.GetOutput().GetCellData().GetArray({1: 'Length', 2: 'Area', 3: 'Volume'}[dimension])
    assert measure.GetRange()[0] > 0, (directory, 'degenerate/reversed geometry', measure.GetRange())
    legacy_error = 0.0
    if (directory / 'legacy.vtk').exists():
        old = vtk.vtkGenericDataObjectReader()
        old.ReadAllScalarsOn()
        old.SetFileName(str(directory / 'legacy.vtk'))
        old.Update()
        legacy = old.GetOutput()
        assert legacy.GetNumberOfCells() == mesh.GetNumberOfCells()
        for field in FIELDS:
            x = vtk_to_numpy(legacy.GetCellData().GetArray(field))
            y = vtk_to_numpy(mesh.GetCellData().GetArray(field))
            if field not in ['Velocity_X', 'Shift_beta_r']:
                assert (x == y).all(), ('legacy field', field)
            error = max(abs(u-v)/max(1.0, abs(u), abs(v)) for u,v in zip(x,y))
            assert error < 1e-11
            legacy_error = max(legacy_error, error)
        # Different point numbering is valid; compare each cell's actual vertices.
        for i in range(mesh.GetNumberOfCells()):
            one, two = legacy.GetCell(i), mesh.GetCell(i)
            a = sorted(set(one.GetPoints().GetPoint(k) for k in range(one.GetNumberOfPoints())))
            b = sorted(set(two.GetPoints().GetPoint(k) for k in range(two.GetNumberOfPoints())))
            if legacy.GetCellType(i) != mesh.GetCellType(i):
                # Legacy partial-phi grids use degenerate polar hexahedra.
                # C04 emits the same six distinct vertices as a proper wedge.
                assert legacy.GetCellType(i) == 12 and mesh.GetCellType(i) == 13 and len(a) == 6
            assert len(a) == len(b)
            assert max(abs(x-y) for p, q in zip(a, b) for x, y in zip(p, q)) < 1e-12
    return dict(checked=True, version=vtk.vtkVersion.GetVTKVersion(), cells=mesh.GetNumberOfCells(),
                minimum_measure=measure.GetRange()[0], legacy_compared=(directory / 'legacy.vtk').exists(),
                legacy_field_error=legacy_error)


def reject_corruptions(source, output):
    results = []
    def mutate(name, action):
        dest = output / name
        shutil.copytree(source, dest)
        action(dest)
        try:
            read_frame(dest)
        except (ValueError, OSError, ET.ParseError, KeyError, struct.error):
            results.append(name)
        else:
            raise AssertionError(('corrupt visualization accepted', name))
    def edit(path, a, b):
        data = path.read_bytes()
        assert a in data
        path.write_bytes(data.replace(a, b, 1))
    piece = 'piece_0000000000.vtu'
    mutate('missing', lambda d: (d/piece).unlink())
    mutate('truncated', lambda d: (d/piece).write_bytes((d/piece).read_bytes()[:-100]))
    mutate('negative-offset', lambda d: edit(d/piece, b'offset="0"', b'offset="-1"'))
    mutate('wrong-count', lambda d: edit(d/piece, b'NumberOfCells="5624"', b'NumberOfCells="5623"'))
    mutate('wrong-shape', lambda d: edit(d/piece, b'>37 19 8</DataArray>', b'>38 19 8</DataArray>'))
    mutate('wrong-local', lambda d: edit(d/piece, b'Name="LocalShape"', b'Name="MissingShape"'))
    mutate('wrong-owner', lambda d: edit(d/piece, b'Name="Rank" NumberOfTuples="1" NumberOfComponents="1" format="ascii">0',
                                        b'Name="Rank" NumberOfTuples="1" NumberOfComponents="1" format="ascii">1'))
    mutate('wrong-time', lambda d: edit(d/'frame.pvd', b'1.2500000000000000E-001', b'2.5000000000000000E-001'))
    mutate('ghosts', lambda d: edit(d/'frame.pvtu', b'GhostLevel="0"', b'GhostLevel="1"'))
    mutate('wrong-index-fields', lambda d: edit(d/'frame.pvtu', b'Name="Density"', b'Name="Other"'))
    mutate('unsafe-source', lambda d: edit(d/'frame.pvtu', piece.encode(), b'../outside.vtu'))
    mutate('duplicate-source', lambda d: edit(d/'frame.pvtu', b'</PUnstructuredGrid>',
             b'<Piece Source="piece_0000000000.vtu"/></PUnstructuredGrid>'))
    mutate('trailing-bytes', lambda d: edit(d/piece, b'</AppendedData>', b'X</AppendedData>'))
    def binary(dest, field, code, value):
        path = dest / piece
        data = bytearray(path.read_bytes())
        marker = b'<AppendedData encoding="raw">_'; start = data.index(marker)
        header = ET.fromstring(data[:start] + b'</VTKFile>')
        node = next(n for n in header.findall('.//Piece//DataArray') if n.attrib['Name'] == field)
        offset = start + len(marker) + int(node.attrib['offset']) + 8
        struct.pack_into(code, data, offset, value)
        path.write_bytes(data)
    mutate('nan-field', lambda d: binary(d, 'Density', '>d', math.nan))
    mutate('wrong-global-id', lambda d: binary(d, 'GlobalCellId', '>q', 42))
    mutate('wrong-cell-owner', lambda d: binary(d, 'RankId', '>q', 3))
    mutate('invalid-connectivity', lambda d: binary(d, 'connectivity', '>q', -1))
    mutate('bad-cell-type', lambda d: binary(d, 'types', '>B', 255))
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--probe', type=Path, required=True)
    parser.add_argument('--driver', type=Path, required=True)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--serial', action='store_true')
    parser.add_argument('--smoke', action='store_true')
    parser.add_argument('--vtk-required', action='store_true')
    args = parser.parse_args()
    args.probe, args.driver = args.probe.resolve(), args.driver.resolve()
    output = args.output.resolve() if args.output else Path(tempfile.mkdtemp(prefix='grial_vtk_matrix_'))
    output.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS='2', OMP_WAIT_POLICY='PASSIVE', OMP_DYNAMIC='FALSE', OPENBLAS_NUM_THREADS='1')
    results, rejections = [], []
    def run(directory, grid, command, success=True):
        directory.mkdir(parents=True, exist_ok=True)
        prefix = [] if args.serial else ['mpirun', '--oversubscribe', '--bind-to', 'none', '-n', str(math.prod(grid))]
        command = [*prefix, *map(str, command)]
        with (directory / f'launch_{len(list(directory.glob("launch_*.log"))):03}.log').open('xb') as stream:
            proc = subprocess.run(command, cwd=directory, env=env, stdout=stream, stderr=subprocess.STDOUT, timeout=120)
        assert (proc.returncode == 0) == success, (directory, proc.returncode, command)
    grids = [(1, 1)] if args.serial else [(1, 1), (2, 1), (1, 2), (4, 1), (1, 4), (2, 2), (4, 2), (2, 4), (8, 1)]
    fixtures = [1, 4, 8, 10, 11] if args.smoke else list(range(1, 16))
    references = {}
    for fixture in fixtures:
        for grid in grids:
            if fixture in [1, 2, 5, 11, 14, 15] and grid[1] > 1:
                continue
            if fixture == 15 and grid != (1, 1):
                continue
            directory = output / f'fixture_{fixture:02}_{grid[0]}x{grid[1]}'
            run(directory, grid, [args.probe, *grid, fixture, 'normal'])
            frame = read_frame(directory / 'good', reconstruct=True)
            coordinates = geometry(directory / 'good')
            if grid == (1, 1):
                references[fixture] = frame, coordinates
                # Analytic density/pressure indexed independently of decomposition.
                for gid, value in enumerate(frame['fields']['Density']):
                    assert value == 1 + 0.001 * gid
                    assert frame['fields']['Pressure'][gid] == 0.1 + 0.0001 * gid
            else:
                baseline, previous = references[fixture]
                compare_fields(frame, baseline)
                assert coordinates == previous, 'geometry changes with decomposition'
            consumer = vtk_check(directory / 'good', args.vtk_required)
            results.append(dict(kind='topology', fixture=fixture, grid=grid, cells=frame['cells'],
                                field_error=compare_fields(frame, references[fixture][0]), geometry_exact=True, vtk=consumer))
            print('TOPOLOGY', fixture, grid, 'PASS', flush=True)
            (output/'PROGRESS.json').write_text(json.dumps(results, indent=2)+'\n')
    for mode in ['duplicate', 'io', 'limit']:
        grid = (1, 1) if args.serial else (2, 2)
        directory = output / ('io_' + mode)
        run(directory, grid, [args.probe, *grid, 8, mode])
        read_frame(directory/'good')
        assert not list((directory/'truncated').glob('frame.pvtu'))
        results.append(dict(kind='io', mode=mode, previous_frame_preserved=True))
    if not args.serial:
        for mode in ['path-disagree', 'epoch-disagree', 'grid-disagree', 'nan', 'spin-disagree']:
            directory = output / ('reject_' + mode)
            run(directory, (2, 2), [args.probe, 2, 2, 8, mode], success=False)
            assert not (directory/'bad').exists() and not (directory/'other').exists()
            read_frame(directory/'good')
            rejections.append(mode)
    corruption = reject_corruptions(output/'fixture_08_1x1/good', output/'corrupt')
    selected = ('wave_0', 'fd_3d_a0.9_0') if args.smoke else (
        'wave_0', 'sod_1', 'khi_0', 'jet_1', 'michel_0',
        'fd_3d_a0.9_0', 'fm_3d_a0_1', 'fm_3d_a0.9_2', 'fd_sag_1', 'fm_sag_0')
    for name, text in cases():
        if name not in selected:
            continue
        parent = output / name
        parent.mkdir()
        par = parent/'case.par'
        par.write_text(text)
        def segment(directory, mode, source, grid, quota, session):
            run(directory, grid, [args.driver, mode, source, *grid, 0, 0, quota, session])
            return receipt(directory, session)
        reference, chain = parent/'reference', parent/'chain'
        segment(reference, 'cold', par, (1, 1), 0, 1)
        first = segment(chain, 'cold', par, (1, 1), 1, 1)
        split = (1, 1) if args.serial else ((2, 2) if 'ny=19\n' in text else (4, 1))
        second = segment(chain, 'resume', first['checkpoint'], split, 2, 2)
        third = segment(chain, 'resume', second['checkpoint'], (1, 1) if args.serial else (2, 1), 0, 3)
        assert third['state'] == 'completed'
        physical = compare(reference, chain, [first, second, third])
        root_a = next((reference/'data').iterdir())/'mpi_fields'
        root_b = next((chain/'data').iterdir())/'mpi_fields'
        frames_a = collect_frames(root_a, root_a/'frames.pvd')
        frames_b = collect_frames(root_b, root_b/'frames.pvd')
        keyed = {row['step']: row for row in frames_b}
        assert {row['step'] for row in frames_a} <= set(keyed)
        frame_error = 0.0
        for row in frames_a:
            other = keyed[row['step']]
            assert row['time'] == other['time']
            a = read_frame(row['directory'], reconstruct=True)
            b = read_frame(other['directory'], reconstruct=True)
            frame_error = max(frame_error, compare_fields(a, b))
        results.append(dict(kind='segments', case=name, common_frames=len(frames_a),
                            actual_frames=len(frames_b), visualization_error=frame_error, **physical))
        print('SEGMENTS', name, 'PASS', flush=True)
    assert len([x for x in results if x['kind'] == 'segments']) == len(selected), selected
    report = dict(passed=True, serial=args.serial, results=results, rejections=rejections,
                  reader_corruptions=corruption, vtk_required=args.vtk_required,
                  topology_comparisons=sum(x['kind']=='topology' for x in results),
                  segment_comparisons=len(selected))
    (output/'VERIFIED.json').write_text(json.dumps(report, indent=2)+'\n')
    print('C04 VISUALIZATION MATRIX PASSED:', output, flush=True)


if __name__ == '__main__':
    main()
