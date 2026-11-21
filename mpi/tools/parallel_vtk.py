#!/usr/bin/env python3
"""Inspect GRIAL C04 raw-appended VTU pieces and build a ParaView PVD collection.

The simulation never gathers fields. read_frame(reconstruct=True) is an explicit
offline operation that allocates the seven global scalar fields in this process.
"""
import argparse
from array import array
from collections import Counter
import json
import math
import os
from pathlib import Path
import struct
import sys
import tempfile
import xml.etree.ElementTree as ET

FIELDS = ('Density', 'Pressure', 'Velocity_X', 'Velocity_Y', 'Velocity_Z', 'Lapse_alpha', 'Shift_beta_r')
TYPES = {'Float64': ('d', 8), 'Int64': ('q', 8), 'UInt8': ('B', 1)}
VERTICES = {3: 2, 9: 4, 12: 8, 13: 6}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_piece(path):
    """Strict reader for the C04 subset; reject truncated/overlapping payloads."""
    path = Path(path)
    data = path.read_bytes()
    marker = b'<AppendedData encoding="raw">_'
    require(data.count(marker) == 1, f'{path}: missing raw appended data')
    start = data.index(marker)
    header = ET.fromstring(data[:start] + b'</VTKFile>')
    tail = b'</AppendedData></VTKFile>\n'
    require(data.endswith(tail), f'{path}: incomplete XML trailer')
    raw = memoryview(data)[start + len(marker):-len(tail)]
    require(header.attrib == dict(type='UnstructuredGrid', version='1.0', byte_order='BigEndian', header_type='UInt64'),
            f'{path}: unsupported VTK encoding')
    grid = header.find('UnstructuredGrid')
    require(grid is not None and len(grid.findall('Piece')) == 1, f'{path}: expected one local piece')
    meta = {}
    for node in grid.findall('FieldData/DataArray'):
        name = node.attrib['Name']
        require(name not in meta and node.attrib['format'] == 'ascii', f'{path}: invalid metadata')
        kind = float if node.attrib['type'] == 'Float64' else int
        values = [kind(x) for x in node.text.split()]
        require(len(values) == int(node.attrib['NumberOfComponents']) and node.attrib['NumberOfTuples'] == '1',
                f'{path}: invalid metadata shape')
        require(all(math.isfinite(x) for x in values), f'{path}: nonfinite metadata')
        meta[name] = values[0] if len(values) == 1 else values
    expected_meta = {'TimeValue', 'Step', 'GlobalShape', 'LocalShape', 'GlobalOffset', 'ProcessGrid',
                     'Rank', 'MappingId', 'LogRadius', 'Spin'}
    require(set(meta) == expected_meta, f'{path}: incomplete metadata')
    require(meta['TimeValue'] >= 0 and meta['Step'] >= 0 and meta['Rank'] >= 0, f'{path}: invalid epoch/rank')
    for key, count in [('GlobalShape', 3), ('LocalShape', 3), ('GlobalOffset', 3), ('ProcessGrid', 2)]:
        require(isinstance(meta[key], list) and len(meta[key]) == count, f'{path}: invalid {key}')
    require(all(x > 0 for x in meta['GlobalShape'] + meta['LocalShape'] + meta['ProcessGrid']), f'{path}: empty grid')
    require(all(0 <= o and o + n <= g for o, n, g in zip(meta['GlobalOffset'], meta['LocalShape'], meta['GlobalShape'])),
            f'{path}: block outside grid')
    piece = grid.find('Piece')
    cells, points = int(piece.attrib['NumberOfCells']), int(piece.attrib['NumberOfPoints'])
    require(cells == math.prod(meta['LocalShape']) and points > 0, f'{path}: wrong local cell count')
    arrays = {}
    next_offset = 0
    for parent in ['CellData', 'Points', 'Cells']:
        for node in piece.findall(parent + '/DataArray'):
            name, kind = node.attrib['Name'], node.attrib['type']
            require(name not in arrays and kind in TYPES and node.attrib['format'] == 'appended', f'{path}: invalid array')
            expected_type = 'UInt8' if name == 'types' else ('Int64' if name in
                            ['GlobalCellId', 'RankId', 'connectivity', 'offsets'] else 'Float64')
            require(kind == expected_type, f'{path}: unexpected data type for {name}')
            offset, components = int(node.attrib['offset']), int(node.attrib['NumberOfComponents'])
            require(offset == next_offset and offset + 8 <= len(raw), f'{path}: invalid binary offset')
            length = struct.unpack_from('>Q', raw, offset)[0]
            code, width = TYPES[kind]
            require(length % width == 0 and offset + 8 + length <= len(raw), f'{path}: truncated array')
            values = array(code)
            values.frombytes(raw[offset + 8:offset + 8 + length])
            if sys.byteorder != 'big' and width > 1:
                values.byteswap()
            require(components > 0 and len(values) % components == 0, f'{path}: invalid tuple shape')
            require(all(math.isfinite(x) for x in values), f'{path}: nonfinite array')
            if parent == 'CellData':
                require(components == 1 and len(values) == cells, f'{path}: wrong field size')
            if parent == 'Points':
                require(name == 'Points' and components == 3 and len(values) == 3 * points, f'{path}: wrong point size')
            arrays[name] = values
            next_offset = offset + 8 + length
    require(next_offset == len(raw), f'{path}: unreferenced binary bytes')
    require(set(arrays) == set(FIELDS) | {'GlobalCellId', 'RankId', 'Points', 'connectivity', 'offsets', 'types'},
            f'{path}: missing/unexpected arrays')
    require(len(arrays['offsets']) == cells and len(arrays['types']) == cells, f'{path}: incomplete cells')
    previous = 0
    for end, kind in zip(arrays['offsets'], arrays['types']):
        require(kind in VERTICES and end - previous == VERTICES[kind], f'{path}: invalid cell topology')
        vertices = arrays['connectivity'][previous:end]
        require(len(set(vertices)) == len(vertices) and all(0 <= x < points for x in vertices), f'{path}: invalid vertices')
        previous = end
    require(previous == len(arrays['connectivity']), f'{path}: trailing connectivity')
    gx, gy, _ = meta['GlobalShape']
    lx, ly, lz = meta['LocalShape']
    ox, oy, oz = meta['GlobalOffset']
    for index, (gid, rank) in enumerate(zip(arrays['GlobalCellId'], arrays['RankId'])):
        i, j, k = index % lx, (index // lx) % ly, index // (lx * ly)
        require(gid == i + ox + gx * (j + oy + gy * (k + oz)) and rank == meta['Rank'], f'{path}: invalid global ID/owner')
    return dict(path=str(path), metadata=meta, points=points, cells=cells, arrays=arrays)


def read_frame(directory, reconstruct=False):
    """Validate exact block coverage. Global arrays are allocated only on request."""
    directory = Path(directory)
    index = ET.parse(directory / 'frame.pvtu').getroot()
    require(index.attrib.get('type') == 'PUnstructuredGrid', 'index is not PUnstructuredGrid')
    grid = index.find('PUnstructuredGrid')
    require(grid is not None and grid.attrib.get('GhostLevel') == '0', 'unexpected ghost cells')
    declared = grid.findall('PCellData/PDataArray')
    require([node.attrib['Name'] for node in declared] == list(FIELDS) + ['GlobalCellId', 'RankId'], 'index fields mismatch')
    for node in declared:
        require(node.attrib.get('NumberOfComponents') == '1' and node.attrib.get('type') ==
                ('Int64' if node.attrib['Name'] in ['GlobalCellId', 'RankId'] else 'Float64'), 'index types mismatch')
    sources = [node.attrib['Source'] for node in grid.findall('Piece')]
    require(sources and len(set(sources)) == len(sources), 'duplicate/missing piece references')
    frame, fields, topology = None, {}, Counter()
    summaries = []
    for rank, name in enumerate(sources):
        require(name == f'piece_{rank:010}.vtu', 'unexpected piece path or rank order')
        piece = read_piece(directory / name)
        meta = piece['metadata']
        require(meta['Rank'] == rank, 'piece rank mismatch')
        if frame is None:
            frame = meta
            require(math.prod(meta['ProcessGrid']) == len(sources), 'process grid does not cover pieces')
            if reconstruct:
                fields = {key: array('d', [0.0]) * math.prod(meta['GlobalShape']) for key in FIELDS}
        for key in ['GlobalShape', 'ProcessGrid', 'TimeValue', 'Step', 'MappingId', 'LogRadius', 'Spin']:
            require(meta[key] == frame[key], f'piece metadata mismatch: {key}')
        expected_shape = list(frame['GlobalShape'])
        expected_offset = [0, 0, 0]
        for axis, coord in enumerate([rank % frame['ProcessGrid'][0], rank // frame['ProcessGrid'][0]]):
            base, remainder = divmod(frame['GlobalShape'][axis], frame['ProcessGrid'][axis])
            expected_shape[axis] = base + (coord < remainder)
            expected_offset[axis] = coord * base + min(coord, remainder)
        require(meta['LocalShape'] == expected_shape and meta['GlobalOffset'] == expected_offset,
                'pieces leave gaps, overlap, or violate declared block layout')
        if reconstruct:
            for key in FIELDS:
                for gid, value in zip(piece['arrays']['GlobalCellId'], piece['arrays'][key]):
                    fields[key][gid] = value
        topology.update(piece['arrays']['types'])
        summaries.append(dict(rank=rank, cells=piece['cells'], points=piece['points'],
                              shape=meta['LocalShape'], offset=meta['GlobalOffset']))
    pvd = ET.parse(directory / 'frame.pvd').getroot()
    datasets = pvd.findall('Collection/DataSet')
    require(pvd.attrib.get('type') == 'Collection' and len(datasets) == 1 and
            datasets[0].attrib.get('file') == 'frame.pvtu' and
            float(datasets[0].attrib['timestep']) == frame['TimeValue'], 'time index mismatch')
    require(sum(x['cells'] for x in summaries) == math.prod(frame['GlobalShape']), 'incomplete global coverage')
    result = dict(directory=str(directory), time=frame['TimeValue'], step=frame['Step'],
                  shape=frame['GlobalShape'], process_grid=frame['ProcessGrid'],
                  cells=math.prod(frame['GlobalShape']), pieces=summaries, topology=dict(topology))
    if reconstruct:
        result['fields'] = fields
    return result


def collect_frames(directory, target):
    directory, target = Path(directory), Path(target)
    candidates = sorted(p for p in directory.iterdir() if p.is_dir() and p.name.startswith('step_'))
    require(candidates, 'no step_* frames found')
    frames = [read_frame(p) for p in candidates]  # An incomplete generation is an error, not hidden.
    require(all(a['time'] < b['time'] and a['step'] < b['step'] for a, b in zip(frames, frames[1:])),
            'duplicate/nonmonotone frames')
    require(all(row['shape'] == frames[0]['shape'] for row in frames), 'global grid changed between frames')
    document = ET.Element('VTKFile', type='Collection', version='0.1', byte_order='BigEndian')
    collection = ET.SubElement(document, 'Collection')
    for row in frames:
        ET.SubElement(collection, 'DataSet', timestep=format(row['time'], '.17g'), group='', part='0',
                      file=Path(os.path.relpath(Path(row['directory']) / 'frame.pvtu', target.parent)).as_posix())
    require(target.suffix == '.pvd', 'collection target must end in .pvd')
    with tempfile.NamedTemporaryFile(dir=target.parent, prefix=target.name + '.', delete=False) as stream:
        temporary = Path(stream.name)
        ET.ElementTree(document).write(stream, encoding='utf-8', xml_declaration=True)
    os.replace(temporary, target)
    return frames


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--collection', type=Path, help='validate all step_* directories and publish this PVD index')
    args = parser.parse_args()
    try:
        result = collect_frames(args.directory, args.collection) if args.collection else read_frame(args.directory)
        print(json.dumps(result, indent=2))
    except (ValueError, OSError, ET.ParseError, KeyError, struct.error) as error:
        parser.exit(1, f'C04 validation failed: {error}\n')


if __name__ == '__main__':
    main()
