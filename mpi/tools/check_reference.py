#!/usr/bin/env python3
"""Verify the immutable hydro reference; optionally compare the editable copy."""
import argparse
import hashlib
import json
from pathlib import Path
import tarfile


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--working-copy', action='store_true')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    provenance = json.loads((root/'PROVENANCE.json').read_text())
    expected = json.loads((root/'BASE_SOURCE_SHA256.json').read_text())
    for name, digest in provenance.get('reference_files_sha256', {}).items():
        assert hashlib.sha256((root/name).read_bytes()).hexdigest() == digest, name
    archive = root/provenance['archive']
    assert hashlib.sha256(archive.read_bytes()).hexdigest() == provenance['archive_sha256'], 'Archive changed'
    with tarfile.open(archive, 'r:gz') as stream:
        members = stream.getmembers()
        assert len(members) == len(expected)
        assert {m.name for m in members} == set(expected)
        for member in members:
            assert member.isfile(), member.name
            data = stream.extractfile(member).read()
            assert hashlib.sha256(data).hexdigest() == expected[member.name], member.name
    if args.working_copy:
        changed = [name for name, digest in expected.items()
                   if not (root/name).is_file() or
                   hashlib.sha256((root/name).read_bytes()).hexdigest() != digest]
        if changed:
            raise SystemExit('Working copy differs from initial base:\n'+'\n'.join(changed))
    print(f'HYDRO REFERENCE VERIFIED: {len(expected)} files; working_copy={args.working_copy}')


if __name__ == '__main__':
    main()
