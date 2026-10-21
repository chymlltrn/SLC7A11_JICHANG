#!/usr/bin/env python3
"""Verify the public code without private data / 不读取私有数据的公开代码检查。"""
import ast
import csv
import hashlib
import json
import re
from pathlib import Path

from refresh_manifest import ROOT, TARGET, package_files


def main():
    files = {str(p.relative_to(ROOT)): p for p in package_files()}
    with TARGET.open(newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert len(rows) == len({r['relative_path'] for r in rows}), 'Duplicate manifest entries'
    expected = {r['relative_path']: r for r in rows}
    assert set(files) == set(expected), 'File inventory differs from FILE_MANIFEST.tsv'
    prohibited = {'Source_data', 'Matched_figures', 'Reference', 'Reconstructed_figures'}
    allowed_suffixes = {'.r', '.py', '.md', '.tsv'}
    patterns = [re.compile(r'(?:gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,})'),
                re.compile(r'BEGIN (?:OPENSSH |RSA |EC )?PRIVATE KEY')]
    for name, path in files.items():
        row = expected[name]
        assert path.stat().st_size == int(row['bytes']), 'Size mismatch: ' + name
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row['sha256'], 'Hash mismatch: ' + name
        assert not path.is_symlink(), 'No symlinks in the public package: ' + name
        assert not prohibited.intersection(path.relative_to(ROOT).parts), 'Private directory: ' + name
        assert (path.suffix.lower() in allowed_suffixes or name in {'.gitignore', 'scripts/panels.json'}), 'Unapproved file type: ' + name
        assert path.stat().st_size < 100_000, 'Unexpectedly large public text file: ' + name
        content = path.read_text()
        assert not any(pattern.search(content) for pattern in patterns), 'Credential pattern: ' + name
        if path.suffix == '.py':
            ast.parse(content, filename=name)
    panels = json.loads((ROOT / 'scripts/panels.json').read_text())
    assert len(panels) == 18 and len({p['panel'] for p in panels}) == 18
    assert {p['panel'] for p in panels if p['kind'] == 'incomplete'} == {'Figure_4B', 'Figure_4E', 'Figure_4F', 'Figure_5G', 'Figure_5H'}
    for panel in panels:
        directory = ROOT / panel['directory']
        assert (directory / 'README.md').is_file()
        assert panel['kind'] in {'r', 'cellchat', 'velocity', 'incomplete'}
        entry = ('Replot_from_source.R' if panel['kind'] == 'r' else
                 'Recovery_replot_from_aggregate_networks.R' if panel['kind'] == 'cellchat' else
                 'replot_velocity_from_source.py' if panel['kind'] == 'velocity' else 'README.md')
        assert (directory / 'Code' / entry).is_file(), panel['panel']
        assert bool(panel['required_source_files']) == (panel['kind'] != 'incomplete')
    with (ROOT / 'CODE_PROVENANCE.tsv').open(newline='') as handle:
        provenance = list(csv.DictReader(handle, delimiter='\t'))
    code_files = {name for name, p in files.items() if 'Code' in p.relative_to(ROOT).parts and p.suffix.lower() in {'.r', '.py'}}
    assert len(provenance) == len(code_files)
    assert {r['relative_path'] for r in provenance} == code_files
    for row in provenance:
        assert expected[row['relative_path']]['sha256'] == row['public_sha256']
    print(f'PASS: {len(files)} public files; 18 panels; code integrity and data exclusion verified.')
    print('通过：公开代码完整性与数据隔离检查；尚未验证私有 Source data 或上游科学分析。')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
