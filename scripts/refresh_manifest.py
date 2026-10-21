#!/usr/bin/env python3
"""Record package hashes after intentional edits; do not run to validate."""
import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT/'FILE_MANIFEST.tsv'


def package_files():
    for path in sorted(ROOT.rglob('*')):
        relative=path.relative_to(ROOT)
        if path.is_file() and path!=TARGET and not any(p in {'.git','validation_outputs','__pycache__'} for p in relative.parts) and path.name!='.DS_Store':
            yield path


def main():
    with TARGET.open('w',encoding='utf-8',newline='') as handle:
        writer=csv.writer(handle,delimiter='\t')
        writer.writerow(['relative_path','bytes','sha256'])
        for path in package_files():
            writer.writerow([str(path.relative_to(ROOT)),path.stat().st_size,hashlib.sha256(path.read_bytes()).hexdigest()])
    print('Recorded',TARGET)


if __name__=='__main__':
    main()
