#!/usr/bin/env python3
"""Redraw available panels with authorized external data / 使用授权的外部数据重绘。"""
import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument('--set', choices=['core', 'cellchat', 'all'], default='core')
    selection.add_argument('--panel', action='append', help='Panel ID; repeatable / 可重复指定子图')
    parser.add_argument('--data-root', type=Path, default=os.environ.get('SLC7A11_SOURCE_DATA_ROOT'),
                        help='Authorized private archive root / 已授权下载的私有资料库根目录')
    parser.add_argument('--rscript', default='Rscript')
    parser.add_argument('--out', type=Path, default=ROOT / 'validation_outputs')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    panels = json.loads((ROOT / 'scripts/panels.json').read_text())
    if args.panel:
        requested = set(args.panel)
        unknown = requested - {p['panel'] for p in panels}
        if unknown:
            parser.error('Unknown panels / 未知子图: ' + ', '.join(sorted(unknown)))
        selected = [p for p in panels if p['panel'] in requested]
        incomplete = [p['panel'] for p in selected if p['kind'] == 'incomplete']
        if incomplete:
            parser.error('Final source/code unavailable / 最终来源未找回: ' + ', '.join(incomplete))
    else:
        kinds = {'r', 'velocity'} if args.set == 'core' else {'cellchat'} if args.set == 'cellchat' else {'r', 'velocity', 'cellchat'}
        selected = [p for p in panels if p['kind'] in kinds]
    if args.data_root is None and not args.dry_run:
        parser.error('Source data are private and excluded from this public repository. '
                     'Supply --data-root with an authorized archive copy. '
                     '绘图数据暂未公开，请通过通讯作者获得授权后指定 --data-root。')
    data_root = args.data_root.expanduser().resolve() if args.data_root is not None else Path('/AUTHORIZED_SOURCE_DATA')
    if not args.dry_run:
        missing = [str(data_root / p['directory'] / 'Source_data' / name)
                   for p in selected for name in p['required_source_files']
                   if not (data_root / p['directory'] / 'Source_data' / name).is_file()]
        if missing:
            parser.error('Missing source files / 缺少输入文件:\n' + '\n'.join(missing))
    rscript = shutil.which(args.rscript)
    if not args.dry_run and any(p['kind'] in {'r', 'cellchat'} for p in selected) and not rscript:
        parser.error('Rscript not found / 找不到 Rscript; supply --rscript /absolute/path/to/Rscript')
    rscript = rscript or args.rscript
    output = args.out.expanduser().resolve()
    # Prevent an explicit output argument from overwriting any authorized input.
    protected = [data_root / p['directory'] / 'Source_data' for p in panels]
    if any(output == p or p in output.parents for p in protected):
        parser.error('Output must be outside Source_data / 输出目录不能位于 Source_data 中。')
    commands = []
    for panel in selected:
        directory = ROOT / panel['directory']
        source = data_root / panel['directory'] / 'Source_data'
        out = output / panel['panel']
        if panel['kind'] == 'r':
            commands.append((panel['panel'], 'R', [rscript, str(directory / 'Code/Replot_from_source.R'), str(out)]))
        elif panel['kind'] == 'cellchat':
            commands.append((panel['panel'], 'CellChat', [rscript, str(directory / 'Code/Recovery_replot_from_aggregate_networks.R'),
                             str(source / panel['required_source_files'][0]), str(out), panel['cellchat_panel']]))
        else:
            for group in ['Control', '1A4_Antibody']:
                commands.append((panel['panel'], group, [sys.executable, str(directory / 'Code/replot_velocity_from_source.py'),
                                 str(source / ('source_' + group) / 'pdf_vector_source.json.gz'), str(out / (group + '.svg'))]))
    if args.dry_run:
        print('SLC7A11_SOURCE_DATA_ROOT=' + shlex.quote(str(data_root)))
        for ident, group, command in commands:
            print(ident, group, shlex.join(command))
        return 0
    output.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    environment['SLC7A11_SOURCE_DATA_ROOT'] = str(data_root)
    records = []
    for ident, group, command in commands:
        directory = output / ident
        directory.mkdir(parents=True, exist_ok=True)
        started = time.monotonic()
        run = subprocess.run(command, cwd=ROOT, env=environment, capture_output=True, text=True)
        (directory / (group + '.stdout.log')).write_text(run.stdout)
        (directory / (group + '.stderr.log')).write_text(run.stderr)
        records.append(dict(panel=ident, group=group, command=command, returncode=run.returncode,
                            elapsed_seconds=round(time.monotonic() - started, 3)))
        print(ident, group, 'PASS' if run.returncode == 0 else 'FAIL', flush=True)
    (output / 'run_status.json').write_text(json.dumps(records, indent=2, ensure_ascii=False) + '\n')
    failed = sum(r['returncode'] != 0 for r in records)
    print(f'{len(records)-failed}/{len(records)} commands passed / 命令通过; logs / 日志: {output}')
    return 1 if failed else 0


if __name__ == '__main__':
    raise SystemExit(main())
