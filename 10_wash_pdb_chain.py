#!/usr/bin/env python3
"""
统一单链到A，这样输出的时候方便看。
Wash a PDB: if a MODEL (or whole file without MODEL cards) contains exactly one
polymer chain (ATOM records), unify its chain ID to 'A'.

By default, only ATOM lines are changed; HETATM and other records are preserved.

单文件:
  python bin/fingeRNAt/10_wash_pdb_chain.py --input Puzzles/original/PZ38/DasFARFAR2NNHolo.pdb

批量遍历目录:
  python bin/fingeRNAt/10_wash_pdb_chain.py --dir Puzzles/original

行为规范:
  - 仅当检测到“单链+配体(HETATM)”时，才对该文件进行洗链（把所有 ATOM 的链统一为 'A'），并就地覆盖原文件。
  - 对真实多链或无配体的文件不做改动。
"""

from __future__ import annotations

import argparse
from pathlib import Path
import math


def change_chain_char(line: str, new_chain: str) -> str:
    if len(line) < 22:
        return line
    return line[:21] + new_chain[0] + line[22:]

def change_icode_char(line: str, new_icode: str | None) -> str:
    # insertion code at column 27 (index 26), write space if None/''
    if len(line) < 27:
        return line
    ch = (new_icode or ' ')
    return line[:26] + ch[0] + line[27:]

def parse_xyz(line: str):
    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        return (x, y, z)
    except Exception:
        return None

def connected_components(indices: list[int], coords: list[tuple[float,float,float]], threshold: float = 1.9) -> list[list[int]]:
    # Simple O(n^2) graph clustering by distance threshold
    n = len(indices)
    adj = [[] for _ in range(n)]
    for i in range(n):
        xi, yi, zi = coords[i]
        for j in range(i+1, n):
            xj, yj, zj = coords[j]
            dx = xi - xj
            dy = yi - yj
            dz = zi - zj
            if dx*dx + dy*dy + dz*dz <= threshold*threshold:
                adj[i].append(j)
                adj[j].append(i)
    comps: list[list[int]] = []
    visited = [False]*n
    for i in range(n):
        if visited[i]:
            continue
        stack = [i]
        visited[i] = True
        comp = []
        while stack:
            k = stack.pop()
            comp.append(k)
            for nb in adj[k]:
                if not visited[nb]:
                    visited[nb] = True
                    stack.append(nb)
        comps.append(comp)
    return comps


def wash_pdb(in_path: Path, out_path: Path, inplace: bool = False, force: bool = False) -> None:
    text = in_path.read_text().splitlines(keepends=True)

    # Determine if file has explicit MODEL/MODULE cards
    has_model = any(l.startswith(('MODEL', 'MODULE')) for l in text)
    out_lines = []

    def process_block(lines: list[str]) -> list[str]:
        # First pass: unify polymer chain to 'A' if single-chain
        polymer_seen = set()
        new_lines = []
        for ln in lines:
            if ln.startswith('ATOM'):
                polymer_seen.add(ln[21].strip())
        unify_polymer = force or (len([c for c in polymer_seen]) <= 1)
        for idx, ln in enumerate(lines):
            if ln.startswith('ATOM') and unify_polymer:
                new_lines.append(change_chain_char(ln, 'A'))
            else:
                new_lines.append(ln)

        # Second pass: split HETATM groups that are actually disconnected into separate insertion codes
        # Group HETATM by (resname, chain, resseq, icode)
        groups: dict[tuple[str,str,str,str], list[int]] = {}
        for i, ln in enumerate(new_lines):
            if not ln.startswith('HETATM'):
                continue
            resname = ln[17:20].strip()
            chain = ln[21].strip()
            resseq = ln[22:26]
            icode = ln[26]
            key = (resname, chain, resseq, icode)
            groups.setdefault(key, []).append(i)

        letters = [chr(c) for c in range(ord('A'), ord('Z')+1)]

        for key, idxs in groups.items():
            # collect coordinates
            coords = []
            pick_idxs = []
            for i in idxs:
                xyz = parse_xyz(new_lines[i])
                if xyz is None:
                    continue
                coords.append(xyz)
                pick_idxs.append(i)
            if len(pick_idxs) <= 1:
                continue
            comps = connected_components(pick_idxs, coords, threshold=1.9)
            if len(comps) <= 1:
                continue
            # Assign icode per component; keep original icode for first component, others get A/B/C/... not equal to original
            orig_icode = key[3]
            avail = [ch for ch in letters if ch != orig_icode.strip()]
            comp_icodes: list[str|None] = []
            comp_icodes.append(orig_icode if orig_icode != '' else None)
            for _ in range(1, len(comps)):
                comp_icodes.append(avail.pop(0) if avail else 'Z')
            # Apply changes to lines
            for ci, comp in enumerate(comps):
                new_ic = comp_icodes[ci]
                for pos in comp:
                    line_idx = pick_idxs[pos] # 将局部索引映射回原始行号
                    new_lines[line_idx] = change_icode_char(new_lines[line_idx], new_ic)

        return new_lines

    if not has_model:
        out_lines.extend(process_block(text))
    else:
        block: list[str] = []
        in_model = False
        for ln in text:
            rec = ln[:6].strip()
            if rec in ('MODEL', 'MODULE'):
                # flush previous block
                if block:
                    out_lines.extend(process_block(block))
                    block = []
                in_model = True
                out_lines.append(ln)
                continue
            if rec == 'ENDMDL':
                block.append(ln)
                out_lines.extend(process_block(block))
                block = []
                in_model = False
                continue
            if in_model:
                block.append(ln)
            else:
                out_lines.append(ln)
        if block:
            out_lines.extend(process_block(block))

    # Write result
    # Always write to out_path; caller decides in-place or not
    if inplace:
        in_path.write_text(''.join(out_lines))
    else:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(''.join(out_lines))


def analyze_pdb(in_path: Path) -> tuple[bool, bool, bool]:
    """返回 (has_model, all_blocks_single_polymer, has_hetatm).
    all_blocks_single_polymer: 若文件中每个模型（或整文件）聚合物链(ATOM)数量<=1 则为 True。
    has_hetatm: 文件内是否存在 HETATM 行（视作存在配体/异原子）。
    """
    text = in_path.read_text().splitlines()
    has_model = any(l.startswith(('MODEL', 'MODULE')) for l in text)
    blocks: list[list[str]] = []
    if not has_model:
        blocks = [text]
    else:
        cur: list[str] = []
        in_m = False
        for ln in text:
            rec = ln[:6].strip()
            if rec in ('MODEL', 'MODULE'):
                if cur:
                    blocks.append(cur)
                    cur = []
                in_m = True
                continue
            if rec == 'ENDMDL':
                cur.append(ln)
                blocks.append(cur)
                cur = []
                in_m = False
                continue
            if in_m:
                cur.append(ln)
        if cur:
            blocks.append(cur)

    all_single = True
    has_hetatm = False
    for blk in blocks:
        seen = set()
        for ln in blk:
            if ln.startswith('ATOM'):
                seen.add((ln[21].strip() or '_'))
            elif ln.startswith('HETATM'):
                has_hetatm = True
        if len(seen) > 1:
            all_single = False
    return has_model, all_single, has_hetatm


def main():
    ap = argparse.ArgumentParser(description='Normalize polymer chain IDs to A for single-chain(+ligand) models; single file or directory traversal. Always overwrites original files when applicable.')
    ap.add_argument('--input', type=Path, help='Input PDB path')
    ap.add_argument('--dir', type=Path, help='Traverse this directory recursively for *.pdb')
    args = ap.parse_args()

    if not args.input and not args.dir:
        raise SystemExit('Please provide --input or --dir')

    if args.input:
        if not args.input.exists():
            raise SystemExit(f'Input not found: {args.input}')
        _, all_single, has_hetatm = analyze_pdb(args.input)
        feature = '单链+配体' if (all_single and has_hetatm) else ('多链' if not all_single else '单链')
        if all_single and has_hetatm:
            wash_pdb(args.input, args.input, inplace=True, force=False)
        print(f"{args.input.name}:{feature}")
        return

    # 目录模式（就地覆盖）
    root = args.dir
    if not root.exists():
        raise SystemExit(f'Directory not found: {root}')
    summary_lines = []
    for pdb in root.rglob('*.pdb'):
        rel = pdb.relative_to(root)
        _, all_single, has_hetatm = analyze_pdb(pdb)
        feature = '单链+配体' if (all_single and has_hetatm) else ('多链' if not all_single else '单链')
        if all_single and has_hetatm:
            wash_pdb(pdb, pdb, inplace=True, force=False)
        summary_lines.append(f"{str(rel)}:{feature}")

    for line in summary_lines:
        print(line)


if __name__ == '__main__':
    main()
