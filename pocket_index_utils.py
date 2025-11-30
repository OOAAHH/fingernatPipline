#!/usr/bin/env python3
"""
Helpers for mapping solution/prediction residues to a canonical pocket index
based on split.txt and xx.index.

Design:
- split.txt defines the solution-side pocket span(s) per PDB stem.
- xx.index (if present) defines the prediction-side span(s); if absent,
  prediction is assumed to share the same span as solution.
- Both spans are expanded into linear sequences of (chain, resseq), which
  must have identical length. This defines a canonical pocket index i=1..N.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple


@dataclass
class PocketIndexMapping:
    length: int
    # For each canonical index i (0-based), the corresponding residue identity
    # on the solution side and prediction side.
    sol_positions: List[Tuple[str, int]]
    pred_positions: List[Tuple[str, int]]


def _parse_index_spec(spec: str) -> List[Tuple[str, int, int]]:
    """
    Parse an index spec like 'A:17:69' or 'A:1:27,A:29:13' into segments.

    Semantics follow your index files:
    - Items are encoded as 'chain:start:length'
    - 'start' is the first residue number on that chain
    - 'length' is the number of residues in this segment
    """
    parts = [s.strip() for s in spec.replace(";", ",").split(",") if s.strip()]
    segs: List[Tuple[str, int, int]] = []
    for p in parts:
        items = p.split(":")
        if len(items) != 3:
            continue
        chain = items[0].strip()
        try:
            start = int(items[1])
            length = int(items[2])
        except Exception:
            continue
        if not chain or length <= 0:
            continue
        segs.append((chain, start, length))
    return segs


def _expand_segments(segs: List[Tuple[str, int, int]]) -> List[Tuple[str, int]]:
    """
    Expand index segments into a flat list of (chain, resseq).
    """
    out: List[Tuple[str, int]] = []
    for chain, start, length in segs:
        for offset in range(length):
            out.append((chain, start + offset))
    return out


def _read_split_spec_for_solution(split_path: Path, solution_stem: str) -> Optional[str]:
    """
    Read split.txt and return the index spec string for the given solution stem.

    Expected line format (first token is stem, separated by tab/space/comma):
        <stem> <tab/space/comma> <index_spec>
    """
    if not split_path.exists():
        return None
    try:
        lines = [ln.strip() for ln in split_path.read_text().splitlines() if ln.strip()]
    except Exception:
        return None
    for ln in lines:
        # accept comma/space/tab separated, use first field as stem
        if "\t" in ln:
            stem, rest = ln.split("\t", 1)
        elif "," in ln:
            stem, rest = ln.split(",", 1)
        else:
            parts = ln.split(maxsplit=1)
            if len(parts) == 1:
                continue
            stem, rest = parts[0], parts[1]
        stem = stem.strip()
        if stem == solution_stem:
            return rest.strip()
    return None


def _read_xx_index_spec(index_dir: Path) -> Optional[str]:
    """
    Read xx.index from an index directory, if present.
    """
    xx = index_dir / "xx.index"
    if not xx.exists():
        return None
    try:
        # Concatenate non-comment, non-empty lines
        parts: List[str] = []
        for ln in xx.read_text().splitlines():
            s = ln.strip()
            if not s or s.startswith("#"):
                continue
            parts.append(s)
        return ",".join(parts) if parts else None
    except Exception:
        return None


def build_pocket_index_mapping_from_paths(
    puzzle: str,
    solution_stem: str,
    *,
    solutions_roots: List[Path],
    index_roots: List[Path],
) -> Optional[PocketIndexMapping]:
    """
    Build a PocketIndexMapping for a given puzzle + solution stem, searching
    split.txt and xx.index under the provided candidate roots.

    solutions_roots:
      list of roots that may contain solution/PZxx/split.txt
    index_roots:
      list of roots that may contain PZxx/index/xx.index
    """
    split_spec: Optional[str] = None
    for s_root in solutions_roots:
        split_path = s_root / puzzle / "split.txt"
        split_spec = _read_split_spec_for_solution(split_path, solution_stem)
        if split_spec:
            break
    if not split_spec:
        return None

    xx_spec: Optional[str] = None
    for i_root in index_roots:
        idx_dir = i_root / puzzle / "index"
        xx_spec = _read_xx_index_spec(idx_dir)
        if xx_spec:
            break

    sol_segs = _parse_index_spec(split_spec)
    if not sol_segs:
        return None
    sol_list = _expand_segments(sol_segs)

    if xx_spec:
        pred_segs = _parse_index_spec(xx_spec)
        if not pred_segs:
            return None
        pred_list = _expand_segments(pred_segs)
        if not pred_list:
            return None
        # If lengths differ (e.g. solution uses a pocket slice while xx.index
        # covers the full chain), align by the minimum shared length.
        if len(sol_list) != len(pred_list):
            L = min(len(sol_list), len(pred_list))
            sol_list = sol_list[:L]
            pred_list = pred_list[:L]
    else:
        # No xx.index provided: use split for both sides
        pred_list = list(sol_list)

    return PocketIndexMapping(
        length=len(sol_list),
        sol_positions=sol_list,
        pred_positions=pred_list,
    )
