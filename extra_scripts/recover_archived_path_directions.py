#!/usr/bin/env python3
"""Recover Table 1 directions with its archived pooled-baseline optimizer.

This audit wrapper is not a production fitting option. The October analytic
baseline gradient changed some estimates, so signed directions from that
implementation must not be paired with September discovery p-values.
"""

import ast
from pathlib import Path
import subprocess

from tealeaf.sc import ec_block_glmm
from extra_scripts.run_paired_path_test import main


ARCHIVED_REVISION = "d31da47"


def restore_archived_baseline():
    source = subprocess.check_output(["git", "show", f"{ARCHIVED_REVISION}:tealeaf/sc/ec_block_glmm.py"], cwd=Path(__file__).resolve().parents[1], text=True)
    nodes = [node for node in ast.parse(source).body if isinstance(node, ast.FunctionDef) and node.name == "pooled_isoform_weights"]
    if len(nodes) != 1:
        raise ValueError("archived baseline function not found uniquely")
    module = ast.Module(body=nodes, type_ignores=[])
    exec(compile(module, f"archived:{ARCHIVED_REVISION}:pooled_isoform_weights", "exec"), ec_block_glmm.__dict__)
    print(f"Recovering archived directions with pooled baseline from {ARCHIVED_REVISION}", flush=True)


if __name__ == "__main__":
    restore_archived_baseline()
    main()
