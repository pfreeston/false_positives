#!/usr/bin/env python3
"""Make AlphaFold-Multimer result pickles readable without jax.

The cluster AFM run (old jax) stored arrays as jax DeviceArrays, which unpickle
via jax._src.device_array.reconstruct_device_array. That module no longer
exists in current jax, and none of the local envs have jax, so AlphaJudge
cannot read iptm/ptm/PAE/distogram from them.

This stubs reconstruct_device_array to return the underlying numpy array
(identical values), then rewrites each result_*.pkl.gz under the same name.
Originals are moved to <run>/original_jax_pkls/ (AlphaJudge only globs the
run directory itself, so the backups are not picked up).

Usage: python convert_afm_pickles.py RUN_DIR [RUN_DIR ...]
"""
import gzip
import pickle
import sys
import types
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np


def reconstruct_device_array(fun, args, arr_state, aval_state):
    arr = fun(*args)
    arr.__setstate__(arr_state)
    return arr


def _install_stub():
    for name in ("jax", "jax._src"):
        sys.modules.setdefault(name, types.ModuleType(name))
    mod = types.ModuleType("jax._src.device_array")
    mod.reconstruct_device_array = reconstruct_device_array
    sys.modules["jax._src.device_array"] = mod


def _no_jax(obj):
    """Fail loudly if anything jax-typed survived."""
    if isinstance(obj, dict):
        for v in obj.values():
            _no_jax(v)
    elif type(obj).__module__.startswith("jax"):
        raise TypeError(f"jax object left in pickle: {type(obj)}")


def convert(pkl: Path) -> str:
    _install_stub()
    backup_dir = pkl.parent / "original_jax_pkls"
    backup = backup_dir / pkl.name
    if backup.exists():
        return f"skip (already converted) {pkl}"
    with gzip.open(pkl, "rb") as f:
        data = pickle.load(f)
    _no_jax(data)
    tmp = pkl.with_name(pkl.name + ".tmp")
    with gzip.open(tmp, "wb", compresslevel=4) as f:
        pickle.dump(data, f, protocol=4)
    # verify round-trip with no stub needed
    with gzip.open(tmp, "rb") as f:
        check = pickle.load(f)
    assert np.allclose(check["predicted_aligned_error"], data["predicted_aligned_error"])
    backup_dir.mkdir(exist_ok=True)
    pkl.rename(backup)
    tmp.rename(pkl)
    return f"ok {pkl}  iptm={float(data['iptm']):.3f} ptm={float(data['ptm']):.3f}"


def main():
    pkls = [p for d in sys.argv[1:] for p in sorted(Path(d).glob("result_*.pkl.gz"))]
    with ProcessPoolExecutor(max_workers=6) as ex:
        for msg in ex.map(convert, pkls):
            print(msg, flush=True)


if __name__ == "__main__":
    main()
