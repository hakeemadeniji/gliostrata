"""The stage 3b subtype analysis scripts: populations, helpers, refusals.

scripts/build_subtype_centroids.py, scripts/subtype_scoring_set_dependence.py
and scripts/compare_subtype_labels.py are loaded by path and driven with
synthetic frames; no real data and no CGGA path are read.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"


def _load(name):
    sys.path.insert(0, str(SCRIPTS))
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


build = _load("build_subtype_centroids")
dep = _load("subtype_scoring_set_dependence")
cmp_ = _load("compare_subtype_labels")


def test_scoring_population_is_cohort_patients_with_expression_in_cohort_order():
    cohort = pd.DataFrame({"patient_id": ["C", "A", "B", "D"]})
    expr = pd.DataFrame(np.ones((2, 3)), columns=["A", "C", "X"])
    assert build.scoring_population(cohort, expr) == ["C", "A"]


def test_store_population_is_an_availability_rule():
    cohort = pd.DataFrame({"patient_id": ["A", "B", "C", "D", "E"],
                           "age": [60, 61, None, 50, 55], "kps": [80, None, 90, 70, 60],
                           "sex": ["M", "F", "M", "F", "M"]})
    joined = pd.DataFrame({"patient_id": ["A", "B", "C", "D"],
                           "mgmt_methylated": pd.array([1, 0, 0, None], dtype="Int8"),
                           "idh_mutated": pd.array([0, 0, 1, 0], dtype="Int8")})
    # B lacks KPS, C lacks age, D lacks a MGMT call, E has no joined row.
    assert dep.store_population(cohort, joined, ["A", "B", "C", "D", "E"]) == ["A"]


def test_changed_counts_only_shared_patients():
    a = pd.Series({"A": "Classical", "B": "Neural", "C": "Proneural"})
    b = pd.Series({"B": "Classical", "C": "Proneural", "Z": "Neural"})
    assert dep.n_changed(a, b) == 1


def test_seed_is_documented_and_not_a_registered_seed():
    assert dep.SEED == 20261127
    assert dep.SEED != 20251126  # the seed of the registered TCGA partitions


def test_exact_interval_and_agreement():
    assert cmp_.exact_ci(10, 10)[1] == 1.0
    lo, hi = cmp_.exact_ci(56, 76)
    assert 0.62 < lo < 0.63 and 0.83 < hi < 0.84
    calls = pd.Series({"A": "Classical", "B": "Neural", "C": "Proneural", "D": "Proneural"})
    labels = pd.Series({"A": "Classical", "B": "Classical", "C": "G-CIMP", "D": "Proneural"})
    r = cmp_.agreement(calls, labels, ["A", "B", "C", "D"])
    assert (r["n"], r["agree"]) == (3, 2)  # G-CIMP is not counted
    conf = cmp_.confusion(calls, labels, ["A", "B", "C", "D"])
    assert conf["G-CIMP"]["Proneural"] == 1 and conf["Classical"]["Neural"] == 1


@pytest.mark.parametrize("mod", [build, dep, cmp_])
def test_every_script_refuses_a_cgga_path(mod, tmp_path):
    with pytest.raises(ValueError, match="refusing"):
        mod.main(["--cohort", str(tmp_path / "absent" / "cgga" / "clinical_cohort.parquet")])
