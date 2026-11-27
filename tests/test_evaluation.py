"""Tests for gliostrata.evaluation -- pre-registration and frozen partitions.

These tests are deliberately adversarial about the two properties the project's
defensibility rests on:

* the split is a pure function of ``(seed, patient_id)`` -- not of row order,
  cohort size, or library version;
* the pre-registered thresholds cannot be changed without the build failing.
"""
from __future__ import annotations

import json

import pandas as pd
import pytest

from gliostrata.evaluation import (
    THRESHOLD_TABLE_SHA256,
    THRESHOLDS,
    CIBound,
    Direction,
    PartitionSealedError,
    Status,
    UnregisteredMetricError,
    assign_split,
    build_manifest,
    canonical_table,
    evaluate,
    evaluate_all,
    freeze_partition,
    hash_unit_interval,
    partition_cohort,
    read_manifest,
    registered_metrics,
    sha256_over_ids,
    stratum_balance,
    threshold_table_checksum,
    verify_partition,
    write_manifest,
)

SEED = 20251126


def make_cohort(n: int, *, prefix: str = "PT", n_strata: int = 4) -> pd.DataFrame:
    """A synthetic cohort with stable ids and a balanced stratification column."""
    strata = ["Classical", "Mesenchymal", "Proneural", "Neural"][:n_strata]
    return pd.DataFrame(
        {
            "patient_id": [f"{prefix}-{i:04d}" for i in range(n)],
            "subtype": [strata[i % len(strata)] for i in range(n)],
            "age": [40 + (i % 40) for i in range(n)],
        }
    )


# ---------------------------------------------------------------------------
# Hash primitives
# ---------------------------------------------------------------------------


class TestHashPrimitives:
    def test_hash_unit_interval_is_in_range_and_stable(self):
        values = [hash_unit_interval(f"k{i}") for i in range(500)]
        assert all(0.0 <= v < 1.0 for v in values)
        # Stable across calls (and, by construction, across processes).
        assert hash_unit_interval("k7") == values[7]

    def test_hash_unit_interval_is_roughly_uniform(self):
        vals = [hash_unit_interval(f"PT-{i:05d}") for i in range(5000)]
        below_30 = sum(v < 0.30 for v in vals) / len(vals)
        assert 0.27 < below_30 < 0.33

    def test_assign_split_depends_only_on_id_and_seed(self):
        assert assign_split("PT-0001", seed=SEED) == assign_split("PT-0001", seed=SEED)
        # A different seed is a different partition.
        flips = sum(
            assign_split(f"PT-{i:04d}", seed=SEED)
            != assign_split(f"PT-{i:04d}", seed=SEED + 1)
            for i in range(500)
        )
        assert flips > 0

    def test_assign_split_rejects_bad_fraction(self):
        with pytest.raises(ValueError, match="test_fraction"):
            assign_split("PT-0001", seed=SEED, test_fraction=0.0)
        # 1.0 is valid (an all-test partition); anything above it is not.
        with pytest.raises(ValueError, match=r"must be in \(0, 1\]"):
            assign_split("PT-0001", seed=SEED, test_fraction=1.5)

    def test_sha256_over_ids_is_order_independent(self):
        assert sha256_over_ids(["b", "a", "c"]) == sha256_over_ids(["c", "b", "a"])
        assert sha256_over_ids(["a", "b"]) != sha256_over_ids(["a", "b", "c"])


# ---------------------------------------------------------------------------
# The determinism property -- the core requirement
# ---------------------------------------------------------------------------


class TestPartitionDeterminism:
    def test_identical_cohort_gives_identical_split(self):
        cohort = make_cohort(300)
        a = partition_cohort(cohort, seed=SEED)
        b = partition_cohort(cohort.copy(), seed=SEED)
        assert a.train_ids == b.train_ids
        assert a.test_ids == b.test_ids

    def test_row_order_does_not_change_the_split(self):
        cohort = make_cohort(300)
        base = partition_cohort(cohort, seed=SEED)
        for frac_seed in (0, 1, 2, 99):
            shuffled = cohort.sample(frac=1.0, random_state=frac_seed).reset_index(drop=True)
            assert partition_cohort(shuffled, seed=SEED).assignments() == base.assignments()
        # Also robust to a reversed / sorted-descending order.
        reversed_df = cohort.iloc[::-1].reset_index(drop=True)
        assert partition_cohort(reversed_df, seed=SEED).assignments() == base.assignments()

    def test_adding_unrelated_patients_does_not_move_existing_ones(self):
        cohort = make_cohort(200)
        base = partition_cohort(cohort, seed=SEED).assignments()
        extended = pd.concat(
            [cohort, make_cohort(50, prefix="NEW")], ignore_index=True
        )
        after = partition_cohort(extended, seed=SEED).assignments()
        for pid, split in base.items():
            assert after[pid] == split, f"{pid} moved from {split} to {after[pid]}"

    def test_removing_unrelated_patients_does_not_move_existing_ones(self):
        cohort = make_cohort(200)
        base = partition_cohort(cohort, seed=SEED).assignments()
        trimmed = cohort[~cohort["patient_id"].isin(base and list(base)[:40])]
        after = partition_cohort(trimmed, seed=SEED).assignments()
        for pid, split in after.items():
            assert base[pid] == split

    def test_split_does_not_depend_on_cohort_size(self):
        """Patient PT-0007's assignment must be the same in an n=10 and an n=5000 cohort."""
        small = partition_cohort(make_cohort(10), seed=SEED).assignments()
        large = partition_cohort(make_cohort(5000), seed=SEED).assignments()
        for pid in small:
            assert small[pid] == large[pid]

    def test_duplicate_patient_rows_are_collapsed(self):
        """A multi-sample patient must never straddle the split."""
        cohort = make_cohort(100)
        doubled = pd.concat([cohort, cohort], ignore_index=True)
        base = partition_cohort(cohort, seed=SEED)
        dup = partition_cohort(doubled, seed=SEED)
        assert dup.n_total == base.n_total == 100
        assert dup.assignments() == base.assignments()

    def test_test_fraction_is_approximately_honoured(self):
        result = partition_cohort(make_cohort(4000), seed=SEED, test_fraction=0.30)
        assert 0.27 < result.n_test / result.n_total < 0.33

    def test_extra_columns_and_dtypes_do_not_matter(self):
        cohort = make_cohort(200)
        base = partition_cohort(cohort, seed=SEED).assignments()
        mutated = cohort.copy()
        mutated["age"] = mutated["age"].astype("float32")
        mutated["unrelated"] = range(len(mutated))
        assert partition_cohort(mutated, seed=SEED).assignments() == base


class TestPartitionValidation:
    def test_missing_id_column_raises(self):
        with pytest.raises(ValueError, match="no 'patient_id' column"):
            partition_cohort(pd.DataFrame({"x": [1, 2]}), seed=SEED)

    def test_empty_cohort_raises(self):
        with pytest.raises(ValueError, match="empty cohort"):
            partition_cohort(pd.DataFrame({"patient_id": []}), seed=SEED)

    def test_null_patient_id_raises(self):
        df = pd.DataFrame({"patient_id": ["a", None]})
        with pytest.raises(ValueError, match="null patient_id"):
            partition_cohort(df, seed=SEED)

    def test_missing_stratify_column_raises(self):
        with pytest.raises(ValueError, match="stratify_by column"):
            partition_cohort(make_cohort(10), seed=SEED, stratify_by="nope")

    def test_inconsistent_stratum_per_patient_raises(self):
        df = pd.DataFrame(
            {"patient_id": ["a", "a"], "subtype": ["Classical", "Mesenchymal"]}
        )
        with pytest.raises(ValueError, match="conflicting"):
            partition_cohort(df, seed=SEED, stratify_by="subtype")

    def test_bad_split_name_raises(self):
        result = partition_cohort(make_cohort(10), seed=SEED)
        with pytest.raises(ValueError, match="split must be one of"):
            result.ids("validation")


# ---------------------------------------------------------------------------
# Stratification
# ---------------------------------------------------------------------------


class TestStratification:
    def test_stratified_split_is_balanced_across_strata(self):
        result = partition_cohort(
            make_cohort(4000), seed=SEED, test_fraction=0.30, stratify_by="subtype"
        )
        balance = stratum_balance(result)
        assert set(balance) == {"Classical", "Mesenchymal", "Proneural", "Neural"}
        for stratum, row in balance.items():
            assert abs(row["test_fraction"] - 0.30) < 0.05, (stratum, row)

    def test_stratified_split_keeps_hash_determinism(self):
        cohort = make_cohort(500)
        base = partition_cohort(cohort, seed=SEED, stratify_by="subtype").assignments()
        shuffled = cohort.sample(frac=1.0, random_state=3).reset_index(drop=True)
        assert (
            partition_cohort(shuffled, seed=SEED, stratify_by="subtype").assignments() == base
        )
        extended = pd.concat([cohort, make_cohort(50, prefix="NEW")], ignore_index=True)
        after = partition_cohort(extended, seed=SEED, stratify_by="subtype").assignments()
        for pid, split in base.items():
            assert after[pid] == split

    def test_stratification_changes_the_split_vs_unstratified(self):
        cohort = make_cohort(500)
        plain = partition_cohort(cohort, seed=SEED).assignments()
        strat = partition_cohort(cohort, seed=SEED, stratify_by="subtype").assignments()
        assert plain != strat  # the stratum salt genuinely enters the hash

    def test_missing_stratum_values_become_one_na_stratum(self):
        df = pd.DataFrame(
            {"patient_id": [f"PT-{i:03d}" for i in range(50)], "subtype": [None] * 50}
        )
        result = partition_cohort(df, seed=SEED, stratify_by="subtype")
        assert set(result.strata.values()) == {"__NA__"}

    def test_unstratified_balance_report_is_empty(self):
        assert stratum_balance(partition_cohort(make_cohort(50), seed=SEED)) == {}


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------


@pytest.fixture
def cohort_file(tmp_path):
    cohort = make_cohort(300)
    path = tmp_path / "cohort.parquet"
    cohort.to_parquet(path, index=False)
    return cohort, path


class TestManifest:
    def test_round_trip_preserves_everything(self, cohort_file, tmp_path):
        cohort, path = cohort_file
        result = partition_cohort(cohort, seed=SEED, stratify_by="subtype")
        manifest = build_manifest(
            result,
            cohort_source=path,
            partition_name="unit_test_v1",
            cohort_rows=len(cohort),
        )
        out = write_manifest(manifest, tmp_path / "manifest.json")
        reloaded = read_manifest(out)

        assert reloaded.train_ids == manifest.train_ids
        assert reloaded.test_ids == manifest.test_ids
        assert reloaded.seed == SEED
        assert reloaded.test_fraction == result.test_fraction
        assert reloaded.stratify_by == "subtype"
        assert reloaded.cohort_sha256 == manifest.cohort_sha256
        assert reloaded.cohort_rows == 300
        assert reloaded.n_train + reloaded.n_test == 300
        assert reloaded.train_ids_sha256 == sha256_over_ids(result.train_ids)
        assert reloaded.all_ids_sha256 == manifest.all_ids_sha256
        assert reloaded.to_dict() == manifest.to_dict()

    def test_manifest_records_required_provenance(self, cohort_file, tmp_path):
        cohort, path = cohort_file
        _, manifest, _ = freeze_partition(
            cohort,
            cohort_source=path,
            manifest_path=tmp_path / "m.json",
            partition_name="prov_v1",
            seed=SEED,
        )
        data = json.loads((tmp_path / "m.json").read_text(encoding="utf-8"))
        for key in (
            "seed",
            "train_ids",
            "test_ids",
            "train_ids_sha256",
            "test_ids_sha256",
            "cohort_source",
            "cohort_sha256",
            "cohort_rows",
            "n_train",
            "n_test",
            "created_utc",
            "sealed",
            "schema_version",
        ):
            assert key in data, key
        # ISO-8601 UTC, parseable, Z-suffixed.
        assert manifest.created_utc.endswith("Z")
        assert len(manifest.created_utc) == 20
        assert manifest.cohort_source.endswith("cohort.parquet")
        assert len(manifest.cohort_sha256) == 64

    def test_manifest_from_dict_rejects_incomplete_payload(self):
        from gliostrata.evaluation import PartitionManifest

        with pytest.raises(ValueError, match="missing required field"):
            PartitionManifest.from_dict({"seed": 1})

    def test_ids_are_sorted_in_the_manifest(self, cohort_file, tmp_path):
        cohort, path = cohort_file
        _, manifest, _ = freeze_partition(
            cohort,
            cohort_source=path,
            manifest_path=tmp_path / "m.json",
            partition_name="sorted_v1",
            seed=SEED,
        )
        assert list(manifest.train_ids) == sorted(manifest.train_ids)
        assert list(manifest.test_ids) == sorted(manifest.test_ids)
        assert not set(manifest.train_ids) & set(manifest.test_ids)


# ---------------------------------------------------------------------------
# Sealing
# ---------------------------------------------------------------------------


class TestSeal:
    def test_sealed_manifest_refuses_overwrite(self, cohort_file, tmp_path):
        cohort, path = cohort_file
        mpath = tmp_path / "sealed.json"
        freeze_partition(
            cohort,
            cohort_source=path,
            manifest_path=mpath,
            partition_name="sealed_v1",
            seed=SEED,
            seal=True,
        )
        with pytest.raises(PartitionSealedError) as exc:
            freeze_partition(
                cohort,
                cohort_source=path,
                manifest_path=mpath,
                partition_name="sealed_v1",
                seed=SEED + 1,
            )
        msg = str(exc.value)
        # The refusal must explain the science, not just the mechanics.
        assert "not a file-permission problem" in msg
        assert "test" in msg and "train" in msg
        assert "pre-registration amendment" in msg
        # And the file on disk must be untouched.
        assert read_manifest(mpath).seed == SEED

    def test_force_allows_overwrite(self, cohort_file, tmp_path):
        cohort, path = cohort_file
        mpath = tmp_path / "sealed.json"
        freeze_partition(
            cohort,
            cohort_source=path,
            manifest_path=mpath,
            partition_name="sealed_v1",
            seed=SEED,
        )
        freeze_partition(
            cohort,
            cohort_source=path,
            manifest_path=mpath,
            partition_name="sealed_v2",
            seed=SEED + 1,
            force=True,
        )
        reloaded = read_manifest(mpath)
        assert reloaded.seed == SEED + 1
        assert reloaded.partition_name == "sealed_v2"

    def test_unsealed_manifest_may_be_overwritten(self, cohort_file, tmp_path):
        cohort, path = cohort_file
        mpath = tmp_path / "draft.json"
        freeze_partition(
            cohort,
            cohort_source=path,
            manifest_path=mpath,
            partition_name="draft_v1",
            seed=SEED,
            seal=False,
        )
        freeze_partition(
            cohort,
            cohort_source=path,
            manifest_path=mpath,
            partition_name="draft_v1",
            seed=SEED + 5,
            seal=False,
        )
        assert read_manifest(mpath).seed == SEED + 5

    def test_seal_check_happens_before_anything_is_written(self, cohort_file, tmp_path):
        cohort, path = cohort_file
        mpath = tmp_path / "sealed.json"
        freeze_partition(
            cohort,
            cohort_source=path,
            manifest_path=mpath,
            partition_name="sealed_v1",
            seed=SEED,
        )
        before = mpath.read_bytes()
        with pytest.raises(PartitionSealedError):
            freeze_partition(
                cohort,
                cohort_source=path,
                manifest_path=mpath,
                partition_name="other",
                seed=1,
            )
        assert mpath.read_bytes() == before


# ---------------------------------------------------------------------------
# Verification / drift detection
# ---------------------------------------------------------------------------


class TestAllTestPartition:
    """``test_fraction=1.0`` holds out every patient: the seal for an external cohort.

    An external cohort has no training data, so its frozen partition is an empty
    training list and a test list of every patient id.
    """

    def test_every_patient_is_held_out(self):
        cohort = make_cohort(300)
        result = partition_cohort(cohort, seed=SEED, test_fraction=1.0)
        assert result.n_train == 0
        assert result.train_ids == ()
        assert result.test_ids == tuple(sorted(cohort["patient_id"]))
        assert result.n_test == 300
        assert result.test_fraction == 1.0
        assert set(result.assignments().values()) == {"test"}

    def test_duplicate_ids_are_collapsed(self):
        cohort = pd.concat([make_cohort(50), make_cohort(50).head(7)], ignore_index=True)
        result = partition_cohort(cohort, seed=SEED, test_fraction=1.0)
        assert result.n_test == 50
        assert result.n_train == 0

    def test_is_independent_of_row_order_and_of_other_patients(self):
        cohort = make_cohort(300)
        full = partition_cohort(cohort, seed=SEED, test_fraction=1.0)
        part = partition_cohort(
            cohort.sample(frac=0.5, random_state=3), seed=SEED, test_fraction=1.0
        )
        assert part.n_train == 0
        assert set(part.test_ids) <= set(full.test_ids)
        assert part.n_test == len(part.test_ids)

    def test_fraction_above_one_is_rejected(self):
        with pytest.raises(ValueError, match=r"must be in \(0, 1\]"):
            partition_cohort(make_cohort(10), seed=SEED, test_fraction=1.01)
        with pytest.raises(ValueError, match=r"must be in \(0, 1\]"):
            assign_split("PT-0001", seed=SEED, test_fraction=1.5)

    def test_assign_split_at_one_is_always_test(self):
        assert all(
            assign_split(f"PT-{i:04d}", seed=SEED, test_fraction=1.0) == "test"
            for i in range(200)
        )

    def test_stratified_all_test_balance_is_one(self):
        result = partition_cohort(make_cohort(80), seed=SEED, stratify_by="subtype",
                                  test_fraction=1.0)
        assert all(row["test_fraction"] == 1.0 for row in stratum_balance(result).values())

    def test_freeze_seals_with_an_empty_training_list(self, cohort_file, tmp_path):
        cohort, path = cohort_file
        result, manifest, mpath = freeze_partition(
            cohort,
            cohort_source=path,
            manifest_path=tmp_path / "ext.json",
            partition_name="ext_v1",
            seed=SEED,
            test_fraction=1.0,
        )
        assert result.n_train == 0
        assert manifest.sealed is True
        assert manifest.test_fraction == 1.0
        assert manifest.n_train == 0
        assert manifest.n_test == 300
        assert manifest.train_ids == ()
        assert manifest.train_ids_sha256 == sha256_over_ids([])
        data = json.loads(mpath.read_text(encoding="utf-8"))
        assert data["train_ids"] == []
        assert len(data["test_ids"]) == 300
        assert read_manifest(mpath).to_dict() == manifest.to_dict()

    def test_verify_passes_and_detects_a_dropped_patient(self, cohort_file, tmp_path):
        cohort, path = cohort_file
        _, manifest, _ = freeze_partition(
            cohort,
            cohort_source=path,
            manifest_path=tmp_path / "ext.json",
            partition_name="ext_v2",
            seed=SEED,
            test_fraction=1.0,
        )
        assert verify_partition(manifest, cohort, cohort_source=path).ok
        dropped = tuple(sorted(cohort["patient_id"]))[:1]
        report = verify_partition(manifest, cohort.iloc[1:], cohort_source=path)
        assert not report.ok
        assert report.removed_ids == dropped
        assert report.added_ids == () and report.moved_ids == ()

    def test_sealed_all_test_manifest_refuses_overwrite(self, cohort_file, tmp_path):
        cohort, path = cohort_file
        mpath = tmp_path / "ext.json"
        freeze_partition(
            cohort,
            cohort_source=path,
            manifest_path=mpath,
            partition_name="ext_v3",
            seed=SEED,
            test_fraction=1.0,
        )
        with pytest.raises(PartitionSealedError):
            freeze_partition(
                cohort,
                cohort_source=path,
                manifest_path=mpath,
                partition_name="ext_v3",
                seed=SEED + 1,
                test_fraction=1.0,
            )
        assert read_manifest(mpath).seed == SEED


class TestVerify:
    @pytest.fixture
    def frozen(self, cohort_file, tmp_path):
        cohort, path = cohort_file
        _, manifest, mpath = freeze_partition(
            cohort,
            cohort_source=path,
            manifest_path=tmp_path / "frozen.json",
            partition_name="frozen_v1",
            seed=SEED,
        )
        return cohort, path, manifest, mpath

    def test_no_drift_is_reported_as_verified(self, frozen):
        cohort, path, manifest, _ = frozen
        report = verify_partition(manifest, cohort, cohort_source=path)
        assert report.ok
        assert report.n_drifted == 0
        assert report.added_ids == () and report.removed_ids == () and report.moved_ids == ()
        assert report.cohort_hash_matches
        assert "VERIFIED" in report.describe()

    def test_no_drift_under_row_shuffle(self, frozen):
        cohort, path, manifest, _ = frozen
        shuffled = cohort.sample(frac=1.0, random_state=11).reset_index(drop=True)
        assert verify_partition(manifest, shuffled, cohort_source=path).ok

    def test_detects_added_ids(self, frozen):
        cohort, _, manifest, _ = frozen
        extended = pd.concat([cohort, make_cohort(3, prefix="EXTRA")], ignore_index=True)
        report = verify_partition(manifest, extended)
        assert not report.ok
        assert set(report.added_ids) == {"EXTRA-0000", "EXTRA-0001", "EXTRA-0002"}
        assert report.removed_ids == ()
        assert report.moved_ids == ()
        assert "not in the frozen partition" in report.describe()

    def test_detects_removed_ids(self, frozen):
        cohort, _, manifest, _ = frozen
        dropped = sorted(manifest.test_ids)[:2]
        trimmed = cohort[~cohort["patient_id"].isin(dropped)]
        report = verify_partition(manifest, trimmed)
        assert not report.ok
        assert set(report.removed_ids) == set(dropped)
        assert report.added_ids == ()
        assert "missing from the" in report.describe()

    def test_detects_moved_ids(self, frozen):
        """A manifest hand-edited to put a train patient in test must be caught."""
        cohort, _, manifest, mpath = frozen
        data = json.loads(mpath.read_text(encoding="utf-8"))
        victim = sorted(data["train_ids"])[0]
        data["train_ids"] = [i for i in data["train_ids"] if i != victim]
        data["test_ids"] = sorted([*data["test_ids"], victim])
        data["n_train"] -= 1
        data["n_test"] += 1
        mpath.write_text(json.dumps(data, indent=2), encoding="utf-8")

        report = verify_partition(mpath, cohort)
        assert not report.ok
        assert report.moved_ids == ((victim, "test", "train"),)
        assert report.added_ids == () and report.removed_ids == ()

    def test_detects_train_to_test_move_as_contamination(self, frozen):
        """Changing the seed moves patients both ways; the train->test direction is flagged."""
        cohort, _, manifest, _ = frozen
        from gliostrata.evaluation import PartitionManifest

        drifted = PartitionManifest.from_dict({**manifest.to_dict(), "seed": SEED + 1})
        report = verify_partition(drifted, cohort)
        assert not report.ok
        assert report.moved_ids
        assert any("CONTAMINATION RISK" in m for m in report.messages)

    def test_detects_source_file_rewrite(self, frozen, tmp_path):
        cohort, path, manifest, _ = frozen
        # Same patient ids, different column values -> ids verify, content hash does not.
        rewritten = cohort.copy()
        rewritten["age"] = rewritten["age"] + 1
        rewritten.to_parquet(path, index=False)
        report = verify_partition(manifest, rewritten, cohort_source=path)
        assert not report.ok
        assert report.added_ids == () and report.moved_ids == ()
        assert not report.cohort_hash_matches
        assert any("content hash changed" in m for m in report.messages)

    def test_missing_source_file_is_reported(self, frozen, tmp_path):
        cohort, _, manifest, _ = frozen
        report = verify_partition(manifest, cohort, cohort_source=tmp_path / "gone.parquet")
        assert not report.ok
        assert any("does not exist" in m for m in report.messages)

    def test_row_count_change_alone_does_not_fail_verification(self, frozen):
        """Duplicating rows changes the row count but must not un-freeze the split."""
        cohort, _, manifest, _ = frozen
        doubled = pd.concat([cohort, cohort], ignore_index=True)
        report = verify_partition(manifest, doubled)
        assert report.ok
        assert any("Row count changed" in m for m in report.messages)

    def test_verify_accepts_a_path(self, frozen):
        cohort, _, _, mpath = frozen
        assert verify_partition(mpath, cohort).ok

    def test_report_is_json_serialisable(self, frozen):
        cohort, path, manifest, _ = frozen
        report = verify_partition(manifest, cohort, cohort_source=path)
        assert json.loads(json.dumps(report.to_dict()))["ok"] is True


# ---------------------------------------------------------------------------
# Threshold table integrity
# ---------------------------------------------------------------------------


class TestThresholdTableIntegrity:
    def test_checksum_matches_the_frozen_constant(self):
        """If this fails, a PRE-REGISTERED NUMBER CHANGED.

        Do not "fix" the test by pasting the new hash unless the change is a
        deliberate pre-registration amendment that is also recorded, with its
        date and reason, in docs/preregistration.md.
        """
        assert threshold_table_checksum() == THRESHOLD_TABLE_SHA256, (
            "The pre-registered threshold table has changed. See the docstring "
            "of this test before touching THRESHOLD_TABLE_SHA256."
        )

    def test_checksum_is_sensitive_to_any_change(self):
        """Mutating one byte of the canonical table must change the hash."""
        import hashlib

        original = canonical_table()
        tampered = original.replace('"value":0.88', '"value":0.87')
        assert tampered != original
        assert (
            hashlib.sha256(tampered.encode()).hexdigest() != THRESHOLD_TABLE_SHA256
        )

    def test_canonical_table_is_stable_and_order_independent(self):
        assert canonical_table() == canonical_table()
        rows = json.loads(canonical_table())
        assert [r["key"] for r in rows] == sorted(r["key"] for r in rows)

    def test_all_eleven_protocol_thresholds_are_registered(self):
        assert len(THRESHOLDS) == 11
        per_test = {}
        for t in THRESHOLDS.values():
            per_test[t.test_id] = per_test.get(t.test_id, 0) + 1
        assert per_test == {"T1": 2, "T2": 3, "T3": 2, "T4": 2, "T5": 2}

    def test_registry_is_read_only(self):
        with pytest.raises(TypeError):
            THRESHOLDS["T1_subtype_accuracy_cgga"] = None  # type: ignore[index]

    def test_thresholds_are_frozen_dataclasses(self):
        import dataclasses

        spec = THRESHOLDS["T2_survival_cindex_cgga"]
        with pytest.raises(dataclasses.FrozenInstanceError):
            spec.value = 0.01  # type: ignore[misc]

    def test_every_threshold_has_a_justification_and_unit(self):
        for key, t in THRESHOLDS.items():
            assert t.justification.strip(), key
            assert t.unit.strip(), key
            assert t.test_id in {"T1", "T2", "T3", "T4", "T5"}, key

    def test_ci_bound_matches_direction(self):
        """'Higher is better' must be judged on the lower bound, and vice versa."""
        for key, t in THRESHOLDS.items():
            expected = CIBound.LOWER if t.direction is Direction.GTE else CIBound.UPPER
            assert t.ci_bound is expected, key

    def test_protocol_values_are_exactly_as_pre_registered(self):
        expected = {
            "T1_subtype_accuracy_cgga": (Direction.GTE, 0.88),
            "T1_state_deconvolution_pearson_r": (Direction.GTE, 0.75),
            "T2_survival_cindex_cgga": (Direction.GTE, 0.65),
            "T2_cindex_improvement_vs_baseline": (Direction.GTE, 0.03),
            "T2_brier_18mo": (Direction.LTE, 0.25),
            "T3_growth_rmse_pct_final_volume": (Direction.LTE, 15.0),
            "T3_nsc_migration_r_in_vivo": (Direction.GTE, 0.65),
            "T4_optimizer_regret_reduction_pct": (Direction.GTE, 25.0),
            "T4_personalization_distinct_configs": (Direction.GTE, 3.0),
            "T5_endtoend_rmst_success_fraction": (Direction.GTE, 0.50),
            "T5_safety_predicted_harm_fraction": (Direction.LTE, 0.15),
        }
        assert set(expected) == set(registered_metrics())
        for key, (direction, value) in expected.items():
            assert THRESHOLDS[key].direction is direction, key
            assert THRESHOLDS[key].value == pytest.approx(value), key
        assert THRESHOLDS["T5_safety_predicted_harm_fraction"].point_estimate_limit == 0.10
        assert (
            THRESHOLDS["T4_optimizer_regret_reduction_pct"].auxiliary_requirement
            == "Wilcoxon signed-rank p < 0.01"
        )


# ---------------------------------------------------------------------------
# Threshold evaluation
# ---------------------------------------------------------------------------


class TestEvaluateUnregistered:
    def test_unregistered_metric_raises(self):
        with pytest.raises(UnregisteredMetricError, match="not a pre-registered metric"):
            evaluate("T2_survival_cindex_my_favourite_subset", 0.9, 0.8, 1.0)

    def test_unregistered_metric_error_lists_the_registered_ones(self):
        with pytest.raises(UnregisteredMetricError) as exc:
            evaluate("accuracy", 1.0, 1.0, 1.0)
        assert "T1_subtype_accuracy_cgga" in str(exc.value)

    def test_unregistered_metric_is_a_keyerror_subclass(self):
        with pytest.raises(KeyError):
            evaluate("nope", 0.5, 0.4, 0.6)

    def test_evaluate_all_raises_on_any_unregistered_key(self):
        with pytest.raises(UnregisteredMetricError):
            evaluate_all({"T2_brier_18mo": (0.2, 0.1, 0.24), "made_up": (1, 1, 1)})


class TestEvaluateGTE:
    """'Higher is better' metrics: the CI LOWER bound must clear the threshold."""

    KEY = "T2_survival_cindex_cgga"  # >= 0.65, lower bound

    def test_pass_when_lower_bound_clears(self):
        v = evaluate(self.KEY, 0.70, 0.66, 0.74)
        assert v.status is Status.PASS
        assert v.passed
        assert v.decisive_value == 0.66
        assert "95% CI lower bound 0.66 >= 0.65" in v.reason

    def test_pass_at_exact_boundary(self):
        assert evaluate(self.KEY, 0.70, 0.65, 0.75).status is Status.PASS

    def test_inconclusive_when_point_clears_but_lower_bound_does_not(self):
        v = evaluate(self.KEY, 0.68, 0.62, 0.74)
        assert v.status is Status.INCONCLUSIVE
        assert not v.passed
        assert "does not clear >= 0.65" in v.reason
        assert "point estimate 0.68 >= 0.65" in v.reason

    def test_fail_when_point_estimate_misses(self):
        v = evaluate(self.KEY, 0.60, 0.55, 0.64)
        assert v.status is Status.FAIL
        assert "point estimate 0.6 does not clear >= 0.65" in v.reason

    def test_a_high_upper_bound_cannot_rescue_a_fail(self):
        """The upper bound is irrelevant for a GTE metric -- no cherry-picking."""
        assert evaluate(self.KEY, 0.60, 0.40, 0.99).status is Status.FAIL

    def test_verdict_carries_the_protocol_justification(self):
        v = evaluate(self.KEY, 0.70, 0.66, 0.74)
        assert "0.65-0.68" in v.justification
        assert v.test_id == "T2"
        assert v.ci_bound is CIBound.LOWER
        assert v.threshold == 0.65

    def test_verdict_is_frozen_and_serialisable(self):
        import dataclasses

        v = evaluate(self.KEY, 0.70, 0.66, 0.74)
        with pytest.raises(dataclasses.FrozenInstanceError):
            v.status = Status.PASS  # type: ignore[misc]
        assert json.loads(json.dumps(v.to_dict()))["status"] == "PASS"


class TestEvaluateLTE:
    """'Lower is better' metrics: the CI UPPER bound must clear the threshold."""

    KEY = "T2_brier_18mo"  # <= 0.25, upper bound

    def test_pass_when_upper_bound_clears(self):
        v = evaluate(self.KEY, 0.20, 0.17, 0.23)
        assert v.status is Status.PASS
        assert v.decisive_value == 0.23
        assert "upper bound 0.23 <= 0.25" in v.reason

    def test_pass_at_exact_boundary(self):
        assert evaluate(self.KEY, 0.20, 0.15, 0.25).status is Status.PASS

    def test_inconclusive_when_point_clears_but_upper_bound_does_not(self):
        v = evaluate(self.KEY, 0.23, 0.19, 0.30)
        assert v.status is Status.INCONCLUSIVE
        assert "point estimate 0.23 <= 0.25" in v.reason

    def test_fail_when_point_estimate_misses(self):
        v = evaluate(self.KEY, 0.31, 0.27, 0.36)
        assert v.status is Status.FAIL

    def test_a_low_lower_bound_cannot_rescue_a_fail(self):
        assert evaluate(self.KEY, 0.40, 0.01, 0.80).status is Status.FAIL

    def test_t3_rmse_is_a_percentage_lte_metric(self):
        key = "T3_growth_rmse_pct_final_volume"  # <= 15% of final volume
        assert evaluate(key, 11.0, 8.0, 14.2).status is Status.PASS
        assert evaluate(key, 11.0, 8.0, 18.0).status is Status.INCONCLUSIVE
        assert evaluate(key, 22.0, 18.0, 26.0).status is Status.FAIL


class TestEvaluateCompositeRequirements:
    SAFETY = "T5_safety_predicted_harm_fraction"
    REGRET = "T4_optimizer_regret_reduction_pct"

    def test_safety_pass_requires_both_numbers(self):
        v = evaluate(self.SAFETY, 0.07, 0.04, 0.12)
        assert v.status is Status.PASS
        assert "satisfies the pre-registered limit <= 0.1" in v.reason

    def test_safety_fails_when_point_estimate_exceeds_ten_percent(self):
        # CI upper bound clears 0.15 but the point estimate breaches 0.10.
        v = evaluate(self.SAFETY, 0.12, 0.09, 0.14)
        assert v.status is Status.FAIL
        assert "point-estimate limit violated" in v.reason

    def test_safety_fails_when_upper_bound_exceeds_fifteen_percent(self):
        v = evaluate(self.SAFETY, 0.09, 0.05, 0.20)
        assert v.status is Status.INCONCLUSIVE  # point ok, CI does not settle it
        v2 = evaluate(self.SAFETY, 0.18, 0.14, 0.24)
        assert v2.status is Status.FAIL

    def test_regret_requires_the_wilcoxon_test(self):
        clearing = dict(point_estimate=32.0, ci_lower=27.0, ci_upper=38.0)
        assert evaluate(self.REGRET, **clearing, auxiliary_satisfied=True).status is Status.PASS
        assert (
            evaluate(self.REGRET, **clearing, auxiliary_satisfied=False).status is Status.FAIL
        )
        unknown = evaluate(self.REGRET, **clearing)
        assert unknown.status is Status.INCONCLUSIVE
        assert "was not reported" in unknown.reason

    def test_regret_records_the_auxiliary_requirement_in_the_verdict(self):
        v = evaluate(self.REGRET, 32.0, 27.0, 38.0, auxiliary_satisfied=True)
        assert v.auxiliary_requirement == "Wilcoxon signed-rank p < 0.01"
        assert v.auxiliary_satisfied is True

    def test_personalization_count_threshold(self):
        key = "T4_personalization_distinct_configs"  # >= 3 distinct configs
        assert evaluate(key, 5.0, 5.0, 5.0).status is Status.PASS  # degenerate CI
        assert evaluate(key, 2.0, 2.0, 2.0).status is Status.FAIL
        assert evaluate(key, 4.0, 2.0, 5.0).status is Status.INCONCLUSIVE


class TestEvaluateInputValidation:
    KEY = "T1_subtype_accuracy_cgga"

    def test_malformed_ci_raises(self):
        with pytest.raises(ValueError, match="Malformed confidence interval"):
            evaluate(self.KEY, 0.9, 0.95, 0.99)  # point below lower bound
        with pytest.raises(ValueError, match="Malformed confidence interval"):
            evaluate(self.KEY, 0.9, 0.85, 0.88)  # point above upper bound
        with pytest.raises(ValueError, match="Malformed confidence interval"):
            evaluate(self.KEY, 0.9, 0.95, 0.85)  # inverted interval

    def test_degenerate_ci_is_allowed(self):
        v = evaluate(self.KEY, 0.9, 0.9, 0.9)
        assert v.status is Status.PASS

    def test_non_finite_inputs_raise(self):
        for bad in (float("nan"), float("inf"), float("-inf")):
            with pytest.raises(ValueError, match="finite number"):
                evaluate(self.KEY, bad, 0.8, 1.0)


class TestEvaluateAll:
    def test_evaluates_a_full_t1_t5_result_set(self):
        results = {
            "T1_subtype_accuracy_cgga": (0.91, 0.885, 0.93),
            "T1_state_deconvolution_pearson_r": (0.79, 0.76, 0.82),
            "T2_survival_cindex_cgga": (0.63, 0.60, 0.66),
            "T2_cindex_improvement_vs_baseline": (0.01, -0.01, 0.03),
            "T2_brier_18mo": (0.21, 0.18, 0.24),
            "T3_growth_rmse_pct_final_volume": (10.0, 7.0, 13.0),
            "T3_nsc_migration_r_in_vivo": (0.70, 0.66, 0.74),
            "T4_optimizer_regret_reduction_pct": (31.0, 26.0, 36.0),
            "T4_personalization_distinct_configs": (4.0, 4.0, 4.0),
            "T5_endtoend_rmst_success_fraction": (0.56, 0.51, 0.61),
            "T5_safety_predicted_harm_fraction": (0.08, 0.05, 0.13),
        }
        verdicts = evaluate_all(
            results, auxiliary={"T4_optimizer_regret_reduction_pct": True}
        )
        assert set(verdicts) == set(registered_metrics())
        assert verdicts["T2_survival_cindex_cgga"].status is Status.FAIL
        assert verdicts["T2_cindex_improvement_vs_baseline"].status is Status.FAIL
        passed = {k for k, v in verdicts.items() if v.passed}
        assert len(passed) == 9

    def test_format_threshold_table_renders_every_row(self):
        from gliostrata.evaluation import format_threshold_table

        table = format_threshold_table()
        assert table.count("\n") == len(THRESHOLDS) + 1  # header + separator + rows - 1
        for t in THRESHOLDS.values():
            assert t.metric in table
        assert "Wilcoxon signed-rank p < 0.01" in table
