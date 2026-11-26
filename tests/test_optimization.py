"""Tests for L6 (configuration space + optimizer) and L7 (handoff)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from gliostrata.handoff.protocol import generate_protocols, protocols_to_json
from gliostrata.optimization.bayesian import (
    evaluate_personalization,
    optimize_configurations,
)
from gliostrata.optimization.config_space import (
    CELL_TYPES,
    Configuration,
    DOSE_MAX,
    DOSE_MIN,
    PAYLOADS,
    SCHEDULES,
    configuration_space_summary,
    sample_random_configuration,
)


class TestConfigurationSpace:
    def test_valid_configuration(self):
        cfg = Configuration(
            cell_type="iPSC-NSC",
            payload="CRAd",
            dose_cells=2_000_000,
            n_injection_sites=4,
            schedule="single",
        )
        assert cfg.cell_type == "iPSC-NSC"

    def test_invalid_dose_too_low(self):
        with pytest.raises(ValueError, match="dose_cells"):
            Configuration("iPSC-NSC", "CRAd", 100, 4, "single")

    def test_invalid_dose_too_high(self):
        with pytest.raises(ValueError, match="dose_cells"):
            Configuration("iPSC-NSC", "CRAd", 100_000_000, 4, "single")

    def test_invalid_injection_sites(self):
        with pytest.raises(ValueError, match="n_injection_sites"):
            Configuration("iPSC-NSC", "CRAd", 2_000_000, 12, "single")

    def test_incompatible_cell_payload(self):
        # HB1_F3_line × CAR is incompatible
        with pytest.raises(ValueError, match="Incompatible"):
            Configuration("HB1_F3_line", "CAR", 2_000_000, 4, "single")
        # MSC × CRAd is incompatible
        with pytest.raises(ValueError, match="Incompatible"):
            Configuration("MSC", "CRAd", 2_000_000, 4, "single")

    def test_random_sampling_valid(self):
        rng = np.random.default_rng(42)
        for _ in range(100):
            cfg = sample_random_configuration(rng)
            # Should not raise
            Configuration(cfg.cell_type, cfg.payload, cfg.dose_cells,
                          cfg.n_injection_sites, cfg.schedule)

    def test_configuration_space_summary(self):
        summary = configuration_space_summary()
        assert summary["n_cell_types"] == len(CELL_TYPES)
        assert summary["n_payloads"] == len(PAYLOADS)
        assert summary["n_schedules"] == len(SCHEDULES)
        assert summary["dose_range"] == [DOSE_MIN, DOSE_MAX]
        # Should have at least 10 valid discrete combinations
        assert summary["n_discrete_combos"] > 10


class TestOptimizer:
    def test_returns_ranked_configs(self, feature_store, survival_pred):
        configs = optimize_configurations(feature_store, survival_pred,
                                           n_top=3, n_candidates=20, seed=42)
        # 30 patients × 3 top configs = 90 rows
        assert len(configs) == 30 * 3
        assert "rank" in configs.columns
        # Each patient has ranks 1, 2, 3
        for _, group in configs.groupby("patient_id"):
            assert sorted(group["rank"].tolist()) == [1, 2, 3]

    def test_required_output_columns(self, feature_store, survival_pred):
        configs = optimize_configurations(feature_store, survival_pred,
                                           n_top=2, n_candidates=10, seed=42)
        for col in ["patient_id", "config_id", "rank", "score",
                    "cell_type", "payload", "dose_cells", "n_injection_sites", "schedule",
                    "predicted_survival_with_config", "predicted_tumor_kill_fraction",
                    "survival_gain_vs_stupp_months"]:
            assert col in configs.columns, f"Missing column: {col}"

    def test_score_descending_within_patient(self, feature_store, survival_pred):
        """Rank 1 should have higher score than rank 2."""
        configs = optimize_configurations(feature_store, survival_pred,
                                           n_top=3, n_candidates=20, seed=42)
        for pid, group in configs.groupby("patient_id"):
            scores = group.sort_values("rank")["score"].values
            assert scores[0] >= scores[1] >= scores[2]

    def test_personalization_variance(self, feature_store, survival_pred):
        """T4 threshold: optimizer should return at least 3 distinct top-5 configs."""
        configs = optimize_configurations(feature_store, survival_pred,
                                           n_top=5, n_candidates=50, seed=42)
        result = evaluate_personalization(configs)
        assert result["n_distinct_top5_cell_payload_combos"] >= 3
        assert result["passes_threshold"] is True

    def test_all_cell_payload_combos_explored(self, feature_store, survival_pred):
        """With enough candidates, all valid cell×payload combos should be reachable."""
        configs = optimize_configurations(feature_store, survival_pred,
                                           n_top=5, n_candidates=200, seed=42)
        combos = configs[["cell_type", "payload"]].drop_duplicates()
        # We have 10 valid cell×payload combos (3 cell types × 4 payloads - 2 incompatible)
        assert len(combos) >= 5  # at least half should be reached

    def test_seed_reproducible(self, feature_store, survival_pred):
        c1 = optimize_configurations(feature_store, survival_pred,
                                      n_top=3, n_candidates=20, seed=42)
        c2 = optimize_configurations(feature_store, survival_pred,
                                      n_top=3, n_candidates=20, seed=42)
        pd.testing.assert_frame_equal(c1, c2)


class TestEndToEndPipeline:
    """Full L0→L7 pipeline integration test."""

    def test_pipeline_runs_on_synthetic_data(self):
        from gliostrata.data.harmonize import harmonize
        from gliostrata.data.ingest import generate_synthetic_cohort
        from gliostrata.features.store import build_feature_store
        from gliostrata.handoff.protocol import generate_protocols
        from gliostrata.optimization.bayesian import optimize_configurations
        from gliostrata.prediction.survival import predict_survival
        from gliostrata.profiling.states import deconvolve_states
        from gliostrata.profiling.subtype import classify_subtypes
        from gliostrata.profiling.tme import deconvolve_tme

        # L0
        cohort = generate_synthetic_cohort(n=20, seed=42)
        # L1
        harmonized = harmonize(cohort)
        # L2
        subtypes = classify_subtypes(harmonized)
        states = deconvolve_states(harmonized)
        tme = deconvolve_tme(harmonized)
        # L3
        fs = build_feature_store(harmonized, subtypes, states, tme)
        # L5
        surv = predict_survival(fs)
        # L6
        configs = optimize_configurations(fs, surv, n_top=2, n_candidates=20, seed=42)
        # L7
        protocols = generate_protocols(configs, fs, surv, top_n=1)

        # Assertions
        assert len(fs) == 20
        assert len(configs) == 20 * 2
        assert len(protocols) == 20
        # Each protocol has the required structure
        for pid, proto in protocols.items():
            assert "recommended_config" in proto
            assert "predicted_outcomes" in proto
            assert "controls" in proto
            assert "completeness_checklist" in proto
            assert len(proto["controls"]) == 3  # 3 arms


class TestEnrichedProtocolSchema:
    """L7 — the protocol JSON must be a complete wet-lab handoff document."""

    REQUIRED_SECTIONS = [
        "recommended_config",
        "patient_summary",
        "cell_line_provenance",
        "payload_specification",
        "dose_specification",
        "schedule_detail",
        "pdo_panel_specification",
        "endpoints",
        "assay_plan",
        "statistical_analysis_plan",
        "regulatory_context",
        "randomization_and_blinding",
        "predicted_outcomes",
        "controls",
        "biosafety",
        "completeness_checklist",
    ]

    @staticmethod
    def _protocols(n=10, seed=42):
        from gliostrata.data.harmonize import harmonize
        from gliostrata.data.ingest import generate_synthetic_cohort
        from gliostrata.features.store import build_feature_store
        from gliostrata.handoff.protocol import generate_protocols
        from gliostrata.optimization.bayesian import optimize_configurations
        from gliostrata.prediction.survival import predict_survival
        from gliostrata.profiling.states import deconvolve_states
        from gliostrata.profiling.subtype import classify_subtypes
        from gliostrata.profiling.tme import deconvolve_tme

        harmonized = harmonize(generate_synthetic_cohort(n=n, seed=seed))
        fs = build_feature_store(
            harmonized,
            classify_subtypes(harmonized),
            deconvolve_states(harmonized),
            deconvolve_tme(harmonized),
        )
        surv = predict_survival(fs, seed=seed)
        configs = optimize_configurations(fs, surv, n_top=2, n_candidates=20, seed=seed)
        return generate_protocols(configs, fs, surv, top_n=1), fs

    def test_all_sections_present(self):
        protocols, _ = self._protocols()
        assert protocols
        for pid, proto in protocols.items():
            for section in self.REQUIRED_SECTIONS:
                assert section in proto, f"{pid} missing section {section!r}"
                assert proto[section], f"{pid} has empty section {section!r}"
            assert proto["ethical_notes"].strip()

    def test_patient_summary_is_faithful_projection(self):
        protocols, fs = self._protocols()
        required = [
            "kps", "extent_of_resection", "sex", "age", "subtype",
            "tme_immune", "tme_stromal", "tme_endothelial",
            "state_npc", "state_opc", "state_ac", "state_mes",
            "mgmt_methylated", "idh_mutated",
        ]
        for pid, proto in protocols.items():
            summary = proto["patient_summary"]
            assert len(summary) >= 15, f"{pid} patient_summary has only {len(summary)} fields"
            for field in required:
                assert field in summary, f"{pid} patient_summary missing {field!r}"
                assert summary[field] is not None, f"{pid} patient_summary[{field!r}] is None"
            # Values must match the feature store row, with one-hots decoded.
            row = fs.set_index("patient_id").loc[pid]
            assert summary["age"] == int(row["age"])
            assert summary["kps"] == int(row["kps"])
            assert summary["sex"] == ("M" if int(row["sex_male"]) == 1 else "F")
            assert summary["tme_immune"] == pytest.approx(float(row["tme_immune"]))
            assert summary["state_mes"] == pytest.approx(float(row["state_mes"]))
            assert summary["extent_of_resection"] in {"GTR", "STR", "biopsy"}

    def test_endpoints_section_is_populated(self):
        protocols, _ = self._protocols()
        for pid, proto in protocols.items():
            endpoints = proto["endpoints"]
            assert endpoints["primary_endpoint"].strip(), f"{pid} empty primary_endpoint"
            assert not endpoints["primary_endpoint"].startswith("placeholder")
            assert len(endpoints["secondary_endpoints"]) >= 2
            assert len(endpoints["measurement_methods"]) >= 1
            assert len(endpoints["timepoints"]) >= 1

    def test_dose_specification_splits_per_site(self):
        protocols, _ = self._protocols()
        for proto in protocols.values():
            dose = proto["dose_specification"]
            cfg = proto["recommended_config"]
            assert dose["total_cells"] == cfg["dose_cells"]
            assert dose["cells_per_site"] == cfg["dose_cells"] // cfg["n_injection_sites"]

    def test_provenance_and_payload_defaults_are_specific(self):
        protocols, _ = self._protocols()
        for proto in protocols.values():
            line = proto["cell_line_provenance"]["line_name"]
            construct = proto["payload_specification"]["construct_name"]
            assert not line.startswith("placeholder"), line
            assert not construct.startswith("placeholder"), construct
            assert proto["payload_specification"]["category"] == proto["recommended_config"]["payload"]

    def test_checklist_reflects_content_not_flags(self):
        from gliostrata.handoff.protocol import _completeness_checklist

        protocols, _ = self._protocols()
        for pid, proto in protocols.items():
            checklist = proto["completeness_checklist"]
            # Sections with real content are True.
            assert checklist["endpoints"] is True
            assert checklist["assay_plan"] is True
            assert checklist["statistical_analysis_plan"] is True
            assert checklist["controls"] is True
            assert checklist["randomization_and_blinding"] is True
            assert checklist["power_analysis"] is True
            # Sections still awaiting wet-lab specifics are not claimed as done.
            assert checklist["dose"] == "placeholder"
            assert checklist["cell_line_source"] == "placeholder"
            assert checklist["all_complete"] is False
            assert "dose" in checklist["unresolved_items"]
            assert "endpoints" not in checklist["unresolved_items"]

        # Removing a section must flip its checklist item to False, proving the
        # checklist verifies content rather than trusting a code path.
        proto = dict(next(iter(protocols.values())))
        assert _completeness_checklist(proto)["endpoints"] is True
        proto.pop("endpoints")
        assert _completeness_checklist(proto)["endpoints"] is False
        proto["endpoints"] = {"primary_endpoint": "", "secondary_endpoints": [],
                              "measurement_methods": [], "timepoints": []}
        assert _completeness_checklist(proto)["endpoints"] is False

    def test_protocol_generation_is_deterministic(self):
        first, _ = self._protocols()
        second, _ = self._protocols()
        assert first == second

    def test_protocol_json_round_trips(self, tmp_path):
        import json

        from gliostrata.handoff.protocol import protocols_to_json
        from gliostrata.utils.io import write_json

        protocols, _ = self._protocols()
        for pid, proto in protocols_to_json(protocols).items():
            path = write_json(proto, tmp_path / f"protocol_{pid}.json")
            text = path.read_text()
            loaded = json.loads(text)
            assert '\n  "patient_id"' in text  # indent=2, human readable
            for section in self.REQUIRED_SECTIONS:
                assert section in loaded

    def test_patient_summary_excludes_outcome_labels(self):
        """Outcome labels are training signal, never prospective-recommendation content."""
        from gliostrata.handoff.protocol import OUTCOME_DENYLIST, _patient_summary

        protocols, fs = self._protocols()
        fs_indexed = fs.set_index("patient_id")
        for pid, proto in protocols.items():
            # The source row really does carry the outcome labels...
            row = fs_indexed.loc[pid]
            assert "overall_survival_months" in row.index
            assert pd.notna(row["overall_survival_months"])
            assert "vital_status" in row.index
            assert "treatment_arm" in row.index
            # ...but the projection must not surface them.
            summary = proto["patient_summary"]
            for leaked in ("observed_overall_survival_months", "overall_survival_months",
                           "vital_status", "treatment_arm"):
                assert leaked not in summary, f"{pid} patient_summary leaked {leaked!r}"
            assert not (OUTCOME_DENYLIST & set(summary))
            # The patient profile fields are still there.
            assert summary["subtype"] == row["subtype"]
            assert summary["age"] == int(row["age"])

        # A denylisted column added to the feature store is still filtered, so the
        # trailing "surface remaining columns" pass can't reintroduce leakage.
        row = fs_indexed.iloc[0].copy()
        row["vital_status"] = 1
        row["treatment_arm"] = "Stupp"
        assert not (OUTCOME_DENYLIST & set(_patient_summary(row)))

    def test_biosafety_derives_from_payload_category(self):
        from gliostrata.handoff.protocol import PAYLOAD_BIOSAFETY, _biosafety

        # CRAd is a replication-competent adenovirus, not AAV or lentiviral.
        crad = _biosafety("CRAd")
        assert "adenovirus" in crad["viral_vector_class"].lower()
        assert "aav" not in crad["viral_vector_class"].lower()
        assert "lentiviral" not in crad["viral_vector_class"].lower()
        assert crad["bsl_level"] == "BSL-2"
        assert crad["payload_category"] == "CRAd"

        # A different category gets its own classification.
        car = _biosafety("CAR")
        assert "lentiviral" in car["viral_vector_class"].lower()
        assert "chimeric antigen receptor" in car["viral_vector_class"].lower()
        assert car["viral_vector_class"] != crad["viral_vector_class"]
        assert car["payload_category"] == "CAR"

        # Every category the optimizer can emit is classified.
        from gliostrata.optimization.config_space import PAYLOADS

        for payload in PAYLOADS:
            assert payload in PAYLOAD_BIOSAFETY, f"unclassified payload category {payload!r}"
            assert not _biosafety(payload)["viral_vector_class"].startswith("placeholder")

        # An unknown category placeholders rather than guessing.
        unknown = _biosafety("mRNA_LNP")
        assert unknown["viral_vector_class"] == "placeholder: vector class for mRNA_LNP required"
        assert unknown["bsl_level"] == "placeholder"

        # Generated protocols stay consistent with their own payload_specification.
        protocols, _ = self._protocols()
        for pid, proto in protocols.items():
            payload = proto["recommended_config"]["payload"]
            assert proto["biosafety"]["payload_category"] == payload
            assert proto["payload_specification"]["category"] == payload
            assert proto["biosafety"]["viral_vector_class"] == \
                PAYLOAD_BIOSAFETY[payload]["viral_vector_class"]
            assert proto["biosafety"]["containment_notes"] == \
                "Follow institutional biosafety committee protocols."

    def test_combo_payload_defers_bsl_level_to_ibc(self):
        """Co-delivery is not plain BSL-2 — the level is an IBC determination."""
        from gliostrata.handoff.protocol import _biosafety, _completeness_checklist

        combo = _biosafety("combo")
        assert combo["bsl_level"].startswith("placeholder:"), combo["bsl_level"]
        assert "adenovirus" in combo["viral_vector_class"].lower()
        assert "lentiviral" in combo["viral_vector_class"].lower()
        # The vector class itself is still a real classification, not a placeholder.
        assert not combo["viral_vector_class"].startswith("placeholder")

        # The other categories are unaffected.
        assert _biosafety("CRAd")["bsl_level"] == "BSL-2"
        assert _biosafety("prodrug_enzyme")["bsl_level"] == "BSL-2"
        assert _biosafety("CAR")["bsl_level"] == "BSL-2"

        # The checklist reports the unresolved level rather than claiming BSL-2.
        proto = {"biosafety": combo}
        assert _completeness_checklist(proto)["biosafety_level"] == "placeholder"
        assert _completeness_checklist(
            {"biosafety": _biosafety("CRAd")}
        )["biosafety_level"] is True

        # And it holds for generated protocols of each payload category.
        protocols, _ = self._protocols(n=20)
        by_payload = {p["recommended_config"]["payload"]: p for p in protocols.values()}
        for payload, proto in by_payload.items():
            checklist_item = proto["completeness_checklist"]["biosafety_level"]
            if payload == "combo":
                assert proto["biosafety"]["bsl_level"].startswith("placeholder:")
                assert checklist_item == "placeholder"
                assert "biosafety_level" in proto["completeness_checklist"]["unresolved_items"]
            else:
                assert proto["biosafety"]["bsl_level"] == "BSL-2"
                assert checklist_item is True
