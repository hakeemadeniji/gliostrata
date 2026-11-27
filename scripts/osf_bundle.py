"""Assemble the OSF deposit attachments from a clean, committed working tree.

Every file listed here is copied into ``--out``, keeping its repository path, and a
``MANIFEST.txt`` records the registered commit and each file's SHA-256. The script refuses
to run on a dirty tree, so the bundle is exactly the registered commit plus three
git-ignored result files, whose checksums are checked against the ones recorded in
``docs/project_record/06_evidence_reports.md``.

    python scripts/osf_bundle.py --out ../gliostrata_osf_bundle
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
import subprocess
from pathlib import Path

COMMITTED = [
    "docs/preregistration.md",
    "docs/osf_final_fields.md",
    "src/gliostrata/evaluation/thresholds.py",
    "data/partitions/tcga_gbm_clinical_v1.json",
    "data/partitions/tcga_gbm_clinical_v2.json",
    "data/partitions/cgga_external_v1.json",
    "docs/cgga_validation.md",
    "docs/t1_deconvolution_spec.md",
    "docs/t2_comparator_shortlist.md",
    "docs/t2_comparator_search/README.md",
    "docs/t2_comparator_search/search.json",
    "docs/t2_comparator_search/records.jsonl",
    "docs/t2_comparator_search/protocols/1_title_abstract.md",
    "docs/t2_comparator_search/protocols/2_full_text.md",
    "docs/t2_comparator_search/protocols/3_second_review.md",
    "docs/t2_comparator_search/protocols/4_gene_sets.md",
    "configs/t2_comparator/primary.json",
    "configs/t2_comparator/chosen.json",
    "configs/t2_comparator/gene_sets.json",
    "data/provenance/cgga_fpkm_sources.json",
    "data/provenance/cgga_terms.json",
    "data/provenance/neftel2019_sources.json",
    "data/access_ledger.jsonl",
    "docs/project_record/commit_map.md",
    "docs/osf_study_design_fields.md",
    "docs/osf_supplements/measured_variables.md",
    "docs/osf_supplements/indices.md",
    "docs/osf_supplements/statistical_models.md",
    "docs/dry_runs/t2_improvement_row_dry_run.json",
]
#: Git-ignored results the registration cites, with their recorded SHA-256.
RESULTS = {
    "data/results/cgga/cox_external_validation.json":
        "2f2c48a5c17197f1ddf29d50fbd4ee2afa81d0b497ef98877c9571766470595d",
    "data/results/cgga/brier18_external.json":
        "aca808f0be72b2b7eba2bf25b3ea93e81cb46a2f4308a7f5adf82c9d76e1e2a1",
    "data/results/t1/neftel_mixture_deconvolution.json":
        "a6921697b1b3af92ca1ba483fb6b74ea5d00dce231ceb862341cf6b8d60dc5e4",
}


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    dirty = subprocess.run(["git", "status", "--porcelain"], capture_output=True, text=True,
                           check=True).stdout.strip()
    if dirty:
        raise SystemExit(f"working tree is not clean; commit first:\n{dirty}")
    commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True,
                            check=True).stdout.strip()
    if args.out.exists() and any(args.out.iterdir()):
        raise SystemExit(f"{args.out} is not empty")
    lines = ["GlioStrata OSF deposit attachments", f"registered commit: {commit}", ""]
    for rel, expected in [(p, None) for p in COMMITTED] + list(RESULTS.items()):
        src = Path(rel)
        digest = _sha(src)
        if expected and digest != expected:
            raise SystemExit(f"{rel}: checksum {digest} differs from the recorded {expected}")
        dst = args.out / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dst)
        lines.append(f"{digest}  {rel}{'  (git-ignored result; checksum as recorded)' if expected else ''}")
    # The three field PDFs, rendered from the committed supplements.
    import importlib.util

    spec = importlib.util.spec_from_file_location("render_osf_pdfs",
                                                  Path(__file__).with_name("render_osf_pdfs.py"))
    render = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(render)
    pdf_dir = args.out / "osf_field_uploads"
    pdf_dir.mkdir(parents=True, exist_ok=True)
    lines.append("")
    for src, name in render.DOCS.items():
        render.render(Path(src), pdf_dir / name)
        lines.append(f"{_sha(pdf_dir / name)}  osf_field_uploads/{name}  (rendered from {src})")
    (args.out / "MANIFEST.txt").write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    print(f"{len(COMMITTED) + len(RESULTS)} files and {len(render.DOCS)} PDFs -> {args.out} "
          f"(commit {commit[:12]})")


if __name__ == "__main__":
    main()
