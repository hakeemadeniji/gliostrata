"""Tests for the reproducibility infrastructure: Snakemake, DVC, CI, .gitignore.

The point of this module is DRIFT DETECTION. The original `dvc.yaml` in this
repo was a stale shell: it invoked `gliostrata ingest --source tcga` (a stub
that prints "not yet implemented" and writes nothing) and declared outputs
under `data/raw/tcga/` that no code has ever produced. Nothing caught that,
because nothing cross-checked the pipeline definition against the code.

So the central test here is `TestDvcStagePathsMatchReality`: every path a DVC
stage names must be something the CLI or scripts/ genuinely produces, and
every command must be one that genuinely exists. If someone renames a CLI
flag, moves an output, or deletes a module, these tests fail instead of the
pipeline definition quietly becoming fiction.
"""
from __future__ import annotations

import ast
import json
import re
from pathlib import Path

import pytest

yaml = pytest.importorskip("yaml")

REPO_ROOT = Path(__file__).resolve().parents[1]

SNAKEFILE = REPO_ROOT / "workflow" / "Snakefile"
RULES_DIR = REPO_ROOT / "workflow" / "rules"
DVC_YAML = REPO_ROOT / "dvc.yaml"
DVCIGNORE = REPO_ROOT / ".dvcignore"
CI_YAML = REPO_ROOT / ".github" / "workflows" / "ci.yml"
GITIGNORE = REPO_ROOT / ".gitignore"
MAKEFILE = REPO_ROOT / "Makefile"
PARTITIONS_DIR = REPO_ROOT / "data" / "partitions"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_dvc() -> dict:
    return yaml.safe_load(DVC_YAML.read_text(encoding="utf-8"))


def _dvc_stages() -> dict:
    return _load_dvc()["stages"]


def _stage_cmds() -> dict[str, str]:
    """Stage name -> command, with YAML line folding collapsed to single spaces."""
    return {
        name: " ".join(stage["cmd"].split())
        for name, stage in _dvc_stages().items()
    }


def _cli_command_names() -> set[str]:
    """Click command names registered on the `gliostrata` group, read from source.

    Parsed from the AST rather than imported so that a broken decorator shows
    up as a test failure here rather than as a collection error.
    """
    cli_src = REPO_ROOT / "src" / "gliostrata" / "cli.py"
    tree = ast.parse(cli_src.read_text(encoding="utf-8"))
    names: set[str] = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.FunctionDef):
            continue
        for dec in node.decorator_list:
            # @cli.command() or @cli.command("explicit-name")
            if not isinstance(dec, ast.Call):
                continue
            func = dec.func
            if not (isinstance(func, ast.Attribute) and func.attr == "command"):
                continue
            if not (isinstance(func.value, ast.Name) and func.value.id == "cli"):
                continue
            if dec.args and isinstance(dec.args[0], ast.Constant):
                names.add(str(dec.args[0].value))
            else:
                names.add(node.name.replace("_", "-"))
    return names


def _cli_options_for(command: str) -> set[str]:
    """Long option strings accepted by one `gliostrata` subcommand."""
    cli_src = REPO_ROOT / "src" / "gliostrata" / "cli.py"
    tree = ast.parse(cli_src.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if not isinstance(node, ast.FunctionDef):
            continue
        explicit = None
        is_command = False
        options: set[str] = set()
        for dec in node.decorator_list:
            if not isinstance(dec, ast.Call):
                continue
            func = dec.func
            if not isinstance(func, ast.Attribute):
                continue
            if func.attr == "command" and isinstance(func.value, ast.Name):
                is_command = True
                if dec.args and isinstance(dec.args[0], ast.Constant):
                    explicit = str(dec.args[0].value)
            elif func.attr in {"option", "argument"}:
                for arg in dec.args:
                    if isinstance(arg, ast.Constant) and str(arg.value).startswith("--"):
                        options.add(str(arg.value))
        if is_command and (explicit or node.name.replace("_", "-")) == command:
            return options
    return set()


def _script_options_for(script: Path) -> set[str]:
    """Long argparse options accepted by a script in scripts/."""
    tree = ast.parse(script.read_text(encoding="utf-8"))
    options: set[str] = set()
    for node in ast.walk(tree):
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "add_argument"
        ):
            for arg in node.args:
                if isinstance(arg, ast.Constant) and str(arg.value).startswith("--"):
                    options.add(str(arg.value))
    return options


def _smk_text() -> str:
    """The Snakefile plus every included rule file, concatenated."""
    parts = [SNAKEFILE.read_text(encoding="utf-8")]
    parts += [p.read_text(encoding="utf-8") for p in sorted(RULES_DIR.glob("*.smk"))]
    return "\n".join(parts)


def _shell_commands() -> dict[str, str]:
    """Rule name -> its `shell:` command, with adjacent string literals joined.

    Snakemake shell commands are written as implicitly-concatenated Python
    string literals across several lines. Reading them as raw text and
    pattern-matching line by line gives false results, so reconstruct the
    actual command string instead.
    """
    commands: dict[str, str] = {}
    for path in [SNAKEFILE, *sorted(RULES_DIR.glob("*.smk"))]:
        text = path.read_text(encoding="utf-8")
        for block in re.split(r"^(?=rule\s+\w+\s*:)", text, flags=re.MULTILINE)[1:]:
            name = re.match(r"rule\s+(\w+)\s*:", block).group(1)
            m = re.search(r"^\s{4}shell:\s*\n((?:\s{8}.*\n?)+)", block,
                          flags=re.MULTILINE)
            if not m:
                continue
            literals = re.findall(r'"((?:[^"\\]|\\.)*)"', m.group(1))
            commands[name] = "".join(literals)
    return commands


def _smk_code() -> str:
    """Workflow source with comments AND docstrings removed.

    The docstrings legitimately discuss the GDC and the unwired modules, so a
    naive comment strip leaves false positives.
    """
    out = []
    for path in [SNAKEFILE, *sorted(RULES_DIR.glob("*.smk"))]:
        text = path.read_text(encoding="utf-8")
        # Drop triple-quoted blocks (module docstrings and rule docstrings).
        text = re.sub(r'"""(?:[^"]|"(?!""))*"""', "", text)
        out += [
            line for line in text.splitlines() if not line.lstrip().startswith("#")
        ]
    return "\n".join(out)


# ---------------------------------------------------------------------------
# Snakemake
# ---------------------------------------------------------------------------

class TestSnakefileParses:
    def test_snakefile_exists(self):
        assert SNAKEFILE.is_file(), "workflow/Snakefile is missing"

    def test_rule_files_exist(self):
        found = {p.name for p in RULES_DIR.glob("*.smk")}
        assert found, "workflow/rules/ contains no .smk files"
        # common.smk defines every path constant the other files use.
        assert "common.smk" in found

    def test_every_included_file_exists(self):
        """An `include:` naming a missing file breaks the workflow at parse time."""
        includes = re.findall(r'^include:\s*"([^"]+)"', SNAKEFILE.read_text(
            encoding="utf-8"), flags=re.MULTILINE)
        assert includes, "Snakefile includes no rule files"
        for rel in includes:
            assert (SNAKEFILE.parent / rel).is_file(), f"include: {rel} does not exist"

    def test_every_rule_file_is_included(self):
        """A rule file nobody includes is dead code that silently rots."""
        includes = set(re.findall(r'^include:\s*"rules/([^"]+)"', SNAKEFILE.read_text(
            encoding="utf-8"), flags=re.MULTILINE))
        on_disk = {p.name for p in RULES_DIR.glob("*.smk")}
        assert on_disk - includes == set(), (
            f"rule files present but never included: {sorted(on_disk - includes)}"
        )

    @pytest.mark.parametrize("path", [SNAKEFILE, *sorted(RULES_DIR.glob("*.smk"))],
                             ids=lambda p: p.name)
    def test_tokenizes_cleanly(self, path: Path):
        """Lexical validation that needs no snakemake install.

        Snakemake's grammar is a Python superset, so a .smk file is not valid
        Python — but it IS valid Python *tokens*. Tokenizing therefore catches
        unbalanced brackets, unterminated strings and inconsistent indentation,
        which is the overwhelming majority of how these files break.
        `test_snakemake_dry_run_resolves_the_dag` does the real parse when
        snakemake is installed.
        """
        import tokenize

        with open(path, "rb") as fh:
            list(tokenize.tokenize(fh.readline))

    def test_only_known_section_keywords_are_used(self):
        """A typo'd section (`ouput:`) is silently ignored by Snakemake."""
        known = {
            "input", "output", "params", "log", "shell", "run", "threads",
            "resources", "wildcard_constraints", "conda", "container",
            "benchmark", "priority", "retries", "message", "group", "cache",
            "default_target", "localrule", "envmodules", "notebook", "script",
            "wrapper", "template_engine", "shadow", "handover",
        }
        for path in [SNAKEFILE, *sorted(RULES_DIR.glob("*.smk"))]:
            lines = path.read_text(encoding="utf-8").splitlines()
            in_rule = False
            for line in lines:
                if re.match(r"^(rule|checkpoint)\s+\w+\s*:", line):
                    in_rule = True
                    continue
                if line and not line[0].isspace():
                    in_rule = False
                    continue
                if not in_rule:
                    continue
                m = re.match(r"^\s{4}(\w+)\s*:", line)
                if m:
                    assert m.group(1) in known, (
                        f"{path.name}: unknown rule section {m.group(1)!r} — "
                        f"Snakemake ignores sections it does not recognise"
                    )

    def test_every_rule_has_an_output_or_is_an_aggregator(self):
        """A rule with neither output nor input does nothing."""
        for path in [SNAKEFILE, *sorted(RULES_DIR.glob("*.smk"))]:
            text = path.read_text(encoding="utf-8")
            for block in re.split(r"^(?=rule\s+\w+\s*:)", text, flags=re.MULTILINE)[1:]:
                name = re.match(r"rule\s+(\w+)\s*:", block).group(1)
                assert re.search(r"^\s{4}(output|input|shell|run)\s*:", block,
                                 flags=re.MULTILINE), (
                    f"rule {name} in {path.name} declares nothing"
                )

    def test_snakemake_dry_run_resolves_the_dag(self, tmp_path):
        """The real thing: let snakemake parse the workflow and build the DAG.

        Skipped when snakemake is absent — it lives in the optional `workflow`
        extra. CI installs that extra, so this runs there.
        """
        pytest.importorskip(
            "snakemake",
            reason="snakemake is in the optional `workflow` extra; "
                   "install with `pip install -e \".[workflow]\"`",
        )
        import subprocess
        import sys

        proc = subprocess.run(
            [
                sys.executable, "-m", "snakemake", "-n", "--cores", "1",
                "--snakefile", str(SNAKEFILE),
                "--config", "source=synthetic", "cohort=pytestfixture",
                "n_synthetic=6", "with_expression=false",
                "verify_partition=false",
                "raw_root=.test/pytest/raw", "interim_root=.test/pytest/interim",
                "processed_root=.test/pytest/processed",
                "results_root=.test/pytest/results",
            ],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
        )
        assert proc.returncode == 0, (
            f"snakemake dry run failed:\n{proc.stdout}\n{proc.stderr}"
        )
        combined = proc.stdout + proc.stderr
        for rule in ("ingest_synthetic", "harmonize", "profile",
                     "build_feature_store", "optimize", "generate_protocol"):
            assert rule in combined, f"the DAG does not contain `{rule}`"

    def test_declares_the_real_l0_to_l7_rules(self):
        """The DAG must cover the layers the code actually implements."""
        text = _smk_text()
        rules = set(re.findall(r"^rule\s+(\w+)\s*:", text, flags=re.MULTILINE))
        for expected in (
            "ingest_tcga_clinical",     # L0 clinical
            "ingest_tcga_expression",   # L0 expression
            "ingest_synthetic",         # L0 offline fixture
            "harmonize",                # L1
            "profile",                  # L2
            "build_feature_store",      # L3
            "optimize",                 # L5 + L6
            "generate_protocol",        # L7
            "verify_partition",         # guardrail
            "all",
        ):
            assert expected in rules, f"workflow is missing `rule {expected}`"

    def test_config_driven_paths(self):
        """Nothing may be hardcoded, or the fixture run is impossible."""
        common = (RULES_DIR / "common.smk").read_text(encoding="utf-8")
        for key in ("source", "cohort", "raw_root", "interim_root",
                    "processed_root", "results_root", "partitions_root"):
            assert f'config.setdefault("{key}"' in common, (
                f"config key {key!r} has no default in common.smk"
            )

    def test_stub_layers_are_marked_as_stubs(self):
        """A stubbed layer must say so, not pass itself off as real science."""
        for name in ("harmonize", "profile", "optimize", "handoff"):
            text = (RULES_DIR / f"{name}.smk").read_text(encoding="utf-8")
            assert "STUB" in text, (
                f"workflow/rules/{name}.smk wraps a stubbed layer but never says so"
            )

    def test_only_the_ingest_rules_touch_the_network(self):
        """CI relies on `source=synthetic` being offline. Guard that property."""
        offline_rules = {
            "harmonize", "profile", "build_feature_store", "optimize",
            "generate_protocol", "verify_partition", "ingest_synthetic", "all",
        }
        for rule, cmd in _shell_commands().items():
            if rule not in offline_rules:
                continue
            for script in ("download_tcga_clinical", "download_tcga_expression",
                           "download_cgga"):
                assert script not in cmd, (
                    f"rule {rule} invokes {script}; network access must stay "
                    f"confined to the ingest_tcga_* rules"
                )

    def test_network_rules_are_confined_to_ingest_smk(self):
        network_rules = {
            rule
            for rule, cmd in _shell_commands().items()
            if "download_tcga" in cmd or "download_cgga" in cmd
        }
        in_ingest = set(
            re.findall(r"^rule\s+(\w+)\s*:",
                       (RULES_DIR / "ingest.smk").read_text(encoding="utf-8"),
                       flags=re.MULTILINE)
        )
        assert network_rules <= in_ingest, (
            f"network-touching rules outside ingest.smk: "
            f"{sorted(network_rules - in_ingest)}"
        )

    def test_no_rule_writes_a_partition_manifest(self):
        """Freezing a partition is a human, pre-registered act — never a build step."""
        commands = _shell_commands()
        freezers = {r: c for r, c in commands.items() if "freeze-partition" in c}
        assert freezers, "the partition guardrail is not wired into the workflow"
        for rule, cmd in freezers.items():
            assert "--verify" in cmd, (
                f"rule {rule} invokes `freeze-partition` without --verify; the "
                f"workflow must never create or overwrite a manifest"
            )
            assert "--force" not in cmd
            assert "--no-seal" not in cmd

    def test_unwired_modules_are_not_presented_as_integrated(self):
        """Real modules built alongside stubs must not appear as pipeline steps.

        src/gliostrata/profiling/verhaak.py, profiling/deconvolve.py,
        prediction/cox.py and optimization/bayes_opt.py are real, but the CLI
        still calls the stubs beside them. A rule invoking them would claim an
        integration that does not exist. (The comments and docstrings SAY they
        are unwired, which is why this checks executable content only.)
        """
        for rule, cmd in _shell_commands().items():
            for module in ("verhaak", "bayes_opt", "deconvolve", "cox"):
                assert module not in cmd, (
                    f"rule {rule} invokes {module}, which is not wired into the "
                    f"CLI pipeline; the workflow must not imply otherwise"
                )

    def test_the_unwired_modules_really_are_unwired(self):
        """The premise of the test above. If the CLI gains them, revisit the DAG."""
        cli = (REPO_ROOT / "src" / "gliostrata" / "cli.py").read_text(encoding="utf-8")
        still_unwired = [
            m for m in ("verhaak", "bayes_opt", "deconvolve", "prediction.cox")
            if m not in cli
        ]
        assert still_unwired, (
            "every module previously documented as unwired is now imported by "
            "cli.py — update workflow/Snakefile's honesty notes and this test"
        )


# ---------------------------------------------------------------------------
# DVC — the drift check that matters
# ---------------------------------------------------------------------------

class TestDvcYamlParses:
    def test_dvc_yaml_is_valid_yaml_with_stages(self):
        doc = _load_dvc()
        assert isinstance(doc, dict)
        assert "stages" in doc
        assert doc["stages"], "dvc.yaml declares no stages"

    def test_every_stage_has_cmd_deps_outs(self):
        for name, stage in _dvc_stages().items():
            assert "cmd" in stage, f"stage {name} has no cmd"
            assert stage.get("deps"), f"stage {name} has no deps"
            assert stage.get("outs"), f"stage {name} has no outs"

    def test_no_duplicate_outs_across_stages(self):
        """Two stages claiming the same out is a DVC error at repro time."""
        seen: dict[str, str] = {}
        for name, stage in _dvc_stages().items():
            for out in stage["outs"]:
                assert out not in seen, (
                    f"{out} is an out of both {seen[out]} and {name}"
                )
                seen[out] = name

    def test_stage_graph_is_connected_and_acyclic(self):
        """Every stage's non-source deps must be produced by another stage."""
        stages = _dvc_stages()
        producer = {
            out: name for name, stage in stages.items() for out in stage["outs"]
        }
        for name, stage in stages.items():
            for dep in stage["deps"]:
                if dep.startswith(("src/", "scripts/", "configs/")):
                    continue  # a source-code dep, not a pipeline artifact
                if dep.startswith("data/partitions/"):
                    continue  # a frozen manifest in git: read, never produced
                assert dep in producer, (
                    f"stage {name} depends on {dep}, which no stage produces and "
                    f"which is not source code — this is exactly the staleness "
                    f"that made the previous dvc.yaml fiction"
                )
                assert producer[dep] != name, f"stage {name} depends on its own out"

    def test_artifacts_agree_with_snakemake(self):
        """DVC and Snakemake must describe the same artifacts.

        Snakemake builds its paths from config, so the literal
        `data/interim/tcga_gbm_harmonized.parquet` never appears there —
        `f"{INTERIM}/{COHORT}_harmonized.parquet"` does. Compare on the
        cohort-independent tail.
        """
        smk = _smk_text()
        common = (RULES_DIR / "common.smk").read_text(encoding="utf-8")
        cohort = re.search(r'config\.setdefault\("cohort",\s*"([^"]+)"\)', common)
        assert cohort, "common.smk has no default cohort name"
        prefix = cohort.group(1)

        for name, stage in _dvc_stages().items():
            for out in stage["outs"]:
                tail = Path(out).name
                if not tail:
                    tail = Path(out.rstrip("/")).name  # directory out
                # Strip the cohort token however it appears in the path.
                tail = tail.replace(f"{prefix}_", "")
                assert tail in smk, (
                    f"dvc.yaml stage {name} produces {out}, but no Snakemake rule "
                    f"mentions {tail!r} — the two pipeline definitions have drifted"
                )

    def test_every_dvc_stage_has_a_snakemake_rule(self):
        """Stage-level agreement, not just artifact-level."""
        rules = set(re.findall(r"^rule\s+(\w+)\s*:", _smk_text(), flags=re.MULTILINE))
        for name in _dvc_stages():
            assert name in rules, (
                f"dvc.yaml declares stage {name!r} with no matching Snakemake "
                f"rule; the two definitions have drifted"
            )


class TestDvcStagePathsMatchReality:
    """Every path and command in dvc.yaml must correspond to real code.

    This is the test that would have caught the original stale dvc.yaml.
    """

    def test_source_code_deps_exist_on_disk(self):
        for name, stage in _dvc_stages().items():
            for dep in stage["deps"]:
                if not dep.startswith(("src/", "scripts/", "configs/", "data/partitions/")):
                    continue
                assert (REPO_ROOT / dep).exists(), (
                    f"stage {name} depends on {dep}, which does not exist"
                )

    def test_every_command_is_a_real_cli_command_or_script(self):
        cli_names = _cli_command_names()
        for name, cmd in _stage_cmds().items():
            tokens = cmd.split()
            if tokens[0] == "gliostrata":
                assert tokens[1] in cli_names, (
                    f"stage {name} runs `gliostrata {tokens[1]}`, which is not a "
                    f"registered CLI command. Registered: {sorted(cli_names)}"
                )
            elif tokens[0] == "python":
                script = tokens[1]
                assert (REPO_ROOT / script).is_file(), (
                    f"stage {name} runs {script}, which does not exist"
                )
            else:
                pytest.fail(f"stage {name} has an unrecognised command: {cmd!r}")

    def test_every_flag_is_accepted_by_its_command(self):
        """A renamed flag must break this test, not `dvc repro` six months later."""
        for name, cmd in _stage_cmds().items():
            tokens = cmd.split()
            used = {t for t in tokens if t.startswith("--")}
            if tokens[0] == "gliostrata":
                accepted = _cli_options_for(tokens[1])
                label = f"gliostrata {tokens[1]}"
            else:
                accepted = _script_options_for(REPO_ROOT / tokens[1])
                label = tokens[1]
            unknown = used - accepted
            assert not unknown, (
                f"stage {name} passes {sorted(unknown)} to {label}, which does not "
                f"accept them. Accepted: {sorted(accepted)}"
            )

    def test_no_stage_invokes_the_stub_ingest_command(self):
        """The precise bug in the original dvc.yaml.

        `gliostrata ingest --source tcga` is a STUB: it prints a "not yet
        implemented" message and writes nothing, so the stage's declared outs
        could never appear. Real TCGA ingestion is scripts/download_tcga_*.py.
        """
        for name, cmd in _stage_cmds().items():
            if "gliostrata ingest" in cmd:
                assert "--source synthetic" in cmd, (
                    f"stage {name} runs `{cmd}`. `gliostrata ingest` only works "
                    f"with --source synthetic; for TCGA use "
                    f"scripts/download_tcga_clinical.py."
                )

    def test_declared_outs_match_what_the_code_writes(self):
        """Cross-check every out filename against the producing code's source."""
        for name, stage in _dvc_stages().items():
            cmd = _stage_cmds()[name]
            tokens = cmd.split()
            if tokens[0] == "python":
                producer_src = (REPO_ROOT / tokens[1]).read_text(encoding="utf-8")
            else:
                producer_src = (
                    REPO_ROOT / "src" / "gliostrata" / "cli.py"
                ).read_text(encoding="utf-8")
            for out in stage["outs"]:
                basename = Path(out).name
                if not basename:
                    continue
                # Either the command line names the path explicitly (the CLI
                # takes it as --out), or the producing script hardcodes the
                # filename it writes.
                assert out in cmd or basename in producer_src, (
                    f"stage {name} declares out {out}, but neither the command "
                    f"nor {tokens[1] if tokens[0] == 'python' else 'cli.py'} "
                    f"mentions {basename}"
                )

    def test_expression_download_cache_is_never_an_out(self):
        """The 1.2 GB per-sample TSV cache must stay out of DVC entirely."""
        for name, stage in _dvc_stages().items():
            for path in list(stage["outs"]) + list(stage["deps"]):
                assert "_cache" not in path, (
                    f"stage {name} references {path}; the download cache is "
                    f"~1.2 GB of regenerable scratch and must not be tracked"
                )

    def test_partition_manifests_are_never_a_dvc_out(self):
        """Frozen manifests live in git so they are diff-reviewable."""
        for name, stage in _dvc_stages().items():
            for out in stage["outs"]:
                assert "data/partitions" not in out, (
                    f"stage {name} declares {out} as an out. Partition manifests "
                    f"are version-controlled in git, and no pipeline stage may "
                    f"be able to regenerate one."
                )

    def test_no_stage_freezes_a_partition(self):
        for name, cmd in _stage_cmds().items():
            assert "freeze-partition" not in cmd or "--verify" in cmd, (
                f"stage {name} runs `freeze-partition` without --verify"
            )


class TestDvcignore:
    def test_exists_and_excludes_the_download_cache(self):
        assert DVCIGNORE.is_file(), ".dvcignore is missing"
        text = DVCIGNORE.read_text(encoding="utf-8")
        assert "_cache" in text, ".dvcignore must exclude the 1.2 GB download cache"

    def test_excludes_partition_manifests(self):
        text = DVCIGNORE.read_text(encoding="utf-8")
        assert "data/partitions" in text, (
            ".dvcignore must exclude data/partitions/ — those manifests belong to git"
        )


# ---------------------------------------------------------------------------
# GitHub Actions CI
# ---------------------------------------------------------------------------

class TestCiWorkflow:
    def test_ci_yaml_is_valid_yaml(self):
        assert CI_YAML.is_file(), ".github/workflows/ci.yml is missing"
        doc = yaml.safe_load(CI_YAML.read_text(encoding="utf-8"))
        assert isinstance(doc, dict)
        assert doc.get("jobs"), "ci.yml declares no jobs"

    def test_runs_on_push_and_pull_request(self):
        doc = yaml.safe_load(CI_YAML.read_text(encoding="utf-8"))
        # PyYAML parses the bare key `on:` as the boolean True.
        triggers = doc.get("on", doc.get(True))
        assert triggers is not None, "ci.yml has no trigger block"
        keys = set(triggers) if isinstance(triggers, dict) else set(triggers)
        assert "push" in keys
        assert "pull_request" in keys

    def test_every_job_pins_an_exact_python_version(self):
        """A floating Python is how a green CI stops matching the dev environment."""
        doc = yaml.safe_load(CI_YAML.read_text(encoding="utf-8"))
        pyproject = (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
        floor = re.search(r'requires-python\s*=\s*">=\s*(\d+\.\d+)"', pyproject)
        assert floor, "could not read requires-python from pyproject.toml"
        by_job: dict[str, set[str]] = {}
        for job_name, job in doc["jobs"].items():
            versions = {
                str(step["with"]["python-version"])
                for step in job.get("steps", [])
                if isinstance(step, dict)
                and "setup-python" in str(step.get("uses", ""))
                and isinstance(step.get("with"), dict)
                and "python-version" in step["with"]
            }
            assert versions, f"job {job_name} never pins a Python version"
            for v in versions:
                assert re.fullmatch(r"\d+\.\d+(\.\d+)?", v), (
                    f"job {job_name} pins python-version {v!r}; it must be an "
                    f"exact version, not a range or expression"
                )
            by_job[job_name] = versions
        # One project, one Python. The single exception is deliberate: the floor
        # job runs the same suite on the oldest Python pyproject accepts, with
        # the lowest dependency versions it accepts (requirements-floor.txt).
        main = {v for name, versions in by_job.items() if name != "test-floor" for v in versions}
        assert len(main) == 1, f"CI pins inconsistent Python versions: {main}"
        assert by_job.get("test-floor") == {floor.group(1)}, (
            f"test-floor must pin exactly the requires-python floor {floor.group(1)}; "
            f"it pins {by_job.get('test-floor')}"
        )

    def test_pinned_python_matches_the_project_floor(self):
        doc = yaml.safe_load(CI_YAML.read_text(encoding="utf-8"))
        version = next(
            str(step["with"]["python-version"])
            for job in doc["jobs"].values()
            for step in job.get("steps", [])
            if isinstance(step, dict)
            and "setup-python" in str(step.get("uses", ""))
            and isinstance(step.get("with"), dict)
            and "python-version" in step["with"]
        )
        pyproject = (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
        floor = re.search(r'requires-python\s*=\s*">=\s*(\d+)\.(\d+)"', pyproject)
        assert floor, "could not read requires-python from pyproject.toml"
        major, minor = (int(p) for p in version.split(".")[:2])
        assert (major, minor) >= (int(floor.group(1)), int(floor.group(2))), (
            f"CI pins Python {version}, below the project floor "
            f"{floor.group(1)}.{floor.group(2)}"
        )

    def test_runs_the_test_suite(self):
        text = CI_YAML.read_text(encoding="utf-8")
        assert "pytest" in text, "CI never runs the test suite"

    def test_has_a_lint_job(self):
        doc = yaml.safe_load(CI_YAML.read_text(encoding="utf-8"))
        text = CI_YAML.read_text(encoding="utf-8")
        assert "ruff" in text, "CI has no lint step"
        assert any("lint" in name for name in doc["jobs"]), "CI has no lint job"

    def test_has_a_blocking_lint_step(self):
        """An advisory-only lint job is theatre. At least one must be blocking."""
        doc = yaml.safe_load(CI_YAML.read_text(encoding="utf-8"))
        blocking = [
            step
            for job in doc["jobs"].values()
            for step in job.get("steps", [])
            if isinstance(step, dict)
            and "ruff check" in str(step.get("run", ""))
            and not step.get("continue-on-error", False)
        ]
        assert blocking, "every ruff step is continue-on-error; nothing gates merges"

    def test_verifies_the_current_partition_manifest(self):
        """Partition drift must fail CI — a real scientific guardrail."""
        text = CI_YAML.read_text(encoding="utf-8")
        assert "freeze-partition" in text, "CI has no partition-drift job"
        assert "--verify" in text
        assert "tcga_gbm_clinical_v2.json" in text, (
            "CI must verify the CURRENT partition manifest (v2), not a "
            "superseded one"
        )

    def test_the_manifest_ci_verifies_actually_exists(self):
        """Catches CI silently pointing at a manifest that was renamed or removed."""
        text = CI_YAML.read_text(encoding="utf-8")
        for manifest in set(re.findall(r"data/partitions/[\w.]+\.json", text)):
            assert (REPO_ROOT / manifest).is_file(), (
                f"ci.yml references {manifest}, which does not exist"
            )

    def test_never_runs_the_network_ingest(self):
        """CI must not reach the GDC or CGGA, or download gigabytes."""
        doc = yaml.safe_load(CI_YAML.read_text(encoding="utf-8"))
        runs = [
            str(step.get("run", ""))
            for job in doc["jobs"].values()
            for step in job.get("steps", [])
            if isinstance(step, dict)
        ]
        for run in runs:
            # Strip shell comments before matching.
            code = "\n".join(
                line for line in run.splitlines() if not line.lstrip().startswith("#")
            )
            for script in ("download_tcga_clinical", "download_tcga_expression",
                           "download_cgga"):
                assert script not in code, (
                    f"a CI step invokes {script}; CI must never hit the network"
                )
            # A real (non-dry-run) snakemake invocation must be the offline one.
            for line in code.splitlines():
                if "snakemake" in line and " -n" not in line and "--dry-run" not in line:
                    assert "source=synthetic" in code, (
                        "CI runs snakemake for real without source=synthetic; that "
                        "would trigger the GDC ingest rules"
                    )


# ---------------------------------------------------------------------------
# .gitignore
# ---------------------------------------------------------------------------

class TestGitignore:
    DATA_DIRS = ("data/raw", "data/interim", "data/processed",
                 "data/external", "data/manual", "data/results")

    def test_excludes_every_large_data_dir(self):
        text = GITIGNORE.read_text(encoding="utf-8")
        patterns = {
            line.strip()
            for line in text.splitlines()
            if line.strip() and not line.strip().startswith("#")
        }
        for d in self.DATA_DIRS:
            assert f"{d}/**" in patterns or f"{d}/" in patterns, (
                f"{d} is no longer gitignored; large data files could be committed"
            )

    def test_models_are_excluded(self):
        patterns = {
            line.strip() for line in GITIGNORE.read_text(encoding="utf-8").splitlines()
        }
        assert "models/**" in patterns

    def test_partitions_are_NOT_excluded(self):
        """Frozen manifests are version-controlled on purpose.

        This is the inverse of every other rule in this class, and it is
        deliberate: a sealed train/test split must be reviewable in a diff and
        must travel with the code. If someone "tidies up" by adding
        data/partitions to .gitignore, this fails.
        """
        patterns = {
            line.strip()
            for line in GITIGNORE.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.strip().startswith("#")
        }
        for bad in ("data/partitions", "data/partitions/", "data/partitions/**",
                    "data/partitions/*.json", "data/**"):
            assert bad not in patterns, (
                f".gitignore contains {bad!r}, which would hide the frozen "
                f"partition manifests. They must stay version-controlled."
            )

    def test_at_least_one_partition_manifest_is_present(self):
        manifests = list(PARTITIONS_DIR.glob("*.json"))
        assert manifests, "data/partitions/ holds no manifests"

    def test_committed_manifests_are_sealed_and_self_consistent(self):
        """A manifest in git must be sealed, or it is not a commitment."""
        for path in sorted(PARTITIONS_DIR.glob("*.json")):
            manifest = json.loads(path.read_text(encoding="utf-8"))
            assert manifest.get("sealed") is True, f"{path.name} is not sealed"
            assert manifest["n_train"] == len(manifest["train_ids"])
            assert manifest["n_test"] == len(manifest["test_ids"])
            assert not set(manifest["train_ids"]) & set(manifest["test_ids"]), (
                f"{path.name} has ids in both train and test"
            )

    def test_snakemake_and_fixture_scratch_are_excluded(self):
        patterns = {
            line.strip() for line in GITIGNORE.read_text(encoding="utf-8").splitlines()
        }
        assert ".snakemake/" in patterns
        assert ".test/" in patterns, (
            "the fixture pipeline output root .test/ must be gitignored"
        )


# ---------------------------------------------------------------------------
# Makefile
# ---------------------------------------------------------------------------

class TestMakefile:
    def test_has_the_documented_targets(self):
        text = MAKEFILE.read_text(encoding="utf-8")
        targets = set(re.findall(r"^([a-zA-Z][\w-]*):", text, flags=re.MULTILINE))
        for expected in ("install", "test", "lint", "pipeline", "clean"):
            assert expected in targets, f"Makefile lost the `{expected}` target"

    def test_preserved_the_original_targets(self):
        """The Makefile was extended, not replaced."""
        text = MAKEFILE.read_text(encoding="utf-8")
        targets = set(re.findall(r"^([a-zA-Z][\w-]*):", text, flags=re.MULTILINE))
        for expected in ("help", "dev-install", "format", "synthetic"):
            assert expected in targets, (
                f"Makefile dropped the pre-existing `{expected}` target"
            )

    def test_every_phony_target_is_defined(self):
        text = MAKEFILE.read_text(encoding="utf-8")
        phony_block = re.search(
            r"^\.PHONY:((?:[^\n\\]|\\\n)*)", text, flags=re.MULTILINE
        )
        assert phony_block, "Makefile has no .PHONY declaration"
        declared = set(phony_block.group(1).replace("\\\n", " ").split())
        targets = set(re.findall(r"^([a-zA-Z][\w-]*):", text, flags=re.MULTILINE))
        assert declared - targets == set(), (
            f".PHONY names undefined targets: {sorted(declared - targets)}"
        )

    def test_verify_partition_points_at_the_current_manifest(self):
        text = MAKEFILE.read_text(encoding="utf-8")
        assert "tcga_gbm_clinical_v2.json" in text, (
            "make verify-partition must target the CURRENT manifest (v2)"
        )

    def test_no_target_freezes_a_partition(self):
        text = MAKEFILE.read_text(encoding="utf-8")
        code = "\n".join(
            line for line in text.splitlines() if not line.lstrip().startswith("#")
        )
        for call in re.findall(r"freeze-partition[^\n]*", code):
            assert "--verify" in call, (
                f"Makefile recipe `{call.strip()}` could overwrite a sealed manifest"
            )


# ---------------------------------------------------------------------------
# Documentation
# ---------------------------------------------------------------------------

class TestReproducibilityDocs:
    DOC = REPO_ROOT / "docs" / "reproducibility.md"

    def test_exists(self):
        assert self.DOC.is_file(), "docs/reproducibility.md is missing"

    def test_states_what_is_not_reproducible(self):
        """The honesty requirement, enforced."""
        text = self.DOC.read_text(encoding="utf-8")
        assert "NOT reproducible" in text or "not reproducible" in text
        # Each stubbed layer must be named.
        for layer in ("L1", "L2", "L5", "L6", "L7"):
            assert layer in text, f"the docs never mention {layer}"

    def test_explains_the_frozen_partitions(self):
        text = self.DOC.read_text(encoding="utf-8")
        assert "freeze-partition" in text
        assert "regenerat" in text.lower(), (
            "the docs must explain why partitions must not be regenerated"
        )

    def test_documents_the_dvc_decision(self):
        text = self.DOC.read_text(encoding="utf-8")
        assert "dvc init" in text, "the docs must state the dvc init decision"
        assert "_cache" in text, "the docs must address the 1.2 GB download cache"
