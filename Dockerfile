# GlioStrata test image. Built from the lock file, so no version is resolved at
# build time: requirements-lock.txt is the environment that produced the results
# on record (Python 3.14.3), the same install the CI `test` job performs.
#
# No data is copied in. .dockerignore drops data/ except the small git-tracked
# manifests the tests read (data/partitions, data/provenance); nothing under
# data/raw or data/processed, and so no cohort table, enters the image. The
# tests that read the real TCGA data skip, as they do in CI.
#
#   docker build -t gliostrata .
#   docker run --rm gliostrata                 # the test suite
#   docker run --rm gliostrata gliostrata --help
FROM python:3.14-slim

ENV PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_NO_CACHE_DIR=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

RUN useradd --create-home --uid 1000 gliostrata
WORKDIR /app

# Dependencies first, so this layer is cached until the lock file changes.
COPY requirements-lock.txt ./
RUN pip install -r requirements-lock.txt

# The package itself, editable and with --no-deps, exactly as CI does. pyproject
# declares readme = "README.md", so it is needed for the build. The other files
# copied here (workflow, Makefile, dvc.yaml, .gitignore, ...) are read by tests.
COPY pyproject.toml README.md LICENSE CHANGELOG.md Makefile dvc.yaml environment.yml \
     .gitignore .dvcignore requirements.txt requirements-dev.txt requirements-floor.txt ./
COPY .github ./.github
COPY src ./src
COPY configs ./configs
COPY scripts ./scripts
COPY tests ./tests
COPY workflow ./workflow
COPY data ./data
COPY docs ./docs
RUN pip install --no-deps -e . && chown -R gliostrata:gliostrata /app

USER gliostrata

CMD ["python", "-m", "pytest", "-q"]
