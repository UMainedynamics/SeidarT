#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DOCS_DIR="${REPO_ROOT}/docs"
SOURCE_DIR="${DOCS_DIR}/source"
PACKAGE_DIR="${REPO_ROOT}/src/seidart"

RUN=()
if [[ "${CONDA_DEFAULT_ENV:-}" != "documentation" ]] && command -v conda >/dev/null 2>&1; then
    RUN=(conda run -n documentation)
fi

cp "${REPO_ROOT}/README.rst" "${SOURCE_DIR}/README.rst"

"${RUN[@]}" sphinx-apidoc -o "${SOURCE_DIR}" "${PACKAGE_DIR}" -e -f
"${RUN[@]}" make -C "${DOCS_DIR}" clean html
