#!/bin/bash
eval "$(conda shell.bash hook)"

set -e

PY_ENV_PATH=$1
CONFIG_PATH=$2

# Run ENACT pipeline
conda activate $PY_ENV_PATH
cd "${0%/*}" && python -m src.enact.pipeline --configs_path "$CONFIG_PATH"
