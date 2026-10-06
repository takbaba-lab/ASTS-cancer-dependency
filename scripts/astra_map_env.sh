#!/bin/bash

# ============================================================
# ASTRA-Drug micromamba environment configuration
# ============================================================
#
# Purpose:
#   Define the micromamba executable, root prefix, and astra-map
#   environment directory used by the ASTRA-Drug job scripts.
#
# These values can be overridden by environment variables before
# job submission when micromamba is installed in a non-default
# location.
# ============================================================

MICROMAMBA_BIN="${MICROMAMBA_BIN:-$(command -v micromamba || true)}"

export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-${HOME}/micromamba}"

ASTRA_ENV_DIR="${ASTRA_ENV_DIR:-${MAMBA_ROOT_PREFIX}/envs/astra-map}"
