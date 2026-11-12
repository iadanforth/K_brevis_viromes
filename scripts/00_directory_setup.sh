#!/bin/bash
set -euo pipefail
# set up folders in data directory
mkdir -p "$HOME/scr10/data/kbrevis"
mkdir -p "$HOME/scr10/data/kbrevis/raw"
mkdir -p "$HOME/scr10/data/kbrevis/clean"

# set up folders in main directory
mkdir -p "$HOME/K_brevis_viromes/logs"
mkdir -p "$HOME/K_brevis_viromes/output"
mkdir -p "$HOME/K_brevis_viromes/scripts"
mkdir -p "$HOME/K_brevis_viromes/metadata"
mkdir -p "$HOME/K_brevis_viromes/envs"
