#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONTAINER_DIR="$(dirname "$SCRIPT_DIR")"

echo "Building Apptainer containers into $CONTAINER_DIR..."

for def_file in "$SCRIPT_DIR"/*.def; do
    base_name="$(basename "$def_file" .def)"
    sif_path="$CONTAINER_DIR/${base_name}.sif"
    
    echo "=== Building ${base_name}.sif ==="
    apptainer build "$sif_path" "$def_file"
done

echo "All containers built successfully."
