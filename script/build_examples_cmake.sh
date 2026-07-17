#!/usr/bin/env bash

set -o pipefail

if [[ "${1:-}" == "--one" ]]; then
    openfpm_dir=$2
    sweep_root=$3
    source_dir=$4
    name=${source_dir#example/}
    safe_name=${name//\//__}
    safe_name=${safe_name// /_}
    build_dir="$sweep_root/$safe_name"
    configure_log="$sweep_root/$safe_name.configure.log"
    build_log="$sweep_root/$safe_name.build.log"

    if cmake -S "$source_dir" -B "$build_dir" \
            -Dopenfpm_DIR="$openfpm_dir" >"$configure_log" 2>&1 && \
       cmake --build "$build_dir" -j "${OPENFPM_EXAMPLE_BUILD_JOBS:-2}" \
            >"$build_log" 2>&1; then
        printf 'PASS %s\n' "$name"
        exit 0
    fi

    printf 'FAIL %s\n' "$name"
    tail -n 25 "$configure_log"
    tail -n 50 "$build_log"
    exit 1
fi

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: $0 OPENFPM_CONFIG_DIR [PARALLEL_BUILDS]" >&2
    exit 2
fi

openfpm_dir=$1
parallel_builds=${2:-1}
script_path=$(cd "$(dirname "$0")" && pwd)/$(basename "$0")
sweep_root=$(mktemp -d "${TMPDIR:-/tmp}/openfpm-examples.XXXXXX")
export OPENFPM_EXAMPLE_BUILD_JOBS=${OPENFPM_EXAMPLE_BUILD_JOBS:-2}

find example -type f \( -name '*.cpp' -o -name '*.cu' \) -print \
    | xargs -n 1 dirname \
    | sort -u \
    | xargs -P "$parallel_builds" -I '{}' \
        bash "$script_path" --one "$openfpm_dir" "$sweep_root" '{}'
result=$?

printf 'Example build logs: %s\n' "$sweep_root"
exit "$result"
