#!/usr/bin/env bash

set -euo pipefail

repo_dir=$(cd "$(dirname "$0")/.." && pwd)
parallel_script="$repo_dir/do-lopnor-parallel.sh"
anisotropic_script="$repo_dir/do-lopnor-anistropic.sh"
big_script="$repo_dir/do-lopnor-big.sh"

assert_contains() {
    local file="$1"
    local text="$2"

    if ! rg --fixed-strings --quiet -- "$text" "$file"; then
        echo "Expected '$text' in $file" >&2
        exit 1
    fi
}

for script in "$parallel_script" "$anisotropic_script" "$big_script"; do
    test -x "$script"
    bash -n "$script"
    assert_contains "$script" 'source scripts/do-fundamentals.sh'
    assert_contains "$script" 'MODIDX=1'
    assert_contains "$script" 'RunSimulation || exit $?'
done

assert_contains "$parallel_script" '--workers=4'
assert_contains "$parallel_script" '--seed=0x5eedc0de12345678'
assert_contains "$parallel_script" 'NUMPHONS=10M'
if rg --fixed-strings --quiet -- '--scatter-horizontal=' "$parallel_script" || \
   rg --fixed-strings --quiet -- '--scatter-vertical=' "$parallel_script"; then
    echo "Parallel recipe must retain isotropic scattering" >&2
    exit 1
fi

assert_contains "$anisotropic_script" '--workers=1'
assert_contains "$anisotropic_script" '--seed=0x5eedc0de87654321'
assert_contains "$anisotropic_script" '--scatter-horizontal=0.25'
assert_contains "$anisotropic_script" '--scatter-vertical=1.25'

assert_contains "$big_script" 'NUMPHONS=10M'
if ! rg --quiet --regexp '^ADDITIONAL="--workers=1 --seed=0x5eedc0de12345678"$' "$big_script"; then
    echo "Expected a top-level serial ADDITIONAL assignment in $big_script" >&2
    exit 1
fi
assert_contains "$big_script" '## ___FIG_GEN_START___'
assert_contains "$big_script" '## ___FIG_GEN_END___'
