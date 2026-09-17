#!/bin/bash

set -euo pipefail

repo_root=$(cd "$(dirname "$0")/.." && pwd)
scratch_dir=$(mktemp -d)
trap 'rm -rf "$scratch_dir"' EXIT

cd "$scratch_dir"
octave --no-gui -qf <<EOF
addpath("$repo_root/vis/seisplot");
figinit(5.0, 3.75, "paperunits", "inches", "visible", "off");
plot(1:3);
colorbar("eastoutside", "linewidth", 0.5, "fontsize", 9);
h = papertext(0.01, 0.01, "plotting-regression", ...
               "horizontalalignment", "left", ...
               "verticalalignment", "bottom", "fontsize", 5);
assert(ishghandle(h));
print("plotting-regression.pdf");
EOF

test -s plotting-regression.pdf
echo "Octave plotting regression passed"
