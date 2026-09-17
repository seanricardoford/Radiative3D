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

arrayimage_output="$scratch_dir/arrayimage-output.txt"
octave --no-gui -qf > "$arrayimage_output" 2>&1 <<EOF
addpath("$repo_root/vis/seisplot");
for i = 0:1
  SEIS.TraceXYZ = [1.0 0.0 0.0;
                   0.5 0.2 0.1;
                   0.2 0.1 0.3;
                   0.1 0.1 0.1];
  SEIS.Location = [i + 1, 0, 0];
  SEIS.EventLoc = [0, 0, 0];
  SEIS.TimeWindow = [0, 4];
  SEIS.NumBins = 4;
  SEIS.Frequency = 1;
  SEIS.EventMT = [1, 0, 0; 0, -1, 0; 0, 0, 0];
  save(sprintf("seis_%03d.octv", i), "SEIS");
end
AR = array("seis_%03d.octv", 0, 1);
arrayimage(AR, [1 1 1], 2.0, [0 4], 0.5, 0);
annotate_array();
set(findobj("tag", "colorbar"), "yticklabel", []);
set(findobj("tag", "colorbar"), "ytick", []);
newstr = get(findobj("tag", "OLCapt2"), "string");
newstr = sprintf("%s | %s", "Channel:  Ex+Ey+Ez", newstr);
set(findobj("tag", "OLCapt2"), "string", newstr);
print("arrayimage-regression.pdf");
EOF

if rg --fixed-strings --quiet -- "invalid command" "$arrayimage_output"; then
  cat "$arrayimage_output" >&2
  exit 1
fi
test -s arrayimage-regression.pdf

for script in do-lopnor.sh do-halfspace.sh do-halfspace-nearsrc25.sh \
              do-halfspace-nearsrc50.sh dot-mltw.sh; do
  if rg --fixed-strings --quiet -- 'sprintf("%s\n%s","Channel:  Ex+Ey+Ez"' "$repo_root/$script"; then
    echo "$script must not inject literal newlines into gnuplot captions" >&2
    exit 1
  fi
done
echo "Octave plotting regression passed"
