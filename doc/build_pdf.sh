#!/usr/bin/env bash
# Builds Realphy_manual.pdf from Realphy_manual.md using pandoc and XeLaTeX.
set -euo pipefail

cd "$(dirname "$0")"

for tool in pandoc xelatex; do
  command -v "$tool" >/dev/null || { echo "Error: '$tool' not found. See doc/README.md for setup." >&2; exit 1; }
done

pandoc Realphy_manual.md \
  --from markdown+smart \
  --pdf-engine=xelatex \
  --include-in-header=header.tex \
  --resource-path=. \
  -V documentclass=article \
  -V papersize=a4 \
  -V fontsize=12pt \
  -V geometry:margin=2.5cm \
  -V mainfont="Times New Roman" \
  -V monofont="Courier New" \
  -V monofontoptions=Scale=0.85 \
  -V colorlinks=true \
  -V linkcolor=blue \
  -V urlcolor=blue \
  -o Realphy_manual.pdf

echo "Done: $PWD/Realphy_manual.pdf"
