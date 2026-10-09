#!/usr/bin/env bash
# Builds a standalone RealPhy.jar using only the JDK and the jars in lib/.
set -euo pipefail

cd "$(dirname "$0")"

for tool in javac jar; do
  command -v "$tool" >/dev/null || { echo "Error: '$tool' not found. Install a JDK." >&2; exit 1; }
done

LIBS=(lib/RCaller-2.5.jar lib/sam-1.80.jar)
for j in "${LIBS[@]}"; do
  [ -f "$j" ] || { echo "Error: missing $j" >&2; exit 1; }
done

BIN=bin
rm -rf "$BIN"
mkdir -p "$BIN"

echo "==> Compiling"
find src -name '*.java' > sources.txt
CP=$(IFS=:; echo "${LIBS[*]}")
javac --release 8 -Xlint:-options -cp "$CP" -d "$BIN" @sources.txt
echo "==> Unpacking libraries"
ABS_LIBS=()
for j in "${LIBS[@]}"; do ABS_LIBS+=("$PWD/$j"); done
(
  cd "$BIN"
  for j in "${ABS_LIBS[@]}"; do jar -xf "$j"; done
  rm -rf META-INF
)

echo "==> Packaging RealPhy.jar"
jar cfe RealPhy.jar buildTreesFromSequenceData.RealPhy -C "$BIN" .

echo "Done: $PWD/RealPhy.jar"
