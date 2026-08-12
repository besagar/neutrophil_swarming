#!/usr/bin/env bash
# Convert the Mathcha PDF export into per-page vector SVGs for the reveal deck.
#
# Workflow:
#   1. In Mathcha: Export → PDF (respects your page breaks).
#   2. Save it as the PDF path below (kept in the gitignored "slides html/" dir —
#      it's a local authoring source; the generated SVGs are what get committed).
#   3. Run this script:  ./slides/build-math.sh
#   4. Refresh deck.html — it lists slides/mathcha/*.svg automatically.
#
# SVG keeps the math crisp at any projector resolution. Pass "png" as the first
# arg for a rasterized fallback (2× density) if any glyph misrenders as SVG.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PDF="${PDF:-$ROOT/slides html/pdf_slides.pdf}"
OUT="$ROOT/slides/mathcha"
FMT="${1:-svg}"

if [[ ! -f "$PDF" ]]; then
  echo "PDF not found: $PDF" >&2
  echo "Export it from Mathcha (Export → PDF) and place it there, or set PDF=/path." >&2
  exit 1
fi

mkdir -p "$OUT"
rm -f "$OUT"/page-*.svg "$OUT"/page-*.png

PAGES=$(pdfinfo "$PDF" | awk '/^Pages:/ {print $2}')
echo "Converting $PAGES pages → $FMT in $OUT"

: > "$OUT/manifest.txt"   # reference list of generated pages (deck order lives in deck.html)
for i in $(seq 1 "$PAGES"); do
  n=$(printf "%02d" "$i")
  if [[ "$FMT" == "png" ]]; then
    # -r 192 ≈ 2× a 96dpi slide; crisp on projectors without huge files.
    pdftocairo -png -r 192 -f "$i" -l "$i" -singlefile "$PDF" "$OUT/page-$n"
  else
    pdftocairo -svg -f "$i" -l "$i" "$PDF" "$OUT/page-$n.svg"
  fi
  echo "page-$n.$FMT" >> "$OUT/manifest.txt"   # built in-loop → space-in-path safe
done

echo "Wrote $OUT/manifest.txt ($PAGES pages)"
