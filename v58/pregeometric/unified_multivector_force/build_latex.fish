#!/usr/bin/env fish
# build_latex.fish
#
# Hybrid build script: Markdown → standalone .tex (via pandoc) → PDF (via LuaLaTeX + STIX Two Math)
#
# This route is more reliable for long, math-heavy documents than direct pandoc → PDF.
#
# Usage:
#   ./build_latex.fish                              # builds the default paper
#   ./build_latex.fish PREGEOMETRIC_EMERGENCE_OF_EFFECTIVE_GEOMETRY.md
#
# Output:
#   <basename>.tex
#   <basename>.pdf
#
# The script runs lualatex twice to resolve cross-references.

if test (count $argv) -gt 0
    set PAPER $argv[1]
else
    set PAPER "PREGEOMETRIC_EMERGENCE_OF_EFFECTIVE_GEOMETRY.md"
end

set -l BASENAME (basename "$PAPER" .md)
set -l TEXFILE "$BASENAME.tex"
set -l PDFFILE "$BASENAME.pdf"
set -l LOGFILE "$BASENAME"_latex_build.log

echo "=== Generating standalone .tex from $PAPER ==="

pandoc "$PAPER" -o "$TEXFILE" --standalone \
    -V mainfont="DejaVu Sans" \
    -V mathfont="STIX Two Math" \
    --variable=header-includes="\usepackage{unicode-math}\setmathfont{STIX Two Math}" \
    --variable=geometry:margin=1in \
    --variable=fontsize:11pt \
    -V papersize=a4 \
    --highlight-style=tango 2>&1 | tee "$LOGFILE"

if not test -f "$TEXFILE"
    echo "❌ Failed to generate .tex file. See $LOGFILE"
    exit 1
end

echo
echo "=== Compiling with LuaLaTeX (pass 1) ==="
lualatex -interaction=nonstopmode "$TEXFILE" 2>&1 | tee -a "$LOGFILE"

echo
echo "=== Compiling with LuaLaTeX (pass 2 - cross-references) ==="
lualatex -interaction=nonstopmode "$TEXFILE" 2>&1 | tee -a "$LOGFILE"

if test -f "$PDFFILE"
    echo
    echo "✅ PDF successfully produced: $PDFFILE"
    ls -lh "$PDFFILE"

    # Optional: clean up auxiliary files (comment out if you want to keep them for debugging)
    rm -f "$BASENAME".aux "$BASENAME".log "$BASENAME".out "$BASENAME".toc "$BASENAME".lof "$BASENAME".lot 2>/dev/null

    echo "Done. Primary output is $PDFFILE"
else
    echo
    echo "❌ PDF was not produced. Check $LOGFILE for errors."
    exit 1
end