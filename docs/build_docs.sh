#!/bin/bash

SCRIPT=$(basename "$0")
PREFIX="[$SCRIPT]"

echo "$PREFIX Building documentation (Sphinx + Breathe + Doxygen)..."
echo
# --------------------------------------

doxygen
make html

echo
echo "$PREFIX Done."
# --------------------------------------
