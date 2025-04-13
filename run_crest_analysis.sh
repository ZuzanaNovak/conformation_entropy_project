#!/bin/bash
set -e  # Ukončí skript, pokud nastane chyba

echo "Extracting Sconf"
python CREST/programs/extract_sconf.py

echo
echo "Computing other molecular information"
python CREST/programs/007_JAK1_other_informations.py

echo
echo "Plotting and correlation analysis"
python CREST/programs/plotting_and_correlation.py

echo
echo "Done!"