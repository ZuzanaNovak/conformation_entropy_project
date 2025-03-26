Konformační entropie ligandů
============================

Vstupní data
------------

Pracujeme s ligandy z protein-ligand komplexů z [PL-REX datasetu](https://doi.org/10.1038/s41467-024-45431-8). 3D struktury molekul ve formátu SDF jsou k dispozici v adresáři `ligands`.


Odhad pomocí počtu rotovatelných vazeb
-----------------------------------------

Jako měřítko konformační energie molekuly se často používá počet volně rotovatelných vazeb. Použijte Python a knihovnu [rdkit](https://www.rdkit.org/docs/index.html) ke spočítání rotovatelných vazeb ve všech ligandech.


2) Hledání konformerů programem crest


