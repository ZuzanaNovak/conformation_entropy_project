Konformační entropie ligandů
============================

Vstupní data
------------

Pracujeme s ligandy z protein-ligand komplexů z [PL-REX datasetu](https://doi.org/10.1038/s41467-024-45431-8). 3D struktury molekul ve formátu SDF jsou k dispozici v adresáři `ligands`.

V datasetu jsou samostatné série ligandů pro deset cílových proteinů. První část jména souboru je číslo série a zkratka názvu cílového proteinu, např. "001-CA2". Druhá část je označení samotného ligandu (např. 5NXG, což je PDB kód struktury ze keteré byl ligand získán).


Úkol 1: Odhad pomocí počtu rotovatelných vazeb
----------------------------------------------

Jako měřítko konformační energie molekuly se často používá počet volně rotovatelných vazeb. Použijte Python a knihovnu [rdkit](https://www.rdkit.org/docs/index.html) ke spočítání rotovatelných vazeb ve všech ligandech.

Výstupem by měla být tabulka výsledků a histogram distribuce počtu rotovatelných vazeb v celém datasetu.


Úkol 2: Hledání konformerů programem crest
------------------------------------------

Vyberete jednu sérii ligandů ve které jsou největší rozdíly počtu rotovatelných vazeb. Na těchto ligandech proveďte hledání konformerů programem [crest](https://crest-lab.github.io/crest-docs/). Program je pro linux, na windows jde použít WSL nebo virtuální stroj s linuxem.

Příklad k hledání konformerů je [zde](https://crest-lab.github.io/crest-docs/page/examples/example_1.html). Až zprovozníte program a vyzkoušíte tento příklad, konzultujte nastavení programu pro další výpočty.

K tomu je třeba přidat i výpočet [konformační entropie](https://crest-lab.github.io/crest-docs/page/examples/entropy.html)

Analýza výsledků
----------------

Srovnejte a korelujte počty rotovatelných vazeb a vypočítané konformační entropie.


