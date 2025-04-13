# První část projektu
Složka : rotatable_bonds
## Řešení Confirmational Entropy
Zde je mé řešení první části zadání projektu Confirmational Entropy, které zkoumá distribuci rotovatelných vazeb v ligandech z protein-ligand komplexů z PL-REX datasetu.
## Program
Program, který slouží ke spočítání počtu rotovatelných vazeb v jednotlivých molekulách, uložení výsledků do tabulek a následnému vytvoření histogramu, se nachází v souboru conformation_entropy_calculation.py.
## Výsledky
Výsledky, tedy histogram distribuce rotovatelných vazeb v ligandech a shrnující tabulka ve formátu .xlsx a .csv jsou k nalezení ve složce rotatable_bonds/results. Dále zde lze nalézt shrnující tabulku s minimem a maximem rotovatelných vazeb v ligandech pro každý target protein.

# Druhá část projektu - analýza CREST
složka: CREST
## Programy
Programy jsou k nalezení ve složce CREST/programs. 
###
Pro snadné spuštění programů ve správném pořadí slouží skript run_crest_analysis.sh v hlavní složce.
#### extract_sconf
Program který extrahuje konformační entropie všech ligandů JAK1.
#### other_information
Program získává další informace o ligandech zájmu. Jde o: molekulovou hmotnost, množství těžkých atomů, množství donorů/akceptorů vodíkových můstků.
####  plotting_and_correlation
Program, který pro každou z výše zmíněných metrik + z předchozí části úkolu již známé množství rotovatelných vazeb provádí korelační analýzy a plotuje vztah mezi těmito metrikami a konformační entropií.


## Výsledky
Výsledky jsou k nalezení ve složce CREST/results. V podsložce plot jsou uloženy grafy závislostí jednotlivých metrik s konformační entropií. V podsložce tables je tabulka všech metrik pro každý ligand a také tabulka výsledků korelačních testů pro každou metriku všetně p-value.

