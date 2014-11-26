TIGER-rates: DAAD_project/main.cpp basic-DNA-RNA-AA-routines.h CFile/CFile2_1.h CSequence_Mol2_1.h CSequences2.h CSplit2.h fast-dynamic-bitset/fast-dynamic-bitset.h fast-realloc-vector.h faststring2.h site_pattern.h
	g++ -O2 -Wall DAAD_project/main.cpp -o TIGER-rates -I .

clean:
	rm TIGER-rates
