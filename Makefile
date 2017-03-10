mlgene: main.c sequence.c genealogy.c
	gcc -Wall -g -o mlgene main.c sequence.c genealogy.c -lm sfmt/SFMT.c -DSFMT_MEXP=19937

debug:  main.c sequence.c genealogy.c
	gcc -Wall -g -o mlgene main.c sequence.c genealogy.c -lm sfmt/SFMT.c -DSFMT_MEXP=19937 -DDEBUG

clean:
	rm mlgene *.o

