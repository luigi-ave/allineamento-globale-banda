#include "BLOSUM62.h"
#include "needleman_wunsch_funct.h"

int main (int argc, char *argv[]) {
	char s1[STRING_SIZE];
	char s2[STRING_SIZE];

	do{
		printf("Sequenza S1: ");
		scanf("%s", s1);
		toUpper(s1);
		if (!check_string(s1))
			printf("Caratteri stringa non validi\n");
	} while (!check_string(s1));

	do{
		printf("Sequenza S2: ");
		scanf("%s", s2);
		toUpper(s2);
		if (!check_string(s2))
			printf("Caratteri stringa non validi\n");
	} while (!check_string(s2));


	assert(strlen(s1) < STRING_SIZE && "Limite stringa superato\n");
	assert(strlen(s2) < STRING_SIZE && "Limite stringa superato\n");

	if (strlen(s2) < strlen(s1)) {
		char tmp[STRING_SIZE];
		strcpy(tmp, s1); 
		strcpy(s1, s2); 
		strcpy(s2, tmp); 
	}

	int score = INT_MIN;
	int pscore = INT_MIN;
	uint64_t l1 = strlen(s1);
	uint64_t l2 = strlen(s2);
	uint64_t base = l2 - l1 + 1;
	uint64_t extra = 1;

#ifdef DEBUG
	printf("Sequenza 1: %s\n", s1);
	printf("Sequenza 2: %s\n", s2);
	printf("Lunghezza istanze: %ld,%ld\n", l1, l2);
#endif

	char* nwab = band_align(base, 0, s1, s2, &score);
	nwab = band_align(base, extra, s1, s2, &pscore);

	while(pscore < score || nwab == NULL){
		nwab = band_align(base, extra, s1, s2, &pscore);
		nwab = band_align(base, extra*2, s1, s2, &score);
		extra *= 2;
	}

	if (extra == 1){
		nwab = band_align(base, 1, s1, s2, &score);
	}else{
		if(pscore == score){
			nwab = band_align(base, extra, s1, s2, &score);
		}else{
			nwab = band_align(base, extra/2, s1, s2, &score);
		}
	}


	// Stampa allineamento e valore di score
	if (nwab != NULL) {
		printf("Allineamento trovato\n");
		int i = strlen(nwab) - 1;
		do{
			printf("%c", nwab[i]);
			i = i - 2;
		} while (i >= 0);
		printf("\n");
		i = strlen(nwab) - 2;
		do{
			printf("%c", nwab[i]);
			i = i - 2;
		} while (i >= 0);
	}
	printf("\nscore: %i\n", score);

	free(nwab);
	return 0;
}

