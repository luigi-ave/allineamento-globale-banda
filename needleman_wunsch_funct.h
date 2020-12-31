#ifndef NEEDLEMAN_WUNSCH_FUNCT
#define NEEDLEMAN_WUNSCH_FUNCT
// #define DEBUG
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>

/**
 * Struttura informativa della cella della matrice
*/
typedef struct cell_s {
	int cell; ///< valore della cella 
	uint64_t prev_x; ///< coordinata precedente x
	uint64_t prev_y; ///< coordinata precedente y
	uint64_t val; ///< valore corrispondente alla provenienza (0-1-2)
} cell_s;

/**
 * Struttura informativa range della banda
*/
typedef struct range_s {
	uint64_t sx; ///< limite sinistro
	uint64_t dx; ///< limite destro
} range_s;

/**
 * Struttura informativa sulla banda
*/
typedef struct band_s {
	uint64_t l1; ///< lunghezza stringa 1
	uint64_t l2; ///< lunghezza stringa 2
	uint64_t base; ///< larghezza banda base
	uint64_t extra; ///< ampliamento banda
} band_s;


/**
 * Restituisce il valore di score tra i due caratteri
 * 
 * @param c1 Primo carattere
 * @param c2 Secondo carattere
 * 
 * @return Il valore di score tra i caratteri
*/
int score_value(char c1, char c2);

/**
 * Restituisce, data la banda, i limiti minimo e massimo in essa
 * 
 * @param band Informazioni sulla banda
 * @param y Coordinata y
*/
range_s range(band_s band, uint64_t y);

/**
 * Converte le coordinate in coordinata lineare
 * 
 * @param band Informazioni sulla banda
 * @param x Coordinata x
 * @param y Coordinata y
 * 
 * @return Coordinata lineare
*/
uint64_t conv(band_s band, uint64_t x, uint64_t y);

/**
 * Allineamento globale di due sequenze
 * 
 * @param base Valore della base (differenza stringhe + 1)
 * @param extra Ampliamento banda aggiuntivo
 * @param s1 Prima stringa
 * @param s2 Seconda stringa
 * @param score Valore di score che verrà sovrascritto
 * 
 * @return Allineamento ottimo delle due stringhe, NULL altrimenti
*/
char *band_align(uint64_t base, uint64_t extra, char *s1, char *s2, int *score);

/**
 * Controlla che la stringa inserita sia formata da 
 * caratteri accettati
 * 
 * @param s1 Stringa da controllare
 * 
 * @return true se la stringa è accettata, false altrimenti
*/
bool check_string(char *s1);

/**
 * Trasforma una stringa da minuscolo in maiuscolo
 * 
 * @param s Stringa da trasformare
*/
void toUpper(char *s);

#endif