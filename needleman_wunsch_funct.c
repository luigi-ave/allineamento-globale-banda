#include "BLOSUM62.h"
#include "needleman_wunsch_funct.h"

int MIN = INT16_MIN;

int score_value(char c1, char c2) {
	unsigned int x, y;
	for (unsigned int i = 0; i < SIZE; ++i){
		if (alfabeto[i] == c1)
			x = i;
	}
	for (unsigned int i = 0; i < SIZE; ++i){
		if (alfabeto[i] == c2)
			y = i;
	}
	return matrix_score[x][y];
}


bool check_string(char *s1){
	bool val = true;
	for (unsigned int i = 0; i < strlen(s1) && val;i++) {
		val = false;
		for (unsigned int k = 0; k < SIZE - 1; k++) {
			if(s1[i] == alfabeto[k])
				val = true;
		}
	}
	return val;
}

void toUpper(char *s){
	for (unsigned int i = 0; i < strlen(s); i++){
		if(s[i] >= 'a' && s[i] <= 'z') {
			s[i] = s[i] - 32;
		}
	}
}

range_s range(band_s band, uint64_t y) {
	range_s r = {0, 0};
	r.sx = y - band.extra;
	if (y < band.extra) r.sx = 0;
	r.dx = y + (band.base - 1) + band.extra;
	if (r.dx > band.l2) r.dx = band.l2;
	return r;
}

uint64_t conv(band_s band, uint64_t x, uint64_t y) {
	assert(y <= band.l1);

	assert(x >= y || y - x <= band.extra);
	assert(x <= y || x - y <= band.base + band.extra);

	uint64_t width = band.base + 2 * band.extra;
	return (width * y) + (x - y + band.extra);
}


char* band_align(uint64_t base, uint64_t extra, char* s1, char* s2, int *score) {
	uint64_t l1 = strlen(s1);
	uint64_t l2 = strlen(s2);
	uint64_t band = base + 2 * extra;

#ifdef DEBUG
	printf("\nBANDA: %lu - BASE: %lu L1:%lu", band, base, l1);
#endif
	cell_s* m = (cell_s*)malloc(band * (l1 + 1) * sizeof(cell_s));
	assert(m != NULL && "Cannot allocate matrix\n");

	m[0].cell = 0;
	
	printf("\n");
	band_s b = {
		.l1 = l1,
		.l2 = l2,
		.base = base,
		.extra = extra,
	};

	for (uint64_t y = 1; y <= l1; y++) {
		range_s r = range(b, y);
		{
			cell_s* tp = m + conv(b, r.sx, y);
			tp->cell = MIN;
			tp->prev_x = 0;
			tp->prev_y = 0;
			tp = m + conv(b, r.dx, y);
			tp->cell = 0;
			tp->prev_x = 0;
			tp->prev_y = 0;
		}
#ifdef DEBUG
		printf("range: %lu %lu\n", r.sx, r.dx);
#endif
		for (uint64_t x = r.sx + 1; x <= r.dx; x++) { 
#ifdef DEBUG
			printf("StepB: %lu, %lu\n", x, y);
#endif
			// Needleman-Wunsch
			cell_s* tp = m + conv(b, x, y);

			int diag, up, left;
			if(x - 1 == 0 && y - 1 == 0){
				diag = 0  + score_value(s1[y-1], s2[x-1]);
				up = -4 - 4;
				left = -4 - 4;
			}else{
				if(x - 1 == 0 ){
					diag = -4 * (y - 1) + score_value(s1[y-1], s2[x-1]);
					up = m[conv(b, x, y - 1)].cell + score_value('*', s2[x-1]);
					left = MIN;
				}else{
					if(y - 1 == 0){
						diag = -4 * (x - 1) + score_value(s1[y-1], s2[x-1]);
						up = -4 * x - 4;
						left = m[conv(b, x - 1, y)].cell + score_value(s1[y-1], '*');
					}else{
						diag = m[conv(b, x - 1, y - 1)].cell + score_value(s1[y-1], s2[x-1]);
						up = m[conv(b, x, y - 1)].cell + score_value('*', s2[x-1]);
						left = m[conv(b, x - 1, y)].cell + score_value(s1[y-1], '*');
					}
				}
			}

			// Cerco il massimo
			if (diag >= up && diag >= left){
				tp->cell = diag;
				tp->prev_x = x - 1;
				tp->prev_y = y - 1;
				tp->val = 0;
			}else{
				if (up >= diag && up >= left){
					tp->cell = up;
					tp->prev_x = x;
					tp->prev_y = y - 1;
					tp->val = 2;
				}else{
					tp->cell = left;
					tp->prev_x = x - 1;
					tp->prev_y = y;
					tp->val = 1;
				}
			}
						
#ifdef DEBUG
			printf("StepB result: %i, %lu, %lu\n", m[conv(b, x, y)].cell, m[conv(b, x, y)].prev_x, m[conv(b, x, y)].prev_y);
#endif
		}
	}
	
	
#ifdef DEBUG
	printf("CONV: base=%lu, extra=%lu\n", b.base, b.extra);
	printf("MATRIX: base=%lu, extra=%lu\n", base, extra);
#endif


/* ricostruzione allineamento, se possibile */
	uint64_t x = l2;
	uint64_t y = l1;
#ifdef DEBUG
	printf("StepC: Ricostruzione allineamento\n");
#endif

	char* allign = (char*)malloc((STRING_SIZE * 4) * sizeof(char));
	allign[STRING_SIZE * 4] = '\0';
	int k = 0;
	*score = m[conv(b, x, y)].cell;

	for (cell_s c = m[conv(b, x, y)]; x > 0 && y > 0; x = c.prev_x, y = c.prev_y, c = m[conv(b, x, y)]) {
#ifdef DEBUG
			printf("StepC: %lu,%lu,%i,%lu,%lu=%c%c\n", x, y, c.cell, c.prev_x, c.prev_y, s1[y-1], s2[x-1]);
#endif

		range_s r = range(b, y);
		if ((x == r.sx && r.sx > 0) || (x == r.dx && r.dx < l2)) 
			return NULL;

		if (c.val == 0){
			allign[k] = s1[y - 1];
			k++;
			allign[k] = s2[x - 1];
			k++;
		}
		if(c.val == 1){
			allign[k] = '-';
			k++;
			allign[k] = s2[x - 1];
			k++;
		}
		if(c.val == 2){
			allign[k] = s1[y - 1];
			k++;
			allign[k] = '-';
			k++;
		}
	}

	while(y != 0){
		allign[k] = s1[y - 1];
		k++;
		allign[k] = '-';
		k++;
		y--;
	}

	while(x != 0){
		allign[k] = '-';
		k++;
		allign[k] = s2[x - 1];
		k++;
		x--;
	}
	
	free(m);
	
	return allign;
}