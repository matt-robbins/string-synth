#define D 2

#ifdef MATLAB
	#include <mex.h>
	#define malloc mxMalloc
	#define free mxFree
	#define calloc mxCalloc
	#define realloc mxRealloc
#endif

#define float_data float
#define int_data int32_t

void sim(float_data *r, float_data *v, float_data *aud, float_data *damp,
			float_data h, int n, int nS, float_data nP, float_data L,
			float_data K, float_data S);
