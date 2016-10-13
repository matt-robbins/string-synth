#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "simR.h"

typedef struct
{
    float_data *r;
	float_data *v; 
	float_data *d;
	float_data h;
	int n; 
	int nS; 
	float_data nP;
	float_data L;
	float_data K; 
	float_data S;
}
stringState;

#define NTEST 200

void string_init(stringState * state){
	
	int k;
	
	float_data damp_f = 10000.;
	
	state->r = (float_data *) calloc(state->n * 2 + 5, sizeof(float_data));
	state->v = (float_data *) calloc(state->n * 2 + 5, sizeof(float_data));
	state->d = (float_data *) calloc(state->n * 2, sizeof(float_data));
	
	state->S = 1.2; state->K = 20;
	
	for (k = 0; k < state->n; k++){
		
		state->r[2*k + 1] = (float_data) k + 1;
		state->d[k] = 0.00001;
	}
	
	state->v[6] = 3.0; 
}



int main(){
	
	int i;
	stringState state;
	
	state.h = 0.001;
	state.n = 20;
	state.nS = 1;
	state.nP = 20;
	state.L = 20.;
	state.K = 3.;
	state.S = 1.3;
	
	float * out = calloc(NTEST, sizeof(float));
	
	string_init(&state);
	
	sim(
		state.r, state.v, out, state.d, state.h, state.n, NTEST, 
		state.nP, state.L, state.K, state.S
	);
	/*
	for (i = 0; i < NTEST; i++){
		printf("%f\n", out[i]);
	}
	
	printf("\n");
	*/
}
