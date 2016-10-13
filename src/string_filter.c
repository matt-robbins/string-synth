#include "string_filter.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

string_filter_t * string_filter_create(float freq, float rate, float damp, float fdamp)
{
	float delta;
	string_filter_t * s = calloc(1,sizeof(string_filter_t));
	float fL = rate/(freq) - 1.1;
	s->L = (int) floorf(fL);
	delta = fL - (float) s->L;
		
	s->comb = calloc(s->L, sizeof(float));
	s->damp = damp * (0.5 - freq/rate);
	s->fdamp = fdamp * (0.5 - freq/rate);
	
	s->nu = (1 - delta);
	s->disp = -0.01;
	
	if (freq < 55)
	{
		//s->disp = -0.5; s->L -= lrintf(s->L * 0.1); 
	}
	
	return s;
}

void string_filter_destroy(string_filter_t * s)
{
	free(s->comb);
	free(s);
}


static inline int circ_index_add(int ix, int d, int N)
{
	int out;
	if (d > N)
		d = d % N;
		
	out = ix + d;
	
	if (out < 0)
		out = N + out;
		
	if (out > N - 1)
		out = out - N;
	
	return out;
}

void string_filter_init_data(string_filter_t * s, float * in, int N)
{
	int k;
	
	for (k = 0; k < N && k < s->L; k++)
	{
		s->comb[k] += in[k];
	}
}

void string_filter_init_impulse(string_filter_t * s, float in)
{
	int k, ix, N;

	if (in < 0.0001)
		return;
	
	float fN = (8) / in;
	float val;
	
	N = (int) lrintf(fN);
	
	if (N > s->L-2)
		N = s->L-2;
	if (N < 3)
		N = 3;
						
	ix = circ_index_add(s->ix, 1, s->L);
	for (k = 0; k < N; k++)
	{
		val = in * 0.4 * (1 - cos(2*M_PI * (float)k / (float)N));
		
		s->comb[ix] += val;
		ix = circ_index_add(ix, 1, s->L);
	}
}

void string_filter_init_rand(string_filter_t * s, float in)
{
	int k;
	
	for (k = 0; k < s->L; k++)
	{
		s->comb[k] += in*(((float) rand() / RAND_MAX) - 0.5);
	}
}

void string_filter_damp(string_filter_t * s, int on)
{
	s->pedal_damp = on ? 0.1 : 0.0;
}



static inline float allpass_update(float * xy, float nu, float new)
{
	float out;
	float nu_loc = nu;
	
	if (nu_loc > 1)
		nu_loc = 1;
	out = xy[1] = (new - xy[1]) * nu + xy[0]; xy[0] = new;
	return out;
}

static void string_filter_compute_step(string_filter_t * s, float * out, float fb_in)
{
	float tuned;
	s->dbuf = s->fdamp*s->dbuf + (1 - s->fdamp)*s->comb[s->ix];
	s->lbuf = s->pedal_damp * s->lbuf - (1 - s->pedal_damp) * s->dbuf;
	tuned = allpass_update(s->abuf, s->nu, s->lbuf);
	s->comb[s->ix] = allpass_update(s->disp_buf, s->disp, tuned);
	
	//s->comb[s->ix] *= (1 - s->damp) * (1 - s->pedal_damp);
	s->comb[s->ix] += fb_in;
	*out += (s->comb[s->ix] - s->comb[circ_index_add(s->ix, -lrintf(s->L/7.0), s->L)]);
		
	s->ix = circ_index_add(s->ix, 1, s->L);
}

void string_filter_compute(string_filter_t * s, float * out, int N)
{
	int last_ix, k;
	float new;
	
	for (k = 0; k < N; k++)
	{
		out[k] = 0.0f;
		string_filter_compute_step(s, out + k, 0.0f);
	}
}

float midi_to_freq(int N)
{
	return 440.0f * powf(2.0f, ((float) N - 69.0f)/12.0f);
}

string_filter_bank_t * string_filter_bank_create(int low, int hi, float rate, float damp, float fdamp)
{
	string_filter_bank_t * b = calloc(1, sizeof(string_filter_bank_t));
	int k;
	float f;
	
	b->N = hi - low + 1;
	b->low = low;
	b->hi = hi;
	b->fb = 0.002;
	b->damp_on = 0;
	
	b->filters = calloc(b->N, sizeof(string_filter_t *));
	
	for (k = 0; k < b->N; k++)
	{
		f = midi_to_freq(low + k);
		b->filters[k] = string_filter_create(
			f, rate, damp, fdamp * (1 - 0.9 * ((float) k/b->N)));
			
		string_filter_damp(b->filters[k], 1);
	}
	
	return b;
}

void string_filter_bank_compute(string_filter_bank_t * b, float * out, int N)
{
	int j, k;
		
	for (j = 0; j < N; j++)
	{
		out[j] = 0;
		
		for (k = 0; k < b->N; k++)
		{
			string_filter_compute_step(b->filters[k], out + j, -b->fb_out * b->fb);
		}
		
		b->fb_out = out[j];
	}
}

void string_filter_bank_noteon(string_filter_bank_t * b, int midinote, int midivel)
{
	int k, n;
	float vel;
	if (midinote < b->low || midinote > b->hi)
		return;
		
	n = midinote - b->low;
	
	vel = (float) midivel / 127.0f;
	
	for (k = 0; k < b->N; k++)
	{
		if (k == n)
		{
			string_filter_init_impulse(b->filters[k], vel);
			b->filters[k]->note_on = 1;
			string_filter_damp(b->filters[n], 0);
		}
		else {}
			//string_filter_init_impulse(b->filters[k], vel*0.01);
	}
}

void string_filter_bank_noteoff(string_filter_bank_t * b, int midinote, int midivel)
{
	int n;
	
	if (midinote < b->low || midinote > b->hi)
		return;
		
	n = midinote - b->low;
	
	b->filters[n]->note_on = 0;
	string_filter_damp(b->filters[n], b->damp_on);
		
	//string_filter_damp(b->filters[n], midivel * 1000);
	
}

void string_filter_bank_dampset(string_filter_bank_t * b, int val)
{
	int k;
	int d_on = val < 64;
	
	b->damp_on = d_on;
	
	for (k = 0; k < b->N; k++)
	{
		string_filter_damp(b->filters[k], !b->filters[k]->note_on & d_on);
		
		if (!d_on)
		{
			string_filter_init_rand(b->filters[k], 0.0002);
		}
	}
}


void string_filter_bank_destroy(string_filter_bank_t * b)
{
	int k;
	
	for (k = 0; k < b->N; k++)
	{
		string_filter_destroy(b->filters[k]);
	}
	
	free(b->filters);
	free(b);
}


#ifdef TEST
int main(int argc, char * argv[])
{
	string_filter_t * f = string_filter_create(1300, 44100, 0.001, 0.2);
	int N = 1000;
	float * out = calloc(N, sizeof(float));
	int k;
	
	string_filter_init_impulse(f, 1.0);
	
	string_filter_compute(f, out, N);
	
	for (k = 0; k < N; k++)
	{
		fprintf(stdout, "%f\n", out[k]);
	}
	
	string_filter_init_impulse(f, 1.0);
	
	string_filter_compute(f, out, N);
	
	for (k = 0; k < N; k++)
	{
		fprintf(stdout, "%f\n", out[k]);
	}
	
	free(out);
	string_filter_destroy(f);
	
	return 0;
}


#endif
