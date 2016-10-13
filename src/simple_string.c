#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "simple_string.h"

simple_string_t * simple_string_create(
int npoints, float mass, float spring_k, float tension, float damp, float gdamp, float step_size)
{
	simple_string_t * s = calloc(1, sizeof(simple_string_t));
	int k;
	
	s->npoints = npoints;
	s->step_size = step_size;
	s->mass = mass;
	s->spring_k = spring_k;
	s->tension = tension;
	s->strike = 0;
	s->damp = damp;
	s->gdamp = gdamp;
	
	s->particles = calloc(s->npoints, sizeof(particle_t));
	
	for (k = 0; k < npoints; k++)
	{
		s->particles[k].damp = damp;
		s->particles[k].m = mass;
		s->particles[k].k = spring_k;
	}
		
	return s;
}

void simple_string_destroy(simple_string_t * s)
{
	free(s->particles);
	free(s);
}

static inline float particle_force(float dyl, float dyr, float v, float spring_k, float damp)
{
	return -dyl * spring_k - dyr * spring_k - v * damp;
}

inline double fastPow(double a, double b) {
  union {
    double d;
    int x[2];
  } u = { a };
  u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
  u.x[0] = 0;
  return u.d;
}

static inline float spring_nonlinearity(float in)
{
	return in + fastPow(in/3, 3);
}

static inline void particle_update_force(particle_t * p, particle_t * l, particle_t * r, simple_string_t * s)
{
	p->f = spring_nonlinearity((l->y - p->y)) * s->spring_k + 
		   spring_nonlinearity((r->y - p->y)) * s->spring_k - 
		   p->damp * (p->v - (l->v + r->v) / 2) - 
		   s->gdamp * p->v;
		   
	
}

static inline void particle_integrate(particle_t * p, float step_size)
{
	p->v += p->f/p->m * step_size;
	p->y += p->v * step_size;
}

void simple_string_init_sin(simple_string_t * s)
{
	int k;
	for (k = 0; k < s->npoints; k++)
	{
		s->particles[k].f = 0;
		s->particles[k].v = 0;
		s->particles[k].y = 0.1*sin(M_PI * 1 * ((float) (k + 0.5) / s->npoints));
	}
}

void simple_string_init_hammer(simple_string_t * s)
{
	int nhammer = s->npoints / 10;
	float hammer_pos = 1.0/5.0;
	int k, ix;
	
	if (s->npoints < 10)
	{
		s->particles[0].f += 30;
		return;
	}
	
	int hn = lrintf(hammer_pos * s->npoints);
	
	if (hn == 0)
		hn = 1;
	
	if (nhammer == 0)
		nhammer = 1;
				
	for (k = 0; k < lrintf(nhammer); k++)
	{
		ix = hn + k - lrintf(nhammer);
		s->particles[ix].f += 50.0*(cos(2 * M_PI * (k / hn)) + 1)/2;
	}
	
}

static inline void set_damp(simple_string_t * s, float damp)
{
	int k;
	
	for (k = 0; k < s->npoints; k++)
	{
		s->particles[k].damp = damp;
	}
}

void simple_string_compute(simple_string_t * s, float * out, int N)
{
	int i,j;
	int e = s->npoints-1;
	
	for (i = 0; i < N; i++)
	{
			
		/* update force*/
		
		particle_update_force(s->particles, &s->lbound, s->particles + 1, s);
			
		for (j = 1; j < e; j++)
		{
			particle_update_force(s->particles+j, s->particles+j-1, s->particles+j+1, s);
		}
		
		particle_update_force(s->particles+e, s->particles+e-1, &s->rbound, s);
		
		if (s->strike != 0)
		{
			set_damp(s, 1);
			s->strike = 0;
			s->strike_ix = 200;
		}
		
		if (s->strike_ix > 0)
		{
			simple_string_init_hammer(s);
			s->strike_ix--;
			
			if (s->strike_ix == 0)
			{
				set_damp(s, s->damp);
			}
		}
		
		/* update velocity, positions: discrete time Newton's method -> v += f/m * step_size */
		
		
		
		for (j = 0; j < s->npoints; j++)
		{
			particle_integrate(s->particles+j, s->step_size);
		}
				
		out[i] = s->particles[0].y;
	}
}

#ifdef TEST
int main(int argc, char * argv[])
{
	int N = 20;
	simple_string_t * s = simple_string_create(N, 0.01, 500, 0, 0.01, 0.00002);
	int i,k;
	float data[100];
	
	s->strike = 1.0;
	
	for (i = 0; i < 100; i++)
	{
		simple_string_compute(s, data, 100);
		
		for (k = 0; k < 100; k++)
		{
			fprintf(stdout, "%f\n", data[k]);
		}
	}
	
	simple_string_destroy(s);
}
#endif
