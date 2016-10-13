#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

typedef struct particle_s {
	float y;
	float v;
	float f;
	float damp;
	float m;
	float k;
} particle_t;

typedef struct sstring_s {
	int npoints;
	float step_size;
	float mass;
	float spring_k;
	float tension;
	float gdamp;
	float damp;
	particle_t * particles;
	particle_t lbound, rbound;
	float strike;
	int strike_ix;
} simple_string_t;

simple_string_t * simple_string_create(
int npoints, float mass, float spring_k, float tension, 
float damp, float gdamp, float step_size);

void simple_string_init_sin(simple_string_t * s);

void simple_string_destroy(simple_string_t * s);

static inline float particle_force(float dyl, float dyr, float v, float spring_k, float damp);

static inline void particle_update_force(particle_t * p, particle_t * l, particle_t * r, simple_string_t * s);

void simple_string_compute(simple_string_t * s, float * out, int N);
