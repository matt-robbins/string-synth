#ifndef STRING_FILTER_H
#define STRING_FILTER_H
typedef struct string_filter_s {
	int L;
	float * comb;
	float dbuf;
	float lbuf;
	float abuf[2];
	float disp_buf[2];
	float damp;
	float fdamp;
	float nu;
	float disp;
	int ix;
	float pedal_damp;
	int damp_on;
	int note_on;
} string_filter_t;

typedef struct string_filter_bank_s {
	string_filter_t ** filters;
	float fb_out;
	float fb;
	int N;
	int low;
	int hi;
	float pedal_damp;
	int damp_on;
} string_filter_bank_t;

float midi_to_freq(int N);

string_filter_bank_t * string_filter_bank_create(int low, int hi, float rate, float damp, float fdamp);

string_filter_t * string_filter_create(
	float freq, float rate, float damp, float fdamp);
	
void string_filter_bank_destroy(string_filter_bank_t * b);

void string_filter_init_data(string_filter_t * s, float * in, int N);
void string_filter_init_impulse(string_filter_t * s, float in);

void string_filter_destroy(string_filter_t * s);

void string_filter_compute(string_filter_t * s, float * out, int N);
void string_filter_bank_compute(string_filter_bank_t * s, float * out, int N);

void string_filter_bank_noteon(string_filter_bank_t * b, int midinote, int vel);
void string_filter_bank_noteoff(string_filter_bank_t * b, int midinote, int vel);
void string_filter_bank_dampset(string_filter_bank_t * b, int val);

#endif
