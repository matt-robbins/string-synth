
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "simR.h"

#include <jack/jack.h>
#include <jack/midiport.h>
#include "simple_string.c"
#include "string_filter.c"

jack_port_t *input_port;
jack_port_t *output_port;

simple_string_t * String;
string_filter_t * Filt;
string_filter_bank_t * Harp;

static inline void string_process_midi_event(string_filter_bank_t * f, jack_midi_event_t * e)
{
	int ix;
	uint8_t st, note, vel;
		
	if (e->size != 3)
		return;
	
	st = e->buffer[0] >> 4;
	note = e->buffer[1];
	vel = e->buffer[2];
		
	switch (st)
	{
		case 0xb:
			if (e->buffer[1] == 0x40)
				string_filter_bank_dampset(f, (int) e->buffer[2]);
			break;
		case 0x9:
			string_filter_bank_noteon(f, (int) note, (int) vel);
			break;
		case 0x8:
			string_filter_bank_noteoff(f, (int) note, (int) vel);
			break;
		default:
			break;
	}
}

int process(jack_nframes_t nframes, void *arg)
{
    string_filter_bank_t * f = (string_filter_bank_t *) arg;
    
    jack_default_audio_sample_t * out = (jack_default_audio_sample_t *) 
		jack_port_get_buffer (output_port, nframes);
		
	void * in = jack_port_get_buffer(input_port, nframes);
	int i;
	jack_midi_event_t e;
	jack_nframes_t ix = 0;
	jack_nframes_t n_events = jack_midi_get_event_count(in);
	
	if (!n_events)
	{
		string_filter_bank_compute(f, out, nframes);
		return 0;
	}
	
	ix = 0;
	for(i = 0; i < n_events; i++)
	{
		jack_midi_event_get(&e, in, i);		
		string_filter_bank_compute(f, out + ix, e.time - ix + 1);
		
		string_process_midi_event(f, &e);
		
		ix = e.time + 1;
	}
	
	string_filter_bank_compute(f, out + ix, nframes - ix);
	    	
    return 0;
}

void jack_shutdown(void *arg)
{
	exit(1);
}

/*******************************************************************/

int main(int argc, char ** argv)
{
    float rate;
    int N = 50;;
    int i, k, err;

	jack_client_t *client;
	
	float_data damp = 0.000001;
	float_data fdamp = 0.0001;
	int h_pos = 12;
	float_data h_vel = 3.;
	float_data stretch = 1.3;
	float_data K = 500.;
	
	if (argc > 1)
		K = atof(argv[1]);
	
	if (argc > 2)
		damp = atof(argv[2]);
	
	if (argc > 3)
		fdamp = atof(argv[3]);
		
	if (argc > 4)
		N = atoi(argv[4]);
		
	fprintf(stderr, "K = %f, damp = %f, N = %d\n", K, damp, N);
		
	String = simple_string_create(N, 0.001, K, stretch, damp, 0.00002, 0.0001);
		            
    if ((client = jack_client_open ("string", JackNullOption, NULL)) == 0)
	{
		fprintf(stderr, "jack server not running?\n");
		return 1;
	}
	
	rate = jack_get_sample_rate (client);
	N = jack_get_buffer_size(client);	
    printf("Play string. SR = %f, BufSize = %d\n", rate, N);
    
    Filt = string_filter_create(K, rate, damp, fdamp);
    
    Harp = string_filter_bank_create(21, 108, rate, damp, fdamp);
	
	jack_set_process_callback (client, process, Harp);
	jack_on_shutdown (client, jack_shutdown, 0);
	
	output_port = jack_port_register (client, "audio_out", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
    input_port = jack_port_register (client, "midi_in", JACK_DEFAULT_MIDI_TYPE, JackPortIsInput, 0);

    const char ** ports  = jack_get_ports(client, NULL, NULL, JackPortIsPhysical|JackPortIsInput);
    	
	if (jack_activate (client))
	{
		fprintf(stderr, "cannot activate client");
		return 1;
	}
	
	for (k = 0; ports[k] != NULL; k++)
	{
		jack_connect(client, jack_port_name(output_port), ports[k]);
	}
	
	free(ports);
	
	while (1)
	{
		if (String->strike == 0)
			String->strike = 1.0;
			
		//simple_string_init_hammer(String);
		
		//string_filter_init_impulse(Filt, 0.9);
		usleep(10000000);
	}
	
	jack_client_close(client);
	
	exit (0);
}
