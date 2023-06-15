#include <stdio.h>
#include <stdlib.h>
#include "debug_info.h"
#include "gf_cal.h"
#include "rnd.h"
#include "math.h"

float eb2n0 = 5.0;

float awgn_gen(float eb2n0_snr)
{
	float val = 0;

	#define E_B		1.0
	float w = 0, r = 0;

	w = (float)rand() / (float)RAND_MAX;
	if (w == 1.0)
	{
		w = 0.999999;
	}

	r = gaussrand();
	/*it is for BPSK*/
	val = r * (sqrt((E_B / ((float)MESSAGE_LEN / (float)CODEWORD_LEN)) / pow(10, eb2n0_snr / 10) / 2));

	return val;
}
