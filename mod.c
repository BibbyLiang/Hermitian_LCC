#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "debug_info.h"
#include "gf_cal.h"
#include "mod.h"
#include "channel.h"
#include "encoding.h"

float **recv_seq;
float recv_rel[CODEWORD_LEN];
long long symbol_num = CODEWORD_LEN * GF_Q * BITS_PER_SYMBOL_BPSK;

int mod_init()
{
	long long i = 0;
	long long symbol_num = CODEWORD_LEN * GF_Q * BITS_PER_SYMBOL_BPSK;

	recv_seq = (float**)malloc(sizeof(float*) * symbol_num);
	for (i = 0; i < symbol_num; i++)
	{
		recv_seq[i] = (float*)malloc(sizeof(float) * 2);
	}

	DEBUG_SYS("mod_init OK\n");

	return 0;
}

int mod_exit()
{
	long long i = 0;
	long long symbol_num = CODEWORD_LEN * GF_Q * BITS_PER_SYMBOL_BPSK;

	for (i = 0; i < symbol_num; i++)
	{
  		free(recv_seq[i]);
		recv_seq[i] = NULL;
  	}
	free(recv_seq);
	recv_seq = NULL;

	return 0;
}

int bpsk_mod(unsigned char *input_seq, 
				 unsigned int input_len,
				 float **output_seq,
				 unsigned int output_len)
{
	unsigned int i = 0, j = 0, k = 0;
	unsigned char tmp_input = 0;

	for(i = 0; i < input_len; i++)
	{
		for(j = 0; j < GF_Q; j++)
		{
			/*LSB first, small-endian*/
			if(0xFF == input_seq[i])
			{
				tmp_input = (power_polynomial_table[0][1] >> j) & 0x1;
			}
			else
			{
				tmp_input = (power_polynomial_table[input_seq[i] + 0x1][1] >> j) & 0x1;
			}
			/*0 -> 1, 1 -> -1*/
			output_seq[k][0] = 1 - (float)(2 * tmp_input);
			output_seq[k][1] = 0.0;

			k = k + 1;
			if(k > output_len)
			{
				DEBUG_NOTICE("mod_len_err.\n");
			}
		}
	}

	if(k < output_len)
	{
		for(i = k; i < output_len; i++)
		{
			output_seq[i][0] = 0.0;
			output_seq[i][1] = 0.0;
		}
	}

#if (1 == CFG_DEBUG_INFO)
	DEBUG_INFO("bpsk modulation:\n");
	for(i = 0; i < output_len; i++)
	{
		DEBUG_INFO("%f %f\n", output_seq[i][0], output_seq[i][1]);
	}
	DEBUG_INFO("\n");
#endif

	return 0;
}

int bpsk_demod(float **input_seq, 
				    unsigned int input_len,
				    unsigned char *output_seq,
				    unsigned int output_len)
{
	unsigned int i = 0, j = 0, k = 0;
	unsigned char tmp_bit_decision = 0, tmp_output = 0;
	unsigned char tmp_output_seq[output_len];
	float d0 = 0, d1 = 0;

	j = 0, k = 0;
	for(i = 0; i < input_len; i++)
	{
		/*LSB first, small-endian*/
		/*hard-decision, don't think about Q*/
		//tmp_bit_decision = ((unsigned char)((1 - input_seq[i][0]) / 2)) & 0x1;
		/*soft-like-decision*/
		d0 = (input_seq[i][0] - (1.0)) * (input_seq[i][0] - (1.0))
				+ (input_seq[i][1] - (0.0)) * (input_seq[i][1] - (0.0));
		d1 = (input_seq[i][0] - (-1.0)) * (input_seq[i][0] - (-1.0))
				+ (input_seq[i][1] - (0.0)) * (input_seq[i][1] - (0.0));
		if(d0 > d1)
		{
			tmp_bit_decision = 1;
		}
		else
		{
			tmp_bit_decision = 0;
		}
		tmp_output = tmp_output | (tmp_bit_decision << j);

		j = j + 1;
		if(GF_Q == j)
		{
			tmp_output_seq[k] = tmp_output;
			tmp_output = 0;
			j = 0;
			k = k + 1;
		}
	}

#if 1
	for(i = 0; i < output_len; i++)
	{
		for(j = 0; j < GF_FIELD; j++)
		{
			if(power_polynomial_table[j][1] == tmp_output_seq[i])
			{
				output_seq[i] = power_polynomial_table[j][0];
				break;
			}
		}
	}
#else
	for(i = 0; i < output_len; i++)
	{
		output_seq[i] = power_polynomial_table[tmp_output_seq[i] + 0x1][0];
	}
#endif	

#if (1 == CFG_DEBUG_INFO)
	DEBUG_INFO("bpsk demodulation:\n");
	for(i = 0; i < output_len; i++)
	{
		DEBUG_INFO("%x\n", output_seq[i]);
	}
	DEBUG_INFO("\n");
#endif

	return 0;
}

int trans_over_chnl()
{
	long long i = 0;

	/*modulation*/
	bpsk_mod(cwd_poly,
			 CODEWORD_LEN,
			 recv_seq,
			 symbol_num);

	/*transmission over channel*/
	DEBUG_NOTICE("Transmission over Channel:\n");
	for(i = 0; i < symbol_num; i++)
	{
		if(0 == (i % GF_Q))
		{
			DEBUG_NOTICE("---------------\n");
		}
		recv_seq[i][0] = recv_seq[i][0] + awgn_gen(eb2n0);
		recv_seq[i][1] = recv_seq[i][1] + awgn_gen(eb2n0);
		DEBUG_NOTICE("%f %f\n", recv_seq[i][0], recv_seq[i][1]);
	}
	DEBUG_NOTICE("\n");

	/*demodulation*/
	bpsk_demod((float **)recv_seq,
				symbol_num,
				recv_poly,
				CODEWORD_LEN);

#if (1 == FAC_FREE_ERR)
	memcpy(recv_poly, cwd_poly, sizeof(unsigned char) * CODEWORD_LEN);
#if 0
	recv_poly[2] = 0xd;
	recv_poly[12] = 0x9;
	recv_poly[22] = 0x1;
	recv_poly[36] = 0xc;
	recv_poly[43] = 0xFF;
	recv_poly[45] = 0xd;
	recv_poly[46] = 0x4;
	recv_poly[49] = 0x9;
	recv_poly[61] = 0x1;
#endif
	recv_poly[2] = 0xFF;
	recv_poly[12] = 0xFF;
	recv_poly[22] = 0xFF;
	recv_poly[36] = 0xFF;
	recv_poly[43] = 0xFF;
	recv_poly[45] = 0xFF;
	recv_poly[46] = 0xFF;
	recv_poly[49] = 0xFF;
	recv_poly[61] = 0x0;
#endif

	return 0;
}
