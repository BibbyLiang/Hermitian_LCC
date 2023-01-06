#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "debug_info.h"
#include "math.h"
#include "gf_cal.h"
#include "encoding.h"
#include "channel.h"
#include "lcc_decoding.h"

#define PI	3.1415926

/*col(locator, 0xff~GF_FIELD-2)-row(mesg, CODEWORD_LEN), same as matlab, contrary to most papers*/
float chnl_rel_matrix[CODEWORD_LEN + 1][CODEWORD_LEN];
float chnl_rel_order[CODEWORD_LEN];
long long chnl_rel_order_idx[CODEWORD_LEN];
long long chnl_rel_max_id[CODEWORD_LEN];
long long chnl_rel_scd_id[CODEWORD_LEN];
long long tst_vct_num = 0;
unsigned char **tst_vct;

long long x_term_degree_table[MAX_POLY_TERM_SIZE];
long long y_term_degree_table[MAX_POLY_TERM_SIZE];
long long z_term_degree_table[MAX_POLY_TERM_SIZE];

unsigned char intp_poly_coef[KOT_INTP_POLY_NUM][MAX_POLY_TERM_SIZE];

int chnl_rel_seq_order()
{
	long long i = 0, j = 0;
	float tmp = 0, max_val = 0, scd_val = 0;
	long long max_idx = 0, scd_idx = 0;
	float chnl_rel[CODEWORD_LEN];

	memset(chnl_rel_order, 0, sizeof(float) * CODEWORD_LEN);
	memset(chnl_rel, 0, sizeof(float) * CODEWORD_LEN);
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		chnl_rel_order_idx[i] = i;
	}

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		max_val = 0;
		scd_val = 0;
	
		for(j = 0; j < (CODEWORD_LEN + 1); j++)
		{
			if(max_val < chnl_rel_matrix[j][i])
			{
				max_val = chnl_rel_matrix[j][i];
				max_idx = j;
			}

			chnl_rel_max_id[i] = max_idx;
		}
		for(j = 0; j < (CODEWORD_LEN + 1); j++)
		{
			if((scd_val < chnl_rel_matrix[j][i])
				&& (max_val > chnl_rel_matrix[j][i]))
			{
				scd_val = chnl_rel_matrix[j][i];
				scd_idx = j;
			}

			chnl_rel_scd_id[i] = scd_idx;
		}

		if(0 == scd_val)
		{
			scd_val = 0.00001;
		}
		chnl_rel[i] = max_val / scd_val;
		DEBUG_NOTICE("chnl_rel_order: %f %f %f %ld %ld\n",
		             max_val,
		             scd_val,
		             chnl_rel[i],
		             chnl_rel_max_id[i],
		             chnl_rel_scd_id[i]);
	}

	memcpy(chnl_rel_order, chnl_rel, sizeof(float) * CODEWORD_LEN);

	BubbleSort4(chnl_rel_order, CODEWORD_LEN, chnl_rel_order_idx);

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		DEBUG_NOTICE("chnl_rel_order_idx: %ld %f %ld\n",
					 i,
					 chnl_rel_order[i],
		             chnl_rel_order_idx[i]);
	}

	return 0;
}

int chnl_rel_cal(float **input_seq,
				    long long input_len)
{
	long long i = 0, j = 0, k = 0;
	float n0, d0 = 0, d1 = 0, temp = 0;
	unsigned char tmp_bit = 0;
	/*for BPSK*/
	float map[input_len][2];

	for(i = 0; i < (CODEWORD_LEN + 1); i++)
	{
		memset(chnl_rel_matrix[i], 0, sizeof(float) * CODEWORD_LEN);
	}

	n0 = 1 / ((float)MESSAGE_LEN / (float)CODEWORD_LEN) / (pow(10, eb2n0 / 10) / 2);

	long long seq_idx = 0;
	float base_symbol = 0.0;
	float d_sum = 0;
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		d_sum = 0;

		for(j = 0; j < GF_FIELD; j++)
		{
			d0 = 0;
			d1 = 0;

			for(k = 0; k < GF_Q; k++)
			{
				seq_idx = i * GF_Q + k;
				base_symbol = (float)((power_polynomial_table[j][1] >> k) & 0x1);
				base_symbol = 1 - 2 * base_symbol;
				
				d0 = (input_seq[seq_idx][0] - base_symbol) * (input_seq[seq_idx][0] - base_symbol)
					  + (input_seq[seq_idx][1] - (0.0)) * (input_seq[seq_idx][1] - (0.0));
				if(d0 >= d1)
				{
					d1 = d0;
				}

			}
			
			if(d_sum < d1)
			{
				d_sum = d1;
			}
		}
		DEBUG_INFO("d0: %ld | %f %f| %f\n", i, d0, d1, d_sum);
	}

	for(i = 0; i < input_len; i++)
	{
		d0 = (input_seq[i][0] - (1.0)) * (input_seq[i][0] - (1.0))
				+ (input_seq[i][1] - (0.0)) * (input_seq[i][1] - (0.0));
		d1 = (input_seq[i][0] - (-1.0)) * (input_seq[i][0] - (-1.0))
				+ (input_seq[i][1] - (0.0)) * (input_seq[i][1] - (0.0));

#if 1/*for float precision problem*/
		d0 = d0 / 1e2;
		d1 = d1 / 1e2;
#endif

		/*for BPSK*/
		map[i][0] = 1 / (PI * n0) * exp((-d0) / n0);
        map[i][1] = 1 / (PI * n0) * exp((-d1) / n0);
        DEBUG_INFO("map: %d %f | %f %f | %f %f\n",
        		   i,
        		   n0,
        		   d0,
        		   d1,
        		   map[i][0],
        		   map[i][1]);
        if((0 == map[i][0])
        	|| (0 == map[i][1]))
        {
        	DEBUG_INFO("map: %d %f | %f %f | %f %f\n",
        			  i,
        			  n0,
        			  map[i][0],
        			  map[i][0],
        			  d0,
        			  d1);
        }
	}

	for(i = 0; i < GF_FIELD; i++)
	{
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			chnl_rel_matrix[i][j] = 1;

			for(k = 0; k < GF_Q; k++)
			{
				/*for BPSK*/
				tmp_bit = ((power_polynomial_table[i][1] >> k) & 0x1);
				DEBUG_INFO("chnl_rel_matrix: %d %d | %d | %d | %f %f\n",
		        		   i,
		        		   j,
		        		   k,
		        		   tmp_bit,
		        		   chnl_rel_matrix[i][j],
		        		   map[j * GF_Q + k][tmp_bit]);

				chnl_rel_matrix[i][j] = chnl_rel_matrix[i][j]
									   * map[j * GF_Q + k][tmp_bit];
			}
		}
	}

	for (i = 0; i < CODEWORD_LEN; i++)
    {
		DEBUG_NOTICE("recv: %x\n", recv_poly[i]);

        temp = 0;
        for (j = 0; j < GF_FIELD; j++)
        {
            temp = temp + chnl_rel_matrix[j][i];
        }        

        for (j = 0; j < GF_FIELD; j++)
        {
            chnl_rel_matrix[j][i] = chnl_rel_matrix[j][i] / temp;
            if(chnl_rel_matrix[j][i] != chnl_rel_matrix[j][i])
            {
            	DEBUG_SYS("NAN: %d %d | %f %f | %f %f %f %f %f %f %f %f\n",
            			   i,
            	           j,
            	           temp,
            	           chnl_rel_matrix[j][i],
            	           input_seq[i * GF_Q][0],
            	           input_seq[i * GF_Q + 1][0],
            	           input_seq[i * GF_Q + 2][0],
            	           input_seq[i * GF_Q + 3][0],
            	           input_seq[i * GF_Q + 4][0],
            	           input_seq[i * GF_Q + 5][0],
            	           input_seq[i * GF_Q + 6][0],
            	           input_seq[i * GF_Q + 7][0]);
            }
			DEBUG_NOTICE("chnl_rel: %ld %ld | %f\n", i, power_polynomial_table[j][0], chnl_rel_matrix[j][i]);
        }
    }

#if 0//(1 == TEST_MODE)
	for (i = 0; i < CODEWORD_LEN; i++)
    {
        for (j = 0; j < (CODEWORD_LEN + 1); j++)
        {
            chnl_rel_matrix[j][i] = 0;
        }
    }
    chnl_rel_matrix[3][5] = 1;
    chnl_rel_matrix[5][1] = 1;
    chnl_rel_matrix[0][3] = 1;

	chnl_rel_matrix[1][6] = 0.9;
	chnl_rel_matrix[7][6] = 0.1;

	chnl_rel_matrix[3][4] = 0.8;
	chnl_rel_matrix[5][4] = 0.2;

	chnl_rel_matrix[7][2] = 0.7;
	chnl_rel_matrix[6][2] = 0.3;

	chnl_rel_matrix[3][0] = 0.6;
	chnl_rel_matrix[7][0] = 0.4;
#endif

	chnl_rel_seq_order();

	return 0;
}

int tst_vct_num_cal()
{
	tst_vct_num = (long long)pow(2, ETA);
	DEBUG_NOTICE("tst_vct_num: %ld\n", tst_vct_num);

	return 0;
}

int tst_vct_init()
{
	long long i = 0;

	tst_vct = (unsigned char**)malloc(sizeof(unsigned char*) * tst_vct_num);
	for (i = 0; i < tst_vct_num; i++)
	{
		tst_vct[i] = (unsigned char*)malloc(sizeof(unsigned char) * CODEWORD_LEN);
	}
	
	for (i = 0; i < tst_vct_num; i++)
	{
		memset(tst_vct[i], 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
	}

	DEBUG_NOTICE("tst_vct_init OK\n");

	return 0;
}

int tst_vct_exit()
{
	long long i = 0;

	for (i = 0; i < tst_vct_num; i++)
	{
  		free(tst_vct[i]);
		tst_vct[i] = NULL;
  	}
	free(tst_vct);
	tst_vct = NULL;

	return 0;
}

int tst_vct_form()
{
	long long i = 0, j = 0, pickout_place = 0;
	long long cwd_idx = 0;
	unsigned char pickout_sym = 0xFF;

	for(i = 0; i < tst_vct_num; i++)
	{
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			cwd_idx = chnl_rel_order_idx[j];

			if(ETA < j)
			{
				pickout_place = 0;
				pickout_sym = chnl_rel_max_id[cwd_idx];
			}
			else
			{
				pickout_place = (i >> j) & 0x1;
				if(0 == pickout_place)
				{
					pickout_sym = chnl_rel_max_id[cwd_idx];
				}
				else
				{
					pickout_sym = chnl_rel_scd_id[cwd_idx];
				}
			}
			DEBUG_NOTICE("tst_vct_cal: %ld %ld | %ld %ld | %ld\n",
			             i,
			             j,
			             cwd_idx,
			             pickout_place,
			             pickout_sym);
			tst_vct[i][cwd_idx] = power_polynomial_table[pickout_sym][0];
		}
	}

#if (1 == CFG_DEBUG_NOTICE)
	for(i = 0; i < tst_vct_num; i++)
	{
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			DEBUG_NOTICE("tst_vct: %ld %ld | %x\n",
			             i,
			             j,
			             tst_vct[i][j]);
		}
	}
#endif	

	return 0;
}

int term_degree_table_init()
{
	long long i = 0;
	
	for(i = 0; i < ((MAX_DEGREE + 1) * (MAX_DEGREE + 1) * (MAX_DEGREE + 1)); i++)
	{
		x_term_degree_table[i] = i / ((MAX_DEGREE + 1) * (MAX_DEGREE + 1));
		y_term_degree_table[i] = (i / (MAX_DEGREE + 1)) % (MAX_DEGREE + 1);
		z_term_degree_table[i] = i % (MAX_DEGREE + 1);
#if 0		
		DEBUG_NOTICE("term_degree_table: %ld | %ld %ld %ld\n",
		             i,
		             x_term_degree_table[i],
		             y_term_degree_table[i],
		             z_term_degree_table[i]);
#endif		             
	}

	return 0;
}

int koetter_interpolation_hermitian()
{
	term_degree_table_init();
	poly_init();
	
	return 0;
}

int poly_init()
{
	long long i = 0, j = 0;
	long long y_degree = 0, z_degree = 0;
	
	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		memset(intp_poly_coef[i], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}
	
	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		y_degree = i % GF_Q;
		z_degree = i / GF_Q;
		
		for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
		{
			if((0 == x_term_degree_table[j])
				&& (y_degree == y_term_degree_table[j])
				&& (z_degree == z_term_degree_table[j]))
			{
				intp_poly_coef[i][j] = 0x0;
				DEBUG_NOTICE("intp_poly_coef: %ld | %ld %ld %ld | %x\n",
				             i,
				             x_term_degree_table[j],
				             y_term_degree_table[j],
				             z_term_degree_table[j],
				             intp_poly_coef[i][j]);
				break;
			}
		}
	}
	
	return 0;
}
