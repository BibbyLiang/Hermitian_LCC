#include <stdio.h>
#include <string.h>
#include "debug_info.h"
#include "rnd.h"
#include "gf_cal.h"
#include "interpolation.h"
#include "encoding.h"

unsigned char af_pnt[CODEWORD_LEN][2];
long long pole_basis_pow[MESSAGE_LEN][2];
long long mono_degree_table[GF_Q + 1][MAX_DEGREE][2];
long long mono_order_table[GF_Q + 1][MAX_DEGREE][2];
unsigned char msg_poly[MESSAGE_LEN];
unsigned char cwd_poly[CODEWORD_LEN];
unsigned char recv_poly[CODEWORD_LEN];
unsigned char est_msg_poly[MESSAGE_LEN];
unsigned char est_cwd_poly[CODEWORD_LEN];

int affine_points_cal()
{
	long long i = 0, j = 0, k = 0;;
	long long af_pnt_cnt = 0;
	unsigned char x_term = 0x0, y_term = 0x0;

	/*clear*/
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		for(j = 0; j < 2; j++)
		{
			af_pnt[i][j] = 0xFF;
		}
	}

	for(i = 0; i < GF_FIELD; i++)
	{
		x_term = 0x0;
		for(k = 0; k < (GF_Q + 1); k++)
		{
			x_term = gf_multp(x_term, power_polynomial_table[i][0]);
		}

		for(j = 0; j < GF_FIELD; j++)
		{
			y_term = 0x0;
			for(k = 0; k < GF_Q; k++)
			{
				y_term = gf_multp(y_term, power_polynomial_table[j][0]);
			}
#if 0			
			DEBUG_NOTICE("tmp_term: %x %x | %x + %x = %x\n",
						 power_polynomial_table[i][0],
						 power_polynomial_table[j][0],
			             x_term,
			             y_term,
			             gf_add(x_term, y_term));
#endif			
			if(power_polynomial_table[j][0]
			   == gf_add(x_term, y_term))
			{
				if(CODEWORD_LEN <= af_pnt_cnt)
				{
					return 0;
				}

				af_pnt[af_pnt_cnt][0] = power_polynomial_table[i][0];
				af_pnt[af_pnt_cnt][1] = power_polynomial_table[j][0];
				
				DEBUG_NOTICE("af_pnt: %ld | %x %x\n",
				             af_pnt_cnt,
				             af_pnt[af_pnt_cnt][0],
				             af_pnt[af_pnt_cnt][1]);
				             
				af_pnt_cnt++;             
			}
		}
	}

	return 0;
}

int pole_basis_cal()
{
	long long i = 0, j = 0, k = 0, l = 0;
	long long cmp_val = 0;
	
	for(i = 0; i < MESSAGE_LEN; i++)
	{
		for(j = 0; j < 2; j++)
		{
			pole_basis_pow[i][j] = 0;
		}
	}
#if 0	
	for(i = 0; i < (MESSAGE_LEN - 1); i++)
	{
		if(0 == pole_basis_pow[i][0])
		{
			pole_basis_pow[i + 1][0] = pole_basis_pow[i][1] + 1;
			pole_basis_pow[i + 1][1] = 0;
		}
		else
		{
			pole_basis_pow[i + 1][0] = pole_basis_pow[i][0] - 1;
			pole_basis_pow[i + 1][1] = pole_basis_pow[i][1] + 1;;
		}
		DEBUG_NOTICE("pole_basis_pow: %ld | %ld %ld | %ld\n",
		             i + 1,
		             pole_basis_pow[i + 1][0],
		             pole_basis_pow[i + 1][1],
		             GF_Q * pole_basis_pow[i + 1][0] + (GF_Q + 1) * pole_basis_pow[i + 1][1]);
	}
#else
	for(i = 0; i < (GF_Q + 1); i++)
	{
		for(j = 0; j < MAX_DEGREE; j++)
		{
			mono_order_table[i][j][0] = -1;
			mono_order_table[i][j][1] = -1;
			
			mono_degree_table[i][j][0] = i * GF_Q + j * (GF_Q + 1);
		}
	}
	
	for(i = 0; i < (GF_Q + 1); i++)
	{
		for(j = 0; j < MAX_DEGREE; j++)
		{
			cmp_val = 0;
			
			for(k = 0; k < (GF_Q + 1); k++)
			{
				for(l = 0; l < MAX_DEGREE; l++)
				{
					if((k == i)
						&& (l == j))
					{
						continue;
					}
					else
					{
						if(mono_degree_table[i][j][0] > mono_degree_table[k][l][0])
						{
							cmp_val++;
						}
					}
				}
			}
			mono_order_table[i][j][0] = cmp_val;
			DEBUG_NOTICE("mono_order_table: %ld %ld | %ld\n",
			             i,
			             j,
			             mono_order_table[i][j][0]);
		}
	}
	
	for(k = 0; k < MESSAGE_LEN; k++)
	{
		for(i = 0; i < (GF_Q + 1); i++)
		{
			for(j = 0; j < MAX_DEGREE; j++)
			{
				if(k == mono_order_table[i][j][0])
				{
					pole_basis_pow[k][0] = i;
					pole_basis_pow[k][1] = j;
					DEBUG_NOTICE("pole_basis_pow: %ld | %ld %ld | %ld\n",
					             k,
					             pole_basis_pow[k][0],
					             pole_basis_pow[k][1],
					             GF_Q * pole_basis_pow[k][0] + (GF_Q + 1) * pole_basis_pow[k][1]);
				}
			}
		}
	}
#endif
	return 0;
}

int rnd_msg_gen()
{
	long long i = 0, j = 0;

	for(i = 0; i < MESSAGE_LEN; i++)
	{
		j = genrand_int32() % GF_FIELD;
		msg_poly[i] = power_polynomial_table[j][0];
		
		DEBUG_NOTICE("msg_poly: %ld | %x\n",
		             i,
		             msg_poly[i]);
	}

#if 0//(1 == TEST_MODE)
	memset(msg_poly, 0xFF, sizeof(unsigned char) * MESSAGE_LEN);
	msg_poly[0] = 0xFF;
	msg_poly[1] = 0x0;
	msg_poly[2] = 0x1;
	msg_poly[3] = 0x1;
	msg_poly[4] = 0x1;
	for(i = 0; i < MESSAGE_LEN; i++)
	{
		DEBUG_NOTICE("test_msg_poly: %ld | %x\n",
		             i,
		             msg_poly[i]);
	}
#endif

	return 0;
}

int her_encoding(unsigned char *msg, unsigned char *cwd)
{
	long long i = 0, j = 0;
	unsigned char tmp_x = 0x0, tmp_y = 0x0, tmp_prod = 0x0;

#if 0	
	for(i = 0; i < MESSAGE_LEN; i++)
	{
		msg_poly[i] = 0xFF;
	}
#endif	
	
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		cwd[i] = 0xFF;
		
		for(j = 0; j < MESSAGE_LEN; j++)
		{
			tmp_x = gf_pow_cal(af_pnt[i][0], pole_basis_pow[j][0]);
			tmp_y = gf_pow_cal(af_pnt[i][1], pole_basis_pow[j][1]);
			tmp_prod = gf_multp(tmp_x, tmp_y);
			cwd[i] = gf_add(cwd[i], gf_multp(tmp_prod, msg[j]));
#if 0			
			DEBUG_NOTICE("cwd_poly_cal: %ld | %x %x | %ld %ld | %x\n",
		                 i,
		                 af_pnt[i][0],
		                 af_pnt[i][1],
		                 pole_basis_pow[j][0],
		                 pole_basis_pow[j][1],
		                 cwd_poly[i]);
#endif		                 
		}
		DEBUG_NOTICE("cwd_poly: %ld | %x\n",
		             i,
		             cwd[i]);
	}
	
	return 0;
}

int her_convert(unsigned char *poly)
{
	long long i = 0, j = 0;
	int degree_over_flag = 0;
	int convert_flag = 0;
	unsigned char poly_tmp[MAX_POLY_TERM_SIZE];

	while(0 != her_degree_check(poly))
	{
		memcpy(poly_tmp, poly, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	
		for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
		{
			degree_over_flag = 0;

			if((0xFF != poly_tmp[i])
				&& (GF_Q < x_term_degree_table[i]))
			{
				degree_over_flag = 1;
			}
			/*x^(w + 1) = y^w + y*/
			if(0 != degree_over_flag)
			{
#if (0 == CFG_QUICK_POLY_SEARCH)
				for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
				{
					if(((x_term_degree_table[i] - GF_Q - 1) == x_term_degree_table[j])
						&& ((GF_Q + y_term_degree_table[i]) == y_term_degree_table[j])
						&& (z_term_degree_table[i] == z_term_degree_table[j]))
					{
						poly[j] = gf_add(poly_tmp[i], poly[j]);
					}
					if(((x_term_degree_table[i] - GF_Q - 1) == x_term_degree_table[j])
						&& ((1 + y_term_degree_table[i]) == y_term_degree_table[j])
						&& (z_term_degree_table[i] == z_term_degree_table[j]))
					{
						poly[j] = gf_add(poly_tmp[i], poly[j]);
					}
				}
#else
				j = term_search((x_term_degree_table[i] - GF_Q - 1), (GF_Q + y_term_degree_table[i]), z_term_degree_table[i]);
				poly[j] = gf_add(poly_tmp[i], poly[j]);

				j = term_search((x_term_degree_table[i] - GF_Q - 1), (1 + y_term_degree_table[i]), z_term_degree_table[i]);
				poly[j] = gf_add(poly_tmp[i], poly[j]);
#endif
				/*clear this over-degree term*/
				poly[i] = 0xFF;
				convert_flag++;
			}
		}
	}
#if 1
	if(0 != convert_flag)
	{
		DEBUG_NOTICE("convert_flag: %ld\n", convert_flag);
		for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
		{
			if(0xFF != poly_tmp[i])
			{
				DEBUG_NOTICE("poly_convert_before: %ld %ld %ld | %x\n",
				             x_term_degree_table[i],
				             y_term_degree_table[i],
				             z_term_degree_table[i],
				             poly_tmp[i]);
			}
		}
		for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
		{
			if(0xFF != poly[i])
			{
				DEBUG_NOTICE("poly_convert: %ld %ld %ld | %x\n",
				             x_term_degree_table[i],
				             y_term_degree_table[i],
				             z_term_degree_table[i],
				             poly[i]);
			}
		}
	}
#endif

	return 0;
}

unsigned char her_degree_check(unsigned char *poly)
{
	long long i = 0;
	unsigned char degree_over_flag = 0;

	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if((0xFF != poly[i])
			&& (GF_Q < x_term_degree_table[i]))
		{
			degree_over_flag = 1;
			break;
		}
	}

	return degree_over_flag;
}
