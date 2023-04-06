#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "time.h"
#include "debug_info.h"
#include "math.h"
#include "gf_cal.h"
#include "encoding.h"
#include "channel.h"
#include "mod.h"
#include "interpolation.h"
#include "factorization.h"
#include "re_encoding.h"
#include "br.h"
#include "lcc_decoding.h"

unsigned char **br_lag_poly;
unsigned char **br_k_poly;
unsigned char gamma_c_poly[MAX_POLY_TERM_SIZE];
unsigned char gamma_c_val[CODEWORD_LEN];
unsigned char ***br_m_poly;

int br_poly_init()
{
	long long i = 0, j = 0;
	
	br_lag_poly = (unsigned char**)malloc(sizeof(unsigned char*) * CODEWORD_LEN);
	for (i = 0; i < CODEWORD_LEN; i++)
	{
		br_lag_poly[i] = (unsigned char*)malloc(sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}
	for (i = 0; i < CODEWORD_LEN; i++)
	{
		memset(br_lag_poly[i], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}

	br_k_poly = (unsigned char**)malloc(sizeof(unsigned char*) * tst_vct_num);
	for (i = 0; i < tst_vct_num; i++)
	{
		br_k_poly[i] = (unsigned char*)malloc(sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}
	for (i = 0; i < tst_vct_num; i++)
	{
		memset(br_k_poly[i], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}
	
	memset(gamma_c_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	memset(gamma_c_val, 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
	
	br_m_poly = (unsigned char***)malloc(sizeof(unsigned char**) * tst_vct_num);
	for (i = 0; i < tst_vct_num; i++)
	{
		br_m_poly[i] = (unsigned char**)malloc(sizeof(unsigned char*) * BR_BASIS_NUM);
		for(j = 0; j < BR_BASIS_NUM; j++)
		{
			br_m_poly[i][j] = (unsigned char*)malloc(sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		}
	}
	for (i = 0; i < tst_vct_num; i++)
	{
		for(j = 0; j < BR_BASIS_NUM; j++)
		{
			memset(br_m_poly[i][j], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		}
	}

	return 0;
}

int br_poly_exit()
{
	long long i = 0, j = 0;

	for (i = 0; i < CODEWORD_LEN; i++)
	{
  		free(br_lag_poly[i]);
		br_lag_poly[i] = NULL;
  	}
	free(br_lag_poly);
	br_lag_poly = NULL;

	for (i = 0; i < tst_vct_num; i++)
	{
  		free(br_k_poly[i]);
		br_k_poly[i] = NULL;
  	}
	free(br_k_poly);
	br_k_poly = NULL;
	
	for (i = 0; i < tst_vct_num; i++)
	{
		for(j = 0; j < BR_BASIS_NUM; j++)
		{
			free(br_m_poly[i][j]);
			br_m_poly[i][j] = NULL;
		}
		free(br_m_poly[i]);
		br_m_poly[i] = NULL;
	}

	return 0;
}

int br_lag_poly_construct()
{
	long long i = 0, j = 0, k = 0, l = 0;
	unsigned char x_poly[MAX_POLY_TERM_SIZE], y_poly[MAX_POLY_TERM_SIZE];
	unsigned char tmp_poly_x[MAX_POLY_TERM_SIZE], tmp_poly_y[MAX_POLY_TERM_SIZE];
	unsigned char y_her_curve_find_flag = 0;
	unsigned char x_div = 0x0, y_div = 0x0;
	unsigned char check_val = 0xFF;

	unsigned ret_x_flag = 0;

	for(i = 0; i < CODEWORD_LEN; i++)//af. point
	{
		if(1 == keep_flag[i])
		{
			DEBUG_NOTICE("gamma skip: %d\n", i);
			continue;
		}

		memset(x_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		x_poly[0] = 0x0;
		memset(y_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		y_poly[0] = 0x0;
		x_div = 0x0;
		y_div = 0x0;

		for(j = 0; j < GF_FIELD; j++)//locator
		{
			ret_x_flag = 0;
			for(k = 0; k < CODEWORD_LEN; k++)
			{
#if 1
				if((0 == keep_flag[k])
					&& (af_pnt[k][0] == power_polynomial_table[j][0]))
#else
				if(af_pnt[k][0] == power_polynomial_table[j][0])
#endif
				{
					ret_x_flag = 1;
					break;
				}
			}

			if(af_pnt[i][0] == power_polynomial_table[j][0])
			{							 
				//continue;
			}
			else if(0 == ret_x_flag)
			{
				DEBUG_NOTICE("ret_x_flag_0: %ld | %x | %x %x\n",
							 i,
							 power_polynomial_table[j][0],
							 af_pnt[i][0],
							 af_pnt[i][1]);
			}
			else
			{
				memset(tmp_poly_x, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

				/*this can be cal. offline*/
				x_div = gf_multp(x_div, gf_add(af_pnt[i][0], power_polynomial_table[j][0]));

				for(k = (MAX_POLY_TERM_SIZE - 1); k >= 0; k--)
				{
					if(0xFF != x_poly[k])
					{
						tmp_poly_x[k] = gf_multp(x_poly[k], power_polynomial_table[j][0]);
						for(l = 0; l < MAX_POLY_TERM_SIZE; l++)
						{
							if((x_term_degree_table[l] == (x_term_degree_table[k] + 1))
								&& (y_term_degree_table[l] == y_term_degree_table[k])
								&& (z_term_degree_table[l] == z_term_degree_table[k]))
							{
								x_poly[l] = x_poly[k];
								x_poly[k] = 0xFF;
								break;
							}
						}
					}
				}
				for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
				{
					if(0xFF == tmp_poly_x[k])
					{
						continue;
					}
					else if((0xFF == x_poly[k])
							&& (0xFF != tmp_poly_x[k]))
					{
						x_poly[k] = tmp_poly_x[k];
					}
					else
					{
						x_poly[k] = gf_add(x_poly[k], tmp_poly_x[k]);
					}
				}
			}

			if(af_pnt[i][1] == power_polynomial_table[j][0])
			{
				//continue;
			}
			else
			{
				memset(tmp_poly_y, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
				y_her_curve_find_flag = 0;
				for(k = 0; k < CODEWORD_LEN; k++)
				{
					if((af_pnt[k][1] == power_polynomial_table[j][0])
						&& (af_pnt[k][0] == af_pnt[i][0]))
					{
						y_her_curve_find_flag = 1;
						break;
					}
				}

				/*this can be cal. offline*/
				y_div = gf_multp(y_div, gf_add(af_pnt[i][1], power_polynomial_table[j][0]));

				for(k = (MAX_POLY_TERM_SIZE - 1); k >= 0; k--)
				{
					if((0xFF != y_poly[k])
						&& (1 == y_her_curve_find_flag))//notice this
					{
						tmp_poly_y[k] = gf_multp(y_poly[k], power_polynomial_table[j][0]);
						for(l = 0; l < MAX_POLY_TERM_SIZE; l++)
						{
							if((x_term_degree_table[l] == x_term_degree_table[k])
								&& (y_term_degree_table[l] == (y_term_degree_table[k] + 1))
								&& (z_term_degree_table[l] == z_term_degree_table[k]))
							{
								y_poly[l] = y_poly[k];
								y_poly[k] = 0xFF;
								break;
							}
						}
					}
				}
				for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
				{
					if(0xFF == tmp_poly_y[k])
					{
						continue;
					}
					else if((0xFF == y_poly[k])
							&& (0xFF != tmp_poly_y[k]))
					{
						y_poly[k] = tmp_poly_y[k];
					}
					else
					{
						y_poly[k] = gf_add(y_poly[k], tmp_poly_y[k]);
					}
				}
			}
		}

		for(j = 0; j < MAX_POLY_TERM_SIZE; j++)//for y
		{
			if(0xFF != y_poly[j])
			{
				for(k = 0; k < MAX_POLY_TERM_SIZE; k++)//for x
				{
					if(0xFF != x_poly[k])
					{
						for(l = 0; l < MAX_POLY_TERM_SIZE; l++)//for new term
						{
							if((x_term_degree_table[l] == x_term_degree_table[k])
								&& (y_term_degree_table[l] == (y_term_degree_table[k] + y_term_degree_table[j]))
								&& (z_term_degree_table[l] == z_term_degree_table[k]))
							{
								br_lag_poly[i][l] = gf_add(br_lag_poly[i][l],
														gf_multp(x_poly[k], y_poly[j]));
								br_lag_poly[i][l] = gf_div(br_lag_poly[i][l], x_div);
								br_lag_poly[i][l] = gf_div(br_lag_poly[i][l], y_div);
							}
						}
					}
				}
			}
		}

		her_convert(br_lag_poly[i]);

#if (1 == TEST_MODE)
		for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
		{
			if(0xFF != br_lag_poly[i][j])
			{
				DEBUG_NOTICE("br_lag_poly: %ld | %ld %ld %ld | %x\n",
							  i,
							  x_term_degree_table[j],
							  y_term_degree_table[j],
							  z_term_degree_table[j],
							  br_lag_poly[i][j]);
			}
		}
#endif		
	}

	return 0;
}

int br_k_poly_construct()
{
	long long i = 0, j = 0, k = 0;
	unsigned char tmp_poly[MAX_POLY_TERM_SIZE];
	memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

	for(i = 0; i < tst_vct_num; i++)
	{
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			if(1 == keep_flag[j])
			{
				DEBUG_NOTICE("gamma skip: %d\n", j);
				continue;
			}

			memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != br_lag_poly[j][k])
				{
					if((0xFF == ret_tst_vct[i][j])
						|| (0xFF == br_lag_poly[j][k]))
					{
						tmp_poly[k] = 0xFF;
					}
					else if(0x0 == ret_tst_vct[i][j])
					{
						tmp_poly[k] = br_lag_poly[j][k];
					}
					else if(0x0 == br_lag_poly[j][k])
					{
						tmp_poly[k] = ret_tst_vct[i][j];
					}
					else
					{
						tmp_poly[k] = gf_multp(ret_tst_vct[i][j], br_lag_poly[j][k]);
					}

					if((0xFF == br_k_poly[i][k])
						&& (0xFF == tmp_poly[k]))
					{
						br_k_poly[i][k] = 0xFF;
					}
					else if(0xFF == br_k_poly[i][k])
					{
						br_k_poly[i][k] = tmp_poly[k];
					}
					else if(0xFF == tmp_poly[k])
					{
						br_k_poly[i][k] = br_k_poly[i][k];
					}
					else
					{
						br_k_poly[i][k] = gf_add(br_k_poly[i][k], tmp_poly[k]);
					}

#if 1
					DEBUG_NOTICE("br_k_poly_cal: %ld | %x %x | %ld %ld %ld | %x\n",
								  i,
								  ret_tst_vct[i][j],
								  br_lag_poly[i][k],
					              x_term_degree_table[k],
					              y_term_degree_table[k],
					              z_term_degree_table[k],
					              br_k_poly[i][k]);
#endif
				}
			}
		}

		her_convert(br_k_poly[i]);

		for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
		{
			if(0xFF != br_k_poly[i][k])
			{
				DEBUG_NOTICE("br_k_poly: %d | %ld %ld %ld | %x\n",
				              i,
				              x_term_degree_table[k],
				              y_term_degree_table[k],
				              z_term_degree_table[k],
				              br_k_poly[i][k]);
			}
		}
	}

	return 0;
}

int br_gamma_c_poly_constrcut()
{
	long long i = 0, j = 0, k = 0;
	unsigned char tmp_poly_x[MAX_POLY_TERM_SIZE];
	unsigned char keep_gf_flag[GF_FIELD];
	memset(keep_gf_flag, 1, sizeof(unsigned char) * GF_FIELD);
	
	for(i = 0; i < GF_FIELD; i++)
	{
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			if((af_pnt[j][0] == power_polynomial_table[i][0])
				&& (1 == keep_flag[j]))
			{
				keep_gf_flag[i] = 0;
				break;
			}
		}
	}

	memset(gamma_c_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	gamma_c_poly[0] = 0x0;
	memset(gamma_c_val, 0xFF, sizeof(unsigned char) * CODEWORD_LEN);

	for(i = 0; i < GF_FIELD; i++)
	{
		if(0 == keep_gf_flag[i])
		{
			continue;
		}
		else
		{
			memset(tmp_poly_x, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			
			for(j = (MAX_POLY_TERM_SIZE - 1); j >= 0; j--)
			{
				if(0xFF != gamma_c_poly[j])
				{
					tmp_poly_x[j] = gf_multp(gamma_c_poly[j], power_polynomial_table[i][0]);

					k = term_search((x_term_degree_table[j] + 1), y_term_degree_table[j], z_term_degree_table[j]);
					gamma_c_poly[k] = gamma_c_poly[j];
					gamma_c_poly[j] = 0xFF;

				}
			}
			
			for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
			{
				if(0xFF == tmp_poly_x[j])
				{
					continue;
				}
				else if((0xFF == gamma_c_poly[j])
						&& (0xFF != tmp_poly_x[j]))
				{
					gamma_c_poly[j] = tmp_poly_x[j];
				}
				else
				{
					gamma_c_poly[j] = gf_add(gamma_c_poly[j], tmp_poly_x[j]);
				}
			}
		}
		
		for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
		{
			if(0xFF != gamma_c_poly[j])
			{
				DEBUG_NOTICE("gamma_c_poly_cal: %ld %ld %ld | %x\n",
							  x_term_degree_table[j],
							  y_term_degree_table[j],
							  z_term_degree_table[j],
							  gamma_c_poly[j]);
			}
		}
	}

	her_convert(gamma_c_poly);

	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != gamma_c_poly[i])
		{
			DEBUG_NOTICE("gamma_c_poly: %ld %ld %ld | %x\n",
						  x_term_degree_table[i],
						  y_term_degree_table[i],
						  z_term_degree_table[i],
						  gamma_c_poly[i]);
		}
	}

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		gamma_c_val[i] = poly_eva_x_y(gamma_c_poly, af_pnt[i][0], af_pnt[i][1]);
		DEBUG_NOTICE("gamma_c_val: %ld | %x %x | %x\n",
					 i,
					 af_pnt[i][0],
					 af_pnt[i][1],
					 gamma_c_val[i]);
	}

	return 0;
}


int br_m_poly_construct()
{
	long long i = 0, j = 0, k = 0, l = 0;

	for(j = 0; j < (BR_BASIS_NUM / 2); j++)
	{
		for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
		{
			if(0xFF != gamma_c_poly[k])
			{
				l = term_search(x_term_degree_table[k], (y_term_degree_table[k] + j), z_term_degree_table[k]);
				br_m_poly[0][j][l] = gamma_c_poly[k];
			}
		}
	}
	for(i = 1; i < tst_vct_num; i++)
	{
		for(j = 0; j < (BR_BASIS_NUM / 2); j++)
		{
			memcpy(br_m_poly[i][j], br_m_poly[0][j], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		}
	}

	for(i = 0; i < tst_vct_num; i++)
	{
		memcpy(br_m_poly[i][(BR_BASIS_NUM / 2)], br_k_poly[i], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		k = term_search(0, 0, 1);
		br_m_poly[i][(BR_BASIS_NUM / 2)][k] = gf_add(br_m_poly[i][(BR_BASIS_NUM / 2)][k], br_k_poly[i][k]);
		for(j = (BR_BASIS_NUM / 2 + 1); j < BR_BASIS_NUM; j++)
		{
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != br_m_poly[i][(BR_BASIS_NUM / 2)][k])
				{
					l = term_search(x_term_degree_table[k], (y_term_degree_table[k] + (j - (BR_BASIS_NUM / 2))), z_term_degree_table[k]);
					br_m_poly[i][j][l] = br_m_poly[i][(BR_BASIS_NUM / 2)][k];
				}
			}
		}
	}

#if (1 == TEST_MODE)
	for(i = 0; i < tst_vct_num; i++)
	{
		for(j = 0; j < BR_BASIS_NUM; j++)
		{
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != br_m_poly[i][j][k])
				{
					DEBUG_NOTICE("br_m_poly: %d %d | %ld %ld %ld | %x\n",
					              i,
					              j,
								  x_term_degree_table[k],
								  y_term_degree_table[k],
								  z_term_degree_table[k],
								  br_m_poly[i][j][k]);
				}
			}
		}		
	}
#endif

	return 0;
}

int br_v_matric_gen()
{
	long long i = 0, j = 0;



	return 0;
}

int br_test()
{
	DEBUG_NOTICE("br test\n");
	br_poly_init();
	br_lag_poly_construct();
	br_k_poly_construct();
	br_gamma_c_poly_constrcut();
	br_m_poly_construct();
	br_poly_exit();

	return 0;
}
