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

unsigned char br_g_poly[MAX_POLY_TERM_SIZE];
unsigned char **br_lag_poly;
unsigned char **br_k_poly;
unsigned char br_k_com_poly[MAX_POLY_TERM_SIZE];
unsigned char **br_zk_poly;
unsigned char gamma_c_poly[MAX_POLY_TERM_SIZE];
unsigned char gamma_c_val[CODEWORD_LEN];
unsigned char ***br_m_poly;
unsigned char ****br_v_matrix_poly;
long long br_v_deg[BR_BASIS_NUM][BR_BASIS_NUM];
long long br_v_lp[BR_BASIS_NUM];
long long br_v_lt[BR_BASIS_NUM];
unsigned char br_v_lc[BR_BASIS_NUM];
long long br_v_lp_hist[BR_BASIS_NUM];
unsigned char br_q_poly[BR_BASIS_NUM][MAX_POLY_TERM_SIZE];

int br_poly_init()
{
	long long i = 0, j = 0, k = 0;
	
	br_lag_poly = (unsigned char**)malloc(sizeof(unsigned char*) * CODEWORD_LEN);
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		br_lag_poly[i] = (unsigned char*)malloc(sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		memset(br_lag_poly[i], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}

	br_k_poly = (unsigned char**)malloc(sizeof(unsigned char*) * tst_vct_num);
	for(i = 0; i < tst_vct_num; i++)
	{
		br_k_poly[i] = (unsigned char*)malloc(sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}
	for(i = 0; i < tst_vct_num; i++)
	{
		memset(br_k_poly[i], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}
	memset(br_k_com_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	
	br_zk_poly = (unsigned char**)malloc(sizeof(unsigned char*) * tst_vct_num);
	for(i = 0; i < tst_vct_num; i++)
	{
		br_zk_poly[i] = (unsigned char*)malloc(sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}
	for(i = 0; i < tst_vct_num; i++)
	{
		memset(br_zk_poly[i], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}
	
	memset(gamma_c_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	memset(gamma_c_val, 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
	
	br_m_poly = (unsigned char***)malloc(sizeof(unsigned char**) * tst_vct_num);
	for(i = 0; i < tst_vct_num; i++)
	{
		br_m_poly[i] = (unsigned char**)malloc(sizeof(unsigned char*) * BR_BASIS_NUM);
		for(j = 0; j < BR_BASIS_NUM; j++)
		{
			br_m_poly[i][j] = (unsigned char*)malloc(sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		}
	}
	for(i = 0; i < tst_vct_num; i++)
	{
		for(j = 0; j < BR_BASIS_NUM; j++)
		{
			memset(br_m_poly[i][j], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		}
	}

	br_v_matrix_poly = (unsigned char****)malloc(sizeof(unsigned char***) * tst_vct_num);
	for(i = 0; i < tst_vct_num; i++)
	{
		br_v_matrix_poly[i] = (unsigned char***)malloc(sizeof(unsigned char**) * BR_BASIS_NUM);
		for(j = 0; j < BR_BASIS_NUM; j++)
		{
			br_v_matrix_poly[i][j] = (unsigned char**)malloc(sizeof(unsigned char*) * BR_BASIS_NUM);
			for(k = 0; k < BR_BASIS_NUM; k++)
			{
				br_v_matrix_poly[i][j][k] = (unsigned char*)malloc(sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			}
		}
	}
	for(i = 0; i < tst_vct_num; i++)
	{
		for(j = 0; j < BR_BASIS_NUM; j++)
		{
			for(k = 0; k < BR_BASIS_NUM; k++)
			{
				memset(br_v_matrix_poly[i][j][k], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			}
		}
	}
	
	for(i = 0; i < BR_BASIS_NUM; i++)
	{
		memset(br_v_deg, 0, sizeof(long long) * BR_BASIS_NUM);
	}
	memset(br_v_lp, 0, sizeof(long long) * BR_BASIS_NUM);
	memset(br_v_lt, 0, sizeof(long long) * BR_BASIS_NUM);
	memset(br_v_lc, 0xFF, sizeof(unsigned char) * BR_BASIS_NUM);
	memset(br_v_lp_hist, 0, sizeof(long long) * BR_BASIS_NUM);

	return 0;
}

int br_poly_exit()
{
	long long i = 0, j = 0, k = 0;

	for(i = 0; i < CODEWORD_LEN; i++)
	{
  		free(br_lag_poly[i]);
		br_lag_poly[i] = NULL;
  	}
	free(br_lag_poly);
	br_lag_poly = NULL;

	for(i = 0; i < tst_vct_num; i++)
	{
  		free(br_k_poly[i]);
		br_k_poly[i] = NULL;
  	}
	free(br_k_poly);
	br_k_poly = NULL;
	
	for(i = 0; i < tst_vct_num; i++)
	{
  		free(br_zk_poly[i]);
		br_zk_poly[i] = NULL;
  	}
	free(br_zk_poly);
	br_zk_poly = NULL;
	
	for(i = 0; i < tst_vct_num; i++)
	{
		for(j = 0; j < BR_BASIS_NUM; j++)
		{
			free(br_m_poly[i][j]);
			br_m_poly[i][j] = NULL;
		}
		free(br_m_poly[i]);
		br_m_poly[i] = NULL;
	}
	
	for(i = 0; i < tst_vct_num; i++)
	{
		for(j = 0; j < BR_BASIS_NUM; j++)
		{
			for(k = 0; k < BR_BASIS_NUM; k++)
			{
				free(br_v_matrix_poly[i][j][k]);
				br_v_matrix_poly[i][j][k] = NULL;
			}
			free(br_v_matrix_poly[i][j]);
			br_v_matrix_poly[i][j] = NULL;
		}
		free(br_v_matrix_poly[i]);
		br_v_matrix_poly[i] = NULL;
	}

	return 0;
}

int br_poly_clear()
{
	long long i = 0, j = 0, k = 0;
#if 0
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		memset(br_lag_poly[i], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}
#endif
	for(i = 0; i < tst_vct_num; i++)
	{
		memset(br_k_poly[i], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}
	
	for(i = 0; i < tst_vct_num; i++)
	{
		memset(br_k_poly[i], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}
	
	memset(br_k_com_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	
	for(i = 0; i < tst_vct_num; i++)
	{
		memset(br_zk_poly[i], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}

	memset(gamma_c_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	memset(gamma_c_val, 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
	
	for(i = 0; i < tst_vct_num; i++)
	{
		for(j = 0; j < BR_BASIS_NUM; j++)
		{
			memset(br_m_poly[i][j], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		}
	}
	
	for(i = 0; i < tst_vct_num; i++)
	{
		for(j = 0; j < BR_BASIS_NUM; j++)
		{
			for(k = 0; k < BR_BASIS_NUM; k++)
			{
				memset(br_v_matrix_poly[i][j][k], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			}
		}
	}
	
	for(i = 0; i < BR_BASIS_NUM; i++)
	{
		memset(br_v_deg, 0, sizeof(long long) * BR_BASIS_NUM);
	}
	memset(br_v_lp, 0, sizeof(long long) * BR_BASIS_NUM);
	memset(br_v_lt, 0, sizeof(long long) * BR_BASIS_NUM);
	memset(br_v_lc, 0xFF, sizeof(unsigned char) * BR_BASIS_NUM);
	memset(br_v_lp_hist, 0, sizeof(long long) * BR_BASIS_NUM);
	
	for(i = 0; i < BR_BASIS_NUM; i++)
	{
		memset(br_q_poly[i], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
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
#if 0	
		if(1 == keep_flag[i])
		{
			DEBUG_NOTICE("gamma skip: %d\n", i);
			continue;
		}
#endif
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
#if 0
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

		//her_convert(br_lag_poly[i]);

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

int br_k_poly_construct(long long tv_idx)
{
	long long i = 0, j = 0, k = 0, l = 0;
	unsigned char tmp_poly[MAX_POLY_TERM_SIZE];
	memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	long long intp_idx = 0;

	//for(i = 0; i < tst_vct_num; i++)
	{
		for(j = 0; j < CODEWORD_LEN; j++)
		{
#if 0		
			if(1 == keep_flag[j])
			{
				DEBUG_NOTICE("gamma skip: %d\n", j);
				continue;
			}
#endif

			intp_idx = chnl_rel_order_idx[CODEWORD_LEN - 1 - j];

			if(0xFF == tst_vct[tv_idx][intp_idx])
			{
				continue;
			}

			memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != br_lag_poly[intp_idx][k])
				{
					if((0xFF == tst_vct[tv_idx][intp_idx])
						|| (0xFF == br_lag_poly[intp_idx][k]))
					{
						tmp_poly[k] = 0xFF;
					}
					else if(0x0 == tst_vct[tv_idx][intp_idx])
					{
						tmp_poly[k] = br_lag_poly[intp_idx][k];
					}
					else if(0x0 == br_lag_poly[intp_idx][k])
					{
						tmp_poly[k] = tst_vct[tv_idx][intp_idx];
					}
					else
					{
						tmp_poly[k] = gf_multp(tst_vct[tv_idx][intp_idx], br_lag_poly[intp_idx][k]);
					}

					if((0xFF == br_k_poly[tv_idx][k])
						&& (0xFF == tmp_poly[k]))
					{
						br_k_poly[tv_idx][k] = 0xFF;
					}
					else if(0xFF == br_k_poly[tv_idx][k])
					{
						br_k_poly[tv_idx][k] = tmp_poly[k];
					}
					else if(0xFF == tmp_poly[k])
					{
						br_k_poly[tv_idx][k] = br_k_poly[tv_idx][k];
					}
					else
					{
						br_k_poly[tv_idx][k] = gf_add(br_k_poly[tv_idx][k], tmp_poly[k]);
					}

#if 1
					DEBUG_NOTICE("br_k_poly_cal: %ld %ld | %x %x | %ld %ld %ld | %x\n",
								  tv_idx,
								  intp_idx,
								  tst_vct[tv_idx][intp_idx],
								  br_lag_poly[intp_idx][k],
					              x_term_degree_table[k],
					              y_term_degree_table[k],
					              z_term_degree_table[k],
					              br_k_poly[tv_idx][k]);
#endif
				}
			}

			if((CODEWORD_LEN - 1 - ETA) == j)
			{
				memcpy(br_k_com_poly, br_k_poly[tv_idx], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			}
		}

		//her_convert(br_k_poly[i]);
#if (1 == TEST_MODE)
		for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
		{
			if(0xFF != br_k_poly[tv_idx][k])
			{
				DEBUG_NOTICE("br_k_poly: %d | %ld %ld %ld | %x\n",
				              tv_idx,
				              x_term_degree_table[k],
				              y_term_degree_table[k],
				              z_term_degree_table[k],
				              br_k_poly[tv_idx][k]);
			}
		}
#endif

		memset(br_zk_poly[tv_idx], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		memcpy(br_zk_poly[tv_idx], br_k_poly[tv_idx], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		l = term_search(0, 0, 1);
		br_zk_poly[tv_idx][l] = 0x0;

#if (1 == TEST_MODE)
		for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
		{
			if(0xFF != br_zk_poly[tv_idx][k])
			{
				DEBUG_NOTICE("br_zk_poly: %d | %ld %ld %ld | %x\n",
				              tv_idx,
				              x_term_degree_table[k],
				              y_term_degree_table[k],
				              z_term_degree_table[k],
				              br_zk_poly[tv_idx][k]);
			}
		}
#endif		
	}

	return 0;
}

int br_k_poly_ded(long long base_tv_idx, long long ded_tv_idx)
{
	long long i = 0, j = 0, k = 0;
	unsigned char tmp_poly[MAX_POLY_TERM_SIZE];
	memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	long long intp_idx = 0;
	unsigned char sym_diff = 0xFF;

	memcpy(br_k_poly[ded_tv_idx], br_k_poly[base_tv_idx], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		intp_idx = chnl_rel_order_idx[CODEWORD_LEN - 1 - i];
		if(tst_vct[base_tv_idx][intp_idx] == tst_vct[ded_tv_idx][intp_idx])
		{
			continue;
		}
		
		sym_diff = gf_add(tst_vct[base_tv_idx][intp_idx], tst_vct[ded_tv_idx][intp_idx]);
		DEBUG_NOTICE("sym_diff: %ld | %x\n", intp_idx, sym_diff);
		for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
		{
			if(0xFF != br_lag_poly[intp_idx][j])
			{
				if(0xFF == br_k_poly[ded_tv_idx][j])
				{
					br_k_poly[ded_tv_idx][j] = gf_multp(br_lag_poly[intp_idx][j], sym_diff);
				}
				else
				{
					br_k_poly[ded_tv_idx][j] = gf_add(br_k_poly[ded_tv_idx][j],
												      gf_multp(br_lag_poly[intp_idx][j], sym_diff));
				}
			}
		}
	}

#if (1 == TEST_MODE)
	for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
	{
		if(0xFF != br_k_poly[ded_tv_idx][k])
		{
			DEBUG_NOTICE("br_k_poly: %d | %ld %ld %ld | %x\n",
			              ded_tv_idx,
			              x_term_degree_table[k],
			              y_term_degree_table[k],
			              z_term_degree_table[k],
			              br_k_poly[ded_tv_idx][k]);
		}
	}
#endif

	memset(br_zk_poly[ded_tv_idx], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	memcpy(br_zk_poly[ded_tv_idx], br_k_poly[ded_tv_idx], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	k = term_search(0, 0, 1);
	br_zk_poly[ded_tv_idx][k] = 0x0;
	
#if (1 == TEST_MODE)
	for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
	{
		if(0xFF != br_zk_poly[ded_tv_idx][k])
		{
			DEBUG_NOTICE("br_zk_poly: %d | %ld %ld %ld | %x\n",
						  ded_tv_idx,
						  x_term_degree_table[k],
						  y_term_degree_table[k],
						  z_term_degree_table[k],
						  br_zk_poly[ded_tv_idx][k]);
		}
	}
#endif

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

int br_y_high_deg_trans(unsigned char *tras_poly)
{
	long long i = 0, j = 0;
	unsigned char y_high_deg_flag = 1;
	unsigned char tmp_poly[MAX_POLY_TERM_SIZE];

	while(1 == y_high_deg_flag)
	{
		memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

		for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
		{
			if((0xFF != tras_poly[i])
				&& (GF_Q <= y_term_degree_table[i]))
			{
				j = term_search(x_term_degree_table[i], (y_term_degree_table[i] + 1 - GF_Q), z_term_degree_table[i]);
				tmp_poly[j] = gf_add(tras_poly[i], tmp_poly[j]);
				j = term_search((x_term_degree_table[i] + GF_Q + 1), (y_term_degree_table[i] - GF_Q), z_term_degree_table[i]);
				tmp_poly[j] = gf_add(tras_poly[i], tmp_poly[j]);
				tras_poly[i] = 0xFF;
			}
		}
		for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
		{
			if((0xFF == tras_poly[i])
				&& (0xFF == tmp_poly[i]))
			{
				tras_poly[i] = 0xFF;
			}
			else if((0xFF != tras_poly[i])
					&& (0xFF == tmp_poly[i]))
			{
				tras_poly[i] = tras_poly[i];
			}
			else if((0xFF == tras_poly[i])
					&& (0xFF != tmp_poly[i]))
			{
				tras_poly[i] = tmp_poly[i];
			}
			else
			{
				tras_poly[i] = gf_add(tras_poly[i], tmp_poly[i]);
			}
		}
	
		y_high_deg_flag = 0;
		for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
		{
			if((0xFF != tras_poly[i])
				&& (GF_Q <= y_term_degree_table[i]))
			{
				y_high_deg_flag= 1;
				break;
			}
		}
	}

	return 0;
}

int br_m_poly_construct(long long tv_idx)
{
	long long i = 0, j = 0, k = 0, l = 0;

	//for(i = 0; i < tst_vct_num; i++)
	{
		for(j = 0; j < (BR_BASIS_NUM / 2); j++)
		{
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != br_g_poly[k])
				{
					l = term_search(x_term_degree_table[k], (y_term_degree_table[k] + j), z_term_degree_table[k]);
					br_m_poly[tv_idx][j][l] = br_g_poly[k];
				}
			}
		}
		
		for(j = (BR_BASIS_NUM / 2); j < BR_BASIS_NUM; j++)
		{
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != br_zk_poly[tv_idx][k])
				{
					l = term_search(x_term_degree_table[k], (y_term_degree_table[k] + j - (BR_BASIS_NUM / 2)), z_term_degree_table[k]);
					br_m_poly[tv_idx][j][l] = br_zk_poly[tv_idx][k];
				}
			}
		}
	}
#if 0
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
#endif

	//for(i = 0; i < tst_vct_num; i++)
	{
		for(j = 0; j < BR_BASIS_NUM; j++)
		{
			br_y_high_deg_trans(br_m_poly[tv_idx][j]);
		}		
	}

#if (1 == TEST_MODE)
	//for(i = 0; i < tst_vct_num; i++)
	{
		for(j = 0; j < BR_BASIS_NUM; j++)
		{
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != br_m_poly[tv_idx][j][k])
				{
					DEBUG_NOTICE("br_m_poly: %d %d | %ld %ld %ld | %x\n",
					              tv_idx,
					              j,
								  x_term_degree_table[k],
								  y_term_degree_table[k],
								  z_term_degree_table[k],
								  br_m_poly[tv_idx][j][k]);
				}
			}
		}		
	}
#endif	

	return 0;
}

int br_v_matric_gen(long long tv_idx)
{
	long long i = 0, j = 0, k = 0, l = 0;
	unsigned char tmp_poly[MAX_POLY_TERM_SIZE];
	long long w_z = GF_Q * pole_basis_pow[MESSAGE_LEN - 1][0]
				  + (GF_Q + 1) * pole_basis_pow[MESSAGE_LEN - 1][1];
	long long up_x_deg = 0;

	for(i = 0; i < BR_BASIS_NUM; i++)//for row
	{
		for(j = 0; j < BR_BASIS_NUM; j++)//for col
		{
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if((0xFF != br_m_poly[tv_idx][i][k])
					&& ((j % GF_Q) == y_term_degree_table[k])
					&& ((j / GF_Q) == z_term_degree_table[k]))
				{
					l = term_search(x_term_degree_table[k], 0, 0);
					br_v_matrix_poly[tv_idx][i][j][l] = br_m_poly[tv_idx][i][k];
				}
			}
		}
	}

#if (1 == TEST_MODE)
	for(i = 0; i < BR_BASIS_NUM; i++)//for row
	{
		for(j = 0; j < BR_BASIS_NUM; j++)//for col
		{
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != br_v_matrix_poly[tv_idx][i][j][k])
				{
					DEBUG_NOTICE("br_v_matrix_poly: %ld | %ld %ld | %ld %ld %ld | %x\n",
					             tv_idx,
					             i,
					             j,
					             x_term_degree_table[k],
								 y_term_degree_table[k],
								 z_term_degree_table[k],
								 br_v_matrix_poly[tv_idx][i][j][k]);
				}
			}
		}
	}
#endif

	for(i = 0; i < BR_BASIS_NUM; i++)//for row
	{
		for(j = 0; j < BR_BASIS_NUM; j++)//for col
		{
			memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

			up_x_deg = (((GF_Q + 1) * (j % GF_Q)) + (j / GF_Q * w_z)) / GF_Q;

			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != br_v_matrix_poly[tv_idx][i][j][k])
				{
					l = term_search((x_term_degree_table[k] + up_x_deg), 0, 0);
					tmp_poly[l] = br_v_matrix_poly[tv_idx][i][j][k];
					br_v_matrix_poly[tv_idx][i][j][k] = 0xFF;
				}
			}
			memcpy(br_v_matrix_poly[tv_idx][i][j], tmp_poly, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
#if (1 == TEST_MODE)
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != br_v_matrix_poly[tv_idx][i][j][k])
				{
					DEBUG_NOTICE("br_v_matrix_poly_up: %ld | %ld %ld | %ld | %ld %ld %ld | %x\n",
					             tv_idx,
					             i,
					             j,
					             up_x_deg,
					             x_term_degree_table[k],
								 y_term_degree_table[k],
								 z_term_degree_table[k],
								 br_v_matrix_poly[tv_idx][i][j][k]);
				}
			}
#endif
		}
	}

	return 0;
}

int br_g_poly_gen()
{
	long long i = 0, j = 0, k = 0;
	unsigned char tmp_poly_x[MAX_POLY_TERM_SIZE];

	memset(br_g_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	br_g_poly[0] = 0x0;

	for(i = 0; i < GF_FIELD; i++)
	{
		memset(tmp_poly_x, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		
		for(j = (MAX_POLY_TERM_SIZE - 1); j >= 0; j--)
		{
			if(0xFF != br_g_poly[j])
			{
				tmp_poly_x[j] = gf_multp(br_g_poly[j], power_polynomial_table[i][0]);

				k = term_search((x_term_degree_table[j] + 1), y_term_degree_table[j], z_term_degree_table[j]);
				br_g_poly[k] = br_g_poly[j];
				br_g_poly[j] = 0xFF;
			}
		}
		
		for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
		{
			if(0xFF == tmp_poly_x[j])
			{
				continue;
			}
			else if((0xFF == br_g_poly[j])
					&& (0xFF != tmp_poly_x[j]))
			{
				br_g_poly[j] = tmp_poly_x[j];
			}
			else
			{
				br_g_poly[j] = gf_add(br_g_poly[j], tmp_poly_x[j]);
			}
		}
	}

	//her_convert(br_g_poly);

	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != br_g_poly[i])
		{
			DEBUG_NOTICE("br_g_poly: %ld %ld %ld | %x\n",
						  x_term_degree_table[i],
						  y_term_degree_table[i],
						  z_term_degree_table[i],
						  br_g_poly[i]);
		}
	}

	return 0;
}

int is_popov_form(long long tv_idx)
{
	long long i = 0, j = 0, k = 0;
	long long tmp_deg = 0, max_row_deg = 0;
	unsigned char tmp_lc = 0xFF;
	memset(br_v_lp_hist, 0 , sizeof(long long) * BR_BASIS_NUM);
	int val = 1;

	for(i = 0; i < BR_BASIS_NUM; i++)
	{
		memset(br_v_deg, 0, sizeof(long long) * BR_BASIS_NUM);
	}
	
	for(i = 0; i < BR_BASIS_NUM; i++)
	{
		max_row_deg = 0;
		for(j = 0; j < BR_BASIS_NUM; j++)
		{
			tmp_deg = 0;
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if((0xFF != br_v_matrix_poly[tv_idx][i][j][k])
					&& (tmp_deg < x_term_degree_table[k]))
				{
					tmp_deg = x_term_degree_table[k];
					tmp_lc = br_v_matrix_poly[tv_idx][i][j][k];
				}
			}
			br_v_deg[i][j] = tmp_deg;
			DEBUG_NOTICE("br_v_deg: %ld %ld | %ld\n", i, j, br_v_deg[i][j]);
			if(max_row_deg <= br_v_deg[i][j])//notice choose the right one because of z-term
			{
				max_row_deg = br_v_deg[i][j];
				br_v_lp[i] = j;
				br_v_lt[i] = br_v_deg[i][j];
				br_v_lc[i] = tmp_lc;
			}
		}
		
		br_v_lp_hist[br_v_lp[i]]++;
	}

#if (1 == TEST_MODE)
	for(i = 0; i < BR_BASIS_NUM; i++)
	{
		DEBUG_NOTICE("br_v_leading: %ld | %ld %ld | %x\n", i, br_v_lp[i], br_v_lt[i], br_v_lc[i]);
	}
	for(i = 0; i < BR_BASIS_NUM; i++)
	{
		DEBUG_NOTICE("tmp_lp_hist: %ld | %ld\n", i, br_v_lp_hist[i]);
	}
#endif

	for(i = 0; i < BR_BASIS_NUM; i++)
	{
		if(1 < br_v_lp_hist[i])
		{
			val = 0;
			break;
		}
	}

	return val;
}

int ms_reduction(long long tv_idx)
{
	long long i = 0, j = 0, k = 0, l = 0;
	long long row_t_up = 0, row_t_down = 0, row_t_tmp[2], tmp_cnt = 0;
	memset(row_t_tmp, 0 , sizeof(long long) * 2);
	unsigned char div_coef = 0xFF;
	long long div_deg = 0;
	unsigned char tmp_poly[MAX_POLY_TERM_SIZE];
	memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

	for(i = 0; i < BR_BASIS_NUM; i++)//for row
	{
		if(1 < br_v_lp_hist[i])
		{
			row_t_up = 0;
			row_t_down = 0;
			tmp_cnt = 0;
			memset(row_t_tmp, 0 , sizeof(long long) * 2);
			for(j = 0; j < BR_BASIS_NUM; j++)
			{
				if(br_v_lp[j] == i)
				{
					row_t_tmp[tmp_cnt] = j;
					tmp_cnt++;
				}
				if(2 <= tmp_cnt)
				{
					break;
				}
			}
			if(br_v_lt[row_t_tmp[0]] > br_v_lt[row_t_tmp[1]])
			{
				row_t_up = row_t_tmp[0];
				row_t_down = row_t_tmp[1];
			}
			else
			{
				row_t_up = row_t_tmp[1];
				row_t_down = row_t_tmp[0];
			}
			div_coef = gf_div(br_v_lc[row_t_up], br_v_lc[row_t_down]);
			div_deg = gf_div(br_v_lt[row_t_up], br_v_lt[row_t_down]);
			DEBUG_NOTICE("row_t: %ld | %ld %ld | %ld | %x\n", i, row_t_up, row_t_down, div_deg, div_coef);
			
			for(j = 0; j < BR_BASIS_NUM; j++)
			{
				memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
				for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
				{
					/*k = l + div_deg*/
					if(div_deg <= x_term_degree_table[k])
					{
						l = term_search((x_term_degree_table[k] - div_deg), y_term_degree_table[k], z_term_degree_table[k]);
					}
					else
					{
						l = -1;
					}

					if(-1 == l)
					{
						tmp_poly[k] = br_v_matrix_poly[tv_idx][row_t_up][j][k];
					}
					else
					{
						if((0xFF == br_v_matrix_poly[tv_idx][row_t_up][j][k])
							&& (0xFF == br_v_matrix_poly[tv_idx][row_t_down][j][l]))
						{
							tmp_poly[k] = 0xFF;
						}
						else if((0xFF != br_v_matrix_poly[tv_idx][row_t_up][j][k])
								&& (0xFF == br_v_matrix_poly[tv_idx][row_t_down][j][l]))
						{
							tmp_poly[k] = br_v_matrix_poly[tv_idx][row_t_up][j][k];
						}
						else if((0xFF == br_v_matrix_poly[tv_idx][row_t_up][j][k])
								&& (0xFF != br_v_matrix_poly[tv_idx][row_t_down][j][l]))
						{
							tmp_poly[k] = gf_multp(div_coef, br_v_matrix_poly[tv_idx][row_t_down][j][l]);
						}
						else
						{
							tmp_poly[k] = gf_add(br_v_matrix_poly[tv_idx][row_t_up][j][k],
							                     gf_multp(div_coef, br_v_matrix_poly[tv_idx][row_t_down][j][l]));
						}
					}
					if(0xFF != tmp_poly[k])
					{
						DEBUG_NOTICE("br_v_matrix_poly_red_cal: %ld | %ld %ld | %ld %ld %ld |%x %x | %x\n",
						             tv_idx,
						             i,
						             j,
						             x_term_degree_table[k],
									 y_term_degree_table[k],
									 z_term_degree_table[k],
									 br_v_matrix_poly[tv_idx][row_t_up][j][k],
									 br_v_matrix_poly[tv_idx][row_t_down][j][k],
									 tmp_poly[k]);
					}			 
				}
				memcpy(br_v_matrix_poly[tv_idx][row_t_up][j], tmp_poly, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			}
		}
	}

#if (1 == TEST_MODE)
	for(i = 0; i < BR_BASIS_NUM; i++)
	{
		for(j = 0; j < BR_BASIS_NUM; j++)
		{
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != br_v_matrix_poly[tv_idx][i][j][k])
				{
					DEBUG_NOTICE("br_v_matrix_poly_red: %ld | %ld %ld | %ld %ld %ld | %x\n",
					             tv_idx,
					             i,
					             j,
					             x_term_degree_table[k],
								 y_term_degree_table[k],
								 z_term_degree_table[k],
								 br_v_matrix_poly[tv_idx][i][j][k]);
				}
			}
		}
	}
#endif

	return 0;
}

int br_g_poly_ret(long long tv_idx)
{
	long long i = 0, j = 0, k = 0, l = 0;
	unsigned char tmp_poly[MAX_POLY_TERM_SIZE];
	long long w_z = GF_Q * pole_basis_pow[MESSAGE_LEN - 1][0]
				  + (GF_Q + 1) * pole_basis_pow[MESSAGE_LEN - 1][1];
	long long down_x_deg = 0;

	for(i = 0; i < BR_BASIS_NUM; i++)//for row
	{
		for(j = 0; j < BR_BASIS_NUM; j++)//for col
		{
			memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

			down_x_deg = (((GF_Q + 1) * (j % GF_Q)) + (j / GF_Q * w_z)) / GF_Q;

			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != br_v_matrix_poly[tv_idx][i][j][k])
				{
					if(0 > (x_term_degree_table[k] - down_x_deg))
					{
						DEBUG_SYS("down_x_deg_err: %ld %ld | %ld %ld\n", i, j, x_term_degree_table[k], down_x_deg);
					}
				
					l = term_search((x_term_degree_table[k] - down_x_deg), 0, 0);
					tmp_poly[l] = br_v_matrix_poly[tv_idx][i][j][k];
					br_v_matrix_poly[tv_idx][i][j][k] = 0xFF;
				}
			}
			memcpy(br_v_matrix_poly[tv_idx][i][j], tmp_poly, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
#if (1 == TEST_MODE)
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != br_v_matrix_poly[tv_idx][i][j][k])
				{
					DEBUG_NOTICE("br_v_matrix_poly_down: %ld | %ld %ld | %ld | %ld %ld %ld | %x\n",
					             tv_idx,
					             i,
					             j,
					             down_x_deg,
					             x_term_degree_table[k],
								 y_term_degree_table[k],
								 z_term_degree_table[k],
								 br_v_matrix_poly[tv_idx][i][j][k]);
				}
			}
#endif
		}
	}

	return 0;
}

int br_q_poly_ret(long long tv_idx)
{
	long long i = 0, j = 0, k = 0, l = 0;
	long long min_q_deg = 65536, tmp_deg = 0;
	long long min_idx = -1;

	for(i = 0; i < BR_BASIS_NUM; i++)
	{
		memset(br_q_poly[i], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		
		for(j = 0; j < BR_BASIS_NUM; j++)
		{
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != br_v_matrix_poly[tv_idx][i][j][k])
				{
					l = term_search(x_term_degree_table[k],
					                y_term_degree_table[k] + (j % GF_Q),
					                z_term_degree_table[k] + (j / GF_Q));
					br_q_poly[i][l] = br_v_matrix_poly[tv_idx][i][j][k];
				}
			}
		}
		
		tmp_deg = poly_degree_cal(br_q_poly[i]);
		if(tmp_deg < min_q_deg)
		{
			min_q_deg = tmp_deg;
			min_idx = i;
		}
		
#if (1 == TEST_MODE)
		for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
		{
			if(0xFF != br_q_poly[i][k])
			{
				DEBUG_NOTICE("br_q_poly: %ld | %ld %ld %ld | %x\n",
				             i,
				             x_term_degree_table[k],
							 y_term_degree_table[k],
							 z_term_degree_table[k],
							 br_q_poly[i][k]);
			}
		}
#endif		
	}

	memcpy(min_intp_poly, br_q_poly[min_idx], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	min_intp_idx = min_idx;
	DEBUG_NOTICE("br_min_idx: %ld\n", min_intp_idx);

	return 0;
}

int br_test()
{
	long long i = 0, j = 0;
	int popov_flag = 1;
	long long genus = GF_Q * (GF_Q - 1) / 2;
	long long radius = (CODEWORD_LEN - MESSAGE_LEN - genus + 1 - 1) / 2;
	long long ms_to_cnt = 0;
	unsigned bug_flag = 0;

	DEBUG_NOTICE("br test\n");	
	br_poly_init();
	br_g_poly_gen();
	br_lag_poly_construct();
	
#if 0//for br debug
	memcpy(tst_vct[0], tst_vct[1], sizeof(unsigned char) * CODEWORD_LEN);
#endif
	
	for(i = 0; i < tst_vct_num; i++)
	{
		//br_poly_clear();
		bug_flag = 0;
		ms_to_cnt = 0;

		if(0 == i)
		{
			br_k_poly_construct(0);
		}
		else
		{
			br_k_poly_ded(0, i);
		}
		br_m_poly_construct(i);
		br_v_matric_gen(i);

		popov_flag = is_popov_form(i);
		while(0 == popov_flag)
		{
			ms_reduction(i);
			popov_flag = is_popov_form(i);
			ms_to_cnt++;
			if(CODEWORD_LEN < ms_to_cnt)
			{
				DEBUG_SYS("br_ms_to: %ld\n", ms_to_cnt);
				break;
			}
		}

		br_g_poly_ret(i);
		br_q_poly_ret(i);
		
		her_fac(min_intp_poly, tv_est_msg[i], tv_est_cwd[i]);
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			DEBUG_NOTICE("br_cwd_check: %ld %ld | %x %x\n", i, j, tv_est_cwd[i][j], cwd_poly[j]);
			if((tv_est_cwd[i][j] != cwd_poly[j])
				&& (radius >= tv_err[i]))
			{
				DEBUG_SYS("br_test_err: %ld | %ld\n", i, tv_err[i]);
				bug_flag = 1;
			}
		}
		if(1 == bug_flag)
		{
			DEBUG_SYS("Transmission over Channel:\n");
			for(i = 0; i < symbol_num; i++)
			{
				DEBUG_SYS("rx_symbol[%ld][0] = %f;\n", i, recv_seq[i][0]);
				DEBUG_SYS("rx_symbol[%ld][1] = %f;\n", i, recv_seq[i][1]);
			}
			for(i = 0; i < MESSAGE_LEN; i++)
			{
				DEBUG_SYS("msg_poly[%ld] = 0x%x;\n",
				          i,
				          msg_poly[i]);
			}
		}
	}

	br_poly_exit();

	return 0;
}
