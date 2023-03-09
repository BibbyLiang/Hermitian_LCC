#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "debug_info.h"
#include "math.h"
#include "gf_cal.h"
#include "encoding.h"
#include "channel.h"
#include "interpolation.h"
#include "factorization.h"
#include "re_encoding.h"
#include "lcc_decoding.h"

unsigned char lag_poly[CODEWORD_LEN][MAX_POLY_TERM_SIZE];
unsigned char ret_poly[MAX_POLY_TERM_SIZE];
unsigned char ret_cwd_poly[CODEWORD_LEN];
unsigned char ret_msg_poly[MESSAGE_LEN];
unsigned char keep_sym[CODEWORD_LEN];
long long keep_flag[CODEWORD_LEN];
long long keep_cnt = 0;
long long ret_keep_sym_cnt = 0;
unsigned char ret_trans_cwd[CODEWORD_LEN];
unsigned char ret_est_msg[MESSAGE_LEN];
unsigned char ret_est_cwd[CODEWORD_LEN];
unsigned char v_poly[MAX_POLY_TERM_SIZE];
unsigned char v_val[CODEWORD_LEN];
unsigned char t_val[GF_FIELD][GF_FIELD][GF_FIELD];
#if (1 == CFG_Y_RET_STORE)
unsigned char y_t_val[CODEWORD_LEN][CODEWORD_LEN];
#endif

int her_lagrange_poly_construct()
{
	long long i = 0, j = 0, k = 0, l = 0;
	unsigned char x_poly[MAX_POLY_TERM_SIZE], y_poly[MAX_POLY_TERM_SIZE];
	unsigned char tmp_poly_x[MAX_POLY_TERM_SIZE], tmp_poly_y[MAX_POLY_TERM_SIZE];
	unsigned char y_her_curve_find_flag = 0;
	unsigned char x_div = 0x0, y_div = 0x0;
	
	unsigned ret_x_flag = 0;

	for(i = 0; i < CODEWORD_LEN; i++)//af. point
	{
#if 1			
		if(0 == keep_flag[i])
		{
			continue;
		}
#endif

		memset(lag_poly[i], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		memset(x_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		x_poly[0] = 0x0;
		memset(y_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		y_poly[0] = 0x0;
		x_div = 0x0;
		y_div = 0x0;
		
		for(j = 0; j < GF_FIELD; j++)//locator
		{
#if 1
			ret_x_flag = 0;
			for(k = 0; k < CODEWORD_LEN; k++)
			{
				if((1 == keep_flag[k])
					&& (af_pnt[k][0] == power_polynomial_table[j][0]))
				{
					ret_x_flag = 1;
					break;
				}
			}
#endif
			if(af_pnt[i][0] == power_polynomial_table[j][0])
			{
#if 1
				DEBUG_NOTICE("x_skip: %ld | %x | %x %x\n",
				             i,
				             power_polynomial_table[j][0],
				             af_pnt[i][0],
				             af_pnt[i][1]);
#endif
				//continue;
			}
#if 1
			else if(0 == ret_x_flag)
			{
				DEBUG_NOTICE("ret_x_flag_0: %ld | %x | %x %x\n",
				             i,
				             power_polynomial_table[j][0],
				             af_pnt[i][0],
				             af_pnt[i][1]);
			}
#endif
			else
			{
				memset(tmp_poly_x, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
				
				x_div = gf_multp(x_div, gf_add(af_pnt[i][0], power_polynomial_table[j][0]));
				
				for(k = (MAX_POLY_TERM_SIZE - 1); k >= 0; k--)
				{
					if(0xFF != x_poly[k])
					{
						tmp_poly_x[k] = gf_multp(x_poly[k], power_polynomial_table[j][0]);
#if 0
						DEBUG_NOTICE("x_poly_cal_process: %ld %ld | %ld %ld %ld | %x %x\n",
									  i,
									  j,
						              x_term_degree_table[k],
						              y_term_degree_table[k],
						              z_term_degree_table[k],
						              x_poly[k],
						              tmp_poly_x[k]);
#endif
#if (0 == CFG_QUICK_POLY_SEARCH)
						for(l = 0; l < MAX_POLY_TERM_SIZE; l++)
						{
#if 0						
							DEBUG_NOTICE("x_poly_cal_find: %ld %ld | %ld %ld %ld | %x %x\n",
										  i,
										  j,
							              x_term_degree_table[l],
							              y_term_degree_table[l],
							              z_term_degree_table[l],
							              x_poly[l],
							              x_poly[k]);
#endif							              
							if((x_term_degree_table[l] == (x_term_degree_table[k] + 1))
								&& (y_term_degree_table[l] == y_term_degree_table[k])
								&& (z_term_degree_table[l] == z_term_degree_table[k]))
							{
								x_poly[l] = x_poly[k];
								x_poly[k] = 0xFF;
#if 0								
								DEBUG_NOTICE("x_poly_cal_set: %ld %ld | %ld %ld %ld | %x %x\n",
											  i,
											  j,
								              x_term_degree_table[l],
								              y_term_degree_table[l],
								              z_term_degree_table[l],
								              x_poly[l],
								              x_poly[k]);
#endif
								break;
							}
						}
#else
						l = term_search((x_term_degree_table[k] + 1), y_term_degree_table[k], z_term_degree_table[k]);
						x_poly[l] = x_poly[k];
						x_poly[k] = 0xFF;
#endif
					}
				}
#if 1
				for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
				{
					if(0xFF != x_poly[k])
					{
						DEBUG_NOTICE("x_poly_cal: %ld %ld | %ld %ld %ld | %x %x\n",
									  i,
									  j,
						              x_term_degree_table[k],
						              y_term_degree_table[k],
						              z_term_degree_table[k],
						              x_poly[k],
						              tmp_poly_x[k]);
					}
				}
#endif				
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
#if 1
				DEBUG_NOTICE("y_skip: %ld | %x | %x %x\n",
				             i,
				             power_polynomial_table[j][0],
				             af_pnt[i][0],
				             af_pnt[i][1]);
#endif				             
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
#if 1
						DEBUG_NOTICE("y_her_curve_find_flag: %ld | %x | %x %x\n",
						             i,
						             power_polynomial_table[j][0],
						             af_pnt[k][0],
						             af_pnt[k][1]);
#endif						
						break;
					}
				}
				
				y_div = gf_multp(y_div, gf_add(af_pnt[i][1], power_polynomial_table[j][0]));
				
				for(k = (MAX_POLY_TERM_SIZE - 1); k >= 0; k--)
				{
					if((0xFF != y_poly[k])
						&& (1 == y_her_curve_find_flag))//notice this
					{
						tmp_poly_y[k] = gf_multp(y_poly[k], power_polynomial_table[j][0]);
#if (0 == CFG_QUICK_POLY_SEARCH)
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
#else
						l = term_search(x_term_degree_table[k], (y_term_degree_table[k] + 1), z_term_degree_table[k]);
						y_poly[l] = y_poly[k];
						y_poly[k] = 0xFF;
#endif
					}
				}
#if 1
				for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
				{
					if(0xFF != y_poly[k])
					{
						DEBUG_NOTICE("y_poly_cal: %ld %ld | %ld %ld %ld | %x %x\n",
									  i,
									  j,
						              x_term_degree_table[k],
						              y_term_degree_table[k],
						              z_term_degree_table[k],
						              y_poly[k],
						              tmp_poly_y[k]);
					}
				}
#endif				
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
#if 1
		for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
		{
			if(0xFF != x_poly[j])
			{
				DEBUG_NOTICE("x_poly: %ld | %ld %ld %ld | %x\n",
							  i,
				              x_term_degree_table[j],
				              y_term_degree_table[j],
				              z_term_degree_table[j],
				              x_poly[j]);
			}
		}
		for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
		{
			if(0xFF != y_poly[j])
			{
				DEBUG_NOTICE("y_poly: %ld | %ld %ld %ld | %x\n",
							  i,
				              x_term_degree_table[j],
				              y_term_degree_table[j],
				              z_term_degree_table[j],
				              y_poly[j]);
			}
		}
#endif
		//memcpy(lag_poly[i], x_poly, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		for(j = 0; j < MAX_POLY_TERM_SIZE; j++)//for y
		{
			if(0xFF != y_poly[j])
			{
				for(k = 0; k < MAX_POLY_TERM_SIZE; k++)//for x
				{
					if(0xFF != x_poly[k])
					{
#if (0 == CFG_QUICK_POLY_SEARCH)					
						for(l = 0; l < MAX_POLY_TERM_SIZE; l++)//for new term
						{
							if((x_term_degree_table[l] == x_term_degree_table[k])
								&& (y_term_degree_table[l] == (y_term_degree_table[k] + y_term_degree_table[j]))
								&& (z_term_degree_table[l] == z_term_degree_table[k]))
							{
								lag_poly[i][l] = gf_add(lag_poly[i][l],
														gf_multp(x_poly[k], y_poly[j]));
								lag_poly[i][l] = gf_div(lag_poly[i][l], x_div);
								lag_poly[i][l] = gf_div(lag_poly[i][l], y_div);
							}
						}
#else
						l = term_search(x_term_degree_table[k], (y_term_degree_table[k] + y_term_degree_table[j]), z_term_degree_table[k]);
						lag_poly[i][l] = gf_add(lag_poly[i][l],
												gf_multp(x_poly[k], y_poly[j]));
						lag_poly[i][l] = gf_div(lag_poly[i][l], x_div);
						lag_poly[i][l] = gf_div(lag_poly[i][l], y_div);
#endif
					}
				}
			}
		}
		
		her_convert(lag_poly[i]);
#if (1 == TEST_MODE)
		for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
		{
			if(0xFF != lag_poly[i][j])
			{
				DEBUG_NOTICE("lag_poly: %ld | %ld %ld %ld | %x\n",
							  i,
				              x_term_degree_table[j],
				              y_term_degree_table[j],
				              z_term_degree_table[j],
				              lag_poly[i][j]);
			}
		}
#endif		
	}

	return 0;
}

int ret_poly_construct()
{
	long long i = 0, j = 0;
	unsigned char tmp_poly[MAX_POLY_TERM_SIZE];
	memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	memset(ret_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

	memcpy(keep_sym, recv_poly, sizeof(unsigned char) * CODEWORD_LEN);

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if(0 == keep_flag[i])
		{
			continue;
		}
		else
		{
			memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
			{
				if(0xFF != lag_poly[i][j])
				{
					tmp_poly[j] = gf_multp(keep_sym[i], lag_poly[i][j]);
					ret_poly[j] = gf_add(ret_poly[j], tmp_poly[j]);
					DEBUG_NOTICE("ret_poly_cal: %ld | %x %x | %ld %ld %ld | %x\n",
								  i,
								  keep_sym[i],
								  lag_poly[i][j],
					              x_term_degree_table[j],
					              y_term_degree_table[j],
					              z_term_degree_table[j],
					              ret_poly[j]);
				}
			}
		}
	}
#if (1 == TEST_MODE)
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != ret_poly[i])
		{
			DEBUG_NOTICE("ret_poly: %ld %ld %ld | %x\n",
			              x_term_degree_table[i],
			              y_term_degree_table[i],
			              z_term_degree_table[i],
			              ret_poly[i]);
		}
	}
#endif
	return 0;
}

int ret_encoding()
{
	long long i = 0, j = 0, k = 0;
	unsigned char x_val = 0xFF, y_val = 0xFF;;
	unsigned char val = 0xFF;

	memset(ret_cwd_poly, 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
	
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if(1 == keep_flag[i])
		{
			ret_cwd_poly[i] = recv_poly[i];
		}
		else
		{
			ret_cwd_poly[i] = poly_eva_x_y(ret_poly, af_pnt[i][0], af_pnt[i][1]);
		}
		DEBUG_NOTICE("ret_cwd_poly: %ld | %x\n",
		             i,
		             ret_cwd_poly[i]);
	}
	
	/*check ret-cwd*/
#if (1 == TEST_MODE)
	for(i = 0; i < GF_FIELD; i++)
	{
		for(j = 0; j < GF_FIELD; j++)
		{
			val = 0xFF;
			for(k = 0; k < CODEWORD_LEN; k++)
			{
				x_val = gf_pow_cal(af_pnt[i][0], i);
				y_val = gf_pow_cal(af_pnt[i][1], j);
				val = gf_add(val,
							 gf_multp(gf_add(recv_poly[i], ret_cwd_poly[i]),
							 		  gf_multp(x_val, y_val)));
			}
			DEBUG_NOTICE("ret_cwd_check: %ld %ld | %x\n",
			             i,
			             j,
			             val);
		}
	}
#endif
	return 0;
}

int ret_trans()
{
	long long i = 0;

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if(1 == keep_flag[i])
		{
			ret_trans_cwd[i] = 0xFF;
		}
		else
		{
			ret_trans_cwd[i] = gf_add(recv_poly[i], ret_cwd_poly[i]);
			ret_trans_cwd[i] = gf_div(ret_trans_cwd[i], v_val[i]);
		}
		DEBUG_NOTICE("ret_trans_cwd: %ld | %x %x | %x\n",
		             i,
		             recv_poly[i],
		             ret_cwd_poly[i],
		             ret_trans_cwd[i]);
	}

	return 0;
}

int ret_tv_trans()
{
	long long i = 0, j = 0;

	for(i = 0; i < tst_vct_num; i++)
	{
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			if(1 == keep_flag[j])
			{
				ret_tst_vct[i][j] = 0xFF;
			}
			else
			{
				ret_tst_vct[i][j] = gf_add(tst_vct[i][j], ret_cwd_poly[j]);
				ret_tst_vct[i][j] = gf_div(ret_tst_vct[i][j], v_val[j]);
			}
			DEBUG_NOTICE("tv_trans: %ld %ld | %x %x | %x\n",
			             i,
			             j,
			             tst_vct[i][j],
			             ret_cwd_poly[j],
			             ret_tst_vct[i][j]);
		}
	}

	return 0;
}

int ret_cwd_recover()
{
	long long i = 0;

#if (0 == CFG_FAC_FREE)
	memset(ret_est_cwd, 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
	her_encoding(ret_est_msg, ret_est_cwd);
#endif

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		est_cwd_poly[i] = gf_add(ret_est_cwd[i], ret_cwd_poly[i]);
		DEBUG_NOTICE("ret_est_cwd: %ld | %x %x | %x\n",
		             i,
		             ret_est_cwd[i],
		             ret_cwd_poly[i],
		             est_cwd_poly[i]);
	}

	return 0;
}

int v_poly_construct()
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
				&& (0 == keep_flag[j]))
			{
				keep_gf_flag[i] = 0;
				break;
			}
		}
	}

	memset(v_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	v_poly[0] = 0x0;
	memset(v_val, 0xFF, sizeof(unsigned char) * CODEWORD_LEN);

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
				if(0xFF != v_poly[j])
				{
					tmp_poly_x[j] = gf_multp(v_poly[j], power_polynomial_table[i][0]);
#if (0 == CFG_QUICK_POLY_SEARCH)
					for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
					{
						if((x_term_degree_table[k] == (x_term_degree_table[j] + 1))
							&& (y_term_degree_table[k] == y_term_degree_table[j])
							&& (z_term_degree_table[k] == z_term_degree_table[j]))
						{
							v_poly[k] = v_poly[j];
							v_poly[j] = 0xFF;
							break;
						}
					}
#else
					k = term_search((x_term_degree_table[j] + 1), y_term_degree_table[j], z_term_degree_table[j]);
					v_poly[k] = v_poly[j];
					v_poly[j] = 0xFF;
#endif
				}
			}
			
			for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
			{
				if(0xFF == tmp_poly_x[j])
				{
					continue;
				}
				else if((0xFF == v_poly[j])
						&& (0xFF != tmp_poly_x[j]))
				{
					v_poly[j] = tmp_poly_x[j];
				}
				else
				{
					v_poly[j] = gf_add(v_poly[j], tmp_poly_x[j]);
				}
			}
		}
	}

	her_convert(v_poly);

	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != v_poly[i])
		{
			DEBUG_NOTICE("v_poly: %ld %ld %ld | %x\n",
			              x_term_degree_table[i],
			              y_term_degree_table[i],
			              z_term_degree_table[i],
			              v_poly[i]);
		}
	}

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		v_val[i] = poly_eva_x_y(v_poly, af_pnt[i][0], af_pnt[i][1]);
		DEBUG_NOTICE("v_val: %ld | %x %x | %x\n",
					 i,
					 af_pnt[i][0],
					 af_pnt[i][1],
		             v_val[i]);
	}

	return 0;
}

int keep_position_set(long long *keep_poition)
{
	long long i = 0, j = 0;
	long long tmp_cnt = 0;
	long long pair_hist[GF_FIELD];
	memset(pair_hist, 0, sizeof(long long) * GF_FIELD);
	long long sym_idx = 0;

	memset(keep_poition, 0, sizeof(long long) * CODEWORD_LEN);
	keep_cnt = 0;
	ret_keep_sym_cnt = 0;

	for(i = 0; i < (CODEWORD_LEN - ETA); i++)
	{
		sym_idx = chnl_rel_order_idx[CODEWORD_LEN - i - 1];
		for(j = 0; j < GF_FIELD; j++)
		{
			if(af_pnt[sym_idx][0] == power_polynomial_table[j][0])
			{
				pair_hist[j]++;
				break;
			}
		}
	}
	for(i = 0; i < GF_FIELD; i++)
	{
		DEBUG_NOTICE("pair_hist: %x | %ld\n", power_polynomial_table[i][0], pair_hist[i]);
	}

	for(i = 0; i < GF_FIELD; i++)
	{
		if(GF_Q <= pair_hist[i])
		{
			for(j = 0; j < CODEWORD_LEN; j++)
			{
				if((af_pnt[j][0] == power_polynomial_table[i][0])
					//&& ((MESSAGE_LEN / GF_Q * GF_Q) > keep_cnt))
					&& (((MESSAGE_LEN - (GF_Q) * (GF_Q - 1) / 2) / GF_Q * GF_Q) > keep_cnt))//?
				{
					keep_poition[j] = 1;
					keep_cnt++;
				}
			}
		}
	}

#if (1 == FAC_FREE_ERR)
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		keep_poition[i] = 0;
	}
	for(i = 0; i < 32; i++)
	{
		keep_poition[i] = 1;
	}
	keep_cnt = 32;
#endif	
#if (1 == TEST_MODE)
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		DEBUG_NOTICE("keep_poition: %ld | %ld %ld\n", i, keep_cnt, keep_poition[i]);
	}
#endif

#if 1/*re-order intp_seq*/
	float *rel_order;
	long long *rel_order_idx;
	float *unrel_order;
	long long *unrel_order_idx;
	rel_order = (float*)malloc(sizeof(float) * keep_cnt);
	rel_order_idx = (long long*)malloc(sizeof(long long) * keep_cnt);
	unrel_order = (float*)malloc(sizeof(float) * (CODEWORD_LEN - keep_cnt));
	unrel_order_idx = (long long*)malloc(sizeof(long long) * (CODEWORD_LEN - keep_cnt));
	tmp_cnt = 0;
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if(1 == keep_flag[chnl_rel_order_idx[i]])
		{
			rel_order[tmp_cnt] = chnl_rel_order[i];
			rel_order_idx[tmp_cnt] = chnl_rel_order_idx[i];
			DEBUG_NOTICE("rel_sel: %ld | %ld | %f %ld\n",
			             tmp_cnt,
			             i,
			             rel_order[tmp_cnt],
			             rel_order_idx[tmp_cnt]);
			tmp_cnt++;
		}
	}
	BubbleSort4(rel_order, (int)keep_cnt, rel_order_idx);
	tmp_cnt = 0;
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if(0 == keep_flag[chnl_rel_order_idx[i]])
		{
			unrel_order[tmp_cnt] = chnl_rel_order[i];
			unrel_order_idx[tmp_cnt] = chnl_rel_order_idx[i];
			DEBUG_NOTICE("unrel_sel: %ld | %ld | %f %ld\n",
			             tmp_cnt,
			             i,
			             unrel_order[tmp_cnt],
			             unrel_order_idx[tmp_cnt]);
			tmp_cnt++;
		}
	}
	BubbleSort4(unrel_order, (int)(CODEWORD_LEN - keep_cnt), unrel_order_idx);
	tmp_cnt = 0;
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if((CODEWORD_LEN - keep_cnt) > tmp_cnt)
		{
			chnl_rel_order[i] = unrel_order[i];
			chnl_rel_order_idx[i] = unrel_order_idx[i];
		}
		else
		{
			chnl_rel_order[i] = rel_order[i - (CODEWORD_LEN - keep_cnt)];
			chnl_rel_order_idx[i] = rel_order_idx[i - (CODEWORD_LEN - keep_cnt)];
		}
		DEBUG_NOTICE("re-order chnl_rel_order_idx: %ld %f %ld\n",
					 i,
					 chnl_rel_order[i],
		             chnl_rel_order_idx[i]);
		tmp_cnt++;
	}
	free(rel_order);
	rel_order = NULL;
	free(rel_order_idx);
	rel_order_idx = NULL;
	free(unrel_order);
	unrel_order = NULL;
	free(unrel_order_idx);
	unrel_order_idx = NULL;
#endif

	ret_keep_sym_cnt = keep_cnt;//notice this

	return 0;
}

unsigned char ret_fac_free_dev(long long sym_idx, unsigned char dev_flag, unsigned char dev_poly_v)
{
	long long dev_cnt = 0;
	unsigned char q0_val = 0xFF, q1_val = 0xFF, v_val = 0xFF;;
	unsigned char val = 0xFF;

	if(0 == dev_flag)
	{
		fac_dev_init();
	}

	q1_val = poly_eva_x_y(store_q1_dev, af_pnt[sym_idx][0], af_pnt[sym_idx][1]);
	DEBUG_NOTICE("dev_q1_val_cal: %x\n", q1_val);

	if(0xFF != q1_val)//notice this
	{
		if(0 == dev_poly_v)
		{
			q0_val = poly_eva_x_y(store_q0_dev, af_pnt[sym_idx][0], af_pnt[sym_idx][1]);
			val = gf_div(q0_val, q1_val);
		}
		else
		{
			v_val = poly_eva_x_y(store_v_dev, af_pnt[sym_idx][0], af_pnt[sym_idx][1]);
			val = gf_div(v_val, q1_val);
		}
	}
	else
	{
		erasure_flag[sym_idx] = 1;
		dev_to_flag[sym_idx] = 1;
		DEBUG_NOTICE("erasure_check: %ld | %x %x\n", sym_idx, recv_poly[sym_idx], cwd_poly[sym_idx]);
	}

	DEBUG_NOTICE("ret_fac_free_dev: %ld  %ld | %d %d | %x %x %x | %x\n",
	             sym_idx,
	             dev_flag,
	             dev_poly_v,
	             erasure_flag[sym_idx],
	             q0_val,
	             q1_val,
	             v_val,
	             val);

	return val;
}

int test_poly_dev()
{
	long long i = 0;
	unsigned char q0_poly_tmp_1[MAX_POLY_TERM_SIZE], q1_poly_tmp_1[MAX_POLY_TERM_SIZE];
	unsigned char q0_poly_tmp_2[MAX_POLY_TERM_SIZE], q1_poly_tmp_2[MAX_POLY_TERM_SIZE];

	memcpy(q0_poly_tmp_1, store_q0_dev, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	memcpy(q1_poly_tmp_1, store_q1_dev, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

	poly_dev_build(q0_poly_tmp_1, q0_poly_tmp_2);
	poly_dev_build(q1_poly_tmp_1, q1_poly_tmp_2);

	return 0;
}

int ret_fac_free(unsigned char *q0_poly, unsigned char *q1_poly, unsigned char *v_poly)
{
	long long i = 0, j = 0;
	unsigned char q0_val = 0xFF, q1_val = 0xFF;
	unsigned char dev_init_flag = 0;
	long long erasure_keep_cnt = 0;
	memset(erasure_flag, 0, sizeof(unsigned char) * CODEWORD_LEN);
	memset(dev_to_flag, 0, sizeof(unsigned char) * CODEWORD_LEN);
	dev_to_cnt = 0;
	keep_cnt = 0;

	unsigned char q0_all_zero = 1;
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != q0_poly_coef[i])
		{
			q0_all_zero = 0;
			break;
		}
	}
	if(1 == q0_all_zero)
	{
		//memcpy(ret_est_cwd, ret_cwd_poly, sizeof(unsigned char) * CODEWORD_LEN);
		memset(ret_est_cwd, 0xFF, sizeof(unsigned char) * CODEWORD_LEN);

		return 2;
	}

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		q1_val = poly_eva_x_y(q1_poly, af_pnt[i][0], af_pnt[i][1]);
		if(0xFF != q1_val)
		{
			if(0xFF == v_val[i])
			{
				ret_est_cwd[i] = 0xFF;
			}
			else
			{
				q0_val = poly_eva_x_y(q0_poly, af_pnt[i][0], af_pnt[i][1]);

				if(0xFF == q0_val)
				{
					ret_est_cwd[i] = 0xFF;
				}
				else
				{
					ret_est_cwd[i] = gf_multp(v_val[i], q0_val);
					ret_est_cwd[i] = gf_div(ret_est_cwd[i], q1_val);
				}
			}
		}
		else
		{
			if(0xFF == v_val[i])
			{
				ret_est_cwd[i] = ret_fac_free_dev(i, dev_init_flag, 1);
				q0_val = poly_eva_x_y(q0_poly, af_pnt[i][0], af_pnt[i][1]);
				ret_est_cwd[i] = gf_multp(ret_est_cwd[i], q0_val);
				dev_init_flag = 1;
			}
			else
			{
				q0_val = poly_eva_x_y(q0_poly, af_pnt[i][0], af_pnt[i][1]);
				if(0xFF != q0_val)
				{
					ret_est_cwd[i] = 0xFF;

					return 1;
				}
				else
				{
					ret_est_cwd[i] = ret_fac_free_dev(i, dev_init_flag, 0);
					ret_est_cwd[i] = gf_multp(ret_est_cwd[i], v_val[i]);
					dev_init_flag = 1;
				}
			}

			if(1 == erasure_flag[i])
			{
				dev_to_cnt++;
			}
		}
#if (1 == TEST_MODE)
		q0_val = poly_eva_x_y(q0_poly, af_pnt[i][0], af_pnt[i][1]);
#endif		
		DEBUG_NOTICE("ret_fac_free_cal: %ld | %d | %x %x %x | %x\n",
				     i,
				     erasure_flag[i],
				     q0_val,
				     q1_val,
				     v_val[i],
				     ret_est_cwd[i]);
#if (1 == TEST_MODE)
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			if(1 == erasure_flag[j])
			{
				DEBUG_IMPORTANT("erasure_flag_check: %d %ld\n", j, dev_to_cnt);
			}
		}
#endif		
	}

	if(0 != dev_to_cnt)
	{
		for(i = 0; i < CODEWORD_LEN; i++)
		{
			if(1 == erasure_flag[i])
			{
				DEBUG_IMPORTANT("erasure_flag_find: %ld\n", i);
				for(j = 0; j < CODEWORD_LEN; j++)
				{
					if(af_pnt[i][0] == af_pnt[j][0])
					{
						erasure_flag[j] = 1;
						DEBUG_IMPORTANT("erasure_flag_set: %ld\n", j);
					}
				}
			}
		}
#if 0//test for erasure decoding
		for(i = 0; i < CODEWORD_LEN; i++)
		{
			if(0xFF == af_pnt[i][0])
			{
				if(1 >= dev_to_cnt)
				{
					erasure_flag[i] = 1;
					dev_to_flag[i] = 1;
				}
			}
			if(0x0 == af_pnt[i][0])
			{
				//erasure_flag[i] = 1;
			}
			if(0x1 == af_pnt[i][0])
			{
				//erasure_flag[i] = 1;
			}
			if(0x2 == af_pnt[i][0])
			{
				//erasure_flag[i] = 1;
			}
			if(0x3 == af_pnt[i][0])
			{
				//erasure_flag[i] = 1;
			}
			if(0x4 == af_pnt[i][0])
			{
				//erasure_flag[i] = 1;
			}
		}
#endif
		DEBUG_IMPORTANT("dev_to_cnt: %ld\n", dev_to_cnt);
#if (0 == CFG_FAST_RET)			
		fac_her_lagrange_poly_construct();
		fac_ret_poly_construct();
		fac_ret_encoding();
#else
		era_cwd_gen();
#endif
	}

	return 0;
}

int ret_t_val_cal()
{
	long long i = 0, j = 0, k = 0;
	unsigned char up_val = 0xFF, down_val = 0xFF, val = 0xFF;

	/*clear*/
	for(i = 0; i < GF_FIELD; i++)
	{
		for(j = 0; j < GF_FIELD; j++)
		{
			memset(t_val[i][j], 0xFF, sizeof(unsigned char) * GF_FIELD);
		}
	}
#if (1 == CFG_Y_RET_STORE)	
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		memset(y_t_val[i], 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
	}
#endif

	for(i = 0; i < GF_FIELD; i++)
	{
		for(j = 0; j < GF_FIELD; j++)
		{
			down_val = gf_add(power_polynomial_table[j][0], power_polynomial_table[i][0]);
			if(0xFF != down_val)
			{
				for(k = 0; k < GF_FIELD; k++)
				{
					up_val = gf_add(power_polynomial_table[i][0], power_polynomial_table[k][0]);
					t_val[i][j][k] = gf_div(up_val, down_val);
					DEBUG_NOTICE("t_val: %ld %ld %ld | %x\n",
					             i,
					             j,
					             k,
					             t_val[i][j][k]);
				}
			}
		}
	}
#if (1 == CFG_Y_RET_STORE)
	for(k = 0; k < CODEWORD_LEN; k++)
	{
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			val = 0x0;
			for(i = 0; i < CODEWORD_LEN; i++)
			{
				if((af_pnt[i][0] == af_pnt[j][0])
					&& (af_pnt[i][1] != af_pnt[j][1]))
				{
					up_val = gf_add(af_pnt[i][1], af_pnt[k][1]);
					down_val = gf_add(af_pnt[i][1], af_pnt[j][1]);
					val = gf_multp(val,
								   gf_div(up_val, down_val));
				}
			}
			y_t_val[j][k] = val;
		}
	}
#endif
	return 0;
}

int ret_cwd_gen()
{
	long long i = 0, j = 0, k = 0, l = 0, m = 0;
	long long k_idx = 0;
	unsigned char l_val = 0x0, r_val = 0xFF;
	unsigned char y_her_curve_find_flag = 0;
	unsigned char x_i_sel_flag = 0, y_i_sel_flag = 0;
	unsigned char x_j_sel_flag = 0, y_j_sel_flag = 0;
	unsigned char tmp_val = 0xFF;

	for(k = 0; k < CODEWORD_LEN; k++)
	{
		if(1 == keep_flag[k])
		{
			ret_cwd_poly[k] = recv_poly[k];
			continue;
		}
		for(k_idx = 0; k_idx < GF_FIELD; k_idx++)
		{
			if(power_polynomial_table[k_idx][0] == af_pnt[k][0])
			{
				break;
			}
		}
		r_val = 0xFF;
		for(m = 0; m < CODEWORD_LEN; m++)
		{
			if(0 == keep_flag[m])
			{
				continue;
			}
			l_val = 0x0;
			x_j_sel_flag = 0;
#if (0 == CFG_Y_RET_STORE)
			y_j_sel_flag = 0;
#else
			y_j_sel_flag = 1;
#endif
#if (0 == CFG_Y_RET_STORE)
			for(j = 0; j < GF_FIELD; j++)
#else			
			j = m / GF_Q;
#endif			
			{
				if((1 == x_j_sel_flag)
					&& (1 == y_j_sel_flag))
				{
#if (0 == CFG_Y_RET_STORE)				
					break;
#endif					
				}

				if(af_pnt[m][0] == power_polynomial_table[j][0])
				{
					x_j_sel_flag = 1;
					for(i = 0; i < GF_FIELD; i++)
					{
						x_i_sel_flag = 0;
#if (0 == CFG_Y_RET_STORE)						
						for(l = 0; l < CODEWORD_LEN; l++)
#else						
						l = i * GF_Q;
#endif						
						{
							if((1 == keep_flag[l])
#if (0 == CFG_Y_RET_STORE)							
								&& (af_pnt[l][0] == power_polynomial_table[i][0])
#endif								
								&& (af_pnt[m][0] != power_polynomial_table[i][0]))
							{
								x_i_sel_flag = 1;
#if (0 == CFG_Y_RET_STORE)
								break;
#endif
							}
						}
						if(1 == x_i_sel_flag)
						{
							//l_val = gf_multp(l_val, t_val[i][j][k_idx]);
							if((0xFF == t_val[i][j][k_idx])
								|| (0xFF == l_val))
							{
								l_val = 0xFF;
							}
							else if(0x0 == t_val[i][j][k_idx])
							{
								l_val = l_val;
							}
							else if(0x0 == l_val)
							{
								l_val = t_val[i][j][k_idx];
							}
							else
							{
								l_val = gf_multp(l_val, t_val[i][j][k_idx]);
							}
							DEBUG_NOTICE("x_r_val_cal: %ld %ld | %x %x | %x %x %x %x\n",
							             k,
							             m,
							             power_polynomial_table[j][0],
							             power_polynomial_table[i][0],
							             t_val[i][j][k_idx],
							             l_val,
							             recv_poly[j],
							             r_val);
						}
					} 
				}

#if (0 == CFG_Y_RET_STORE)
				if(af_pnt[m][1] == power_polynomial_table[j][0])
				{
					y_j_sel_flag = 1;
					for(i = 0; i < GF_FIELD; i++)
					{
						y_i_sel_flag = 0;
						for(l = 0; l < CODEWORD_LEN; l++)
						{
							if((1 == keep_flag[l])
								&& (af_pnt[l][1] == power_polynomial_table[i][0])
								&& (af_pnt[l][0] == af_pnt[m][0])
								&& (af_pnt[l][1] != af_pnt[m][1]))
							{
								y_i_sel_flag = 1;
								break;
							}
						}
						if(1 == y_i_sel_flag)
						{
							l_val = gf_multp(l_val, t_val[i][j][k_idx]);
							DEBUG_NOTICE("y_r_val_cal: %ld %ld | %x %x | %x %x %x %x\n",
							             k,
							             m,
							             power_polynomial_table[j][0],
							             power_polynomial_table[i][0],
							             t_val[i][j][k_idx],
							             l_val,
							             recv_poly[j],
							             r_val);
						}
					} 
				}		
#endif
			}

#if (1 == CFG_Y_RET_STORE)
			if((0xFF == y_t_val[m][k])
				|| (0xFF == l_val))
			{
				l_val = 0xFF;
			}
			else if(0x0 == y_t_val[m][k])
			{
				l_val = l_val;
			}
			else if(0x0 == l_val)
			{
				l_val = y_t_val[m][k];
			}
			else
			{
				l_val = gf_multp(l_val, y_t_val[m][k]);
			}
#endif

			if((0xFF == recv_poly[m])
				|| (0xFF == l_val))
			{
				tmp_val = 0xFF;
			}
			else if(0x0 == recv_poly[m])
			{
				tmp_val = l_val;
			}
			else if(0x0 == l_val)
			{
				tmp_val = recv_poly[m];
			}
			else
			{
				tmp_val = gf_multp(recv_poly[m], l_val);
			}

			if((0xFF == r_val)
				&& (0xFF == tmp_val))
			{
				r_val = 0xFF;
			}
			else if((0xFF != r_val)
					&& (0xFF == tmp_val))
			{
				r_val = r_val;
			}
			else if((0xFF == r_val)
					&& (0xFF != tmp_val))
			{
				r_val = tmp_val;
			}
			else
			{
				r_val = gf_add(r_val, tmp_val);
			}
		}
		DEBUG_NOTICE("ret_cwd_gen: %ld | %x %x\n",
		             k,
		             r_val,
		             ret_cwd_poly[k]);

		ret_cwd_poly[k] = r_val;
	}

	return 0;
}

int era_cwd_gen()
{
	long long i = 0, j = 0, k = 0, l = 0, m = 0;
	long long k_idx = 0;
	unsigned char l_val = 0x0, r_val = 0xFF;
	unsigned char y_her_curve_find_flag = 0;
	unsigned char x_i_sel_flag = 0, y_i_sel_flag = 0;
	unsigned char x_j_sel_flag = 0, y_j_sel_flag = 0;
	unsigned char tmp_val = 0xFF;

	memcpy(keep_sym, ret_est_cwd, sizeof(unsigned char) * CODEWORD_LEN);

	for(k = 0; k < CODEWORD_LEN; k++)
	{
		if((0 == erasure_flag[k])
			|| (0 == dev_to_flag[k]))
		{
			ret_est_cwd[k] = ret_est_cwd[k];//useless
			continue;
		}
		for(k_idx = 0; k_idx < GF_FIELD; k_idx++)
		{
			if(power_polynomial_table[k_idx][0] == af_pnt[k][0])
			{
				break;
			}
		}
		r_val = 0xFF;
		for(m = 0; m < CODEWORD_LEN; m++)
		{
			if(1 == erasure_flag[m])
			{
				continue;
			}
			l_val = 0x0;
			x_j_sel_flag = 0;
			y_j_sel_flag = 0;
#if (0 == CFG_Y_RET_STORE)			
			for(j = 0; j < GF_FIELD; j++)
#else			
			j = m / GF_Q;
#endif			
			{
				if((1 == x_j_sel_flag)
					&& (1 == y_j_sel_flag))
				{
#if (0 == CFG_Y_RET_STORE)				
					break;
#endif					
				}

				if(af_pnt[m][0] == power_polynomial_table[j][0])
				{
					x_j_sel_flag = 1;
					for(i = 0; i < GF_FIELD; i++)
					{
						x_i_sel_flag = 0;
#if (0 == CFG_Y_RET_STORE)						
						for(l = 0; l < CODEWORD_LEN; l++)
#else						
						l = i * GF_Q;
#endif						
						{
							if((0 == erasure_flag[l])
#if (0 == CFG_Y_RET_STORE)							
								&& (af_pnt[l][0] == power_polynomial_table[i][0])
#endif								
								&& (af_pnt[m][0] != power_polynomial_table[i][0]))
							{
								x_i_sel_flag = 1;
#if (0 == CFG_Y_RET_STORE)								
								break;
#endif								
							}
						}
						if(1 == x_i_sel_flag)
						{
							l_val = gf_multp(l_val, t_val[i][j][k_idx]);
							DEBUG_NOTICE("x_r_val_cal: %ld %ld | %x %x | %x %x %x %x\n",
							             k,
							             m,
							             power_polynomial_table[j][0],
							             power_polynomial_table[i][0],
							             t_val[i][j][k_idx],
							             l_val,
							             keep_sym[j],
							             r_val);
						}
					} 
				}
#if (0 == CFG_Y_RET_STORE)
				if(af_pnt[m][1] == power_polynomial_table[j][0])
				{
					y_j_sel_flag = 1;
					for(i = 0; i < GF_FIELD; i++)
					{
						y_i_sel_flag = 0;
						for(l = 0; l < CODEWORD_LEN; l++)
						{
							if((0 == erasure_flag[l])
								&& (af_pnt[l][1] == power_polynomial_table[i][0])
								&& (af_pnt[l][0] == af_pnt[m][0])
								&& (af_pnt[l][1] != af_pnt[m][1]))
							{
								y_i_sel_flag = 1;
								break;
							}
						}
						if(1 == y_i_sel_flag)
						{
							l_val = gf_multp(l_val, t_val[i][j][k_idx]);
							DEBUG_NOTICE("y_r_val_cal: %ld %ld | %x %x | %x %x %x %x\n",
							             k,
							             m,
							             power_polynomial_table[j][0],
							             power_polynomial_table[i][0],
							             t_val[i][j][k_idx],
							             l_val,
							             keep_sym[j],
							             r_val);
						}
					} 
				}
#endif				
			}
#if (1 == CFG_Y_RET_STORE)
			if((0xFF == y_t_val[m][k])
				|| (0xFF == l_val))
			{
				l_val = 0xFF;
			}
			else if(0x0 == y_t_val[m][k])
			{
				l_val = l_val;
			}
			else if(0x0 == l_val)
			{
				l_val = y_t_val[m][k];
			}
			else
			{
				l_val = gf_multp(l_val, y_t_val[m][k]);
			}
#endif
			if((0xFF == keep_sym[m])
				|| (0xFF == l_val))
			{
				tmp_val = 0xFF;
			}
			else if(0x0 == keep_sym[m])
			{
				tmp_val = l_val;
			}
			else if(0x0 == l_val)
			{
				tmp_val = keep_sym[m];
			}
			else
			{
				tmp_val = gf_multp(keep_sym[m], l_val);
			}
			if((0xFF == r_val)
				&& (0xFF == tmp_val))
			{
				r_val = 0xFF;
			}
			else if((0xFF != r_val)
					&& (0xFF == tmp_val))
			{
				r_val = r_val;
			}
			else if((0xFF == r_val)
					&& (0xFF != tmp_val))
			{
				r_val = tmp_val;
			}
			else
			{
				r_val = gf_add(r_val, tmp_val);
			}
		}
		DEBUG_NOTICE("ret_cwd_gen: %ld | %x %x\n",
		             k,
		             r_val,
		             ret_est_cwd[k]);

		ret_est_cwd[k] = r_val;
	}

	return 0;
}

int re_encoding_transform()
{
	keep_position_set(keep_flag);
	v_poly_construct();

#if (0 == CFG_FAST_RET)	
	her_lagrange_poly_construct();
	ret_poly_construct();
	ret_encoding();
#else
	memset(ret_cwd_poly, 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
	ret_cwd_gen();
#endif
	ret_trans();
	ret_tv_trans();

	return 0;
}
