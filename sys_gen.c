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
#include "sys_gen.h"
#include "lcc_decoding.h"

#define K_PLUS_G	CODEWORD_LEN//(MESSAGE_LEN + (GF_Q) * (GF_Q - 1) / 2)

unsigned char sys_intp_point[K_PLUS_G][2];
unsigned char sys_intp_poly[CFG_SYS_GEN_POLY][MAX_POLY_TERM_SIZE];
unsigned char sys_poly_tmp[CFG_SYS_GEN_POLY][MAX_POLY_TERM_SIZE];
unsigned char sys_q0[MAX_POLY_TERM_SIZE], sys_q1[MAX_POLY_TERM_SIZE];

unsigned char sys_proc_col_flag[CODEWORD_LEN];
unsigned char **sys_proc_mat;
unsigned char **sys_sqr_mat;
unsigned char sys_gen_mat[MESSAGE_LEN][CODEWORD_LEN];

int sys_intp_point_gen(long long *rel_order_idx)
{
	long long i = 0;
	long long idx = 0;
	long long x_cod_cnt = 0, x_cod_hist[CODEWORD_LEN / GF_Q];
	memset(x_cod_hist, 0, sizeof(long long) * (CODEWORD_LEN / GF_Q));

	for(i = 0; i < K_PLUS_G; i++)
	{
		memset(sys_intp_point[i], 0xFF, sizeof(unsigned char) * 2);
	}

#if 0
	for(i = 0; i < K_PLUS_G; i++)
	{
		idx = rel_order_idx[CODEWORD_LEN - 1 - i];
		sys_intp_point[i][0] = af_pnt[idx][0];
		sys_intp_point[i][1] = af_pnt[idx][1];
		DEBUG_NOTICE("sys_intp_point: %ld | %ld | %x %x\n",
		             i,
		             idx,
		             sys_intp_point[i][0],
		             sys_intp_point[i][1]);

		if(MESSAGE_LEN > i)
		{
			sys_proc_col_flag[idx] = 1;
		}
	}
#else//for information set selection test
	for(i = 0; i < K_PLUS_G; i++)
	{
		sys_intp_point[i][0] = af_pnt[i][0];
		sys_intp_point[i][1] = af_pnt[i][1];

		if(0xFF == sys_intp_point[i][0])
		{
			x_cod_hist[0]++;
		}
		else
		{
			x_cod_hist[sys_intp_point[i][0] + 1]++;
		}

		DEBUG_NOTICE("sys_intp_point: %ld | %ld | %x %x\n",
		             i,
		             idx,
		             sys_intp_point[i][0],
		             sys_intp_point[i][1]);
	
		if(MESSAGE_LEN > i)
		{
			sys_proc_col_flag[i] = 1;
		}
	}
#endif

#if 0
	for(i = 0; i < MESSAGE_LEN; i++)
	{
		if(0xFF == sys_intp_point[i][0])
		{
			x_cod_hist[0]++;
		}
		else
		{
			x_cod_hist[sys_intp_point[i][0] + 1]++;
		}
	}
#endif

#if 0
	for(i = 0; i < (CODEWORD_LEN / GF_Q); i++)
	{
		if(0 != x_cod_hist[i])
		{
			x_cod_cnt++;
		}
	}
	DEBUG_NOTICE("x_cod_cnt: %ld\n", x_cod_cnt);
#endif

	return 0;
}

int sys_intp_poly_init()
{
	long long i = 0, j = 0;

	for(i = 0; i < CFG_SYS_GEN_POLY; i++)
	{
		memset(sys_intp_poly[i], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}

	for(i = 0; i < CFG_SYS_GEN_POLY; i++)
	{
		j = term_search(0, i, 0);
		sys_intp_poly[i][j] = 0x0;

#if (1 == TEST_MODE)
		for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
		{
			if(0xFF != sys_intp_poly[i][j])
			{
				DEBUG_NOTICE("sys_intp_poly: %ld %ld %ld | %x\n",
				             x_term_degree_table[j],
					         y_term_degree_table[j],
					         z_term_degree_table[j],
					         sys_intp_poly[i][j]);
			}
		}
#endif		
	}

	return 0;
}

int sys_q0q1_get(unsigned char *poly)
{
	long long i = 0, j = 0;
	
	memset(sys_q0, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	memset(sys_q1, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if((0 != y_term_degree_table[i])
			&& (0 == z_term_degree_table[i]))
		{
			j = term_search(x_term_degree_table[i], y_term_degree_table[i], 0);
			sys_q1[j] = poly[i];
		}
		else
		{
			sys_q0[i] = poly[i];
		}
	}

#if (1 == TEST_MODE)	
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != sys_q0[i])
		{
			DEBUG_NOTICE("sys_q0_poly: %ld %ld %ld | %x\n",
			             x_term_degree_table[i],
				         y_term_degree_table[i],
				         z_term_degree_table[i],
				         sys_q0[i]);
		}
	}
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != sys_q1[i])
		{
			DEBUG_NOTICE("sys_q1_poly: %ld %ld %ld | %x\n",
			             x_term_degree_table[i],
				         y_term_degree_table[i],
				         z_term_degree_table[i],
				         sys_q1[i]);
		}
	}
#endif

	return 0;
}

unsigned char sys_hasse_dev_cal(unsigned char *poly, unsigned char *input_point)
{
	unsigned char val = 0xFF;

	val = poly_eva_x_y(poly, input_point[0], input_point[1]);

	return val;
}

int sys_poly_update(unsigned char *poly_tmp,//output in this iter
                         unsigned char *poly_update,//to be update
                         unsigned char *poly_min,//min poly
                         unsigned char hs_dev_min,
                         unsigned char hs_dev_self,
                         unsigned char x_point_val,
                         unsigned char hs_flag)
{
	long long i = 0, j = 0;
	/*warn that there may be stack overflow*/
	unsigned char x_up_poly_min[MAX_POLY_TERM_SIZE];
	memset(x_up_poly_min, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	unsigned char poly_update_val = 0xFF, poly_min_val = 0xFF;

	if(1 == hs_flag)
	{
		for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
		{
			if(0xFF != poly_min[i])
			{
				j = term_search((x_term_degree_table[i] + 1), y_term_degree_table[i], z_term_degree_table[i]);
				x_up_poly_min[j] = poly_min[i];
#if 0				
				DEBUG_NOTICE("x_up_poly_min: %ld %ld %ld | %x\n",
				             x_term_degree_table[j],
				             y_term_degree_table[j],
				             z_term_degree_table[j],
				             x_up_poly_min[j]);
#endif
			}
		}
		
		for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
		{
			if(0xFF != poly_min[i])
			{
				poly_min_val = gf_multp(poly_min[i], x_point_val);
			}
			else
			{
				/*notice that this val should be clear*/
				poly_min_val = 0xFF;
			}
			
			if((0xFF == x_up_poly_min[i])
				&& (0xFF == poly_min_val))
			{
				poly_tmp[i] = 0xFF;
			}
			else if((0xFF == x_up_poly_min[i])
				   && (0xFF != poly_min_val))
			{
				poly_tmp[i] = poly_min_val;
			}
			else if((0xFF != x_up_poly_min[i])
				   && (0xFF == poly_min_val))
			{
				poly_tmp[i] = x_up_poly_min[i];
			}
			else
			{
				poly_tmp[i] = gf_add(x_up_poly_min[i], poly_min_val);
			}
#if 0//(1 == TEST_MODE)		
			if(0xFF != poly_tmp[i])
			{
				/*is it different form the TCOM papaer?*/
				//poly_min[i] = gf_multp(poly_min[i], hs_dev_min);
				DEBUG_NOTICE("poly_min: %ld %ld %ld | %x %x | %x\n",
				             x_term_degree_table[i],
				             y_term_degree_table[i],
				             z_term_degree_table[i],
				             poly_min_val,
				             x_up_poly_min[i],
				             poly_tmp[i]);
			}
#endif		
		}
	}
	else
	{
		for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
		{
			/*notice that these two vals should be clear*/
			poly_update_val = 0xFF, poly_min_val = 0xFF;

			/*is it different form the TCOM papaer?*/
			if(0xFF != poly_update[i])
			{
				//poly_update_val = gf_multp(hs_dev_min, poly_update[i]);
				poly_update_val = poly_update[i];
			}
			if(0xFF != poly_min[i])
			{
				poly_min_val = gf_multp(hs_dev_self, poly_min[i]);
				poly_min_val = gf_div(poly_min_val, hs_dev_min);
			}
			if(0xFF == poly_update_val)
			{
				poly_tmp[i] = poly_min_val;
			}
			else if(0xFF == poly_min_val)
			{
				poly_tmp[i] = poly_update_val;
			}
			else
			{
				poly_tmp[i] = gf_add(poly_update_val, poly_min_val);
			}
		}
	}

#if 0//(1 == TEST_MODE)	
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != poly_tmp[i])
		{
			DEBUG_NOTICE("poly_tmp: %ld %ld %ld | %x\n",
			             x_term_degree_table[i],
			             y_term_degree_table[i],
			             z_term_degree_table[i],
			             poly_tmp[i]);
		}
	}
#endif

	return 0;
}

int sys_has_dev_test()
{
	long long i = 0, j = 0;
	unsigned dev_val = 0xFF;

	for(i = 0; i < MESSAGE_LEN; i++)
	{
		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			dev_val = sys_hasse_dev_cal(sys_intp_poly[j], sys_intp_point[i]);
			DEBUG_NOTICE("sys_has_dev_test: %ld | %ld | %x %x | %x\n",
			             i,
			             j,
			             sys_intp_point[i][0],
			             sys_intp_point[i][1],
			             dev_val);
		}
	}

	return 0;
}

int sys_mat_init()
{
	long long i = 0, j = 0, k = 0;
	unsigned char val = 0xFF;
	unsigned char term_tmp[MAX_POLY_TERM_SIZE];

	memset(sys_proc_col_flag, 0, sizeof(unsigned char) * CODEWORD_LEN);
	for(i = 0; i < MESSAGE_LEN; i++)
	{
		memset(sys_gen_mat[i], 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
	}

	for(i = 0; i < MESSAGE_LEN; i++)
	{
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			memset(term_tmp, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			k = term_search(pole_basis_pow[i][0], pole_basis_pow[i][1], 0);
			term_tmp[k] = 0x0;
			val = poly_eva_x_y(term_tmp, af_pnt[j][0], af_pnt[j][1]);
			sys_gen_mat[i][j] = val;
			DEBUG_NOTICE("sys_gen_mat: %ld %ld | %x\n", i, j, sys_gen_mat[i][j]);
		}
	}

	return 0;
}

int sys_par_col_init(long long row, long long col)
{
	long long i = 0;

	sys_proc_mat = (unsigned char**)malloc(sizeof(unsigned char*) * row);
	for (i = 0; i < row; i++)
	{
		sys_proc_mat[i] = (unsigned char*)malloc(sizeof(unsigned char) * col);
	}

	sys_sqr_mat = (unsigned char**)malloc(sizeof(unsigned char*) * row);
	for (i = 0; i < row; i++)
	{
		sys_sqr_mat[i] = (unsigned char*)malloc(sizeof(unsigned char) * row);
	}

	return 0;
}

int sys_par_col_exit(long long row, long long col)
{
	long long i = 0;

	for (i = 0; i < row; i++)
	{
  		free(sys_proc_mat[i]);
		sys_proc_mat[i] = NULL;
  	}
	free(sys_proc_mat);
	sys_proc_mat = NULL;

	for (i = 0; i < row; i++)
	{
  		free(sys_sqr_mat[i]);
		sys_sqr_mat[i] = NULL;
  	}
	free(sys_sqr_mat);
	sys_sqr_mat = NULL;

	return 0;
}

int sys_par_col_gen(long long row, long long col)
{
	long long i = 0, j = 0;
	long long row_cnt = 0, col_cnt = 0;

	for (i = 0; i < row; i++)
	{
		memset(sys_proc_mat[i], 0xFF, sizeof(unsigned char) * col);
	}

	col_cnt = 0;
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if(0 != sys_proc_col_flag[i])
		{
			row_cnt = 0;
			for(j = (MESSAGE_LEN - row); j < MESSAGE_LEN; j++)
			{
				sys_proc_mat[row_cnt][col_cnt] = sys_gen_mat[j][i];
				//DEBUG_NOTICE("sys_proc_mat: %ld %ld | %x\n", col_cnt, row_cnt, sys_proc_mat[row_cnt][col_cnt]);
				row_cnt++;
			}
			col_cnt++;
			if(col <= col_cnt)
			{
				break;
			}
		}
	}

#if (1 == TEST_MODE)
	for(i = 0; i < row; i++)
	{
		for(j = 0; j < col; j++)
		{
			DEBUG_NOTICE("sys_proc_mat: %ld %ld | %x\n", i, j, sys_proc_mat[i][j]);
		}
	}
#endif

	return 0;
}

int sys_par_mat_trans(long long row, long long col)
{
	long long i = 0, j = 0, k = 0;
	unsigned char val = 0xFF;

	for (i = 0; i < row; i++)
	{
		memset(sys_sqr_mat[i], 0xFF, sizeof(unsigned char) * row);
	}

	for(i = 0; i < row; i++)
	{
		for(j = 0; j < row; j++)
		{
			val = 0xFF;
			for(k = 0; k < col; k++)
			{
				val = gf_add(val,
				             gf_multp(sys_proc_mat[i][k], sys_proc_mat[j][k]));
				//DEBUG_NOTICE("sys_sqr_mat_cal: %ld %ld %ld | %x %x\n", i, j, k, sys_proc_mat[i][k], sys_proc_mat[j][k]);
			}
			sys_sqr_mat[i][j] = val;
			DEBUG_NOTICE("sys_sqr_mat: %ld %ld | %x\n", i, j, sys_sqr_mat[i][j]);
		}
	}

	return 0;
}

int sys_par_mat_gau_elm(long long row)
{
	long long i = 0, j = 0, k = 0;
	unsigned char div_factor = 0xFF;

	for(i = 0; i < row; i++)
	{
		for(j = 0; j < row; j++)
		{
			if(j == i)
			{
				continue;
			}

			div_factor = gf_div(sys_sqr_mat[i][i], sys_sqr_mat[j][i]);
			for(k = 0; k < row; k++)
			{
				sys_sqr_mat[j][k] = gf_add(sys_sqr_mat[j][k],
				                           gf_div(sys_sqr_mat[i][k], div_factor));
				DEBUG_NOTICE("sys_sqr_mat_gau_elm_cal: %ld %ld %ld | %x\n", i, j, k, sys_sqr_mat[i][j]);
			}
		}
	}

#if (1 == TEST_MODE)
	for(i = 0; i < row; i++)
	{
		for(j = 0; j < row; j++)
		{
			DEBUG_NOTICE("sys_sqr_mat_gau_elm: %ld %ld | %x\n", i, j, sys_sqr_mat[i][j]);
		}
	}
#endif

	return 0;
}

int sys_kot_intp()
{
}

int sys_gen_test()
{
	long long i = 0, j = 0, k = 0;
	unsigned char dev_val = 0xFF;
	long long min_idx = -1, min_deg = 65536, deg_tmp = 0;
	unsigned char dev_store[CFG_SYS_GEN_POLY];
	memset(dev_store, 0xFF, sizeof(unsigned char) * CFG_SYS_GEN_POLY);
	long long poly_deg[CFG_SYS_GEN_POLY];
	memset(poly_deg, 0, sizeof(long long) * CFG_SYS_GEN_POLY);

	sys_mat_init();

	for(i = 0; i < CFG_SYS_GEN_POLY; i++)
	{
		memset(sys_poly_tmp[i], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}

	sys_intp_point_gen(chnl_rel_order_idx);
	sys_intp_poly_init();

	sys_par_col_init(2, 3);
	sys_par_col_gen(2, 3);
	sys_par_mat_trans(2, 3);
	sys_par_mat_gau_elm(2);
	sys_par_col_exit(2, 3);

	for(i = 0; i < MESSAGE_LEN; i++)
	{
		min_idx = -1, min_deg = 65536, deg_tmp = 0;
		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			//sys_q0q1_get(sys_intp_poly[j]);
			dev_val = sys_hasse_dev_cal(sys_intp_poly[j], sys_intp_point[i]);
			dev_store[j] = dev_val;
			DEBUG_NOTICE("dev_val: %ld %ld | %x %x | %x\n",
			             i,
			             j,
			             sys_intp_point[i][0],
			             sys_intp_point[i][1],
			             dev_val);
			if(0xFF != dev_val)
			{
				deg_tmp = poly_degree_cal(sys_intp_poly[j]);
				if(min_deg > deg_tmp)
				{
					min_deg = deg_tmp;
					min_idx = j;
				}
			}
		}
		
		DEBUG_NOTICE("min_idx: %ld\n", min_idx);
		
		if(-1 != min_idx)
		{
			for(j = 0; j < CFG_SYS_GEN_POLY; j++)
			{
				memset(sys_poly_tmp[j], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

				if(0xFF == dev_store[j])
				{
					continue;
				}

				if(j == min_idx)
				{
					sys_poly_update(sys_poly_tmp[j],
					                sys_intp_poly[j],
					                sys_intp_poly[min_idx],
					                dev_store[min_idx],
					                dev_store[j],
					                sys_intp_point[i][0],
					                1);
				}
				else
				{
					sys_poly_update(sys_poly_tmp[j],
					                sys_intp_poly[j],
					                sys_intp_poly[min_idx],
					                dev_store[min_idx],
					                dev_store[j],
					                sys_intp_point[i][0],
					                0);
				}
			}
			
			for(j = 0; j < CFG_SYS_GEN_POLY; j++)
			{
				if(0xFF == dev_store[j])
				{
					continue;
				}

				memcpy(sys_intp_poly[j], sys_poly_tmp[j], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);			
			}
		}

#if (1 == TEST_MODE)
		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			poly_deg[j] = poly_degree_cal(sys_intp_poly[j]);
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != sys_intp_poly[j][k])
				{
					DEBUG_NOTICE("sys_intp_poly: %ld | %ld | %ld %ld %ld | %x\n",
								 j,
								 poly_deg[j],
					             x_term_degree_table[k],
					             y_term_degree_table[k],
					             z_term_degree_table[k],
					             sys_intp_poly[j][k]);
				}
			}
		}
		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			DEBUG_NOTICE("sys_poly_deg: %ld | %ld\n", j, poly_deg[j]);
		}
#endif
	}

	sys_has_dev_test();

	return 0;
}
