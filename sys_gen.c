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

long long sys_intp_seq[K_PLUS_G];
unsigned char sys_intp_point[K_PLUS_G][2];
unsigned char sys_intp_poly[CFG_SYS_GEN_POLY][MAX_POLY_TERM_SIZE];
unsigned char sys_poly_tmp[CFG_SYS_GEN_POLY][MAX_POLY_TERM_SIZE];
unsigned char sys_q0[MAX_POLY_TERM_SIZE], sys_q1[MAX_POLY_TERM_SIZE];
unsigned char sys_prev_poly[CFG_SYS_GEN_POLY][MAX_POLY_TERM_SIZE];

unsigned char sys_proc_col_flag[CODEWORD_LEN];
unsigned char **sys_proc_mat;
unsigned char **sys_sqr_mat;
unsigned char sys_gen_mat[MESSAGE_LEN][CODEWORD_LEN];
unsigned char sys_tmp_gau_mat[MESSAGE_LEN][CODEWORD_LEN];
long long this_new_col_idx = 0;

unsigned char aj_factor_flag[GF_FIELD][CFG_SYS_GEN_POLY];
unsigned char aj_factor_a_set_flag[CFG_SYS_GEN_POLY];

unsigned char sys_elm_poly[CODEWORD_LEN][MAX_POLY_TERM_SIZE];
unsigned char sys_act_poly[CFG_SYS_GEN_POLY][MAX_POLY_TERM_SIZE];

unsigned char sys_ret_poly[MAX_POLY_TERM_SIZE];

unsigned char sys_ret_cwd[CODEWORD_LEN];

unsigned char sys_intp_min_poly[MAX_POLY_TERM_SIZE];
unsigned char sys_mul_poly[KOT_INTP_POLY_NUM][MAX_POLY_TERM_SIZE];
unsigned char sys_break_mat_poly[KOT_INTP_POLY_NUM][KOT_INTP_POLY_NUM][MAX_POLY_TERM_SIZE];
unsigned char sys_break_poly[KOT_INTP_POLY_NUM][MAX_POLY_TERM_SIZE];

int sys_intp_point_gen(long long *rel_order_idx)
{
	long long i = 0;
	long long idx = 0;

	for(i = 0; i < K_PLUS_G; i++)
	{
		memset(sys_intp_point[i], 0xFF, sizeof(unsigned char) * 2);
	}
	memset(sys_intp_seq, 0, sizeof(long long) * K_PLUS_G);

#if 1
	for(i = 0; i < K_PLUS_G; i++)
	{
		idx = rel_order_idx[CODEWORD_LEN - 1 - i];
		sys_intp_point[i][0] = af_pnt[idx][0];
		sys_intp_point[i][1] = af_pnt[idx][1];
		
		sys_intp_seq[i] = idx;
		
		DEBUG_NOTICE("sys_intp_point: %ld | %ld | %x %x\n",
		             i,
		             idx,
		             sys_intp_point[i][0],
		             sys_intp_point[i][1]);

		if(MESSAGE_LEN > i)
		{
			//sys_proc_col_flag[idx] = 1;
		}
	}
#else//for information set selection test
	for(i = 0; i < K_PLUS_G; i++)
	{
		idx = i;
	
		sys_intp_point[i][0] = af_pnt[idx][0];
		sys_intp_point[i][1] = af_pnt[idx][1];

		sys_intp_seq[i] = idx;

		DEBUG_NOTICE("sys_intp_point: %ld | %ld | %x %x\n",
		             i,
		             idx,
		             sys_intp_point[i][0],
		             sys_intp_point[i][1]);

		if(MESSAGE_LEN > i)
		{
			//sys_proc_col_flag[K_PLUS_G - 1 - i] = 1;
		}
	}
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
		memset(sys_prev_poly[i], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
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

	for(i = 0; i < K_PLUS_G; i++)
	{
		if(0 == sys_proc_col_flag[i])
		{
			continue;
		}

		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			dev_val = sys_hasse_dev_cal(sys_intp_poly[j], af_pnt[i]);
			DEBUG_NOTICE("sys_has_dev_test: %ld | %ld | %x %x | %x\n",
			             i,
			             j,
			             af_pnt[i][0],
			             af_pnt[i][1],
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
		memset(sys_tmp_gau_mat[i], 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
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

#if (0 == CFG_FAST_FULL_GAU)
	col_cnt = 0;
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if(0 != sys_proc_col_flag[i])
		{
			row_cnt = 0;
			for(j = (MESSAGE_LEN - row); j < MESSAGE_LEN; j++)
			{
				sys_proc_mat[row_cnt][col_cnt] = sys_gen_mat[j][i];
				//DEBUG_NOTICE("sys_proc_mat_cal: %ld %ld | %ld %ld | %x\n", i, j, col_cnt, row_cnt, sys_proc_mat[row_cnt][col_cnt]);
				row_cnt++;
			}
			col_cnt++;
			if(col <= col_cnt)
			{
				break;
			}
		}
	}
#else
	if((col <= ((MESSAGE_LEN - (GF_Q) * (GF_Q - 1) / 2) + 1))
		&& (row == MESSAGE_LEN))
	{
		for(i = 0; i < MESSAGE_LEN; i++)
		{
			memcpy(sys_proc_mat[i], sys_tmp_gau_mat[i], sizeof(unsigned char) * (col - 1));
		}
		for(i = 0; i < MESSAGE_LEN; i++)
		{
			sys_proc_mat[i][col - 1] = sys_gen_mat[i][this_new_col_idx];
		}
	}
	else
	{
		col_cnt = 0;
		for(i = 0; i < CODEWORD_LEN; i++)
		{
			if(0 != sys_proc_col_flag[i])
			{
				row_cnt = 0;
				for(j = (MESSAGE_LEN - row); j < MESSAGE_LEN; j++)
				{
					sys_proc_mat[row_cnt][col_cnt] = sys_gen_mat[j][i];
					row_cnt++;
				}
				col_cnt++;
				if(col <= col_cnt)
				{
					break;
				}
			}
		}
	}
#endif

#if 0//(1 == TEST_MODE)
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
	int val = 0, rank = 0;

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

	for(j = 0; j < row; j++)
	{
		for(i = 0; i < row; i++)
		{
			if(0xFF != sys_sqr_mat[i][j])
			{
				rank++;
				break;
			}
		}
	}

	if(row <= rank)
	{
		val = 0;
	}
	else
	{
		val = 1;
	}

	return val;
}

int sys_proc_mat_gau_elm(long long row, long long col)
{
	long long i = 0, j = 0, k = 0;
	long long sel_row_idx = 0;
	unsigned char div_factor = 0xFF;
	int val = 0, rank = 0;

	for(i = 0; i < col; i++)
	{
		for(j = 0; j < row; j++)
		{
			if(0xFF != sys_proc_mat[j][i])
			{
				sel_row_idx = j;
				break;
			}
		}

		if(j >= row)
		{
			continue;
		}

		for(j = 0; j < col; j++)
		{
			if((i == j)
				|| (0xFF == sys_proc_mat[sel_row_idx][j]))
			{
				continue;
			}

			div_factor = gf_div(sys_proc_mat[sel_row_idx][i], sys_proc_mat[sel_row_idx][j]);
			for(k = 0; k < row; k++)
			{
				sys_proc_mat[k][j] = gf_add(sys_proc_mat[k][j],
				                            gf_div(sys_proc_mat[k][i], div_factor));
				DEBUG_NOTICE("sys_sqr_mat_gau_elm_cal: %ld %ld %ld | %ld | %x | %x\n", i, j, k, sel_row_idx, div_factor, sys_proc_mat[k][j]);
			}
		}
	}

#if (1 == TEST_MODE)
	for(i = 0; i < row; i++)
	{
		for(j = 0; j < col; j++)
		{
			DEBUG_NOTICE("sys_proc_mat_gau_elm: %ld %ld | %x\n", i, j, sys_proc_mat[i][j]);
		}
	}
#endif

	for(j = 0; j < col; j++)
	{
		for(i = 0; i < row; i++)
		{
			if(0xFF != sys_proc_mat[i][j])
			{
				rank++;
				break;
			}
		}
	}

	if(col <= rank)
	{
		val = 0;
	}
	else
	{
		val = 1;
	}

	return val;
}

/*this function apply to sys_proc_mat*/
int sys_mat_gau_elm(long long row, long long col)
{
	long long i = 0, j = 0, k = 0;
	long long sel_row_idx = 0;
	unsigned char div_factor = 0xFF;
	int val = 0, rank = 0;

	for(i = 0; i < col; i++)
	{
		for(j = 0; j < row; j++)
		{
			if(0xFF != sys_proc_mat[j][i])
			{
				sel_row_idx = j;
				break;
			}
		}

		if(j >= row)
		{
			continue;
		}

		for(j = 0; j < col; j++)
		{
			if((i == j)
				|| (0xFF == sys_proc_mat[sel_row_idx][j]))
			{
				continue;
			}

			div_factor = gf_div(sys_proc_mat[sel_row_idx][i], sys_proc_mat[sel_row_idx][j]);
			for(k = 0; k < row; k++)
			{
				sys_proc_mat[k][j] = gf_add(sys_proc_mat[k][j],
				                     gf_div(sys_proc_mat[k][i], div_factor));
				//DEBUG_NOTICE("sys_mat_gau_elm_cal: %ld %ld %ld | %ld | %x | %x\n", i, j, k, sel_row_idx, div_factor, sys_proc_mat[k][j]);
			}
		}
	}

#if (1 == TEST_MODE)
	for(i = 0; i < row; i++)
	{
		for(j = 0; j < col; j++)
		{
			//DEBUG_NOTICE("sys_mat_gau_elm: %ld %ld | %x\n", i, j, sys_proc_mat[i][j]);
		}
	}
#endif

	for(j = 0; j < col; j++)
	{
		for(i = 0; i < row; i++)
		{
			if(0xFF != sys_proc_mat[i][j])
			{
				//DEBUG_NOTICE("sel_col_non_zero_idx: %ld %ld\n", i, j);
				rank++;
				break;
			}
		}
	}

	if(col <= rank)
	{
		val = 0;
	}
	else
	{
		val = 1;
	}
	DEBUG_NOTICE("rank: %ld %ld\n", col, rank);

	return val;
}

int sys_check_poly_w(long long *poly_deg)
{
	long long i = 0;
	long long max_deg = 0, min_deg = 65536;
	int val = 0;

	for(i = 0; i < CFG_SYS_GEN_POLY; i++)
	{
		if(max_deg < poly_deg[i])
		{
			max_deg = poly_deg[i];
		}
		if(min_deg > poly_deg[i])
		{
			min_deg = poly_deg[i];
		}
	}

	if(GF_Q >= (max_deg - min_deg))
	{
		val = 0;
	}
	else
	{
		val = 1;
	}

	return val;
}

int sys_check_part_dep(long long col)
{
	long long i = 0;
	int val = 0;

	sys_par_col_init(MESSAGE_LEN - col, col);
	sys_par_col_gen(MESSAGE_LEN - col, col);
	//sys_par_mat_trans(MESSAGE_LEN - col, col);
	//val = sys_par_mat_gau_elm(MESSAGE_LEN - col);
	//val = sys_proc_mat_gau_elm(MESSAGE_LEN - col, col);
	val = sys_mat_gau_elm(MESSAGE_LEN - col, col);
	sys_par_col_exit(MESSAGE_LEN - col, col);

	return val;
}

int sys_check_full_dep(long long col)
{
	long long i = 0;
	int val = 0;

	sys_par_col_init(MESSAGE_LEN, col);
	sys_par_col_gen(MESSAGE_LEN, col);
	val = sys_mat_gau_elm(MESSAGE_LEN, col);
#if (1 == CFG_FAST_FULL_GAU)
	if(0 == val)
	{
		for(i = 0; i < MESSAGE_LEN; i++)
		{
			memcpy(sys_tmp_gau_mat[i], sys_proc_mat[i], sizeof(unsigned char) * col);
		}
	}
#endif	
	sys_par_col_exit(MESSAGE_LEN, col);

	return val;
}

int check_aj_factor(long long locator_idx)
{
	long long i = 0, j = 0, k = 0, l = 0;
	unsigned char tmp_poly[MAX_POLY_TERM_SIZE];
	memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	unsigned char val = 0xFF;

	memset(aj_factor_flag[locator_idx], 0, sizeof(unsigned char) * CFG_SYS_GEN_POLY);

	for(i = 0; i < CFG_SYS_GEN_POLY; i++)
	{	
		aj_factor_flag[locator_idx][i] = 1;
		for(j = 1; j < GF_Q; j++)
		{
			memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if((0xFF != sys_intp_poly[i][k])
					&& (j == y_term_degree_table[k])
					&& (0 == z_term_degree_table[k]))
				{
					l = term_search(x_term_degree_table[k], 0, 0);
					//tmp_poly[l] = gf_add(tmp_poly[l], sys_intp_poly[i][k]);
					tmp_poly[l] = sys_intp_poly[i][k];
				}
			}

#if (1 == TEST_MODE)
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != tmp_poly[k])
				{
					DEBUG_NOTICE("tmp_poly: %ld %ld %ld | %x\n",
								 x_term_degree_table[k],
								 y_term_degree_table[k],
								 z_term_degree_table[k],
								 tmp_poly[k]);
				}
			}
#endif

			val = poly_eva_x_y(tmp_poly, power_polynomial_table[locator_idx][0], 0xFF);
			if(0xFF != val)
			{
				aj_factor_flag[locator_idx][i] = 0;
				break;
			}
		}

		DEBUG_NOTICE("aj_factor_flag: %x %d | %d\n", power_polynomial_table[locator_idx][0], i, aj_factor_flag[locator_idx][i]);
	}

	return 0;
}

int poly_aj_factor_div(unsigned char poly_seq[][MAX_POLY_TERM_SIZE], long long poly_idx, long long locator_idx)
{
	long long i = 0, j = 0;
	unsigned char q0_tmp[MAX_POLY_TERM_SIZE];
	memset(q0_tmp, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	unsigned char res_tmp[MAX_POLY_TERM_SIZE];
	memset(res_tmp, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	unsigned char all_clean_flag = 0;

	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if((0xFF != poly_seq[poly_idx][i])
			&& (0 == z_term_degree_table[i]))
		{
			j = term_search(x_term_degree_table[i], y_term_degree_table[i], 0);
			q0_tmp[j] = poly_seq[poly_idx][i];
		}
	}
	
#if (1 == TEST_MODE)
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != q0_tmp[i])
		{
			DEBUG_NOTICE("q0_tmp: %ld %ld %ld | %x\n",
						 x_term_degree_table[i],
						 y_term_degree_table[i],
						 z_term_degree_table[i],
						 q0_tmp[i]);
		}
	}
#endif

	all_clean_flag = 0;
	memset(res_tmp, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	while(0 == all_clean_flag)
	{
		for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
		{
			if(0xFF != q0_tmp[i])
			{
#if 1
				DEBUG_NOTICE("q0_tmp_res: %ld %ld %ld | %x\n",
							 x_term_degree_table[i],
							 y_term_degree_table[i],
							 z_term_degree_table[i],
							 q0_tmp[i]);
#endif				
				all_clean_flag = 0;
				break;
			}
		}
		if(MAX_POLY_TERM_SIZE <= i)
		{
			all_clean_flag = 1;
			break;
		}

		for(i = (MAX_POLY_TERM_SIZE - 1); i >= 0; i--)
		{
			if(0xFF != q0_tmp[i])
			{
				j = term_search(x_term_degree_table[i] - 1, y_term_degree_table[i], 0);
				res_tmp[j] = gf_add(q0_tmp[i], res_tmp[j]);

				q0_tmp[i] = 0xFF;
				q0_tmp[j] = gf_add(q0_tmp[j],
								   gf_multp(res_tmp[j], power_polynomial_table[locator_idx][0]));
#if 1
				DEBUG_NOTICE("q0_tmp_res_cal: %ld %ld %ld | %x || %ld %ld %ld | %x || %ld %ld %ld | %x\n",
							 x_term_degree_table[i],
							 y_term_degree_table[i],
							 z_term_degree_table[i],
							 q0_tmp[i],
							 x_term_degree_table[j],
							 y_term_degree_table[j],
							 z_term_degree_table[j],
							 q0_tmp[j],
							 x_term_degree_table[j],
							 y_term_degree_table[j],
							 z_term_degree_table[j],
							 res_tmp[j]);
#endif			
			}
		}
	}
	memset(q0_tmp, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	memcpy(q0_tmp, res_tmp, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

#if (1 == TEST_MODE)
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != q0_tmp[i])
		{
			DEBUG_NOTICE("q0_tmp_div: %x | %ld %ld %ld | %x\n",
						 power_polynomial_table[locator_idx][0],
						 x_term_degree_table[i],
						 y_term_degree_table[i],
						 z_term_degree_table[i],
						 q0_tmp[i]);
		}
	}
#endif

	memset(sys_act_poly[poly_idx], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	memcpy(sys_act_poly[poly_idx], q0_tmp, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

#if (1 == TEST_MODE)
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != sys_act_poly[poly_idx][i])
		{
			DEBUG_NOTICE("sys_act_poly_div: %x | %d | %ld %ld %ld | %x\n",
						 power_polynomial_table[locator_idx][0],
						 poly_idx,
						 x_term_degree_table[i],
						 y_term_degree_table[i],
						 z_term_degree_table[i],
						 sys_act_poly[poly_idx][i]);
		}
	}
#endif

	return 0;
}

int poly_aj_factor_trans(long long locator_idx)
{
	long long i = 0, j = 0, k = 0, l = 0;
	unsigned char aj_tmp_flag[CFG_SYS_GEN_POLY];
	memset(aj_tmp_flag, 0, sizeof(unsigned char) * CFG_SYS_GEN_POLY);
	unsigned char aj_tmp_b_flag = 0;
	unsigned char tmp_poly[MAX_POLY_TERM_SIZE];
	memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	unsigned char tmp_val[CFG_SYS_GEN_POLY];
	memset(tmp_val, 0xFF, sizeof(unsigned char) * CFG_SYS_GEN_POLY);
	long long l_star_idx = 0;
	long long deg_star = 65536;
	long long deg_tmp[CFG_SYS_GEN_POLY];
	memset(deg_tmp, 0, sizeof(long long) * CFG_SYS_GEN_POLY);

	memset(aj_factor_a_set_flag, 0, sizeof(unsigned char) * CFG_SYS_GEN_POLY);

	for(i = 0; i < CFG_SYS_GEN_POLY; i++)
	{
		memset(aj_tmp_flag, 0, sizeof(unsigned char) * CFG_SYS_GEN_POLY);
		aj_tmp_b_flag = 0;
		l_star_idx = 0;
	    deg_star = 65536;
	    memset(deg_tmp, 0, sizeof(long long) * CFG_SYS_GEN_POLY);
	    memset(tmp_val, 0xFF, sizeof(unsigned char) * CFG_SYS_GEN_POLY);

		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			if(1 == aj_factor_a_set_flag[j])
			{
				/*this aj has been trans.*/
				continue;
			}

			memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if((0xFF != sys_act_poly[j][k])
					&& ((CFG_SYS_GEN_POLY - 1 - i) == y_term_degree_table[k])
					&& (0 == z_term_degree_table[k]))
				{
					l = term_search(x_term_degree_table[k], 0, 0);
					tmp_poly[l] = sys_act_poly[j][k];
				}
			}

			tmp_val[j] = poly_eva_x_y(tmp_poly, power_polynomial_table[locator_idx][0], 0xFF);
			if(0xFF != tmp_val[j])
			{
				aj_tmp_flag[j] = 1;
				aj_tmp_b_flag = 1;
			}
			//deg_tmp[j] = poly_degree_cal(tmp_poly);
			deg_tmp[j] = poly_degree_cal(sys_act_poly[j]);
			if((deg_star > deg_tmp[j])
				&& (0 < deg_tmp[j])
				&& (0xFF != tmp_val[j]))
			{
				DEBUG_NOTICE("l_star: %ld %ld | %ld | %d\n", deg_star, deg_tmp[j], j, aj_tmp_flag[j]);
				deg_star = deg_tmp[j];
				l_star_idx = j;
			}
		}

		//memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);	

#if (1 == TEST_MODE)
		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			DEBUG_NOTICE("b_trans: %ld | %ld | %ld | %d | %x\n",
						  j,
						  l_star_idx,
						  aj_tmp_b_flag,
						  aj_tmp_flag[j],
						  tmp_val[j]);
		}
#endif

		if(0 != aj_tmp_b_flag)
		{
			aj_factor_a_set_flag[l_star_idx] = 1;

			for(j = 0; j < CFG_SYS_GEN_POLY; j++)
			{
				if(0 == aj_tmp_flag[j])
				{
					continue;
				}
				if(j == l_star_idx)
				{
					continue;
				}
				
				for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
				{
#if 0				
					sys_intp_poly[j][k] = gf_add(gf_multp(tmp_val[l_star_idx], sys_intp_poly[j][k]),
					                             gf_multp(tmp_val[j], sys_intp_poly[l_star_idx][k]));
#endif
					sys_act_poly[j][k] = gf_add(gf_multp(tmp_val[l_star_idx], sys_act_poly[j][k]),
					                     gf_multp(tmp_val[j], sys_act_poly[l_star_idx][k]));
				}
				//aj_factor_flag[locator_idx][j] = 1;
			}
		}

#if (1 == TEST_MODE)
		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			deg_tmp[j] = poly_degree_cal(sys_act_poly[j]);
		
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != sys_act_poly[j][k])
				{
					DEBUG_NOTICE("sys_intp_poly_trans: %d | %ld | %ld | %ld %ld %ld | %x\n",
								 locator_idx,
								 j,
								 deg_tmp[j],
					             x_term_degree_table[k],
						         y_term_degree_table[k],
						         z_term_degree_table[k],
						         sys_act_poly[j][k]);
				}
			}
		}
#endif
	}

	for(i = 0; i < CFG_SYS_GEN_POLY; i++)
	{
		aj_factor_flag[locator_idx][i] = 1;
		aj_tmp_b_flag = 0;

		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if((0xFF != sys_act_poly[i][k])
					&& (j == y_term_degree_table[k])
					&& (0 == z_term_degree_table[k]))
				{
					l = term_search(x_term_degree_table[k], 0, 0);
					tmp_poly[l] = sys_act_poly[i][k];
				}
			}

			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != tmp_poly[k])
				{
					DEBUG_NOTICE("tmp_poly: %ld %ld %ld | %x\n",
					             x_term_degree_table[k],
						         y_term_degree_table[k],
						         z_term_degree_table[k],
						         tmp_poly[k]);
				}
			}

			tmp_val[j] = poly_eva_x_y(tmp_poly, power_polynomial_table[locator_idx][0], 0xFF);
			DEBUG_NOTICE("tmp_val: %x %x\n", power_polynomial_table[locator_idx][0], tmp_val[j]);
			
			if(0xFF != tmp_val[j])
			{
				aj_tmp_b_flag = 1;
			}
			
			if(1 == aj_tmp_b_flag)
			{
				break;
			}
		}
		
		if(1 == aj_tmp_b_flag)
		{
			aj_factor_flag[locator_idx][i] = 0;
		}
		
		deg_tmp[i] = poly_degree_cal(sys_act_poly[i]);
		
		DEBUG_NOTICE("sys_intp_poly_trans_result: %ld %ld | %ld | %ld\n",
					 locator_idx,
					 i,
					 deg_tmp[i],
					 aj_factor_flag[locator_idx][i]);
	}

	return 0;
}

int sys_elm_intp_poly_proc(unsigned char poly_seq[][MAX_POLY_TERM_SIZE], long long pnt_idx)
{
	long long i = 0, j = 0;
	long long sel_poly_idx = 0;
	unsigned char eva_val = 0xFF, min_eva_val = 0xFF;;
	long long poly_deg[CFG_SYS_GEN_POLY], min_deg = 65536;
	memset(poly_deg, 0, sizeof(long long) * CFG_SYS_GEN_POLY);

	min_deg = 65536;
	memset(poly_deg, 0, sizeof(long long) * CFG_SYS_GEN_POLY);
	for(j = 0; j < CFG_SYS_GEN_POLY; j++)
	{
		poly_deg[j] = poly_degree_cal(poly_seq[j]);
		if(poly_deg[j] <= min_deg)
		{
			eva_val = sys_hasse_dev_cal(poly_seq[j], af_pnt[pnt_idx]);
			if(0xFF != eva_val)
			{
				min_deg = poly_deg[j];
				sel_poly_idx = j;
				min_eva_val = eva_val;
			}
		}
	}

	memcpy(sys_elm_poly[pnt_idx], poly_seq[sel_poly_idx], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
	{
		if(0xFF != sys_elm_poly[pnt_idx][j])
		{
			sys_elm_poly[pnt_idx][j] = gf_div(sys_elm_poly[pnt_idx][j], min_eva_val);
			DEBUG_NOTICE("sys_elm_poly: %d | %ld | %x | %ld %ld %ld | %x\n",
						 pnt_idx,
						 sel_poly_idx,
						 min_eva_val,
			             x_term_degree_table[j],
				         y_term_degree_table[j],
				         z_term_degree_table[j],
				         sys_elm_poly[pnt_idx][j]);
		}
	}

	eva_val = sys_hasse_dev_cal(sys_elm_poly[pnt_idx], af_pnt[pnt_idx]);
	DEBUG_NOTICE("sys_elm_poly_sel_pos: %ld | %x | %ld\n", pnt_idx, eva_val, min_deg);
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if((0 == sys_proc_col_flag[i])
			|| (pnt_idx == i))
		{
			continue;
		}

		eva_val = sys_hasse_dev_cal(sys_elm_poly[pnt_idx], af_pnt[i]);
		DEBUG_NOTICE("sys_elm_poly_other_pos: %ld | %x\n", i, eva_val);
	}

	return 0;
}

int sys_elm_poly_gen()
{
	long long i = 0, j = 0, k = 0;
	long long pnt_idx = 0, locator_idx = 0;
	long long sel_poly_idx = 0;
	long long val = 0;
	long long deg_tmp[CFG_SYS_GEN_POLY];
	memset(deg_tmp, 0, sizeof(long long) * CFG_SYS_GEN_POLY);

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		memset(sys_elm_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}

	for(i = 0; i < GF_FIELD; i++)
	{
		val = 0;
		for(j = 0; j < GF_Q; j++)
		{
			if(1 == sys_proc_col_flag[i * GF_Q + j])
			{
				val = 1;
				break;
			}
		}
		if(1 == val)
		{
			check_aj_factor(i);
		}
	}

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if(0 == sys_proc_col_flag[i])
		{
			continue;
		}

		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			memcpy(sys_act_poly[j], sys_intp_poly[j], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		}

		for(j = 0; j < GF_FIELD; j++)
		{
			if(power_polynomial_table[j][0] == af_pnt[i][0])
			{
				locator_idx = j;
				break;
			}
		}

		val = 0;
		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			if(1 == aj_factor_flag[locator_idx][j])
			{
				val = 1;
				break;
			}
		}

		//if(0 == val)
		{
			for(j = 0; j < CFG_SYS_GEN_POLY; j++)
			{
				if(0 == aj_factor_flag[locator_idx][j])
				{
					poly_aj_factor_trans(locator_idx);
					break;
				}
			}
		}

#if 0
		val = 65536;
		sel_poly_idx = 0;
		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			if(1 == aj_factor_flag[locator_idx][j])
			{
				deg_tmp[j] = poly_degree_cal(sys_act_poly[j]);
				if((val > deg_tmp[j])
					&& (0 < deg_tmp[j]))
				{
					val = deg_tmp[j];
					sel_poly_idx = j;
				}
			}
		}
#endif		

#if 1
		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			if(1 == aj_factor_flag[locator_idx][j])
			{
				poly_aj_factor_div(sys_act_poly, j, locator_idx);
			}
		}
#else
		poly_aj_factor_div(sys_act_poly, sel_poly_idx, locator_idx);
#endif

		for(k = 0; k < CODEWORD_LEN; k++)
		{
			if((k == i)
				|| (af_pnt[k][0] != af_pnt[i][0])
				|| (0 == sys_proc_col_flag[k]))
			{
				continue;
			}

			sys_kot_intp(sys_act_poly, af_pnt[k]);
		}

#if (1 == TEST_MODE)
		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != sys_act_poly[j][k])
				{
					DEBUG_NOTICE("sys_intp_poly_final_get: %d %d | %ld %ld %ld | %x\n",
								 i,
								 j,
					             x_term_degree_table[k],
						         y_term_degree_table[k],
						         z_term_degree_table[k],
						         sys_act_poly[j][k]);
				}
			}
		}
#endif

		sys_elm_intp_poly_proc(sys_act_poly, i);
	}

	return 0;
}

int sys_ret_poly_gen()
{
	long long i = 0, j = 0;
	unsigned char eva_val = 0xFF;
	int val = 0;
	long long cwd_dim = CODEWORD_LEN;

	memset(sys_ret_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if(0 == sys_proc_col_flag[i])
		{
			continue;
		}

		for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
		{
			if(0xFF != sys_elm_poly[i][j])
			{
				sys_ret_poly[j] = gf_add(sys_ret_poly[j], 
										 gf_multp(sys_elm_poly[i][j], recv_poly[i]));
			}
		}
	}

#if 1//(1 == TEST_MODE)
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != sys_ret_poly[i])
		{
			DEBUG_NOTICE("sys_ret_poly: %ld %ld %ld | %x\n",
						 x_term_degree_table[i],
						 y_term_degree_table[i],
						 z_term_degree_table[i],
						 sys_ret_poly[i]);
		}
	}

	memset(sys_ret_cwd, 0xFF ,sizeof(unsigned char) * CODEWORD_LEN);
	for(i = 0; i < CODEWORD_LEN; i++)
	{
#if 0	
		if(0 == sys_proc_col_flag[i])
		{
			continue;
		}
#endif

		eva_val = sys_hasse_dev_cal(sys_ret_poly, af_pnt[i]);
		DEBUG_NOTICE("sys_ret_pos_check: %ld | %x %x %x\n", i, eva_val, cwd_poly[i], recv_poly[i]);
		sys_ret_cwd[i] = eva_val;
		
		if((eva_val != recv_poly[i])
			&& (1 == sys_proc_col_flag[i]))
		{
			DEBUG_SYS("sys_ret_cwd_err: %ld | %x %x\n", i, eva_val, recv_poly[i]);
			val = 1;
		}
	}
	
	cwd_dim = poly_degree_cal(sys_ret_poly) + 1 - (GF_Q) * (GF_Q - 1) / 2;
	//cwd_dim = check_cwd_dimension(sys_ret_cwd);
	
	if((0 == val)
		&& (MESSAGE_LEN >= cwd_dim))
	{
		DEBUG_NOTICE("sys_ret_ok: %d\n", cwd_dim);
	}
	else
	{
		DEBUG_SYS("sys_ret_err: %d\n", cwd_dim);
	}
#endif

	return 0;
}

int sys_kot_intp(unsigned char poly_seq[][MAX_POLY_TERM_SIZE], unsigned char *input_point)
{
	long long j = 0, k = 0;
	unsigned char dev_val = 0xFF;
	long long min_idx = -1, min_deg = 65536, deg_tmp = 0;
	unsigned char dev_store[CFG_SYS_GEN_POLY];

	for(j = 0; j < CFG_SYS_GEN_POLY; j++)
	{
		dev_val = sys_hasse_dev_cal(poly_seq[j], input_point);
		dev_store[j] = dev_val;
		DEBUG_NOTICE("dev_val: %ld | %x %x | %x\n",
		             j,
		             input_point[0],
		             input_point[1],
		             dev_val);

		if(0xFF != dev_val)
		{
			deg_tmp = poly_degree_cal(poly_seq[j]);
			if(min_deg > deg_tmp)
			{
				min_deg = deg_tmp;
				min_idx = j;
			}
		}
	}

	DEBUG_NOTICE("min_idx: %x %x | %ld\n", input_point[0], input_point[1], min_idx);
	
	/*update poly.*/
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
				                poly_seq[j],
				                poly_seq[min_idx],
				                dev_store[min_idx],
				                dev_store[j],
				                input_point[0],
				                1);
			}
			else
			{
				sys_poly_update(sys_poly_tmp[j],
				                poly_seq[j],
				                poly_seq[min_idx],
				                dev_store[min_idx],
				                dev_store[j],
				                input_point[0],
				                0);
			}
		}

		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			if(0xFF == dev_store[j])
			{
				continue;
			}

			memcpy(poly_seq[j], sys_poly_tmp[j], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);			
		
#if (1 == TEST_MODE)
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != poly_seq[j][k])
				{
					DEBUG_NOTICE("sys_act_poly_add: %d | %ld %ld %ld | %x\n",
								 j,
					             x_term_degree_table[k],
						         y_term_degree_table[k],
						         z_term_degree_table[k],
						         poly_seq[j][k]);
				}
			}
#endif
		}
	}

	return 0;
}

int sys_gen_test()
{
	long long i = 0, j = 0, k = 0;
	long long point_idx = 0, locator_idx = 0;
	unsigned char dev_val = 0xFF;
	long long min_idx = -1, min_deg = 65536, deg_tmp = 0;
	unsigned char dev_store[CFG_SYS_GEN_POLY];
	memset(dev_store, 0xFF, sizeof(unsigned char) * CFG_SYS_GEN_POLY);
	long long poly_deg[CFG_SYS_GEN_POLY];
	memset(poly_deg, 0, sizeof(long long) * CFG_SYS_GEN_POLY);
	int val = 0;

	sys_mat_init();

	for(i = 0; i < CFG_SYS_GEN_POLY; i++)
	{
		memset(sys_poly_tmp[i], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}

	sys_intp_point_gen(chnl_rel_order_idx);
	sys_intp_poly_init();

	point_idx = 0;

	for(i = 0; i < MESSAGE_LEN; i++)
	{
#if 0	
		min_idx = -1, min_deg = 65536, deg_tmp = 0;
		
		/*cal. hasse dev.*/
		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			dev_val = sys_hasse_dev_cal(sys_intp_poly[j], sys_intp_point[point_idx]);
			dev_store[j] = dev_val;
			DEBUG_NOTICE("dev_val: %ld %ld | %x %x | %x\n",
			             i,
			             j,
			             sys_intp_point[point_idx][0],
			             sys_intp_point[point_idx][1],
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
		
		DEBUG_NOTICE("min_idx: %ld %ld | %ld\n", i, point_idx, min_idx);
		
		/*update poly.*/
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
					                sys_intp_point[point_idx][0],
					                1);
				}
				else
				{
					sys_poly_update(sys_poly_tmp[j],
					                sys_intp_poly[j],
					                sys_intp_poly[min_idx],
					                dev_store[min_idx],
					                dev_store[j],
					                sys_intp_point[point_idx][0],
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
#else
		sys_kot_intp(sys_intp_poly, sys_intp_point[point_idx]);
#endif

		/*set intp. point flag*/
		sys_proc_col_flag[sys_intp_seq[point_idx]] = 1;
		this_new_col_idx = sys_intp_seq[point_idx];

		/*cal. deg. for Inf. set check*/
		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			poly_deg[j] = poly_degree_cal(sys_intp_poly[j]);
		}

#if 1/*check and determine the next intp. point*/
		if((MESSAGE_LEN - (GF_Q) * (GF_Q - 1) / 2) <= i)
		{
			val = sys_check_poly_w(poly_deg);

			if(0 != val)
			{
				val = sys_check_part_dep(i + 1);
			}

			if(0 != val)
			{
				val = sys_check_full_dep(i + 1);
			}

			if(0 != val)
			{
				for(j = 0; j < CFG_SYS_GEN_POLY; j++)
				{
					memcpy(sys_intp_poly[j], sys_prev_poly[j], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
				}
				sys_proc_col_flag[sys_intp_seq[point_idx]] = 0;
				i = i - 1;
				
				/**notice this to avoid lock loop*/
				point_idx++;
				if(K_PLUS_G <= point_idx)
				{
					break;
				}

				continue;
			}
			else
			{
				for(j = 0; j < CFG_SYS_GEN_POLY; j++)
				{
					memcpy(sys_prev_poly[j], sys_intp_poly[j], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
				}
			}
		}
		else
		{
			for(j = 0; j < CFG_SYS_GEN_POLY; j++)
			{
				memcpy(sys_prev_poly[j], sys_intp_poly[j], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			}
		}
#endif

		/*avoid overflow*/
		point_idx++;
		if(K_PLUS_G <= point_idx)
		{
			break;
		}

#if (1 == TEST_MODE)
		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			//poly_deg[j] = poly_degree_cal(sys_intp_poly[j]);
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

		//sys_has_dev_test();
#endif
	}

	min_idx = -1;
	min_deg = 65536;
	for(i = 0; i < CFG_SYS_GEN_POLY; i++)
	{
		poly_deg[i] = poly_degree_cal(sys_intp_poly[i]);
		if(poly_deg[i] < min_deg)
		{
			min_deg = poly_deg[i];
			min_idx = i;
		}
	}
	memcpy(sys_intp_min_poly, sys_intp_poly[min_idx], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

	/*check if dev. are zero*/
#if (1 == TEST_MODE)	
	sys_has_dev_test();
#endif

	/*check if there is aj_factor*/
#if 0	
	for(i = 0; i < GF_FIELD; i++)
	{
		check_aj_factor(i);
	}
#endif	
	
	/*process the poly. according to aj_factor*/
#if 0	
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		for(j = 0; j < GF_FIELD; j++)
		{
			if(power_polynomial_table[j][0] == af_pnt[i][0])
			{
				locator_idx = j;
			}
		}

		if((1 == aj_factor_flag[locator_idx][0])
			&& (1 == sys_proc_col_flag[i]))
		{
			poly_aj_factor_div(0, locator_idx);
		}
		
		if((0 == aj_factor_flag[locator_idx][0])
			&& (1 == sys_proc_col_flag[i]))
		{
			poly_aj_factor_trans(locator_idx);
		}
	}
#else
	sys_elm_poly_gen();
	
	sys_ret_poly_gen();
#endif

	return 0;
}


/*the poly_array for br dec. should be borrowed.*/
int sys_br_g_poly_gen()
{
	long long i = 0, j = 0, k = 0;
	unsigned char tmp_poly_x[MAX_POLY_TERM_SIZE];
	int val = 0;

	memset(br_g_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	br_g_poly[0] = 0x0;

	for(i = 0; i < GF_FIELD; i++)
	{
#if 0	
		val = 0;
		for(j = 0; j < GF_Q; j++)
		{
			if(0 == sys_proc_col_flag[i * GF_Q + j])
			{
				val = 1;
				break;
			}
		}
		if(0 == val)
		{
			continue;
		}
#endif		
	
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
			DEBUG_NOTICE("sys_br_g_poly: %ld %ld %ld | %x\n",
						  x_term_degree_table[i],
						  y_term_degree_table[i],
						  z_term_degree_table[i],
						  br_g_poly[i]);
		}
	}

	return 0;
}

int sys_br_k_poly_construct(long long tv_idx)
{
	long long i = 0, j = 0, k = 0, l = 0;
	unsigned char tmp_poly[MAX_POLY_TERM_SIZE];
	memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	long long intp_idx = 0;

	//for(i = 0; i < tst_vct_num; i++)
	{
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			intp_idx = chnl_rel_order_idx[CODEWORD_LEN - 1 - j];

			if(0xFF == gf_add(tst_vct[tv_idx][intp_idx], sys_ret_cwd[intp_idx]))
			{
				continue;
			}
#if 0
			if(1 == sys_proc_col_flag[intp_idx])
			{
				continue;
			}
#endif

			memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != br_lag_poly[intp_idx][k])
				{
					if((0xFF == gf_add(tst_vct[tv_idx][intp_idx], sys_ret_cwd[intp_idx]))
						|| (0xFF == br_lag_poly[intp_idx][k]))
					{
						tmp_poly[k] = 0xFF;
					}
					else if(0x0 == gf_add(tst_vct[tv_idx][intp_idx], sys_ret_cwd[intp_idx]))
					{
						tmp_poly[k] = br_lag_poly[intp_idx][k];
					}
					else if(0x0 == br_lag_poly[intp_idx][k])
					{
						tmp_poly[k] = gf_add(tst_vct[tv_idx][intp_idx], sys_ret_cwd[intp_idx]);
					}
					else
					{
						tmp_poly[k] = gf_multp(gf_add(tst_vct[tv_idx][intp_idx], sys_ret_cwd[intp_idx]), br_lag_poly[intp_idx][k]);
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
					DEBUG_NOTICE("sys_br_k_poly_cal: %ld %ld | %x %x | %ld %ld %ld | %x\n",
								  tv_idx,
								  intp_idx,
								  gf_add(tst_vct[tv_idx][intp_idx], sys_ret_cwd[intp_idx]),
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

#if (1 == TEST_MODE)
		for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
		{
			if(0xFF != br_k_poly[tv_idx][k])
			{
				DEBUG_NOTICE("sys_br_k_poly: %d | %ld %ld %ld | %x\n",
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
				DEBUG_NOTICE("sys_br_zk_poly: %d | %ld %ld %ld | %x\n",
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

int sys_br_m_poly_construct(long long tv_idx)
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
					DEBUG_NOTICE("sys_br_m_poly: %d %d | %ld %ld %ld | %x\n",
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

int sys_br_v_matric_gen(long long tv_idx)
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
					DEBUG_NOTICE("sys_br_v_matrix_poly: %ld | %ld %ld | %ld %ld %ld | %x\n",
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
					DEBUG_NOTICE("sys_br_v_matrix_poly_up: %ld | %ld %ld | %ld | %ld %ld %ld | %x\n",
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

int sys_br_test()
{
	long long i = 0, j = 0, k = 0;
	int popov_flag = 1;
	long long ms_to_cnt = 0;
	unsigned char test_poly[MAX_POLY_TERM_SIZE];
	memset(test_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

	sys_br_g_poly_gen();
	sys_br_k_poly_construct(0);
	sys_br_m_poly_construct(0);
	sys_br_v_matric_gen(0);
	
	popov_flag = is_popov_form(i);
	while(0 == popov_flag)
	{
		ms_reduction(0);
		popov_flag = is_popov_form(0);
		ms_to_cnt++;
		if(CODEWORD_LEN < ms_to_cnt)
		{
			DEBUG_SYS("br_ms_to: %ld\n", ms_to_cnt);
			break;
		}
	}

	br_g_poly_ret(0);
	br_q_poly_ret(0);

	her_fac(min_intp_poly, tv_est_msg[0], tv_est_cwd[0]);

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		tv_est_cwd[0][i] = gf_add(tv_est_cwd[0][i], sys_ret_cwd[i]);
		DEBUG_NOTICE("br_sys_tv_est_cwd: %d | %ld | %x %x %x\n", tv_err[0], i, tv_est_cwd[0][i], cwd_poly[i], recv_poly[i]);
	}

	return 0;
}

int sys_poly_q0q1_get(unsigned char *poly_seq)
{
	long long i = 0, j = 0;
	
	memset(q0_poly_coef, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	memset(q1_poly_coef, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		//if(0 != z_term_degree_table[i])
		if(1 == z_term_degree_table[i])
		{
			for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
			{
				if((x_term_degree_table[j] == x_term_degree_table[i])
					&& (y_term_degree_table[j] == y_term_degree_table[i])
					&& (0 == z_term_degree_table[j]))
				{
					break;
				}
			}
			q1_poly_coef[j] = poly_seq[i];
		}
		else
		{
			q0_poly_coef[i] = poly_seq[i];
		}
	}
#if (1 == TEST_MODE)	
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != q0_poly_coef[i])
		{
			DEBUG_NOTICE("q0_poly: %ld %ld %ld | %x\n",
			             x_term_degree_table[i],
				         y_term_degree_table[i],
				         z_term_degree_table[i],
				         q0_poly_coef[i]);
		}
	}
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != q1_poly_coef[i])
		{
			DEBUG_NOTICE("q1_poly: %ld %ld %ld | %x\n",
			             x_term_degree_table[i],
				         y_term_degree_table[i],
				         z_term_degree_table[i],
				         q1_poly_coef[i]);
		}
	}
#endif
	return 0;
}

unsigned char sys_hasse_full_dev_cal(unsigned char *poly_seq, unsigned char *input_point, unsigned char test_sym)
{
	unsigned char val = 0xFF, q0_val = 0xFF, q1_val = 0xFF;

	sys_poly_q0q1_get(poly_seq);
	q0_val = poly_eva_x_y(q0_poly_coef, input_point[0], input_point[1]);
	q1_val = poly_eva_x_y(q1_poly_coef, input_point[0], input_point[1]);
	q1_val = gf_multp(q1_val, test_sym);
	val = gf_add(q0_val, q1_val);

	return val;
}

int sys_kot_poly_init()
{
	long long i = 0, j = 0, k = 0;

	poly_init();
#if 0	
	memcpy(intp_poly_coef[0], sys_intp_poly[0], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	memcpy(intp_poly_coef[1], sys_intp_poly[1], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
#endif

	return 0;
}

int sys_poly_break(unsigned char *poly)
{
	long long i = 0, j = 0;
	long long y_degree = 0, z_degree = 0, break_idx = 0;

	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		memset(sys_break_poly[i], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}
	
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != poly[i])
		{
			y_degree = y_term_degree_table[i];
			z_degree = z_term_degree_table[i];
			
			break_idx = z_degree * GF_Q + y_degree;
			
			j = term_search(x_term_degree_table[i], 0, 0);
			sys_break_poly[break_idx][j] = gf_add(poly[i], sys_break_poly[break_idx][j]);
		}
	}

#if (1 == TEST_MODE)
	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
		{
			if(0xFF != sys_break_poly[i][j])
			{
				DEBUG_NOTICE("sys_break_poly: %d | %ld %ld %ld | %x\n",
							 i,
				             x_term_degree_table[j],
					         y_term_degree_table[j],
					         z_term_degree_table[j],
					         sys_break_poly[i][j]);
			}
		}
	}
#endif

	return 0;
}

int sys_br_v_matric_break()
{
	long long i = 0, j = 0, k = 0, l = 0;

	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		for(j = 0; j < KOT_INTP_POLY_NUM; j++)
		{
			memset(sys_break_mat_poly[i][j], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		}
	}

	for(i = 0; i < (KOT_INTP_POLY_NUM / 2); i++)//for row
	{
		for(j = 0; j < KOT_INTP_POLY_NUM; j++)//for col
		{
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if((0xFF != sys_intp_poly[i][k])
					&& ((j % GF_Q) == y_term_degree_table[k])
					&& ((j / GF_Q) == z_term_degree_table[k]))
				{
					l = term_search(x_term_degree_table[k], 0, 0);
					sys_break_mat_poly[i][j][l] = gf_add(sys_intp_poly[i][k], sys_break_mat_poly[i][j][l]);
				}
			}
		}
	}
	
	for(i = (KOT_INTP_POLY_NUM / 2); i < KOT_INTP_POLY_NUM; i++)
	{
		l = term_search(0, 0, 0);
		sys_break_mat_poly[i][i][l] = 0x0;
	}

#if (1 == TEST_MODE)
	for(i = 0; i < KOT_INTP_POLY_NUM; i++)//for row
	{
		for(j = 0; j < KOT_INTP_POLY_NUM; j++)//for col
		{
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != sys_break_mat_poly[i][j][k])
				{
					DEBUG_NOTICE("sys_break_mat_poly: %ld %ld | %ld %ld %ld | %x\n",
					             i,
					             j,
					             x_term_degree_table[k],
								 y_term_degree_table[k],
								 z_term_degree_table[k],
								 sys_break_mat_poly[i][j][k]);
				}
			}
		}
	}
#endif

	return 0;
}

int sys_poly_mat_mul()
{
	long long i = 0, j = 0, k = 0, l = 0, s = 0;

	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		memset(sys_mul_poly[i], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}
	
	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		for(j = 0; j < KOT_INTP_POLY_NUM; j++)
		{
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != sys_break_poly[i][k])
				{
					for(l = 0; l < MAX_POLY_TERM_SIZE; l++)
					{
						if(0xFF != sys_break_mat_poly[i][j][l])
						{
							s = term_search(x_term_degree_table[k] + x_term_degree_table[l], 0, 0);
							sys_mul_poly[j][s] = gf_add(sys_mul_poly[j][s],
														gf_multp(sys_break_poly[i][k], sys_break_mat_poly[i][j][l]));
						
							DEBUG_NOTICE("sys_mul_poly_cal: %ld %ld || %ld | %x || %ld | %x || %ld | %x\n",
										 i,
										 j,
										 x_term_degree_table[k],
										 sys_break_poly[i][k],
										 x_term_degree_table[l],
										 sys_break_mat_poly[i][j][l],
										 x_term_degree_table[s],
										 sys_mul_poly[j][s]);
						}
					}
				}
			}
		}
	}

#if (1 == TEST_MODE)
	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
		{
			if(0xFF != sys_mul_poly[i][k])
			{
				DEBUG_NOTICE("sys_mul_poly: %ld | %ld %ld %ld | %x\n",
				             i,
				             x_term_degree_table[k],
							 y_term_degree_table[k],
							 z_term_degree_table[k],
							 sys_mul_poly[i][k]);
			}
		}
	}
#endif

	return 0;
}

int sys_mul_poly_recover()
{
	long long i = 0, j = 0, k = 0;
	long long y_deg = 0, z_deg = 0;

	memset(min_intp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		y_deg = i % GF_Q;
		z_deg = i / GF_Q;
		for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
		{
			if(0xFF != sys_mul_poly[i][j])
			{
				k = term_search(x_term_degree_table[j], y_deg, z_deg);
				min_intp_poly[k] = gf_add(min_intp_poly[k], sys_mul_poly[i][j]);
			}
		}
	}

#if (1 == TEST_MODE)
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != min_intp_poly[i])
		{
			DEBUG_NOTICE("sys_min_intp_poly: %ld %ld %ld | %x\n",
			             x_term_degree_table[i],
						 y_term_degree_table[i],
						 z_term_degree_table[i],
						 min_intp_poly[i]);
		}
	}
#endif

	return 0;
}

long long sys_min_mul_poly_idx_determine()
{
	long long i = 0, j = 0, k = 0;
	long long min_idx = -1, min_deg = 65536, deg_tmp = 0, deg_this_ele = 0;
	long long deg_mul[KOT_INTP_POLY_NUM][KOT_INTP_POLY_NUM];
	long long deg_mat[KOT_INTP_POLY_NUM][KOT_INTP_POLY_NUM];
	long long deg_this_row[KOT_INTP_POLY_NUM];
	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		for(j = 0; j < KOT_INTP_POLY_NUM; j++)
		{
			deg_mul[i][j] = -65536;
			deg_mat[i][j] = -65536;
		}
	}
	memset(deg_this_row, 0, sizeof(long long) * KOT_INTP_POLY_NUM);
	long long y_deg = 0, z_deg = 0;
	long long w_z = GF_Q * pole_basis_pow[MESSAGE_LEN - 1][0]
				  + (GF_Q + 1) * pole_basis_pow[MESSAGE_LEN - 1][1];
	
	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		for(j = 0; j < KOT_INTP_POLY_NUM; j++)
		{
			y_deg = j % GF_Q;
			z_deg = j / GF_Q;
			
			deg_this_ele = -65536;
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if((0xFF != intp_poly_coef[i][k])
					&& (y_deg == y_term_degree_table[k])
					&& (z_deg == z_term_degree_table[k]))
				{
					deg_tmp = x_term_degree_table[k] * GF_Q
							+ y_deg * (GF_Q + 1)
							+ z_deg * w_z;
					
					if(deg_tmp > deg_this_ele)
					{
						deg_this_ele = deg_tmp;
					}
				}
			}
			
			deg_mul[i][j] = deg_this_ele;
			DEBUG_NOTICE("deg_mul: %ld %ld | %ld\n", i, j, deg_mul[i][j]);
		}
	}
	
	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		for(j = 0; j < KOT_INTP_POLY_NUM; j++)
		{
			deg_this_ele = -65536;
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != sys_break_mat_poly[i][j][k])
				{
					deg_tmp = x_term_degree_table[k] * GF_Q;
					
					if(deg_tmp > deg_this_ele)
					{
						deg_this_ele = deg_tmp;
					}
				}
			}
			
			deg_mat[i][j] = deg_this_ele;
			DEBUG_NOTICE("deg_mat: %ld %ld | %ld\n", i, j, deg_mat[i][j]);
		}
	}
	
	for(i = 0; i < KOT_INTP_POLY_NUM; i++)//row in mul
	{
		deg_this_ele = 0;
		for(j = 0; j < KOT_INTP_POLY_NUM; j++)//col in mat
		{
			for(k = 0; k < KOT_INTP_POLY_NUM; k++)//row in mat and col in mul
			{
				deg_tmp = deg_mul[i][k] + deg_mat[k][j];
				
				if(deg_tmp > deg_this_ele)
				{
					deg_this_ele = deg_tmp;
				}
			}
		}
		deg_this_row[i] = deg_this_ele;
		
		if(deg_this_row[i] <= min_deg)
		{
			min_deg = deg_this_row[i];
			min_idx = i;
		}
	}

	DEBUG_NOTICE("min_idx: %ld | %ld\n", min_deg, min_idx);

	return min_idx;
}

int sys_kot_test()
{
	long long i = 0, j = 0, k = 0;
	unsigned char dev_val = 0xFF;
	long long min_idx = -1, min_deg = 65536, deg_tmp = 0;
	unsigned char dev_store[KOT_INTP_POLY_NUM];
	memset(dev_store, 0xFF, sizeof(unsigned char) * KOT_INTP_POLY_NUM);

	sys_kot_poly_init();

	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		memcpy(intp_poly_tmp[i], intp_poly_coef[i], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if(1 == sys_proc_col_flag[i])
		{
			continue;
		}
		
		for(j = 0; j < KOT_INTP_POLY_NUM; j++)
		{
			dev_val = sys_hasse_full_dev_cal(intp_poly_coef[j], af_pnt[i], gf_add(tst_vct[0][i], sys_ret_cwd[i]));
			dev_store[j] = dev_val;
			DEBUG_NOTICE("dev_val: %ld | %x %x | %x %x | %x | %x | %x\n",
						 j,
						 af_pnt[i][0],
						 af_pnt[i][1],
						 tst_vct[0][i],
						 sys_ret_cwd[i],
						 recv_poly[i],
						 gf_add(tst_vct[0][i], sys_ret_cwd[i]),
						 dev_val);

			if(0xFF != dev_val)
			{
				deg_tmp = poly_degree_cal(intp_poly_coef[j]);
				if(min_deg > deg_tmp)
				{
					min_deg = deg_tmp;
					min_idx = j;
				}
			}
		}
		
		DEBUG_NOTICE("min_idx: %x %x | %ld\n", af_pnt[i][0], af_pnt[i][1], min_idx);
		
		if(-1 != min_idx)
		{
			for(j = 0; j < KOT_INTP_POLY_NUM; j++)
			{
				memset(intp_poly_tmp[j], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

				if(0xFF == dev_store[j])
				{
					continue;
				}

				if(min_idx == j)
				{
					poly_min_update(intp_poly_tmp[j], intp_poly_coef[j], af_pnt[i][0], dev_store[j]);
				}
				else
				{
					poly_normal_update(intp_poly_tmp[j], intp_poly_coef[j], intp_poly_coef[min_idx], dev_store[min_idx], dev_store[j]);
				}
			}

			for(j = 0; j < KOT_INTP_POLY_NUM; j++)
			{			
				if(0xFF == dev_store[j])
				{
					continue;
				}			

				memcpy(intp_poly_coef[j], intp_poly_tmp[j], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);			
			
#if (1 == TEST_MODE)
				for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
				{
					if(0xFF != intp_poly_coef[j][k])
					{
						DEBUG_NOTICE("sys_act_poly_full: %d | %ld %ld %ld | %x\n",
									 j,
						             x_term_degree_table[k],
							         y_term_degree_table[k],
							         z_term_degree_table[k],
							         intp_poly_coef[j][k]);
					}
				}
#endif
			}
		}
	}

	min_deg = 65536;
	for(j = 0; j < KOT_INTP_POLY_NUM; j++)
	{
		deg_tmp = poly_degree_cal(intp_poly_coef[j]);
		dev_store[j] = deg_tmp;
		if(min_deg > deg_tmp)
		{
			min_deg = deg_tmp;
			min_idx = j;
		}

#if (1 == TEST_MODE)		
		for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
		{
			if(0xFF != intp_poly_coef[j][k])
			{
				DEBUG_NOTICE("sys_res_poly_full: %d | %ld %ld %ld | %x\n",
							 j,
				             x_term_degree_table[k],
					         y_term_degree_table[k],
					         z_term_degree_table[k],
					         intp_poly_coef[j][k]);
			}
		}
#endif		
	}
	//memcpy(min_intp_poly, intp_poly_coef[min_idx], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

	sys_br_v_matric_break();
	sys_min_mul_poly_idx_determine();
	for(j = 0; j < KOT_INTP_POLY_NUM; j++)
	{
		sys_poly_break(intp_poly_coef[j]);
		sys_poly_mat_mul();
		sys_mul_poly_recover();
		
#if 1
		her_fac(min_intp_poly, tv_est_msg[0], tv_est_cwd[0]);
		
		for(i = 0; i < CODEWORD_LEN; i++)
		{
			tv_est_cwd[0][i] = gf_add(tv_est_cwd[0][i], sys_ret_cwd[i]);
			DEBUG_NOTICE("kot_sys_tv_est_cwd: %d | %d | %ld | %x %x %x\n", j, tv_err[0], i, tv_est_cwd[0][i], cwd_poly[i], recv_poly[i]);
		}
#endif
	}

	return 0;
}
