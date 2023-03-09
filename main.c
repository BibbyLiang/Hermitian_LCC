#include <stdlib.h>
#include <string.h>
#include "cfg_decoding.h"
#include "debug_info.h"
#include "gf_cal.h"
#include "encoding.h"
#include "rnd.h"
#include "mod.h"
#include "time.h"
#include "interpolation.h"
#include "factorization.h"
#include "re_encoding.h"
#include "lcc_decoding.h"

void init_simulation()
{
	srand(time(NULL));
	init_genrand((long)(time(NULL)));

	runtime = 0;
	err_cnt = 0;

	DEBUG_SYS("init_simulation OK\n");

	return;
}

void main()
{
	long long i = 0;
	long long sim_cnt = 10, monitor_cnt = 1;
	int err_msg = 0, err_cwd = 0;

	init_simulation();

	mod_init();
	tst_vct_num_cal();
	tst_vct_init();
	kot_node_init();
	
	affine_points_cal();
	term_degree_table_init();
	pole_basis_cal();
#if (1 == CFG_FAST_RET)
	ret_t_val_cal();
#endif
	recover_lag_poly_init();

	for(i = 0; i < sim_cnt; i++)
	{
		if((0 == (i % monitor_cnt))
			&& (0 != i))
		{
			DEBUG_SYS("sim: %ld / %ld, err: %ld\n", i, sim_cnt, err_cnt);
			DEBUG_SYS("Decoding Complexity: %f %f\n", (float)(add_cnt / i), (float)((mul_cnt + div_cnt) / i));
			DEBUG_SYS("Decoding Measured Time: %f\n", runtime / ((float)i));
		}

		rnd_msg_gen();
		her_encoding(msg_poly, cwd_poly);
		trans_over_chnl();
		chnl_rel_cal(recv_seq, symbol_num);

		cnt_switch = 1;
		start = clock();

		tst_vct_form();

#if (1 == CFG_RET)
		re_encoding_transform();
#endif

		//koetter_interpolation_hermitian();
		her_lcc_dec();

		stop = clock();
		runtime = runtime + (stop - start) / 1000.0000;
		cnt_switch = 0;

#if 0
#if (1 == CFG_FAC_FREE)
		err_cwd = check_result_cwd(cwd_poly, est_cwd_poly);
		err_msg = check_result_msg(msg_poly, est_msg_poly);
#else

#if (0 == CFG_RET)
		err_msg = check_result_msg(msg_poly, est_msg_poly);
#else
		err_cwd = check_result_cwd(cwd_poly, est_cwd_poly);
#endif

#endif

#else/*for GS decoding check*/

#if 0
	err_cwd = check_result_cwd(cwd_poly, est_cwd_poly);
	err_msg = check_result_msg(msg_poly, est_msg_poly);
	if((0 != err_cwd)
		|| (0 != err_msg))
	{
		err_cnt++;
	}
#endif	

#endif
	}
	DEBUG_SYS("sim: %ld, err: %ld\n", sim_cnt, err_cnt);
	DEBUG_SYS("Decoding Complexity: %f %f\n", (float)(add_cnt / sim_cnt), (float)((mul_cnt + div_cnt) / sim_cnt));
	DEBUG_SYS("Decoding Measured Time: %f\n", runtime / ((float)sim_cnt));

	tst_vct_exit();
	mod_exit();
	
	return;
}
