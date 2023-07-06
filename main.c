#include <stdlib.h>
#include <string.h>
#include "cfg_decoding.h"
#include "debug_info.h"
#include "gf_cal.h"
#include "encoding.h"
#include "rnd.h"
#include "channel.h"
#include "mod.h"
#include "time.h"
#include "interpolation.h"
#include "factorization.h"
#include "re_encoding.h"
#include "br.h"
#include "sys_gen.h"
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

void clear_sim()
{
	runtime = 0;
	err_cnt = 0;
	gf_count_reset();
#if (1 == DEV_RECORD)	
	dev_1_avg_cnt = 0;
	dev_2_avg_cnt = 0;
	dev_out_cnt = 0;
#endif	
}

int tv_test_mode()
{
	long long i = 0;
	long long err_val = 1;

	long long genus = GF_Q * (GF_Q - 1) / 2;
	long long radius = (CODEWORD_LEN - MESSAGE_LEN - genus + 1 - 1) / 2;


	for(i = 0; i < tst_vct_num; i++)
	{
		DEBUG_NOTICE("tst_vct_num: %ld\n", tst_vct_num);
		if(radius >= tv_err[i])
		{
			err_val = 0;
			break;
		}
	}
	if(0 != err_val)
	{
		err_cnt++;
	}

	return 0;
}

void main()
{
	long long i = 0;
	long long sim_cnt = 1, monitor_cnt = 1;
	float eb2n0_start = eb2n0, eb2n0_stop = eb2n0, eb2n0_step = 1;
	int val = 0;

	int err_msg = 0, err_cwd = 0;

#if (0 == TEST_MODE)
	printf("Please Input Eb/N0 Start: ");
	scanf("%f", &eb2n0_start);
	printf("Please Input Eb/N0 Stop: ");
	scanf("%f", &eb2n0_stop);
	printf("Please Input Eb/N0 Step: ");
	scanf("%f", &eb2n0_step);
	printf("Please Input Simulation Times: ");
	scanf("%ld", &sim_cnt);
	printf("Please Monitor Times: ");
	scanf("%ld", &monitor_cnt);
#endif

#if (1 == OUTPUT_LOG)
	/*init file log*/
	sprintf(log_name, "n_%d-k_%d-eta_%d-snr_%f_%f_%f-cnt_%ld_%ld.txt",
					  CODEWORD_LEN,
					  MESSAGE_LEN,
					  ETA,
					  eb2n0_start,
					  eb2n0_step,
					  eb2n0_stop,
					  sim_cnt,
					  monitor_cnt);
#endif

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
	
#if (1 == CFG_BR)
	br_poly_init();
	br_g_poly_gen();
	br_lag_poly_construct();
#endif

	for(eb2n0 = eb2n0_start; eb2n0 <= eb2n0_stop; eb2n0 = eb2n0 + eb2n0_step)
	{
		clear_sim();

		for(i = 0; i < sim_cnt; i++)
		{
#if (1 == EARLY_TERMINATION)
			if(ET_NUM <= err_cnt)
			{
				break;
			}
#endif
		
			if((0 == (i % monitor_cnt))
				&& (0 != i))
			{
				DEBUG_SYS("Eb/N0: %f dB\n", eb2n0);
				DEBUG_SYS("sim: %ld / %ld, err: %ld\n", i, sim_cnt, err_cnt);
				DEBUG_SYS("Decoding Complexity: %f %f\n", (float)(add_cnt / i), (float)((mul_cnt + div_cnt) / i));
				DEBUG_SYS("Decoding Measured Time: %f\n", runtime * 1000 / ((float)i));

#if (1 == DEV_RECORD)
				DEBUG_SYS("dev_1_avg_cnt: %f\n", (float)(dev_1_avg_cnt) / i);
				DEBUG_SYS("dev_2_avg_cnt: %f\n", (float)(dev_2_avg_cnt) / i);
				DEBUG_SYS("dev_out_cnt: %ld\n", dev_out_cnt);			
#endif

#if (1 == OUTPUT_LOG)
				frc = fopen(log_name, "a+");
				fprintf(frc, "---------------------\n");
				fprintf(frc, "Eb/N0: %f dB\n", eb2n0);
				fprintf(frc, "sim: %ld / %ld, err: %ld\n", i, sim_cnt, err_cnt);
				fprintf(frc, "Decoding Complexity: %f %f\n", (float)(add_cnt / i), (float)((mul_cnt + div_cnt) / i));
				fprintf(frc, "Decoding Measured Time: %f\n", runtime * 1000 / ((float)i));

#if (1 == DEV_RECORD)
				fprintf(frc, "dev_1_avg_cnt: %f\n", (float)(dev_1_avg_cnt) / i);
				fprintf(frc, "dev_2_avg_cnt: %f\n", (float)(dev_2_avg_cnt) / i);
				fprintf(frc, "dev_out_cnt: %ld\n", dev_out_cnt);
#endif

				fclose(frc);
				frc = NULL;
#endif
			}

			rnd_msg_gen();
			her_encoding(msg_poly, cwd_poly);
			trans_over_chnl();
			chnl_rel_cal(recv_seq, symbol_num);

			cnt_switch = 1;
			start = clock();

			tst_vct_form();
			
			//stop = clock();
			//runtime = runtime + (stop - start) / 1000.0000;

#if (1 == TV_TEST)
			tv_test_mode();

			stop = clock();
			runtime = runtime + (stop - start) / 1000.0000;
			
			continue;
#endif

#if (1 == CFG_RET)
			val = re_encoding_transform();
#endif

#if 0//(1 == CFG_BR)
			br_test();
#endif
#if (1 == CFG_SYS_GEN)
			sys_gen_test();
			//sys_br_test();
			sys_kot_test();
#endif

			//koetter_interpolation_hermitian();
			if(0 == val)
			{
				her_lcc_dec();
			}

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

		DEBUG_SYS("*********************************\n");
		DEBUG_SYS("Eb/N0: %f dB\n", eb2n0);
		DEBUG_SYS("sim: %ld, err: %ld\n", i, err_cnt);
		DEBUG_SYS("Decoding Complexity: %f %f\n", (float)(add_cnt / i), (float)((mul_cnt + div_cnt) / i));
		DEBUG_SYS("Decoding Measured Time: %f\n", runtime * 1000 / ((float)i));
		DEBUG_SYS("*********************************\n");
#if (1 == OUTPUT_LOG)
		frc = fopen(log_name, "a+");
		fprintf(frc, "*********************************\n");
		fprintf(frc, "Eb/N0: %f dB\n", eb2n0);
		fprintf(frc, "sim: %ld, err: %ld\n", i, err_cnt);
		fprintf(frc, "Decoding Complexity: %f %f\n", (float)(add_cnt / i), (float)((mul_cnt + div_cnt) / i));
		fprintf(frc, "Decoding Measured Time: %f\n", runtime * 1000 / ((float)i));
		fprintf(frc, "*********************************\n");
		fclose(frc);
		frc = NULL;
#endif
	}

	tst_vct_exit();
	mod_exit();

#if (1  == CFG_BR)
	br_poly_exit();
#endif

	return;
}
