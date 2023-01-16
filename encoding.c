#include <stdio.h>
#include <string.h>
#include "debug_info.h"
#include "rnd.h"
#include "gf_cal.h"
#include "encoding.h"

unsigned char af_pnt[CODEWORD_LEN][2];
long long pole_basis_pow[MESSAGE_LEN][2];
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
	long long i = 0, j = 0;
	
	for(i = 0; i < MESSAGE_LEN; i++)
	{
		for(j = 0; j < 2; j++)
		{
			pole_basis_pow[i][j] = 0;
		}
	}
	
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
		DEBUG_NOTICE("pole_basis_pow: %ld | %ld %ld\n",
		             i + 1,
		             pole_basis_pow[i + 1][0],
		             pole_basis_pow[i + 1][1]);
	}
	
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

	return 0;
}

int her_encoding(unsigned char *msg)
{
	long long i = 0, j = 0;
	unsigned char tmp_x = 0x0, tmp_y = 0x0, tmp_prod = 0x0;

#if 0	
	for(i = 0; i < MESSAGE_LEN; i++)
	{
		msg_poly[i] = 0xFF;
	}
#endif
#if 0//(1 == TEST_MODE)
	msg[0] = 0xFF;
	msg[1] = 0x0;
	msg[2] = 0xFF;
	msg[3] = 0x2;
	msg[4] = 0xFF;
	for(i = 0; i < MESSAGE_LEN; i++)
	{
		DEBUG_NOTICE("test_msg_poly: %ld | %x\n",
		             i,
		             msg_poly[i]);
	}
#endif	
	
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		cwd_poly[i] = 0xFF;
		
		for(j = 0; j < MESSAGE_LEN; j++)
		{
			tmp_x = gf_pow_cal(af_pnt[i][0], pole_basis_pow[j][0]);
			tmp_y = gf_pow_cal(af_pnt[i][1], pole_basis_pow[j][1]);
			tmp_prod = gf_multp(tmp_x, tmp_y);
			cwd_poly[i] = gf_add(cwd_poly[i], gf_multp(tmp_prod, msg[j]));
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
		             cwd_poly[i]);
	}
	
	return 0;
}
