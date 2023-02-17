#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "debug_info.h"
#include "gf_cal.h"

#if (2 == GF_Q)
unsigned char power_polynomial_table[GF_FIELD][2] = 
{
	/*power <---> polynomial*/
	{0xFF, 0x0},
	{0x0, 0x1},
	{0x1, 0x2},
	{0x2, 0x3}
};
#endif

#if (3 == GF_Q)
unsigned char power_polynomial_table[GF_FIELD][2] = 
{
	/*power <---> polynomial*/
	/*These are also the coefficients of generator polynomials.*/
	{0xFF, 0x0},
	{0x0, 0x1},
	{0x1, 0x2},
	{0x2, 0x4},
	{0x3, 0x3},
	{0x4, 0x6},
	{0x5, 0x7},
	{0x6, 0x5}
};
#endif

#if (4 == GF_Q)
unsigned char power_polynomial_table[GF_FIELD][2] = 
{
	/*power <---> polynomial*/
	/*These are also the coefficients of generator polynomials.*/
	{0xFF, 0x0},
	{0x0, 0x1},
	{0x1, 0x2},
	{0x2, 0x4},
	{0x3, 0x8},
	{0x4, 0x3},
	{0x5, 0x6},
	{0x6, 0xc},
	{0x7, 0xb},
	{0x8, 0x5},
	{0x9, 0xa},
	{0xa, 0x7},
	{0xb, 0xe},
	{0xc, 0xf},
	{0xd, 0xd},
	{0xe, 0x9}
};
#endif

#if (6 == GF_Q)
unsigned char power_polynomial_table[GF_FIELD][2] = 
{
	/*power <---> polynomial*/
	/*These are also the coefficients of generator polynomials.*/
	{0xFF, 0},
	{0, 1},
	{1, 2},
	{2, 4},
	{3, 8},
	{4, 16},
	{5, 32},
	{6, 3},
	{7, 6},
	{8, 12},
	{9, 24},
	{10, 48},
	{11, 35},
	{12, 5},
	{13, 10},
	{14, 20},
	{15, 40},
	{16, 19},
	{17, 38},
	{18, 15},
	{19, 30},
	{20, 60},
	{21, 59},
	{22, 53},
	{23, 41},
	{24, 17},
	{25, 34},
	{26, 7},
	{27, 14},
	{28, 28},
	{29, 56},
	{30, 51},
	{31, 37},
	{32, 9},
	{33, 18},
	{34, 36},
	{35, 11},
	{36, 22},
	{37, 44},
	{38, 27},
	{39, 54},
	{40, 47},
	{41, 29},
	{42, 58},
	{43, 55},
	{44, 45},
	{45, 25},
	{46, 50},
	{47, 39},
	{48, 13},
	{49, 26},
	{50, 52},
	{51, 43},
	{52, 21},
	{53, 42},
	{54, 23},
	{55, 46},
	{56, 31},
	{57, 62},
	{58, 63},
	{59, 61},
	{60, 57},
	{61, 49},
	{62, 33}
};
#endif

#if (8 == GF_Q)
unsigned char power_polynomial_table[GF_FIELD][2] = 
{
	/*power <---> polynomial*/
	/*These are also the coefficients of generator polynomials.*/
	{0xFF, 0},
	{0, 1},
	{1, 2},
	{2, 4},
	{3, 8},
	{4, 16},
	{5, 32},
	{6, 64},
	{7, 128},
	{8, 29},
	{9, 58},
	{10, 116},
	{11, 232},
	{12, 205},
	{13, 135},
	{14, 19},
	{15, 38},
	{16, 76},
	{17, 152},
	{18, 45},
	{19, 90},
	{20, 180},
	{21, 117},
	{22, 234},
	{23, 201},
	{24, 143},
	{25, 3},
	{26, 6},
	{27, 12},
	{28, 24},
	{29, 48},
	{30, 96},
	{31, 192},
	{32, 157},
	{33, 39},
	{34, 78},
	{35, 156},
	{36, 37},
	{37, 74},
	{38, 148},
	{39, 53},
	{40, 106},
	{41, 212},
	{42, 181},
	{43, 119},
	{44, 238},
	{45, 193},
	{46, 159},
	{47, 35},
	{48, 70},
	{49, 140},
	{50, 5},
	{51, 10},
	{52, 20},
	{53, 40},
	{54, 80},
	{55, 160},
	{56, 93},
	{57, 186},
	{58, 105},
	{59, 210},
	{60, 185},
	{61, 111},
	{62, 222},
	{63, 161},
	{64, 95},
	{65, 190},
	{66, 97},
	{67, 194},
	{68, 153},
	{69, 47},
	{70, 94},
	{71, 188},
	{72, 101},
	{73, 202},
	{74, 137},
	{75, 15},
	{76, 30},
	{77, 60},
	{78, 120},
	{79, 240},
	{80, 253},
	{81, 231},
	{82, 211},
	{83, 187},
	{84, 107},
	{85, 214},
	{86, 177},
	{87, 127},
	{88, 254},
	{89, 225},
	{90, 223},
	{91, 163},
	{92, 91},
	{93, 182},
	{94, 113},
	{95, 226},
	{96, 217},
	{97, 175},
	{98, 67},
	{99, 134},
	{100, 17},
	{101, 34},
	{102, 68},
	{103, 136},
	{104, 13},
	{105, 26},
	{106, 52},
	{107, 104},
	{108, 208},
	{109, 189},
	{110, 103},
	{111, 206},
	{112, 129},
	{113, 31},
	{114, 62},
	{115, 124},
	{116, 248},
	{117, 237},
	{118, 199},
	{119, 147},
	{120, 59},
	{121, 118},
	{122, 236},
	{123, 197},
	{124, 151},
	{125, 51},
	{126, 102},
	{127, 204},
	{128, 133},
	{129, 23},
	{130, 46},
	{131, 92},
	{132, 184},
	{133, 109},
	{134, 218},
	{135, 169},
	{136, 79},
	{137, 158},
	{138, 33},
	{139, 66},
	{140, 132},
	{141, 21},
	{142, 42},
	{143, 84},
	{144, 168},
	{145, 77},
	{146, 154},
	{147, 41},
	{148, 82},
	{149, 164},
	{150, 85},
	{151, 170},
	{152, 73},
	{153, 146},
	{154, 57},
	{155, 114},
	{156, 228},
	{157, 213},
	{158, 183},
	{159, 115},
	{160, 230},
	{161, 209},
	{162, 191},
	{163, 99},
	{164, 198},
	{165, 145},
	{166, 63},
	{167, 126},
	{168, 252},
	{169, 229},
	{170, 215},
	{171, 179},
	{172, 123},
	{173, 246},
	{174, 241},
	{175, 255},
	{176, 227},
	{177, 219},
	{178, 171},
	{179, 75},
	{180, 150},
	{181, 49},
	{182, 98},
	{183, 196},
	{184, 149},
	{185, 55},
	{186, 110},
	{187, 220},
	{188, 165},
	{189, 87},
	{190, 174},
	{191, 65},
	{192, 130},
	{193, 25},
	{194, 50},
	{195, 100},
	{196, 200},
	{197, 141},
	{198, 7},
	{199, 14},
	{200, 28},
	{201, 56},
	{202, 112},
	{203, 224},
	{204, 221},
	{205, 167},
	{206, 83},
	{207, 166},
	{208, 81},
	{209, 162},
	{210, 89},
	{211, 178},
	{212, 121},
	{213, 242},
	{214, 249},
	{215, 239},
	{216, 195},
	{217, 155},
	{218, 43},
	{219, 86},
	{220, 172},
	{221, 69},
	{222, 138},
	{223, 9},
	{224, 18},
	{225, 36},
	{226, 72},
	{227, 144},
	{228, 61},
	{229, 122},
	{230, 244},
	{231, 245},
	{232, 247},
	{233, 243},
	{234, 251},
	{235, 235},
	{236, 203},
	{237, 139},
	{238, 11},
	{239, 22},
	{240, 44},
	{241, 88},
	{242, 176},
	{243, 125},
	{244, 250},
	{245, 233},
	{246, 207},
	{247, 131},
	{248, 27},
	{249, 54},
	{250, 108},
	{251, 216},
	{252, 173},
	{253, 71},
	{254, 142}
};
#endif

#if (1 == GF_CAL_COUNT)
unsigned char cnt_switch = 0;

long long add_cnt = 0;
long long mul_cnt = 0;
long long div_cnt = 0;
long long real_cbm_cnt = 0;
long long real_mul_ff_cnt = 0;
long long pow_cnt = 0;
long long add_cnt_prev = 0;
long long mul_cnt_prev = 0;
long long div_cnt_prev = 0;
long long real_cbm_cnt_prev = 0;
long long real_mul_ff_cnt_prev = 0;
long long pow_cnt_prev = 0;
#endif


unsigned char gf_pow2poly(unsigned char val_in_pow)
{
	unsigned char val_in_poly = 0;
	if(0xFF == val_in_pow)
	{
		return power_polynomial_table[0][1];
	}
	else
	{
		return power_polynomial_table[val_in_pow + 1][1];
	}
}

unsigned char gf_poly2pow(unsigned char val_in_poly)
{
	unsigned char i = 0;
	unsigned char val_in_pow = 0;

	for(i = 0; i < GF_FIELD; i++)
	{
		if(power_polynomial_table[i][1] == val_in_poly)
		{
			val_in_pow = power_polynomial_table[i][0];
			break;
		}
	}

	return val_in_pow;
}

unsigned char gf_location(unsigned char val)
{
	unsigned char val_location = 0xFF;

	if(0xFF == val)
	{
		val_location = power_polynomial_table[0][0];
	}
	else
	{
		val_location = power_polynomial_table[val + 1][0];
	}

	return val_location;
}

unsigned char gf_add(unsigned char a, unsigned char b)
{
#if (1 == GF_CAL_COUNT)
	if(1 == cnt_switch)
	{
		add_cnt++;
	}	
#endif

	if(0xFF == a)
	{
		return b;
	}
	if(0xFF == b)
	{
		return a;
	}

	unsigned char i = 0;
	unsigned char sum_in_pow = 0;
	
	unsigned char sum_in_poly = gf_pow2poly(a) ^ gf_pow2poly(b);

	sum_in_pow = gf_poly2pow(sum_in_poly);

	return sum_in_pow;
}

unsigned char gf_multp(unsigned char a, unsigned char b)
{
#if (1 == GF_CAL_COUNT)
	if(1 == cnt_switch)
	{
		mul_cnt++;
	}
#endif

	if((0xFF == a) || (0xFF == b))
	{
		return 0xFF;
	}

	if(0x0 == a)
	{
		return b;
	}
	
	if(0x0 == b)
	{
		return a;
	}

	unsigned char product_in_pow = (a + b) % (GF_FIELD - 1);

	return product_in_pow;
}

unsigned char gf_div(unsigned char a, unsigned char b)
{	
	if(0xFF == a)
	{
		return 0xFF;
	}
	if(0xFF == b)
	{
		DEBUG_NOTICE("div err.\n");
		return 0xFF;
	}
#if (1 == GF_CAL_COUNT)
	if(1 == cnt_switch)
	{
		div_cnt++;
	}
#endif
	//DEBUG_NOTICE("div: %x %x\n", a, b);
	unsigned char quotient_in_pow = 0;
	if(a >= b)
	{
		quotient_in_pow = (a - b) % (GF_FIELD - 1);
	}
	else
	{
		quotient_in_pow = ((b / (GF_FIELD - 1) + 1) * (GF_FIELD - 1) + a - b) % (GF_FIELD - 1);
	}

	return quotient_in_pow;
}

unsigned char gf_mod_single_term(unsigned char a, unsigned char b)
{
	unsigned char remainder_in_pow = 0xFF;
	if(a >= b)
	{
		/*remove small term*/
		remainder_in_pow = 0xFF;
	}
	else
	{
		/*keep big term*/
		remainder_in_pow = 0;
	}

	return remainder_in_pow;
}

long long gf_degree(unsigned char* a, long long len_a)
{
	long long i = 0;

	for(i = len_a - 1; i >= 0; i--)
	{
		if(0xFF != a[i])
		{
			break;
		}
	}
	if(0 > i)
	{
		i = 0;
	}

	return i;
}

unsigned char gf_div_q_r(unsigned char* dividend, long long len_dividend,
							unsigned char* divisor, long long len_divisor,
							unsigned char* quotien, long long len_quotien,
							unsigned char* remainder, long long len_remainder)
{
	long long i = 0, j = 0, k = 0;
	long long locator = 0, locator_rmd = 0;
	unsigned char factor = 0, factor_rmd = 0;
	unsigned char dividend_tmp[len_dividend], remainder_tmp[len_remainder];
	memset(dividend_tmp, 0xFF, sizeof(unsigned char) * len_dividend);
	memset(remainder_tmp, 0xFF, sizeof(unsigned char) * len_remainder);

	if(gf_degree(divisor, len_divisor) > gf_degree(dividend, len_dividend))
	{
		for(i = 0; i < len_remainder; i++)
		{
			remainder[i] = dividend[i];
		}
		DEBUG_NOTICE("quotien is zero: %d %d\n", gf_degree(divisor, len_divisor), gf_degree(dividend, len_dividend));

		return 0;
	}

	memcpy(dividend_tmp, dividend, sizeof(unsigned char) * len_dividend);
	for(i = 0; i < len_dividend; i++)
	{
		locator = gf_degree(dividend_tmp, len_dividend) - gf_degree(divisor, len_divisor);
		factor = gf_div(dividend_tmp[gf_degree(dividend_tmp, len_dividend)], divisor[gf_degree(divisor, len_divisor)]);

		quotien[locator] = factor;
		DEBUG_INFO("quotien: %d %x %d %x\n", i, quotien[locator], locator, factor);

		memcpy(remainder_tmp, dividend_tmp, sizeof(unsigned char) * len_remainder);
		for(j = 0; j < len_divisor; j++)
		{
			if(0xFF == factor)
			{
				factor_rmd = 0xFF;
			}
			else if(0x0 == factor)
			{
				factor_rmd = divisor[gf_degree(divisor, len_divisor) - j];
			}
			else if(0xFF == divisor[gf_degree(divisor, len_divisor) - j])
			{
				factor_rmd = 0xFF;
			}
			else if(0x0 == divisor[gf_degree(divisor, len_divisor) - j])
			{
				factor_rmd = factor;
			}
			else
			{
				factor_rmd = gf_multp(factor, divisor[gf_degree(divisor, len_divisor) - j]);
			}
			locator_rmd = locator + gf_degree(divisor, len_divisor) - j;
			remainder_tmp[locator_rmd] = gf_add(dividend_tmp[locator_rmd], factor_rmd);
			DEBUG_NOTICE("remainder_tmp: %x %d %x\n", factor_rmd, locator_rmd, remainder_tmp[locator_rmd]);
		}
#if 1
		DEBUG_NOTICE("remainder_tmp:\n");
		for(k = 0; k < len_remainder; k++)
		{
			DEBUG_NOTICE("%x ", remainder_tmp[k]);
		}
		DEBUG_NOTICE("\n");
#endif
		DEBUG_NOTICE("div_degree: %d %d\n",
		             gf_degree(divisor, len_divisor),
		             gf_degree(remainder_tmp, len_remainder));
		if(gf_degree(divisor, len_divisor) > gf_degree(remainder_tmp, len_remainder))
		{
			for(k = 0; k < len_remainder; k++)
			{
				remainder[k] = remainder_tmp[k];
			}
			break;
		}
		else
		{
			memcpy(dividend_tmp, remainder_tmp, sizeof(unsigned char) * len_dividend);
			memset(remainder_tmp, 0xFF, sizeof(unsigned char) * len_remainder);
		}
	}

	return 0;
}

unsigned char gf_multp_poly(unsigned char* a, long long len_a,
								unsigned char* b, long long len_b,
								unsigned char* product, long long len_product)
{
	unsigned char i = 0, j = 0, idx = 0;
	unsigned char tmp_val = 0xFF;

	for(i = 0; i < len_a; i++)
	{
		for(j = 0; j < len_b; j++)
		{
			idx = i + j;
			if(len_product <= idx)
			{
				//DEBUG_NOTICE("product len err: %d\n", idx);
				continue;
			}
			if((0xFF == a[i])
				|| (0xFF == b[j]))
			{
				tmp_val = 0xFF;
			}
			else if(0x0 == a[i])
			{
				tmp_val = b[j];
			}
			else if(0x0 == b[j])
			{
				tmp_val = a[j];
			}
			else
			{
				tmp_val = gf_multp(a[i], b[j]);
			}
			
			if(0xFF == product[idx])
			{
				product[idx] = tmp_val;
			}
			else if(0xFF == tmp_val)
			{
				product[idx] = product[idx];
			}
			else
			{
				product[idx] = gf_add(product[idx], tmp_val);
			}
			//DEBUG_NOTICE("product: %d %d %d | %x\n", i, j, idx, product[idx]);
		}
	}
}

int gf_multp_poly_hw(unsigned char* a, unsigned char len_a,
				 		   unsigned char* b, unsigned char len_b,
				 		   unsigned char* product, unsigned char len_product)
{
	unsigned char i = 0, j = 0, idx = len_product - 1;
	unsigned char reg[len_a - 1];
	memset(reg, 0xFF, sizeof(unsigned char) * (len_a - 1));
	unsigned char pd_tmp = 0xFF;

	/*high -> low*/
	for(i = 0; i < len_b; i++)
	{
		pd_tmp = gf_multp(*(a + len_a - 1), *(b + (len_b - 1 - i)));
		product[idx] = gf_add(reg[len_a - 2], pd_tmp);
		//DEBUG_NOTICE("unrel_group_seq: %x\n", a[0]);
#if 0
		DEBUG_NOTICE("%x %x %x %x %x %x\n", 
									idx,
									product[idx],
					  	   	 		reg[len_a - 2],
					       	 		pd_tmp,
					       	 		*(a + len_a - 1),
					       	 		*(b + (len_b - 1 - i)));
#endif
		if(0 >= idx)
		{
			//DEBUG_NOTICE("product len err: %d\n", idx);
			break;
		}
		idx = idx - 1;

		for(j = 1; j < (len_a - 1); j++)
		{
			pd_tmp = gf_multp(*(a + len_a - 1 - j), *(b + (len_b - 1 - i)));
			reg[len_a - 1 - j] = gf_add(reg[len_a - 1 - j - 1], pd_tmp);
#if 0
			DEBUG_NOTICE("idx: %d\n", len_a - 1 - j);
			DEBUG_NOTICE("%x %x %x %x %x\n", reg[len_a - 1 - j], reg[len_a - 1 - j - 1], pd_tmp,
								 	   *(a + len_a - 1 - j), *(b + (len_b - 1 - i)));
#endif
		}
		reg[0] = gf_multp(*(a + 0), *(b + (len_b - 1 - i)));
		//DEBUG_NOTICE("%x %x %x\n", reg[0], *(a + 0), *(b + (len_b - 1 - i)));
	}

	for(i = 0; i < len_a; i++)
	{
		pd_tmp = gf_multp(*(a + len_a - 1), 0xFF);
		product[idx] = gf_add(reg[len_a - 2], pd_tmp);
#if 0
		DEBUG_NOTICE("%x %x %x %x %x %x\n", 
									idx,
									product[idx],
					  	   	 		reg[len_a - 2],
					       	 		pd_tmp,
					       	 		*(a + len_a - 1),
					       	 		0xFF);
#endif
		if(0 >= idx)
		{
			//DEBUG_NOTICE("product len err: %d\n", idx);
			break;
		}
		idx = idx - 1;

		for(j = 1; j < (len_a - 1); j++)
		{
			pd_tmp = gf_multp(*(a + len_a - 1 - j - 1), 0xFF);
			reg[len_a - 1 - j] = gf_add(reg[len_a - 1 - j - 1], pd_tmp);
		}
		reg[0] = gf_multp(*(a + 0), 0xFF);
	}
	
	return 0;
}

long long real_combine(long long n, long long k)
{
#if (1 == GF_CAL_COUNT)	
	if(1 == cnt_switch)
	{
		real_cbm_cnt++;
	}
#endif
	long long combine_num = 0;

#if 0//it is useless when values are too large
	int i = 0;
	long tmp_n = 1, tmp_k = 1, tmp_n_k = 1;

	for(i = 1; i < (n + 1); i++)
	{
		tmp_n = tmp_n * i;
	}

	for(i = 1; i < (k + 1); i++)
	{
		tmp_k = tmp_k * i;
	}

	for(i = 1; i < (n - k + 1); i++)
	{
		tmp_n_k = tmp_n_k * i;
	}

	combine_num = tmp_n / tmp_k / tmp_n_k;
#else//fast calculation for finite field
	if(k == (n & k))
	{
		combine_num = 1;
	}
	else
	{
		combine_num = 2;
	}
#endif

	return combine_num;
}

unsigned char gf_real_mutp_ff(long long n, unsigned char ff)
{
#if (1 == GF_CAL_COUNT)
	if(1 == cnt_switch)
	{
		real_mul_ff_cnt++;
	}
#endif
	unsigned char val = 0xFF;

	if(0 != (n % 2))
	{
		val = ff;
	}
	else
	{
		val = 0xFF;
	}

	return val;
}

unsigned char gf_pow_cal(unsigned char ff, long long n)
{
#if (1 == GF_CAL_COUNT)
	if(1 == cnt_switch)
	{
		pow_cnt++;
	}
#endif	
	unsigned char val = 0xFF;
	if(0xFF == ff)
	{
		if(0 != n)
		{
			return 0xFF;
		}
		else
		{
			return 0x0;
		}
	}

	if(0 <= n)
	{
		val = (ff * n) % (GF_FIELD - 1);
	}
	else
	{
		val = (power_polynomial_table[-n + 1][0]) % (GF_FIELD - 1);
		val = (ff * val) % (GF_FIELD - 1);
	}

	return val;
}

unsigned char phase_trans(unsigned char phase)
{
	unsigned char val = 0;

	if(0 == phase)
	{
		val = 1;
	}
	else
	{
		val = 0;
	}

	return val;
}

#if (1 == GF_CAL_COUNT)
void gf_count_reset()
{
	add_cnt = 0;
	mul_cnt = 0;
	div_cnt = 0;
	real_cbm_cnt = 0;
	real_mul_ff_cnt = 0;
	pow_cnt = 0;
	add_cnt_prev = 0;
	mul_cnt_prev = 0;
	div_cnt_prev = 0;
	real_cbm_cnt_prev = 0;
	real_mul_ff_cnt_prev = 0;
	pow_cnt_prev = 0;

	return;
}
#endif

void gf_count_switch(unsigned char count_switch)
{
#if (1 == GF_CAL_COUNT)	
	cnt_switch = count_switch;
#endif
}

/*升序*/
void BubbleSort4(float *A, int len, long long *A_idx)
{
    int low = 0, high = len - 1;
    int i = 0;
    float tmp = 0;

    while(low < high)
    {
        for(i = low; i < high; i++)  // 正向冒泡,找到最大者
        {
            if(A[i] > A[i + 1])
            {
                tmp = A[i];
                A[i] = A[i + 1];
                A[i + 1] = tmp;

                tmp = A_idx[i];
                A_idx[i] = A_idx[i + 1];
                A_idx[i + 1] = tmp;
            }
        }
        high--;            // 修改high值, 前移一位 
        for(i = high; i > low; i--)     // 反向冒泡,找到最小者 
        {
            if(A[i] < A[i - 1])
            {
                tmp = A[i];
                A[i] = A[i - 1];
                A[i - 1] = tmp;

                tmp = A_idx[i];
                A_idx[i] = A_idx[i - 1];
                A_idx[i - 1] = tmp;
            }
        }
        low++;            // 修改low值,后移一位
    }
}

/*降序*/
void BubbleSort5(float *A, int len, long long *A_idx)
{
    int low = 0, high = len - 1;
    int i = 0;
    float tmp = 0;

    while(low < high)
    {
        for(i = low; i < high; i++)  // 正向冒泡,找到最小者
        {
            if(A[i] < A[i + 1])
            {
                tmp = A[i];
                A[i] = A[i + 1];
                A[i + 1] = tmp;

                tmp = A_idx[i];
                A_idx[i] = A_idx[i + 1];
                A_idx[i + 1] = tmp;
            }
        }
        high--;            // 修改high值, 前移一位 
        for(i = high; i > low; i--)     // 反向冒泡,找到最大者 
        {
            if(A[i] > A[i - 1])
            {
                tmp = A[i];
                A[i] = A[i - 1];
                A[i - 1] = tmp;

                tmp = A_idx[i];
                A_idx[i] = A_idx[i - 1];
                A_idx[i - 1] = tmp;
            }
        }
        low++;            // 修改low值,后移一位
    }
}

long long term_search(long long x_degree, long long y_degree, long long z_degree)
{
	long long i = 0;
	long long val = 0;

	val = x_degree * (Y_MAX_SIZE * Z_MAX_SIZE)
	    + y_degree * Z_MAX_SIZE
	    + z_degree;

	return val;
}
