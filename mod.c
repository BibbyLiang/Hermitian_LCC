#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "debug_info.h"
#include "gf_cal.h"
#include "mod.h"
#include "channel.h"
#include "encoding.h"

float **recv_seq;
float recv_rel[CODEWORD_LEN];
long long symbol_num = CODEWORD_LEN * GF_Q * BITS_PER_SYMBOL_BPSK;

int mod_init()
{
	long long i = 0;
	long long symbol_num = CODEWORD_LEN * GF_Q * BITS_PER_SYMBOL_BPSK;

	recv_seq = (float**)malloc(sizeof(float*) * symbol_num);
	for (i = 0; i < symbol_num; i++)
	{
		recv_seq[i] = (float*)malloc(sizeof(float) * 2);
	}

	DEBUG_SYS("mod_init OK\n");

	return 0;
}

int mod_exit()
{
	long long i = 0;
	long long symbol_num = CODEWORD_LEN * GF_Q * BITS_PER_SYMBOL_BPSK;

	for (i = 0; i < symbol_num; i++)
	{
  		free(recv_seq[i]);
		recv_seq[i] = NULL;
  	}
	free(recv_seq);
	recv_seq = NULL;

	return 0;
}

int bpsk_mod(unsigned char *input_seq, 
				 unsigned int input_len,
				 float **output_seq,
				 unsigned int output_len)
{
	unsigned int i = 0, j = 0, k = 0;
	unsigned char tmp_input = 0;

	for(i = 0; i < input_len; i++)
	{
		for(j = 0; j < GF_Q; j++)
		{
			/*LSB first, small-endian*/
			if(0xFF == input_seq[i])
			{
				tmp_input = (power_polynomial_table[0][1] >> j) & 0x1;
			}
			else
			{
				tmp_input = (power_polynomial_table[input_seq[i] + 0x1][1] >> j) & 0x1;
			}
			/*0 -> 1, 1 -> -1*/
			output_seq[k][0] = 1 - (float)(2 * tmp_input);
			output_seq[k][1] = 0.0;

			k = k + 1;
			if(k > output_len)
			{
				DEBUG_NOTICE("mod_len_err.\n");
			}
		}
	}

	if(k < output_len)
	{
		for(i = k; i < output_len; i++)
		{
			output_seq[i][0] = 0.0;
			output_seq[i][1] = 0.0;
		}
	}

#if (1 == CFG_DEBUG_INFO)
	DEBUG_INFO("bpsk modulation:\n");
	for(i = 0; i < output_len; i++)
	{
		DEBUG_INFO("%f %f\n", output_seq[i][0], output_seq[i][1]);
	}
	DEBUG_INFO("\n");
#endif

	return 0;
}

int bpsk_demod(float **input_seq, 
				    unsigned int input_len,
				    unsigned char *output_seq,
				    unsigned int output_len)
{
	unsigned int i = 0, j = 0, k = 0;
	unsigned char tmp_bit_decision = 0, tmp_output = 0;
	unsigned char tmp_output_seq[output_len];
	float d0 = 0, d1 = 0;

	j = 0, k = 0;
	for(i = 0; i < input_len; i++)
	{
		/*LSB first, small-endian*/
		/*hard-decision, don't think about Q*/
		//tmp_bit_decision = ((unsigned char)((1 - input_seq[i][0]) / 2)) & 0x1;
		/*soft-like-decision*/
		d0 = (input_seq[i][0] - (1.0)) * (input_seq[i][0] - (1.0))
				+ (input_seq[i][1] - (0.0)) * (input_seq[i][1] - (0.0));
		d1 = (input_seq[i][0] - (-1.0)) * (input_seq[i][0] - (-1.0))
				+ (input_seq[i][1] - (0.0)) * (input_seq[i][1] - (0.0));
		if(d0 > d1)
		{
			tmp_bit_decision = 1;
		}
		else
		{
			tmp_bit_decision = 0;
		}
		tmp_output = tmp_output | (tmp_bit_decision << j);

		j = j + 1;
		if(GF_Q == j)
		{
			tmp_output_seq[k] = tmp_output;
			tmp_output = 0;
			j = 0;
			k = k + 1;
		}
	}

#if 1
	for(i = 0; i < output_len; i++)
	{
		for(j = 0; j < GF_FIELD; j++)
		{
			if(power_polynomial_table[j][1] == tmp_output_seq[i])
			{
				output_seq[i] = power_polynomial_table[j][0];
				break;
			}
		}
	}
#else
	for(i = 0; i < output_len; i++)
	{
		output_seq[i] = power_polynomial_table[tmp_output_seq[i] + 0x1][0];
	}
#endif	

#if (1 == CFG_DEBUG_INFO)
	DEBUG_INFO("bpsk demodulation:\n");
	for(i = 0; i < output_len; i++)
	{
		DEBUG_INFO("%x %x\n", output_seq[i], cwd_poly[i]);
	}
	DEBUG_INFO("\n");
#endif

	return 0;
}

int trans_over_chnl()
{
	long long i = 0;

	/*modulation*/
	bpsk_mod(cwd_poly,
			 CODEWORD_LEN,
			 recv_seq,
			 symbol_num);

	/*transmission over channel*/
	DEBUG_NOTICE("Transmission over Channel:\n");
	for(i = 0; i < symbol_num; i++)
	{
		if(0 == (i % GF_Q))
		{
			DEBUG_NOTICE("---------------\n");
		}
		recv_seq[i][0] = recv_seq[i][0] + awgn_gen(eb2n0);
		recv_seq[i][1] = recv_seq[i][1] + awgn_gen(eb2n0);
		DEBUG_NOTICE("%f %f\n", recv_seq[i][0], recv_seq[i][1]);
	}
	DEBUG_NOTICE("\n");

#if (1 == FAC_FREE_ERR)
	recv_seq[0][0] = 1.223182;
	recv_seq[0][1] = -0.617231;
	recv_seq[1][0] = 0.289580;
	recv_seq[1][1] = 0.198964;
	recv_seq[2][0] = -1.528407;
	recv_seq[2][1] = 0.365135;
	recv_seq[3][0] = -0.472032;
	recv_seq[3][1] = -0.291959;
	recv_seq[4][0] = 0.802720;
	recv_seq[4][1] = -0.992985;
	recv_seq[5][0] = 1.333641;
	recv_seq[5][1] = -0.431431;
	recv_seq[6][0] = -0.951676;
	recv_seq[6][1] = 0.281279;
	recv_seq[7][0] = 1.662807;
	recv_seq[7][1] = 0.066072;
	recv_seq[8][0] = -0.651184;
	recv_seq[8][1] = 0.259285;
	recv_seq[9][0] = 1.475712;
	recv_seq[9][1] = 0.069349;
	recv_seq[10][0] = 0.192992;
	recv_seq[10][1] = 0.266447;
	recv_seq[11][0] = -0.677493;
	recv_seq[11][1] = -0.464497;
	recv_seq[12][0] = -1.381971;
	recv_seq[12][1] = 0.005600;
	recv_seq[13][0] = 1.540407;
	recv_seq[13][1] = -0.545868;
	recv_seq[14][0] = 0.567315;
	recv_seq[14][1] = 0.054252;
	recv_seq[15][0] = 0.759291;
	recv_seq[15][1] = -0.591471;
	recv_seq[16][0] = 1.436373;
	recv_seq[16][1] = 0.835564;
	recv_seq[17][0] = -0.511762;
	recv_seq[17][1] = -0.447861;
	recv_seq[18][0] = 2.493391;
	recv_seq[18][1] = 0.246297;
	recv_seq[19][0] = -1.384593;
	recv_seq[19][1] = 0.372882;
	recv_seq[20][0] = 0.684008;
	recv_seq[20][1] = 0.020674;
	recv_seq[21][0] = -0.523958;
	recv_seq[21][1] = 0.913362;
	recv_seq[22][0] = 0.345900;
	recv_seq[22][1] = -0.712151;
	recv_seq[23][0] = 0.400717;
	recv_seq[23][1] = 0.279110;
	recv_seq[24][0] = 1.419154;
	recv_seq[24][1] = 0.211216;
	recv_seq[25][0] = -0.416032;
	recv_seq[25][1] = 1.254062;
	recv_seq[26][0] = 1.541175;
	recv_seq[26][1] = -0.722990;
	recv_seq[27][0] = 1.338704;
	recv_seq[27][1] = 0.751498;
	recv_seq[28][0] = 0.661762;
	recv_seq[28][1] = 0.403584;
	recv_seq[29][0] = -0.470503;
	recv_seq[29][1] = -0.520841;
	recv_seq[30][0] = -0.498623;
	recv_seq[30][1] = 1.252401;
	recv_seq[31][0] = 1.179268;
	recv_seq[31][1] = -0.894205;
	recv_seq[32][0] = 1.185005;
	recv_seq[32][1] = -0.634609;
	recv_seq[33][0] = 1.589969;
	recv_seq[33][1] = 0.253451;
	recv_seq[34][0] = -1.696170;
	recv_seq[34][1] = 0.748536;
	recv_seq[35][0] = -0.863389;
	recv_seq[35][1] = -0.213783;
	recv_seq[36][0] = -1.587110;
	recv_seq[36][1] = 0.483907;
	recv_seq[37][0] = -1.363611;
	recv_seq[37][1] = 0.281917;
	recv_seq[38][0] = -0.876465;
	recv_seq[38][1] = 0.300299;
	recv_seq[39][0] = -0.672587;
	recv_seq[39][1] = 0.609814;
	recv_seq[40][0] = -1.757189;
	recv_seq[40][1] = -0.051767;
	recv_seq[41][0] = 0.654777;
	recv_seq[41][1] = 0.661581;
	recv_seq[42][0] = -0.535075;
	recv_seq[42][1] = 0.366644;
	recv_seq[43][0] = -0.911527;
	recv_seq[43][1] = -0.636879;
	recv_seq[44][0] = -0.775678;
	recv_seq[44][1] = 0.224839;
	recv_seq[45][0] = -0.496652;
	recv_seq[45][1] = -0.144938;
	recv_seq[46][0] = 1.039247;
	recv_seq[46][1] = 0.198974;
	recv_seq[47][0] = 0.883920;
	recv_seq[47][1] = 0.344826;
	recv_seq[48][0] = 1.374237;
	recv_seq[48][1] = -0.196803;
	recv_seq[49][0] = -0.761407;
	recv_seq[49][1] = -0.050654;
	recv_seq[50][0] = 0.865824;
	recv_seq[50][1] = 0.200923;
	recv_seq[51][0] = -0.834665;
	recv_seq[51][1] = -0.639207;
	recv_seq[52][0] = 1.640691;
	recv_seq[52][1] = 0.777729;
	recv_seq[53][0] = 1.425614;
	recv_seq[53][1] = -0.586578;
	recv_seq[54][0] = -1.021004;
	recv_seq[54][1] = -1.238273;
	recv_seq[55][0] = -0.499409;
	recv_seq[55][1] = 0.316808;
	recv_seq[56][0] = -0.378879;
	recv_seq[56][1] = 0.329659;
	recv_seq[57][0] = 1.767763;
	recv_seq[57][1] = -0.014929;
	recv_seq[58][0] = 1.634663;
	recv_seq[58][1] = -0.039856;
	recv_seq[59][0] = 1.447458;
	recv_seq[59][1] = 0.534436;
	recv_seq[60][0] = 1.010618;
	recv_seq[60][1] = -0.094954;
	recv_seq[61][0] = 1.107899;
	recv_seq[61][1] = -0.849894;
	recv_seq[62][0] = -0.686804;
	recv_seq[62][1] = 0.256609;
	recv_seq[63][0] = 1.705817;
	recv_seq[63][1] = -0.057623;
	recv_seq[64][0] = -0.721922;
	recv_seq[64][1] = 0.736440;
	recv_seq[65][0] = -0.849979;
	recv_seq[65][1] = -0.767309;
	recv_seq[66][0] = -1.436715;
	recv_seq[66][1] = 0.305519;
	recv_seq[67][0] = 0.349034;
	recv_seq[67][1] = -0.150997;
	recv_seq[68][0] = -1.143129;
	recv_seq[68][1] = -0.221678;
	recv_seq[69][0] = -0.929083;
	recv_seq[69][1] = -0.531951;
	recv_seq[70][0] = -1.075025;
	recv_seq[70][1] = 0.261210;
	recv_seq[71][0] = 1.301120;
	recv_seq[71][1] = -0.446650;
	recv_seq[72][0] = -1.161756;
	recv_seq[72][1] = 0.913113;
	recv_seq[73][0] = 0.206586;
	recv_seq[73][1] = -0.157823;
	recv_seq[74][0] = 0.931418;
	recv_seq[74][1] = 0.032009;
	recv_seq[75][0] = -1.540186;
	recv_seq[75][1] = 0.406192;
	recv_seq[76][0] = -0.148125;
	recv_seq[76][1] = 0.059937;
	recv_seq[77][0] = 0.665442;
	recv_seq[77][1] = -0.024901;
	recv_seq[78][0] = 0.798488;
	recv_seq[78][1] = -0.439726;
	recv_seq[79][0] = -0.254850;
	recv_seq[79][1] = -0.318682;
	recv_seq[80][0] = -0.371480;
	recv_seq[80][1] = -0.879882;
	recv_seq[81][0] = -0.947772;
	recv_seq[81][1] = -0.216643;
	recv_seq[82][0] = 1.485631;
	recv_seq[82][1] = 0.003698;
	recv_seq[83][0] = -0.863925;
	recv_seq[83][1] = 0.122629;
	recv_seq[84][0] = 0.566869;
	recv_seq[84][1] = 0.378134;
	recv_seq[85][0] = -0.838829;
	recv_seq[85][1] = 0.227345;
	recv_seq[86][0] = -0.135093;
	recv_seq[86][1] = -0.070226;
	recv_seq[87][0] = -1.251499;
	recv_seq[87][1] = -0.241566;
	recv_seq[88][0] = 0.797896;
	recv_seq[88][1] = -0.207417;
	recv_seq[89][0] = -1.127127;
	recv_seq[89][1] = -0.039909;
	recv_seq[90][0] = 0.589821;
	recv_seq[90][1] = -0.295312;
	recv_seq[91][0] = 0.530800;
	recv_seq[91][1] = 0.355583;
	recv_seq[92][0] = 0.138188;
	recv_seq[92][1] = 0.424758;
	recv_seq[93][0] = -0.980875;
	recv_seq[93][1] = -0.508447;
	recv_seq[94][0] = -1.005655;
	recv_seq[94][1] = -0.548046;
	recv_seq[95][0] = 1.952927;
	recv_seq[95][1] = -0.775135;
	recv_seq[96][0] = -0.408304;
	recv_seq[96][1] = -0.258052;
	recv_seq[97][0] = 0.944512;
	recv_seq[97][1] = -0.014509;
	recv_seq[98][0] = 0.393298;
	recv_seq[98][1] = -0.775329;
	recv_seq[99][0] = -0.698337;
	recv_seq[99][1] = -0.140452;
	recv_seq[100][0] = -1.147305;
	recv_seq[100][1] = 0.414115;
	recv_seq[101][0] = 1.845977;
	recv_seq[101][1] = 0.867027;
	recv_seq[102][0] = 0.198284;
	recv_seq[102][1] = 0.501903;
	recv_seq[103][0] = -0.827014;
	recv_seq[103][1] = 1.678582;
	recv_seq[104][0] = -0.340250;
	recv_seq[104][1] = 0.432682;
	recv_seq[105][0] = 1.347089;
	recv_seq[105][1] = 0.133094;
	recv_seq[106][0] = 1.589723;
	recv_seq[106][1] = 0.322011;
	recv_seq[107][0] = -0.843361;
	recv_seq[107][1] = -0.482485;
	recv_seq[108][0] = 1.068632;
	recv_seq[108][1] = -0.132104;
	recv_seq[109][0] = -1.267274;
	recv_seq[109][1] = 0.511562;
	recv_seq[110][0] = 0.567061;
	recv_seq[110][1] = 0.121765;
	recv_seq[111][0] = 1.221685;
	recv_seq[111][1] = 0.122934;
	recv_seq[112][0] = -1.599399;
	recv_seq[112][1] = -0.405573;
	recv_seq[113][0] = 0.003402;
	recv_seq[113][1] = 0.129455;
	recv_seq[114][0] = -0.110604;
	recv_seq[114][1] = 0.854538;
	recv_seq[115][0] = 0.289294;
	recv_seq[115][1] = -0.005816;
	recv_seq[116][0] = 0.351644;
	recv_seq[116][1] = -0.417856;
	recv_seq[117][0] = 1.003374;
	recv_seq[117][1] = 0.339170;
	recv_seq[118][0] = 1.417272;
	recv_seq[118][1] = -0.026039;
	recv_seq[119][0] = -0.951955;
	recv_seq[119][1] = -0.201115;
	recv_seq[120][0] = -1.171430;
	recv_seq[120][1] = 0.222627;
	recv_seq[121][0] = -1.260744;
	recv_seq[121][1] = -0.372628;
	recv_seq[122][0] = 1.353389;
	recv_seq[122][1] = 0.270731;
	recv_seq[123][0] = -1.165569;
	recv_seq[123][1] = -0.050727;
	recv_seq[124][0] = -0.190949;
	recv_seq[124][1] = -0.119575;
	recv_seq[125][0] = -0.910489;
	recv_seq[125][1] = -0.000360;
	recv_seq[126][0] = 0.776410;
	recv_seq[126][1] = -0.297800;
	recv_seq[127][0] = 1.452141;
	recv_seq[127][1] = 0.245885;
	recv_seq[128][0] = -0.740130;
	recv_seq[128][1] = -0.072316;
	recv_seq[129][0] = -0.490001;
	recv_seq[129][1] = -0.149365;
	recv_seq[130][0] = 0.919264;
	recv_seq[130][1] = 0.619385;
	recv_seq[131][0] = -0.681473;
	recv_seq[131][1] = 0.211930;
	recv_seq[132][0] = 0.845584;
	recv_seq[132][1] = 0.483132;
	recv_seq[133][0] = 0.205868;
	recv_seq[133][1] = -0.438433;
	recv_seq[134][0] = 0.811491;
	recv_seq[134][1] = -0.471551;
	recv_seq[135][0] = 0.408564;
	recv_seq[135][1] = -0.084568;
	recv_seq[136][0] = 1.255632;
	recv_seq[136][1] = -0.200727;
	recv_seq[137][0] = -0.405975;
	recv_seq[137][1] = -0.090135;
	recv_seq[138][0] = -1.380686;
	recv_seq[138][1] = 0.107264;
	recv_seq[139][0] = -1.498891;
	recv_seq[139][1] = 0.004822;
	recv_seq[140][0] = -1.511555;
	recv_seq[140][1] = -0.354548;
	recv_seq[141][0] = -0.109031;
	recv_seq[141][1] = -0.044451;
	recv_seq[142][0] = -1.418559;
	recv_seq[142][1] = 0.556679;
	recv_seq[143][0] = -1.105240;
	recv_seq[143][1] = 0.494688;
	recv_seq[144][0] = -0.620497;
	recv_seq[144][1] = -0.889405;
	recv_seq[145][0] = 1.179033;
	recv_seq[145][1] = 0.257889;
	recv_seq[146][0] = 0.589107;
	recv_seq[146][1] = -0.613332;
	recv_seq[147][0] = 0.886306;
	recv_seq[147][1] = 0.011780;
	recv_seq[148][0] = 0.789603;
	recv_seq[148][1] = -0.265471;
	recv_seq[149][0] = 0.739893;
	recv_seq[149][1] = -0.330016;
	recv_seq[150][0] = -1.290530;
	recv_seq[150][1] = 0.226661;
	recv_seq[151][0] = -1.594738;
	recv_seq[151][1] = -0.316906;
	recv_seq[152][0] = -0.121591;
	recv_seq[152][1] = 0.709074;
	recv_seq[153][0] = -0.783322;
	recv_seq[153][1] = -0.466263;
	recv_seq[154][0] = -0.860603;
	recv_seq[154][1] = -0.063415;
	recv_seq[155][0] = -0.953303;
	recv_seq[155][1] = 1.330576;
	recv_seq[156][0] = -0.327413;
	recv_seq[156][1] = -0.248817;
	recv_seq[157][0] = 0.771021;
	recv_seq[157][1] = 0.111511;
	recv_seq[158][0] = -0.288825;
	recv_seq[158][1] = -0.172358;
	recv_seq[159][0] = -0.146463;
	recv_seq[159][1] = -0.511625;
	recv_seq[160][0] = 1.105475;
	recv_seq[160][1] = 0.644843;
	recv_seq[161][0] = -0.818426;
	recv_seq[161][1] = 0.235383;
	recv_seq[162][0] = -1.789085;
	recv_seq[162][1] = 0.642197;
	recv_seq[163][0] = 0.182733;
	recv_seq[163][1] = -0.430360;
	recv_seq[164][0] = -0.654271;
	recv_seq[164][1] = -0.576451;
	recv_seq[165][0] = 1.473393;
	recv_seq[165][1] = -0.130025;
	recv_seq[166][0] = -1.546137;
	recv_seq[166][1] = -0.119377;
	recv_seq[167][0] = -1.098724;
	recv_seq[167][1] = 0.166199;
	recv_seq[168][0] = -1.182958;
	recv_seq[168][1] = -0.478985;
	recv_seq[169][0] = 0.440355;
	recv_seq[169][1] = 0.153238;
	recv_seq[170][0] = -1.149320;
	recv_seq[170][1] = -1.115355;
	recv_seq[171][0] = 1.759917;
	recv_seq[171][1] = 0.730710;
	recv_seq[172][0] = -1.354483;
	recv_seq[172][1] = -0.947050;
	recv_seq[173][0] = 1.473485;
	recv_seq[173][1] = -0.227563;
	recv_seq[174][0] = 0.516806;
	recv_seq[174][1] = -0.008235;
	recv_seq[175][0] = 1.500356;
	recv_seq[175][1] = 0.516366;
	recv_seq[176][0] = 0.501815;
	recv_seq[176][1] = 0.199455;
	recv_seq[177][0] = 1.977858;
	recv_seq[177][1] = -0.302449;
	recv_seq[178][0] = 1.454348;
	recv_seq[178][1] = 1.053999;
	recv_seq[179][0] = 1.268303;
	recv_seq[179][1] = 0.073739;
	recv_seq[180][0] = -1.009692;
	recv_seq[180][1] = 0.145854;
	recv_seq[181][0] = -1.323052;
	recv_seq[181][1] = 0.447657;
	recv_seq[182][0] = -1.216761;
	recv_seq[182][1] = -0.696714;
	recv_seq[183][0] = 1.374593;
	recv_seq[183][1] = -0.053687;
	recv_seq[184][0] = -0.371766;
	recv_seq[184][1] = 0.225775;
	recv_seq[185][0] = 0.319254;
	recv_seq[185][1] = 0.561966;
	recv_seq[186][0] = -0.660184;
	recv_seq[186][1] = -0.209443;
	recv_seq[187][0] = 1.853436;
	recv_seq[187][1] = 0.002311;
	recv_seq[188][0] = -1.073101;
	recv_seq[188][1] = 0.576557;
	recv_seq[189][0] = -1.615438;
	recv_seq[189][1] = 0.911785;
	recv_seq[190][0] = 1.547783;
	recv_seq[190][1] = 0.157932;
	recv_seq[191][0] = -0.274094;
	recv_seq[191][1] = -0.812976;
	recv_seq[192][0] = -1.846771;
	recv_seq[192][1] = -0.071138;
	recv_seq[193][0] = 1.154477;
	recv_seq[193][1] = -0.194417;
	recv_seq[194][0] = -0.715917;
	recv_seq[194][1] = 0.434872;
	recv_seq[195][0] = 1.257221;
	recv_seq[195][1] = -0.476764;
	recv_seq[196][0] = -0.342300;
	recv_seq[196][1] = -0.408150;
	recv_seq[197][0] = -1.370703;
	recv_seq[197][1] = 0.116520;
	recv_seq[198][0] = 0.543902;
	recv_seq[198][1] = 0.542876;
	recv_seq[199][0] = -1.104142;
	recv_seq[199][1] = 0.439927;
	recv_seq[200][0] = 1.441155;
	recv_seq[200][1] = -0.348315;
	recv_seq[201][0] = 1.370994;
	recv_seq[201][1] = -0.161057;
	recv_seq[202][0] = -1.336851;
	recv_seq[202][1] = -0.076346;
	recv_seq[203][0] = 0.444972;
	recv_seq[203][1] = 0.788213;
	recv_seq[204][0] = -1.423921;
	recv_seq[204][1] = 0.622129;
	recv_seq[205][0] = -0.174315;
	recv_seq[205][1] = 0.625916;
	recv_seq[206][0] = 0.924853;
	recv_seq[206][1] = -0.751415;
	recv_seq[207][0] = 0.338128;
	recv_seq[207][1] = -0.203695;
	recv_seq[208][0] = -0.008305;
	recv_seq[208][1] = -0.244734;
	recv_seq[209][0] = 0.546377;
	recv_seq[209][1] = -0.462562;
	recv_seq[210][0] = -0.944519;
	recv_seq[210][1] = -0.575435;
	recv_seq[211][0] = 0.602400;
	recv_seq[211][1] = 0.977649;
	recv_seq[212][0] = -1.684682;
	recv_seq[212][1] = -0.471498;
	recv_seq[213][0] = -0.562798;
	recv_seq[213][1] = 0.616111;
	recv_seq[214][0] = 1.166983;
	recv_seq[214][1] = -0.352753;
	recv_seq[215][0] = -1.663735;
	recv_seq[215][1] = -0.233305;
	recv_seq[216][0] = 1.105075;
	recv_seq[216][1] = 0.138165;
	recv_seq[217][0] = 0.172044;
	recv_seq[217][1] = 0.406637;
	recv_seq[218][0] = 1.736266;
	recv_seq[218][1] = 0.225728;
	recv_seq[219][0] = -0.883650;
	recv_seq[219][1] = -0.776135;
	recv_seq[220][0] = 0.620711;
	recv_seq[220][1] = -0.356952;
	recv_seq[221][0] = -0.255797;
	recv_seq[221][1] = -0.043432;
	recv_seq[222][0] = 0.716091;
	recv_seq[222][1] = 0.087897;
	recv_seq[223][0] = 0.918219;
	recv_seq[223][1] = -0.203549;
	recv_seq[224][0] = -0.403385;
	recv_seq[224][1] = -0.300514;
	recv_seq[225][0] = -1.764171;
	recv_seq[225][1] = -0.152641;
	recv_seq[226][0] = 0.722750;
	recv_seq[226][1] = -0.233085;
	recv_seq[227][0] = -1.824811;
	recv_seq[227][1] = -0.722161;
	recv_seq[228][0] = -1.818634;
	recv_seq[228][1] = -0.081352;
	recv_seq[229][0] = -0.488786;
	recv_seq[229][1] = 0.491407;
	recv_seq[230][0] = 1.089487;
	recv_seq[230][1] = 0.365430;
	recv_seq[231][0] = -0.630902;
	recv_seq[231][1] = 0.290781;
	recv_seq[232][0] = 2.247470;
	recv_seq[232][1] = -0.634318;
	recv_seq[233][0] = -1.517969;
	recv_seq[233][1] = 0.177747;
	recv_seq[234][0] = -0.037446;
	recv_seq[234][1] = 0.598214;
	recv_seq[235][0] = -1.300966;
	recv_seq[235][1] = -0.096489;
	recv_seq[236][0] = -1.180854;
	recv_seq[236][1] = -0.049246;
	recv_seq[237][0] = 0.345273;
	recv_seq[237][1] = 0.337547;
	recv_seq[238][0] = 0.508105;
	recv_seq[238][1] = 0.023447;
	recv_seq[239][0] = -0.013967;
	recv_seq[239][1] = 0.338262;
	recv_seq[240][0] = -1.637301;
	recv_seq[240][1] = 0.627288;
	recv_seq[241][0] = -0.426396;
	recv_seq[241][1] = 0.826171;
	recv_seq[242][0] = 0.903485;
	recv_seq[242][1] = -0.734606;
	recv_seq[243][0] = 1.008204;
	recv_seq[243][1] = 0.243225;
	recv_seq[244][0] = 2.058619;
	recv_seq[244][1] = -0.755621;
	recv_seq[245][0] = -0.452009;
	recv_seq[245][1] = 1.091722;
	recv_seq[246][0] = -0.690703;
	recv_seq[246][1] = -0.226176;
	recv_seq[247][0] = 1.261381;
	recv_seq[247][1] = -0.057279;
	recv_seq[248][0] = 0.478669;
	recv_seq[248][1] = -0.369351;
	recv_seq[249][0] = -1.897455;
	recv_seq[249][1] = -0.254226;
	recv_seq[250][0] = -1.209250;
	recv_seq[250][1] = 0.493786;
	recv_seq[251][0] = 0.617716;
	recv_seq[251][1] = -0.045099;
	recv_seq[252][0] = -0.766410;
	recv_seq[252][1] = -0.210499;
	recv_seq[253][0] = -1.541449;
	recv_seq[253][1] = -0.373458;
	recv_seq[254][0] = 1.672555;
	recv_seq[254][1] = 1.265977;
	recv_seq[255][0] = -1.229820;
	recv_seq[255][1] = 0.284234;
#endif

	/*demodulation*/
	bpsk_demod((float **)recv_seq,
				symbol_num,
				recv_poly,
				CODEWORD_LEN);

#if 0//(1 == TEST_MODE)
	recv_poly[3] = 0x7;
	recv_poly[7] = 0x7;
	recv_poly[25] = 0x7;
	recv_poly[31] = 0xff;
	recv_poly[35] = 0x0;
	recv_poly[38] = 0xd;
	recv_poly[46] = 0x2;
	recv_poly[47] = 0xc;
	recv_poly[49] = 0x4;
#endif

	return 0;
}
