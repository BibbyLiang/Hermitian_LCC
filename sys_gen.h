#ifndef SYS_GEN_H
#define SYS_GEN_H

#include "cfg_decoding.h"

extern int sys_kot_intp(unsigned char poly_seq[][MAX_POLY_TERM_SIZE], unsigned char *input_point);
extern int sys_gen_test();
extern int sys_br_test();
extern int sys_kot_test();

#endif
