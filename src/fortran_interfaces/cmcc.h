#ifndef CMCC_H
#define CMCC_H

/* Include the configuration macros describing the linkage status to other libraries */
#include "cmc_config.h"

#include "stdlib.h"
#include "stdio.h"
#include "stdarg.h"
#include "string.h"

void cmcc_debug_msg(const char* _msg_format, ...);

_Noreturn void cmcc_err_msg(const char* _msg_format, ...);


#endif /* CMCC_H */