#include "cmcc.h"

_Noreturn void cmcc_err_msg(const char* _msg_format, ...)
{
    va_list args;
    fprintf(stderr, "[cmc] ERROR: ");
    va_start(args, _msg_format);
    vfprintf(stderr, _msg_format, args);
    va_end(args);
    exit(EXIT_FAILURE);
}

void cmcc_debug_msg(const char* _msg_format, ...)
{
    va_list args;
    fprintf(stdout, "[cmc] DEBUG: ");
    va_start(args, _msg_format);
    vfprintf(stderr, _msg_format, args);
    va_end(args);
}
