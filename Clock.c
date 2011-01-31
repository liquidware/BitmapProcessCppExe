/*
 * clock.c
 *
 *  Created on: Jan 7, 2011
 *      Author: chris
 */
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <sys/time.h>

FILE *clock_log;
typedef struct ClockStack_Typedef {
	char    function_name[128];
	int64_t clock_in;
} cst;

static cst      clock_stack[2048];
//static int64_t  t0, t1;
static int      stack_depth = 0;

void clock_start(const char fname[64]);
void clock_end();


// Testbed function implementation
void clock_start(const char fname[64]) {
    int i, j;
    char indent[32];
    struct    timeval    tv;
    int64_t    now;

    gettimeofday(&tv,NULL);
    now = tv.tv_sec * (int64_t)1000000 + tv.tv_usec;

    strcpy(clock_stack[stack_depth].function_name, fname);

    j=0;
    for (i=0; i<stack_depth; i++) {
        indent[j++] = ' ';
        indent[j++] = ' ';
    }
    indent[j] = 0;

#   ifdef LOGFILE_ENABLED
        fprintf(clock_log, "%s%s - START\n", indent, clock_stack[stack_depth].function_name);
#   else
        printf("%s%s - START\n", indent, clock_stack[stack_depth].function_name);
#   endif

    clock_stack[stack_depth].clock_in = now;

    stack_depth++;
}


void clock_end() {
	int64_t clock_out;
    long int clock_count;
    int i, j;
    char indent[32];

    struct    timeval    tv;
    int64_t    now;

    gettimeofday(&tv,NULL);
    now = tv.tv_sec * (int64_t)1000000 + tv.tv_usec;

    stack_depth--;
    clock_out = now;
    clock_count = (long int)clock_out - (long int)(clock_stack[stack_depth].clock_in);

    j=0;
    for (i=0; i<stack_depth; i++) {
        indent[j++] = ' ';
        indent[j++] = ' ';
    }
    indent[j] = 0;

#   ifdef LOGFILE_ENABLED
        fprintf(clock_log, "%s%s - END %ld\n",  indent, \
                                            clock_stack[stack_depth].function_name, \
                                            clock_count);
#   else
        printf("%s%s - END %ld\n",  indent, clock_stack[stack_depth].function_name, clock_count);
#   endif
}
