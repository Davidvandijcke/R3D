/*
 * simple_test.c - Simplest possible Stata plugin
 */

#include "stplugin.h"

STDLL stata_call(int argc, char *argv[])
{
    // Just return success
    return 0;
}