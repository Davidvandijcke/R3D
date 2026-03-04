/*
 * test_minimal.c - Minimal Stata plugin for testing
 */

#include "stplugin.h"

STDLL stata_call(int argc, char *argv[])
{
    SF_display("Hello from plugin!\n");
    return 0;
}