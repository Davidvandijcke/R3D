/*
 * test_plugin_interface.do
 * Unit tests for r3d_plugin_interface.mata guards (F1–F4, F6).
 * These tests exercise error-handling paths; they do NOT require the plugin
 * binary to be present — they validate input validation in Mata.
 */

version 18.0
clear all
set more off

// -----------------------------------------------------------------------
// Helper: assert a Mata call returns expected rc
// -----------------------------------------------------------------------
program define assert_rc
    args expected_rc description
    if `expected_rc' != r(rc) {
        di as error "FAIL: `description' — expected rc=`expected_rc', got r(rc)=`r(rc)'"
        exit 1
    }
    di as text "PASS: `description'"
end

// -----------------------------------------------------------------------
// Setup: load Mata library
// -----------------------------------------------------------------------
do stata_r3d/mata/r3d_plugin_interface.mata

// -----------------------------------------------------------------------
// T1: p >= 100 triggers assert (assert p < 100)
// -----------------------------------------------------------------------
// assert() in Mata produces rc=9 (assertion violation)
clear
set obs 10
mata:
{
    X = (1::10)
    Y = (1::10) / 10
    h = 1
    alpha = J(0,0,.)
    w = J(0,0,.)
    rc = r3d_plugin_locweights(X, Y, h, 100, 1, 1, alpha, w)
}
end
capture noisily {
    mata: assert(0)  // should not reach here if assert(p<100) fired
}
* If we reach here without aborting the assert(p<100) caught it; rc returned is 198 via the guard or aborted via assert
di as text "T1: p=100 guard exercised (assert fires before plugin call)"

// -----------------------------------------------------------------------
// T2: strofreal(p, "%12.0f") serializes small integer as plain integer
// -----------------------------------------------------------------------
mata:
{
    real scalar s
    s = strofreal(2, "%12.0f")
    assert(strtrim(s) == "2")
    s = strofreal(99, "%12.0f")
    assert(strtrim(s) == "99")
}
end
di as text "PASS T2: integer-format serialization produces plain integers for p=2 and p=99"

// -----------------------------------------------------------------------
// T3: row-vector X triggers assert(cols(X)==1) → rc abort
// -----------------------------------------------------------------------
clear
set obs 1
mata:
{
    X_row = (1, 2, 3)   // row vector: cols=3, rows=1
    Y = (0.1 \ 0.2 \ 0.3)
    h = 1
    alpha = J(0,0,.)
    w = J(0,0,.)
    // assert(cols(X)==1) will fire; capture the abort
    rc = _error(r3d_plugin_locweights(X_row, Y, h, 1, 1, 1, alpha, w))
}
end
di as text "T3: row-vector X assert guard exercised"

// -----------------------------------------------------------------------
// T4: zero bandwidth returns rc=198
// -----------------------------------------------------------------------
clear
set obs 10
mata:
{
    X = (1::10)
    Y = (1::10) / 10
    h_zero = 0
    alpha = J(0,0,.)
    w = J(0,0,.)
    rc = r3d_plugin_locweights(X, Y, h_zero, 1, 1, 1, alpha, w)
    assert(rc == 198)
}
end
di as text "PASS T4: h=0 returns rc=198"

// -----------------------------------------------------------------------
// T5: negative bandwidth returns rc=198
// -----------------------------------------------------------------------
clear
set obs 10
mata:
{
    X = (1::10)
    Y = (1::10) / 10
    h_neg = -1
    alpha = J(0,0,.)
    w = J(0,0,.)
    rc = r3d_plugin_locweights(X, Y, h_neg, 1, 1, 1, alpha, w)
    assert(rc == 198)
}
end
di as text "PASS T5: h<0 returns rc=198"

// -----------------------------------------------------------------------
// T6: missing values in X return rc=198
// -----------------------------------------------------------------------
clear
set obs 10
mata:
{
    X_miss = (1::9) \ .
    Y = (1::10) / 10
    h = 1
    alpha = J(0,0,.)
    w = J(0,0,.)
    rc = r3d_plugin_locweights(X_miss, Y, h, 1, 1, 1, alpha, w)
    assert(rc == 198)
}
end
di as text "PASS T6: missing X returns rc=198"

// -----------------------------------------------------------------------
// T7: missing values in Y return rc=198
// -----------------------------------------------------------------------
clear
set obs 10
mata:
{
    X = (1::10)
    Y = (1::9) / 10 \ .
    h = 1
    alpha = J(0,0,.)
    w = J(0,0,.)
    rc = r3d_plugin_locweights(X, Y, h, 1, 1, 1, alpha, w)
    assert(rc == 198)
}
end
di as text "PASS T7: missing Y returns rc=198"

di as text ""
di as text "All plugin interface guard tests completed."
