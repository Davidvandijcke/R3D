/*
 * verify_package.do - Verify package structure before installation
 */

version 14.0
clear all

di as text _n "{hline 60}"
di as text "R3D Package Structure Verification"
di as text "{hline 60}"

local pwd = c(pwd)
di as text _n "Current directory: " as result "`pwd'"

// Check if we're in the right directory
capture confirm file "stata.toc"
if _rc == 0 {
    di as result "✓ stata.toc found"
}
else {
    di as error "✗ stata.toc not found"
    di as error "You need to be in the stata_r3d directory!"
    exit 1
}

capture confirm file "r3d.pkg"
if _rc == 0 {
    di as result "✓ r3d.pkg found"
}
else {
    di as error "✗ r3d.pkg not found"
    exit 1
}

// Check main ado files
local adofiles "r3d r3d_bootstrap r3d_bwselect r3d_plugin_load r3d_setup"
foreach file of local adofiles {
    capture confirm file "ado/`file'.ado"
    if _rc == 0 {
        di as result "✓ ado/`file'.ado"
    }
    else {
        di as error "✗ ado/`file'.ado missing"
    }
}

// Check mata files
capture confirm file "mata/r3d_mata.mata"
if _rc == 0 {
    di as result "✓ mata/r3d_mata.mata"
}
else {
    di as error "✗ mata/r3d_mata.mata missing"
}

capture confirm file "mata/r3d_plugin_interface.mata"
if _rc == 0 {
    di as result "✓ mata/r3d_plugin_interface.mata"
}
else {
    di as error "✗ mata/r3d_plugin_interface.mata missing"
}

// Check plugin
capture confirm file "ado/r3d_plugin.plugin"
if _rc == 0 {
    di as result "✓ ado/r3d_plugin.plugin"
}
else {
    di as text "• ado/r3d_plugin.plugin missing (optional)"
}

// Check help files
local helpfiles "r3d r3d_bootstrap r3d_bwselect"
foreach file of local helpfiles {
    capture confirm file "doc/`file'.sthlp"
    if _rc == 0 {
        di as result "✓ doc/`file'.sthlp"
    }
    else {
        di as error "✗ doc/`file'.sthlp missing"
    }
}

di as text _n "{hline 60}"
di as result "Package structure verification complete!"

di as text _n "To install, run:"
di as result "  net install r3d, from(.) replace"

di as text _n "To test after installation:"
di as result "  do test_simple.do"

di as text _n "{hline 60}"