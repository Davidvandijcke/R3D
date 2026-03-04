/*
 * install_r3d.do - Installation script for R3D package
 */

version 14.0

di as text _n "Installing R3D package for STATA..."
di as text "{hline 60}"

// Get the current directory
local pwd = c(pwd)

// Find personal ado directory
adopath
local personaldir : sysdir PERSONAL

di as text _n "Installing to: " as result "`personaldir'"

// Copy ado files
di as text _n "Copying ado files..."
foreach file in r3d.ado r3d_bootstrap.ado r3d_bwselect.ado r3d_plugin_load.ado {
    copy "ado/`file'" "`personaldir'`file'", replace
    di as text "  - `file' " as result "installed"
}

// Copy mata file
di as text _n "Copying Mata functions..."
copy "mata/r3d_mata.mata" "`personaldir'r3d_mata.mata", replace
di as text "  - r3d_mata.mata " as result "installed"

// Copy help file
di as text _n "Copying documentation..."
copy "doc/r3d.sthlp" "`personaldir'r3d.sthlp", replace
di as text "  - r3d.sthlp " as result "installed"

// Check for plugin
di as text _n "Checking for plugin..."
capture confirm file "plugin/r3d_plugin.plugin"
if _rc == 0 {
    // Plugin exists
    copy "plugin/r3d_plugin.plugin" "`personaldir'r3d_plugin.plugin", replace
    di as text "  - r3d_plugin.plugin " as result "installed"
    
    // Test if plugin loads
    capture program drop r3d_plugin
    capture program r3d_plugin, plugin
    if _rc == 0 {
        di as text _n "Plugin " as result "successfully loaded"
        di as text "The package will use the Fortran plugin for better performance."
    }
    else {
        di as text _n "Plugin found but " as error "could not be loaded"
        di as text "The package will use the pure Mata implementation."
    }
}
else {
    di as text "  - Plugin not found"
    di as text "The package will use the pure Mata implementation."
    di as text "To use the plugin, build it with 'make' in the plugin directory."
}

// Compile Mata functions
di as text _n "Compiling Mata functions..."
quietly do "`personaldir'r3d_mata.mata"
di as text "  - Mata functions " as result "compiled"

// Test installation
di as text _n "Testing installation..."
help r3d

di as text _n "{hline 60}"
di as text "Installation " as result "complete!"
di as text _n "To get started, type: " as result "help r3d"
di as text "To run tests, type: " as result "do test/test_r3d.do"