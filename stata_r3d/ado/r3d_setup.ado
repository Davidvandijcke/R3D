*! r3d_setup - Post-installation setup for R3D package
*! Version 1.0.0

program define r3d_setup
    version 18.0
    
    di as text _n "Setting up R3D package..."
    
    // Compile Mata functions
    di as text "Compiling Mata functions..."
    
    // Find the mata directory
    local matapath = ""
    capture findfile r3d_mata.mata
    if _rc == 0 {
        local matapath = r(fn)
    }
    else {
        // Try relative path from ado directory
        capture findfile ../mata/r3d_mata.mata
        if _rc == 0 {
            local matapath = r(fn)
        }
    }
    
    if "`matapath'" != "" {
        quietly capture do "`matapath'"
        if _rc == 0 {
            di as text "  - r3d_mata.mata " as result "compiled successfully"
        }
        else {
            di as error "  - Error compiling r3d_mata.mata"
        }
        
        // Also compile plugin interface
        local matadir = substr("`matapath'", 1, length("`matapath'") - 13)
        quietly capture do "`matadir'r3d_plugin_interface.mata"
        if _rc == 0 {
            di as text "  - r3d_plugin_interface.mata " as result "compiled successfully"
        }
        else {
            di as error "  - Error compiling r3d_plugin_interface.mata"
        }
    }
    else {
        di as error "  - Cannot find Mata files"
        di as text "    You may need to manually run the .mata files"
    }
    
    // Test plugin loading
    di as text "Testing plugin availability..."
    quietly capture r3d_plugin_load
    
    if "$R3D_USE_PLUGIN" == "1" {
        di as text "  - Plugin " as result "loaded successfully"
        di as text "    R3D will use the Fortran plugin for better performance"
    }
    else {
        di as text "  - Plugin not available, using pure Mata implementation"
        di as text "    For better performance, build the plugin with 'make' in the plugin directory"
    }
    
    // Test basic functionality
    di as text "Testing basic functionality..."
    quietly capture help r3d
    if _rc == 0 {
        di as text "  - Help file " as result "accessible"
    }
    else {
        di as error "  - Error accessing help file"
    }
    
    di as text _n "Setup complete!"
    di as text "Type " as result "help r3d" as text " to get started"
    
end