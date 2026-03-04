*! r3d_plugin_load - Load the R3D plugin for the current platform
*! Version 1.0.0

program define r3d_plugin_load
    version 18.0

    // Default assumption: plugin unavailable until proven otherwise
    global R3D_USE_PLUGIN 0
    global R3D_PLUGIN_NAME ""

    // First try the generic plugin name shipped with the package
    capture program r3d_plugin, plugin
    if _rc == 0 {
        global R3D_PLUGIN_NAME r3d_plugin
        global R3D_USE_PLUGIN 1
        exit 0
    }

    // If the generic name is missing, fall back to platform-specific builds
    local os = c(os)
    local bit = c(bit)
    local plugin_file ""

    if "`os'" == "Windows" {
        if `bit' == 64 local plugin_file "r3d_plugin_win64"
        else           local plugin_file "r3d_plugin_win32"
    }
    else if "`os'" == "MacOSX" {
        local machine = c(machine_type)
        if strpos("`machine'", "Intel") > 0 local plugin_file "r3d_plugin_mac_intel"
        else                                  local plugin_file "r3d_plugin"
    }
    else if "`os'" == "Unix" {
        if `bit' == 64 local plugin_file "r3d_plugin_linux64"
        else           local plugin_file "r3d_plugin_linux32"
    }
    else {
        di as error "Unsupported operating system: `os'"
        exit 198
    }

    if "`plugin_file'" == "" {
        exit 0
    }

    capture program `plugin_file', plugin
    if _rc == 0 {
        global R3D_PLUGIN_NAME `plugin_file'
        global R3D_USE_PLUGIN 1
        exit 0
    }

    // If we get here no plugin was found; stay on Mata backend but inform user once.
    di as text "Note: R3D plugin not found, reverting to Mata implementation"
    di as text "      Run 'make' in the plugin/ directory to build r3d_plugin.plugin"
end
