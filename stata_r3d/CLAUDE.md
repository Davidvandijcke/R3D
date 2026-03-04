# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build and Development Commands

### Building the Plugin (Optional - for performance)
```bash
cd plugin/
make clean
make
make install  # Copies plugin to ado/ directory
```

### Running Tests
```stata
# Test main functionality
/Applications/Stata/StataSE.app/Contents/MacOS/stataSE -b do example.do

# Test plugin loading
/Applications/Stata/StataSE.app/Contents/MacOS/stataSE -b do diagnose_plugin.do

# Verify package structure
/Applications/Stata/StataSE.app/Contents/MacOS/stataSE -b do verify_package.do

# Run syntax tests
/Applications/Stata/StataSE.app/Contents/MacOS/stataSE -b do syntax_test.do
```

### Installation Commands
```stata
# From local directory
net install r3d, from("/path/to/stata_r3d/") replace

# Or manually
do install_r3d.do

# Post-installation setup
r3d_setup
```

## Architecture Overview

This is a Stata package implementing Regression Discontinuity with Distribution-valued Outcomes (R3D). The package has a three-layer architecture:

1. **Stata Interface Layer** (`ado/*.ado`): User-facing commands that handle syntax parsing, data preparation, and result presentation
2. **Mata Computational Layer** (`mata/*.mata`): Core statistical algorithms implemented in Mata for performance
3. **Optional Plugin Layer** (`plugin/`): C/Fortran plugin for accelerated computation with large datasets

### Key Components

- **Main Command** (`ado/r3d.ado`): Entry point that handles sharp/fuzzy RD designs, automatically compiles Mata functions on first run, and attempts to load the optional plugin
- **Bootstrap Inference** (`ado/r3d_bootstrap.ado`): Post-estimation bootstrap for uniform confidence bands
- **Bandwidth Selection** (`ado/r3d_bwselect.ado`): Data-driven MSE/IMSE-optimal bandwidth selection
- **Plugin Loader** (`ado/r3d_plugin_load.ado`): Handles platform-specific plugin loading with fallback to pure Mata
- **Mata Functions** (`mata/r3d_mata.mata`): Core algorithms including quantile computation, local polynomial regression, and Frechet bounds
- **Plugin Interface** (`mata/r3d_plugin_interface.mata`): Bridge between Mata and the C/Fortran plugin

### Plugin Architecture

The plugin system uses a C wrapper (`r3d_plugin.c`) that interfaces between Stata and Fortran computational routines (`r3d_fortran.f90`). The plugin is optional - the package falls back to pure Mata if unavailable. On macOS, the plugin uses static linking of Fortran libraries to avoid runtime dependencies.

## Important Implementation Details

- The package supports Stata 14.0+ but requires Stata 18.0 for full functionality
- Mata functions are auto-compiled on first use via the main `r3d` command
- The plugin uses the global `R3D_USE_PLUGIN` to track availability
- Multiple outcome variables (y1, y2, ...) represent repeated measurements from each unit's distribution
- The package normalizes option aliases (e.g., `nquantiles` to `nq`) before parsing

## Package Distribution

The package uses standard Stata distribution format:
- `stata.toc`: Package index file
- `r3d.pkg`: Package description with file manifest
- Files marked with `F` in the .pkg file (like the plugin) are optional

## Testing Approach

When testing changes, verify:
1. Package structure with `verify_package.do`
2. Plugin loading with `diagnose_plugin.do`
3. Basic functionality with `example.do`
4. Edge cases and syntax variations with `syntax_test.do`



# Using Gemini CLI for Large Codebase Analysis

When analyzing large codebases or multiple files that might exceed context limits, use the Gemini CLI with its massive
context window. Use `gemini -p` to leverage Google Gemini's large context capacity.

Before calling gemini export this as a variable, GOOGLE_CLOUD_PROJECT="gen-lang-client-0115392130". 

You can give Gemini ALL the relevant files it has a huge context window you just need to give it a long enough timeout to respond 

## File and Directory Inclusion Syntax

Use the `@` syntax to include files and directories in your Gemini prompts. The paths should be relative to WHERE you run the
  gemini command:

### Examples:

**Single file analysis:**
gemini -p "@src/main.py Explain this file's purpose and structure"

Multiple files:
gemini -p "@package.json @src/index.js Analyze the dependencies used in the code"

Entire directory:
gemini -p "@src/ Summarize the architecture of this codebase"

Multiple directories:
gemini -p "@src/ @tests/ Analyze test coverage for the source code"

Current directory and subdirectories:
gemini -p "@./ Give me an overview of this entire project"

# Or use --all_files flag:
gemini --all_files -p "Analyze the project structure and dependencies"

Implementation Verification Examples

Check if a feature is implemented:
gemini -p "@src/ @lib/ Has dark mode been implemented in this codebase? Show me the relevant files and functions"

Verify authentication implementation:
gemini -p "@src/ @middleware/ Is JWT authentication implemented? List all auth-related endpoints and middleware"

Check for specific patterns:
gemini -p "@src/ Are there any React hooks that handle WebSocket connections? List them with file paths"

Verify error handling:
gemini -p "@src/ @api/ Is proper error handling implemented for all API endpoints? Show examples of try-catch blocks"

Check for rate limiting:
gemini -p "@backend/ @middleware/ Is rate limiting implemented for the API? Show the implementation details"

Verify caching strategy:
gemini -p "@src/ @lib/ @services/ Is Redis caching implemented? List all cache-related functions and their usage"

Check for specific security measures:
gemini -p "@src/ @api/ Are SQL injection protections implemented? Show how user inputs are sanitized"

Verify test coverage for features:
gemini -p "@src/payment/ @tests/ Is the payment processing module fully tested? List all test cases"

When to Use Gemini CLI

Use gemini -p when:
- Analyzing entire codebases or large directories
- Comparing multiple large files
- Need to understand project-wide patterns or architecture
- Current context window is insufficient for the task
- Working with files totaling more than 100KB
- Verifying if specific features, patterns, or security measures are implemented
- Checking for the presence of certain coding patterns across the entire codebase

Important Notes

- Paths in @ syntax are relative to your current working directory when invoking gemini
- The CLI will include file contents directly in the context
- No need for --yolo flag for read-only analysis
- Gemini's context window can handle entire codebases that would overflow Claude's context
- When checking implementations, be specific about what you're looking for to get accurate results
