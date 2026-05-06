#!/bin/bash

# Default values
MAIN_MACRO=""
ARGS=""
CLEAN=false
LINK=false
DEPS=("particle_tree.C" "pion.C" "event.C" "pool.C" "utils.C" "kt_parser.C") # List your dependencies here

usage() {
    echo "Usage: $0 [-m macro.C] [-a 'args'] [-c] [-l]"
    echo "  -m : Main macro file"
    echo "  -a : Arguments for the macro (comma-separated)"
    echo "  -c : Clean up build artifacts"
    echo "  -l : Link/Compile"
    exit 1
}

while getopts "m:a:cl" opt; do
    case ${opt} in
        m ) MAIN_MACRO=$OPTARG ;;
        a ) ARGS=$OPTARG ;;
        c ) CLEAN=true ;;
        l ) LINK=true ;;
        * ) usage ;;
    esac
done

# 1. Cleanup
if [ "$CLEAN" = true ]; then
    echo "Cleaning build artifacts..."
    rm -f *.d *.so *.pcm *.pcm.bin
    if [ -z "$MAIN_MACRO" ]; then exit 0; fi
fi

if [ -z "$MAIN_MACRO" ] && [ "$LINK" = false ]; then
    echo "Error: No macro specified."
    usage
fi

# 2. Prepare the command string
LINK_CMDS=""
for dep in "${DEPS[@]}"; do
    if [ -f "$dep" ]; then
        LINK_CMDS+="gInterpreter->LoadMacro(\"$dep+\"); "
    fi
done

# 3. Execution Logic
if [ "$LINK" = true ] && [ -z "$MAIN_MACRO" ]; then
    echo "Linking and Compiling only: $MAIN_MACRO"

    root -l -b -q -e "$LINK_CMDS"

    echo "Compilation successful."
elif [ "$LINK" = true ] && [ -n "$MAIN_MACRO" ]; then
    echo "Linking and Running $MAIN_MACRO with args: ($ARGS)"

    root -l -b -q -e "$LINK_CMDS"
    echo "Compilation succesful."

    root -l -b -q "${MAIN_MACRO}+(${ARGS})"
elif [ "$LINK" = false ] && [ -n "$MAIN_MACRO" ]; then
    echo "Running $MAIN_MACRO with args: ($ARGS)"

    ALL_CMDS="gROOT->ProcessLine(\".x ${MAIN_MACRO}+(${ARGS}\");"

    root -l -b -q "${MAIN_MACRO}+(${ARGS})"
fi