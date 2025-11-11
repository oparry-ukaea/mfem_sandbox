#!/bin/bash

#--------------------------------------------------------------------------------------------------
# Helper functions
echo_usage() {
    echo "Usage:"
    echo "    $0 [solver_name] <-a solver_args> <-b build_dir/install_dir> <-m mesh_file> <-n num_MPI>"
}

execute() {
    local run_cmd=$1
    local run_dir=$2
    echo "----------------------------------------------------------------------------------------------------"
    echo "Executing [$run_cmd] in [$run_dir]"
    cd "$run_dir" || exit 
    eval "$run_cmd"
    cd - || exit
    echo "----------------------------------------------------------------------------------------------------"
}

generate_run_dir() {
    local solver_name="$1"
    
    if [ -z "$ARCHMFEM" ]; then
        run_dir_root="$REPO_ROOT/runs"
    else
        run_dir_root="$ARCHMFEM"
    fi

    run_dir="$run_dir_root/$solver_name"
    if [ -e "$run_dir" ]; then
        read -p "Overwrite existing run directory at $run_dir? (Y/N): " choice && [[ $choice == [yY] || $choice == [yY][eE][sS] ]] || exit 4
        \rm -rf "$run_dir"
    fi
    mkdir -p "$run_dir"
}


find_mesh() {
    local mesh_name=$1


    # If no mesh name was provided, read stored default
    if [ -z "$mesh_name" ]; then
        default_mesh_names_file="$REPO_ROOT/meshes/defaults.txt"
        mesh_name=$(grep "^$solver_name" "$default_mesh_names_file" | cut -d" " -f 2)
    fi

    # If no file extension was provided, append .mesh
    ext="${mesh_name##*.}"
    if [ "$ext" == "$mesh_name" ]; then
        ext="mesh"
    fi
    mesh_name="${mesh_name%.$ext}.$ext"
    
    mfem_mesh_path="$REPO_ROOT/mfem/data/$mesh_name"
    local_mesh_path="$REPO_ROOT/meshes/$mesh_name"

    if [ -f "$local_mesh_path" ]; then
        mesh_path="$local_mesh_path"
    elif [ -f "$mfem_mesh_path" ]; then
        mesh_path="$mfem_mesh_path"
    else
        echo "No mesh file found at either [$mfem_mesh_path] or [$local_mesh_path]"
        exit 6
    fi
}

parse_args() {
    if [ $# -lt 1 ]; then
        echo_usage
        exit 1
    fi
    POSITIONAL_ARGS=()
    while [[ $# -gt 0 ]]; do
    case $1 in
        -a|--args)
        solver_args="$2"
        shift 2
        ;;
        -b|--build-dir)
        exec_loc=$(realpath "$2")
        shift 2
        ;;
        -h|--help)
        echo_usage
        exit 0
        ;;
        -m|--mesh)
        mesh_name="$2"
        shift 2
        ;;
        -n|--num_mpi)
        nmpi="$2"
        shift 2
        ;;
        -*)
        echo "Unknown option $1"
        exit 2
        ;;
        *)
        # Save positional args in an array
        POSITIONAL_ARGS+=("$1")
        shift
        ;;
    esac
    done

    # Restore and extract positional args
    set -- "${POSITIONAL_ARGS[@]}"
    solver_name=$1
}

report_options() {
    echo "Options:"
    echo "    Solver : $solver_name"
    echo "     n MPI : $nmpi"
    echo ""
}

set_default_exec_loc() {
    solver_exec=$(find -L "$REPO_ROOT/views" -maxdepth 3 -name "$solver_name" -printf "%TY-%Tm-%Td %TT %p\n" | sort -n|tail -1|cut -d " " -f 3)
}

set_paths_and_validate() {
    local exec_loc=$1
   
    if [ -z "$exec_loc" ]; then
        # If no executable location was specified, set solver_exec to the most
        # recently modified executable in the views
        set_default_exec_loc
        if [ -z "$solver_exec" ]; then
            echo "No installed solver found in ./views; run spack install or pass -b <build_directory>"
            exit 3
        fi
    else
        # Else check build and install locations relative to the specified $exec_loc
        solver_build_exec="$exec_loc/$solver_name"
        solver_install_exec="$exec_loc/bin/$solver_name"
        if [ -f "$solver_build_exec" ]; then
            solver_exec="$solver_build_exec"
        elif [ -f "$solver_install_exec" ]; then
            solver_exec="$solver_install_exec"
        else
            echo "No solver found at [$solver_build_exec] or [$solver_install_exec]"
            exit 3
        fi
    fi
}

REPO_ROOT=$( cd -- "$(realpath $( dirname -- "${BASH_SOURCE[0]}" )/..)" &> /dev/null && pwd )

# Default options
solver_name='Not set'
mesh_name=''
mesh_path='Not set'
nmpi='4'
exec_loc=''
solver_args=''

# Parse command line args and report resulting options
parse_args "$@"
report_options

solver_exec='Not set'
# Find the executable inside $exec_loc and validate the examples paths
set_paths_and_validate "$exec_loc"

# Set up run directory, confirming overwrite if it already exists
run_dir='Not set'
generate_run_dir "$solver_name"

# Find mesh file
find_mesh "$mesh_name"

# Execute run_cmd in run_dir
run_cmd="mpirun -np $nmpi $solver_exec $solver_args -m $mesh_path"
execute "$run_cmd" "$run_dir"
echo "See $run_dir for output"
