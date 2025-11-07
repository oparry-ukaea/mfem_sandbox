# mfem_sandbox
Sandbox repo to experiment with mfem.

## Installation

To install sandbox library and solvers via spack:

    git submodule update --init
    . activate_mfemsb
    spack install -j8

(Library and solvers are built in `./builds/build-*`)

## Run

Run solvers (from the repo root) with:

    ./scripts/run_solver.sh [solver_name]

By default, this will run a solver executable from the latest spack build. To choose a different build directory pass `-b [build_path]`.
Run the script without any arguments to see all of the available options.

    ./scripts/run_solver.sh --help
