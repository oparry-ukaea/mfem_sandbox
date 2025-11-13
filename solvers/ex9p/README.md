## ex9p refactored

The ex9p example from mfem, refactored and marked up for clarity.

# Run

    ./scripts/run_solver.sh ex9p

... runs with the default configuration, which is:

|                                | Default             |
| ------------------------------ | ------------------- |
| Mesh                           | 3x3 periodic square |
| Raw mesh dimension             | 3x3                 |
| Num. (serial) mesh refinements | 3                   |
| Mesh dimension                 | 24x24               |
| Element order                  | 3                   |
| Timestep                       | 0.005               |
| Solver ID                      | 4 (RK4)             |

To run with a larger timestep and an implicit (Second-order SDIRK) solver:

    ./scripts/run_solver.sh ex9p -a "-dt 0.01 -s 22"


