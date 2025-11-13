#ifndef EX9P_EX9P_HPP
#define EX9P_EX9P_HPP

#include <mfem.hpp>

#include "DG_Solver.hpp"

namespace ex9p {
class Sim {

  static mfem::Vector domain_min, domain_max; // bounding box
  static real_t init_blob_scale;
  // Inflow boundary condition (always zero for this e.g.)
  static real_t inflow_function(const mfem::Vector& x) { return 0.0; }

  static real_t u0_function(const mfem::Vector& x) {
    int dim = x.Size();

    // map to the reference [-1,1] domain
    mfem::Vector X(dim);
    for (int i = 0; i < dim; i++) {
      real_t center = (Sim::domain_min[i] + Sim::domain_max[i]) * 0.5;
      X(i) = 2 * (x(i) - center) / (Sim::domain_max[i] - Sim::domain_min[i]);
    }

    switch (dim) {
    case 1:
      return exp((-X(0) * X(0)) / (init_blob_scale * init_blob_scale));
      break;
    case 2:
      return exp((-X(0) * X(0) - X(1) * X(1)) / (init_blob_scale * init_blob_scale));
      break;
    case 3: {
      return exp((-X(0) * X(0) - X(1) * X(1) - X(2) * X(2))
                 / (init_blob_scale * init_blob_scale));
      break;
    }
    }
    return 0.0;
  }

  static void velocity_function(const mfem::Vector& x, mfem::Vector& v) {
    int dim = x.Size();

    // Map to the reference [-1,1] domain using the bounding box
    mfem::Vector X(dim);
    for (int i = 0; i < dim; i++) {
      real_t center = (Sim::domain_min[i] + Sim::domain_max[i]) * 0.5;
      X(i) = 2 * (x(i) - center) / (Sim::domain_max[i] - Sim::domain_min[i]);
    }

    // Translations in 1D, 2D, and 3D
    switch (dim) {
    case 1:
      v(0) = 1.0;
      break;
    case 2:
      v(0) = sqrt(2. / 3.);
      v(1) = sqrt(1. / 3.);
      break;
    case 3:
      v(0) = sqrt(3. / 6.);
      v(1) = sqrt(2. / 6.);
      v(2) = sqrt(1. / 6.);
      break;
    }
  }

public:
  Sim(int argc, char* argv[]);
  int run();

private:
  mfem::OptionsParser args;

  // Command-line options
  const char* device_config;
  bool do_element_assembly;
  bool do_full_assembly;
  bool do_partial_assembly;
  real_t dt;
  const char* mesh_file = "../data/periodic-hexagon.mesh";
  int num_mesh_refinements_parallel;
  int num_mesh_refinements_serial;
  int ode_solver_type;
  int order;
  bool paraview_output_enabled;
  PrecType prec_type;
  real_t t_final;
  int vis_step_freq;

  // stdout precision
  int precision;

  // MPI
  bool is_root;
  int nranks;
  int rank;

  // Mesh
  mfem::ParMesh* pmesh = nullptr;
  int mesh_dim;

  // FE space
  mfem::DG_FECollection* fe_collection = nullptr;
  mfem::ParFiniteElementSpace* fe_space = nullptr;

  // Function coeffs
  mfem::FunctionCoefficient* inflow;
  mfem::VectorFunctionCoefficient* velocity;

  // Bilinear and linear forms
  mfem::ParBilinearForm* adv_matrix = nullptr;
  mfem::ParLinearForm* bdy_flow = nullptr;
  mfem::HypreParVector* bdy_flow_par = nullptr;
  mfem::ParBilinearForm* mass_matrix = nullptr;

  // Fields / variables
  mfem::ParGridFunction* fld_u = nullptr;
  mfem::HypreParVector* par_fld_u = nullptr;

  // Output data collection
  mfem::ParaViewDataCollection* paraview_handler = nullptr;

  void create_fe_space();
  void do_assembly();
  void free();
  void init_paraview_output();
  void output_to_paraview(const int& step, const real_t time);
  int parse_cl_args();
  void read_and_refine_mesh();
  void setup_forms_and_integrators();
};

} // namespace ex9p

#endif // EX9P_EX9P_HPP