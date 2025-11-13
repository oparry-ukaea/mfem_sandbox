#include "ex9p.hpp"

#include <fstream>
#include <iostream>
#include <mfem.hpp>

#include "DG_Solver.hpp"
#include "FE_Evolution.hpp"

namespace ex9p {
Sim::Sim(int argc, char* argv[]) : args(argc, argv) {
  // Initialize MPI and HYPRE.
  mfem::Mpi::Init();
  this->nranks = mfem::Mpi::WorldSize();
  this->rank = mfem::Mpi::WorldRank();
  this->is_root = mfem::Mpi::Root();
  mfem::Hypre::Init();

  // Set stdout precision
  this->precision = 8;
  std::cout.precision(this->precision);
}

/**
 * @brief Define the parallel DG finite element space on the parallel refined mesh with
 * the given polynomial order.
 *
 */
void Sim::create_fe_space() {
  this->fe_collection =
      new mfem::DG_FECollection(order, this->mesh_dim, mfem::BasisType::GaussLobatto);
  this->fe_space = new mfem::ParFiniteElementSpace(this->pmesh, this->fe_collection);

  HYPRE_BigInt global_vSize = fe_space->GlobalTrueVSize();
  if (this->is_root) {
    std::cout << "Number of unknowns: " << global_vSize << std::endl;
  }
}

void Sim::do_assembly() {
  int skip_zeros = 0;
  this->mass_matrix->Assemble();
  this->adv_matrix->Assemble(skip_zeros);
  this->bdy_flow->Assemble();
  this->mass_matrix->Finalize();
  this->adv_matrix->Finalize(skip_zeros);
  this->bdy_flow_par = this->bdy_flow->ParallelAssemble();
}

void Sim::free() {
  delete this->adv_matrix;
  delete this->bdy_flow;
  delete this->bdy_flow_par;
  delete this->fe_space;
  delete this->fld_u;
  delete this->mass_matrix;
  delete this->par_fld_u;
  delete this->paraview_handler;
  delete this->pmesh;
}

void Sim::init_paraview_output() {
  this->paraview_handler = new mfem::ParaViewDataCollection("ex9p", this->pmesh);
  this->paraview_handler->SetPrefixPath("ParaView");
  this->paraview_handler->SetLevelsOfDetail(this->order);
  this->paraview_handler->SetDataFormat(mfem::VTKFormat::BINARY);
  this->paraview_handler->SetHighOrderOutput(true);
}

void Sim::output_to_paraview(const int& step, const real_t time) {
  this->paraview_handler->SetCycle(step);
  this->paraview_handler->SetTime(time);
  this->paraview_handler->Save();
}

int Sim::parse_cl_args() {
  // Command-line option defaults
  this->device_config = "cpu";
  this->do_element_assembly = false;
  this->do_full_assembly = false;
  this->do_partial_assembly = false;
  this->dt = 0.005;
  this->mesh_file = "./mfem/data/periodic-square.mesh";
  this->num_mesh_refinements_parallel = 0;
  this->num_mesh_refinements_serial = 3;
  this->ode_solver_type = 4;
  this->order = 3;
  this->paraview_output_enabled = true;
#if MFEM_HYPRE_VERSION >= 21800
  this->prec_type = PrecType::AIR;
#else
  this->prec_type = PrecType::ILU;
#endif
  this->t_final = 10.0;
  this->vis_step_freq = 5;

  // static, unfortunately...
  init_blob_scale = 0.1;

  // Set up options
  this->args.AddOption(&this->mesh_file, "-m", "--mesh", "Mesh file to use.");
  this->args.AddOption(&this->num_mesh_refinements_serial, "-rs", "--refine-serial",
                       "Number of times to refine the mesh uniformly in serial.");
  this->args.AddOption(&this->num_mesh_refinements_parallel, "-rp", "--refine-parallel",
                       "Number of times to refine the mesh uniformly in parallel.");
  this->args.AddOption(&this->order, "-o", "--order",
                       "Order (degree) of the finite elements.");
  this->args.AddOption(&this->do_partial_assembly, "-pa", "--partial-assembly", "-no-pa",
                       "--no-partial-assembly", "Enable Partial Assembly.");
  this->args.AddOption(&this->do_element_assembly, "-ea", "--element-assembly", "-no-ea",
                       "--no-element-assembly", "Enable Element Assembly.");
  this->args.AddOption(&this->do_full_assembly, "-fa", "--full-assembly", "-no-fa",
                       "--no-full-assembly", "Enable Full Assembly.");
  this->args.AddOption(&this->device_config, "-d", "--device",
                       "Device configuration string, see Device::Configure().");
  this->args.AddOption(&this->ode_solver_type, "-s", "--ode-solver",
                       mfem::ODESolver::Types.c_str());
  this->args.AddOption(&this->t_final, "-tf", "--t-final",
                       "Final time; start time is 0.");
  this->args.AddOption(&this->dt, "-dt", "--time-step", "Time step.");
  this->args.AddOption((int*)&prec_type, "-pt", "--prec-type",
                       "Preconditioner for "
                       "implicit solves. 0 for ILU, 1 for pAIR-AMG.");
  this->args.AddOption(&this->paraview_output_enabled, "-paraview",
                       "--paraview-datafiles", "-no-paraview", "--no-paraview-datafiles",
                       "Save data files for ParaView (paraview.org) visualization.");
  this->args.AddOption(&this->vis_step_freq, "-vs", "--visualization-steps",
                       "Visualize every n-th timestep.");
  this->args.AddOption(&init_blob_scale, "-i", "--init_blob_scale",
                       "Exponential scale length of the initial blob.");

  // Parse args and return arg status. Print options on success.
  this->args.Parse();
  if (this->args.Good()) {
    if (this->is_root) {
      this->args.PrintOptions(std::cout);
    }
    return 0;
  } else {
    if (this->is_root) {
      this->args.PrintUsage(std::cout);
    }
    return 1;
  }
}

void Sim::read_and_refine_mesh() {

  // Read the serial mesh from the given mesh file on all processors. We can
  // handle geometrically periodic meshes in this code.
  mfem::Mesh* mesh = new mfem::Mesh(this->mesh_file, 1, 1);
  this->mesh_dim = mesh->Dimension();

  /* Refine the mesh in serial to increase the resolution. If the mesh is of NURBS type,
  we convert it to a (piecewise-polynomial) high-order mesh. */
  for (int lev = 0; lev < num_mesh_refinements_serial; lev++) {
    mesh->UniformRefinement();
  }
  if (mesh->NURBSext) {
    mesh->SetCurvature(std::max(order, 1));
  }
  mesh->GetBoundingBox(Sim::domain_min, Sim::domain_max, std::max(order, 1));

  /* Define the parallel mesh by a partitioning of the serial mesh. Refine this mesh
   * further in parallel to increase the resolution. Once the parallel mesh is defined,
   * the serial mesh is deleted.
   */
  this->pmesh = new mfem::ParMesh(MPI_COMM_WORLD, *mesh);
  delete mesh;
  for (int lev = 0; lev < num_mesh_refinements_parallel; lev++) {
    this->pmesh->UniformRefinement();
  }
}

int Sim::run() {

  // Parse command line arguments
  int args_state = parse_cl_args();
  if (args_state) {
    return args_state;
  }

  // Configure device and report config
  mfem::Device device(device_config);
  if (this->is_root) {
    device.Print();
  }

  // Read the mesh and refine according to CL options
  read_and_refine_mesh();

  // Define the ODE solver used for time integration.
  std::unique_ptr<mfem::ODESolver> ode_solver =
      mfem::ODESolver::Select(this->ode_solver_type);

  // Set up the FE space
  create_fe_space();

  // Set up linear, bilinear forms; add integrators
  setup_forms_and_integrators();

  // Do assembly
  do_assembly();

  // Set up a grid function for field 'u' and set the ICs
  this->fld_u = new mfem::ParGridFunction(this->fe_space);
  mfem::FunctionCoefficient u_init(u0_function);
  this->fld_u->ProjectCoefficient(u_init);
  // Set up parallel (HYPRE) vector for field 'u'
  this->par_fld_u = this->fld_u->GetTrueDofs();

  // Set up paraview handler and write the first output
  if (this->paraview_output_enabled) {
    init_paraview_output();
    this->paraview_handler->RegisterField("u", this->fld_u);
    output_to_paraview(0, 0.0);
  }

  // Time variable
  real_t t = 0.0;

  // Set up the time-dependent evolution operator describing the ODE RHS
  FE_Evolution adv(*(this->mass_matrix), *(this->adv_matrix), *(this->bdy_flow_par),
                   prec_type);
  adv.SetTime(t);
  ode_solver->Init(adv);

  // Main loop
  bool done = false;
  for (int step = 0; !done;) {
    // Set timetep for the ODE solver, ensuring no overstep at the end
    real_t dt_real = std::min(dt, t_final - t);

    // Do step
    ode_solver->Step(*(this->par_fld_u), t, dt_real);
    step++;

    // Stop at t_final
    done = (t >= t_final - 1e-8 * dt);

    // Write data on output steps, or on last step
    if (done || step % vis_step_freq == 0) {
      if (this->is_root) {
        std::cout << "time step: " << step << ", time: " << t << std::endl;
      }

      // Extract the parallel grid function corresponding to the finite
      // element approximation U (the local solution on each processor).
      *(this->fld_u) = *(this->par_fld_u);

      if (this->paraview_output_enabled) {
        output_to_paraview(step, t);
      }
    }
  }

  // Free all allocated memory
  free();

  return 0;
}

void Sim::setup_forms_and_integrators() { /* Set up and assemble the parallel bilinear and
     linear forms (and the parallel hypre matrices) corresponding to the DG
     discretization. The DGTraceIntegrator involves integrals over mesh interior faces.
  */
  this->velocity = new mfem::VectorFunctionCoefficient(this->mesh_dim, velocity_function);
  this->inflow = new mfem::FunctionCoefficient(inflow_function);

  this->mass_matrix = new mfem::ParBilinearForm(this->fe_space);
  this->adv_matrix = new mfem::ParBilinearForm(this->fe_space);
  if (this->do_partial_assembly) {
    this->mass_matrix->SetAssemblyLevel(mfem::AssemblyLevel::PARTIAL);
    this->adv_matrix->SetAssemblyLevel(mfem::AssemblyLevel::PARTIAL);
  } else if (this->do_element_assembly) {
    this->mass_matrix->SetAssemblyLevel(mfem::AssemblyLevel::ELEMENT);
    this->adv_matrix->SetAssemblyLevel(mfem::AssemblyLevel::ELEMENT);
  } else if (this->do_full_assembly) {
    this->mass_matrix->SetAssemblyLevel(mfem::AssemblyLevel::FULL);
    this->adv_matrix->SetAssemblyLevel(mfem::AssemblyLevel::FULL);
  }

  // Add integrators
  this->mass_matrix->AddDomainIntegrator(new mfem::MassIntegrator);

  constexpr real_t adv_direction = -1.0;
  this->adv_matrix->AddDomainIntegrator(
      new mfem::ConvectionIntegrator(*(this->velocity), adv_direction));
  this->adv_matrix->AddInteriorFaceIntegrator(
      new mfem::NonconservativeDGTraceIntegrator(*(this->velocity), adv_direction));
  this->adv_matrix->AddBdrFaceIntegrator(
      new mfem::NonconservativeDGTraceIntegrator(*(this->velocity), adv_direction));

  this->bdy_flow = new mfem::ParLinearForm(this->fe_space);
  this->bdy_flow->AddBdrFaceIntegrator(new mfem::BoundaryFlowIntegrator(
      *(this->inflow), *(this->velocity), adv_direction));
}

mfem::Vector Sim::domain_min, Sim::domain_max;
real_t Sim::init_blob_scale;
} // namespace ex9p

/**
 * @brief Run the simulation
 */
int main(int argc, char* argv[]) {
  ex9p::Sim sim = ex9p::Sim(argc, argv);
  return sim.run();
}