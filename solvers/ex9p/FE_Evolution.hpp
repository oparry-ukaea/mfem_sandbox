#ifndef EX9P_FE_EVOLUTION_HPP
#define EX9P_FE_EVOLUTION_HPP

#include <mfem.hpp>

#include "DG_Solver.hpp"

namespace ex9p {

/** A time-dependent operator for the right-hand side of the ODE. The DG weak
    form of du/dt = -v.grad(u) is M du/dt = K u + b, where M and K are the mass
    and advection matrices, and b describes the flow on the boundary. This can
    be written as a general ODE, du/dt = M^{-1} (K u + b), and this class is
    used to evaluate the right-hand side. */
class FE_Evolution : public mfem::TimeDependentOperator {
private:
  mfem::OperatorHandle M, K;
  const mfem::Vector& b;
  mfem::Solver* M_prec;
  mfem::CGSolver M_solver;
  DG_Solver* dg_solver;

  mutable mfem::Vector z;

public:
  FE_Evolution(mfem::ParBilinearForm& M_, mfem::ParBilinearForm& K_,
               const mfem::Vector& b_, PrecType prec_type);

  void Mult(const mfem::Vector& x, mfem::Vector& y) const override;
  void ImplicitSolve(const real_t dt, const mfem::Vector& x, mfem::Vector& k) override;

  ~FE_Evolution() override;
};

// Implementation of class FE_Evolution
FE_Evolution::FE_Evolution(mfem::ParBilinearForm& M_, mfem::ParBilinearForm& K_,
                           const mfem::Vector& b_, PrecType prec_type)
    : mfem::TimeDependentOperator(M_.ParFESpace()->GetTrueVSize()), b(b_),
      M_solver(M_.ParFESpace()->GetComm()), z(height) {
  if (M_.GetAssemblyLevel() == mfem::AssemblyLevel::LEGACY) {
    M.Reset(M_.ParallelAssemble(), true);
    K.Reset(K_.ParallelAssemble(), true);
  } else {
    M.Reset(&M_, false);
    K.Reset(&K_, false);
  }

  M_solver.SetOperator(*M);

  mfem::Array<int> ess_tdof_list;
  if (M_.GetAssemblyLevel() == mfem::AssemblyLevel::LEGACY) {
    mfem::HypreParMatrix& M_mat = *M.As<mfem::HypreParMatrix>();
    mfem::HypreParMatrix& K_mat = *K.As<mfem::HypreParMatrix>();
    mfem::HypreSmoother* hypre_prec =
        new mfem::HypreSmoother(M_mat, mfem::HypreSmoother::Jacobi);
    M_prec = hypre_prec;

    dg_solver = new DG_Solver(M_mat, K_mat, *M_.FESpace(), prec_type);
  } else {
    M_prec = new mfem::OperatorJacobiSmoother(M_, ess_tdof_list);
    dg_solver = NULL;
  }

  M_solver.SetPreconditioner(*M_prec);
  M_solver.iterative_mode = false;
  M_solver.SetRelTol(1e-9);
  M_solver.SetAbsTol(0.0);
  M_solver.SetMaxIter(100);
  M_solver.SetPrintLevel(0);
}

// Solve the equation:
//    u_t = M^{-1}(Ku + b),
// by solving associated linear system
//    (M - dt*K) d = K*u + b
void FE_Evolution::ImplicitSolve(const real_t dt, const mfem::Vector& x,
                                 mfem::Vector& k) {
  K->Mult(x, z);
  z += b;
  dg_solver->SetTimeStep(dt);
  dg_solver->Mult(z, k);
}

void FE_Evolution::Mult(const mfem::Vector& x, mfem::Vector& y) const {
  // y = M^{-1} (K x + b)
  K->Mult(x, z);
  z += b;
  M_solver.Mult(z, y);
}

FE_Evolution::~FE_Evolution() {
  delete M_prec;
  delete dg_solver;
}

} // namespace ex9p

#endif // EX9P_FE_EVOLUTION_HPP