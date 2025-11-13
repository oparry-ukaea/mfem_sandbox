#ifndef EX9P_DG_SOLVER_HPP
#define EX9P_DG_SOLVER_HPP

#include <mfem.hpp>

#include "AIR_Prec.hpp"

namespace ex9p {

// Type of preconditioner for implicit time integrator
enum class PrecType : int { ILU = 0, AIR = 1 };

using real_t = mfem::real_t;

class DG_Solver : public mfem::Solver {
private:
  mfem::HypreParMatrix &M, &K;
  mfem::SparseMatrix M_diag;
  mfem::HypreParMatrix* A;
  mfem::GMRESSolver linear_solver;
  mfem::Solver* prec;
  real_t dt;

public:
  DG_Solver(mfem::HypreParMatrix& M_, mfem::HypreParMatrix& K_,
            const mfem::FiniteElementSpace& fes, PrecType prec_type)
      : M(M_), K(K_), A(NULL), linear_solver(M.GetComm()), dt(-1.0) {
    int block_size = fes.GetTypicalFE()->GetDof();
    if (prec_type == PrecType::ILU) {
      prec = new mfem::BlockILU(block_size,
                                mfem::BlockILU::Reordering::MINIMUM_DISCARDED_FILL);
    } else if (prec_type == PrecType::AIR) {
#if MFEM_HYPRE_VERSION >= 21800
      prec = new AIR_prec(block_size);
#else
      MFEM_ABORT("Must have MFEM_HYPRE_VERSION >= 21800 to use AIR.\n");
#endif
    }
    linear_solver.iterative_mode = false;
    linear_solver.SetRelTol(1e-9);
    linear_solver.SetAbsTol(0.0);
    linear_solver.SetMaxIter(100);
    linear_solver.SetPrintLevel(0);
    linear_solver.SetPreconditioner(*prec);

    M.GetDiag(M_diag);
  }

  void SetTimeStep(real_t dt_) {
    if (dt_ != dt) {
      dt = dt_;
      // Form operator A = M - dt*K
      delete A;
      A = Add(-dt, K, 0.0, K);
      mfem::SparseMatrix A_diag;
      A->GetDiag(A_diag);
      A_diag.Add(1.0, M_diag);
      // this will also call SetOperator on the preconditioner
      linear_solver.SetOperator(*A);
    }
  }

  void SetOperator(const Operator& op) override { linear_solver.SetOperator(op); }

  void Mult(const mfem::Vector& x, mfem::Vector& y) const override {
    linear_solver.Mult(x, y);
  }

  ~DG_Solver() override {
    delete prec;
    delete A;
  }
};

} // namespace ex9p

#endif // EX9P_DG_SOLVER_HPP