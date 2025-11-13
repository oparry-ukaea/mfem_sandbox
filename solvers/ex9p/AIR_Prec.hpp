#ifndef EX9P_AIR_PREC_HPP
#define EX9P_AIR_PREC_HPP

#include <mfem.hpp>

namespace ex9p {

#if MFEM_HYPRE_VERSION >= 21800
// Algebraic multigrid preconditioner for advective problems based on
// approximate ideal restriction (AIR). Most effective when matrix is
// first scaled by DG block inverse, and AIR applied to scaled matrix.
// See https://doi.org/10.1137/17M1144350.
class AIR_prec : public mfem::Solver {
private:
  const mfem::HypreParMatrix* A;
  // Copy of A scaled by block-diagonal inverse
  mfem::HypreParMatrix A_s;

  mfem::HypreBoomerAMG* AIR_solver;
  int blocksize;

public:
  AIR_prec(int blocksize_) : AIR_solver(NULL), blocksize(blocksize_) {}

  void SetOperator(const mfem::Operator& op) override {
    this->width = op.Width();
    this->height = op.Height();

    A = dynamic_cast<const mfem::HypreParMatrix*>(&op);
    MFEM_VERIFY(A != NULL, "AIR_prec requires a mfem::HypreParMatrix.")

    // Scale A by block-diagonal inverse
    BlockInverseScale(A, &A_s, NULL, NULL, blocksize,
                      mfem::BlockInverseScaleJob::MATRIX_ONLY);
    delete AIR_solver;
    AIR_solver = new mfem::HypreBoomerAMG(A_s);
    AIR_solver->SetAdvectiveOptions(1, "", "FA");
    AIR_solver->SetPrintLevel(0);
    AIR_solver->SetMaxLevels(50);
  }

  void Mult(const mfem::Vector& x, mfem::Vector& y) const override {
    // Scale the rhs by block inverse and solve system
    mfem::HypreParVector z_s;
    BlockInverseScale(A, NULL, &x, &z_s, blocksize, mfem::BlockInverseScaleJob::RHS_ONLY);
    AIR_solver->Mult(z_s, y);
  }

  ~AIR_prec() override { delete AIR_solver; }
};
#endif // MFEM_HYPRE_VERSION >= 21800

} // namespace ex9p

#endif // EX9P_AIR_PREC_HPP