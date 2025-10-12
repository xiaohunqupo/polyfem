#pragma once

#include "AugmentedLagrangianForm.hpp"

namespace polyfem::solver
{
	/// @brief Form of the lagrangian in augmented lagrangian
	class MatrixLagrangianForm : public AugmentedLagrangianForm
	{
	public:
		/// @brief Construct a new MatrixLagrangianForm object for the constraints Ax = b, where A is sparse
		/// @param A Constraints matrix
		/// @param b Constraints value
		/// @param A_proj Projection matrix
		/// @param b_proj Projection value
		MatrixLagrangianForm(const StiffnessMatrix &A,
							 const Eigen::MatrixXd &b,
							 const StiffnessMatrix &A_proj = {},
							 const Eigen::MatrixXd &b_pro = {});

		std::string name() const override { return "generic-lagrangian"; }

		/// @brief Compute the value of the form
		/// @param x Current solution
		/// @return Computed value
		double value_unweighted(const Eigen::VectorXd &x) const override;

		/// @brief Compute the first derivative of the value wrt x
		/// @param[in] x Current solution
		/// @param[out] gradv Output gradient of the value wrt x
		void first_derivative_unweighted(const Eigen::VectorXd &x, Eigen::VectorXd &gradv) const override;

		/// @brief Compute the second derivative of the value wrt x
		/// @param[in] x Current solution
		/// @param[out] hessian Output Hessian of the value wrt x
		void second_derivative_unweighted(const Eigen::VectorXd &x, StiffnessMatrix &hessian) const override;

	public:
		void update_lagrangian(const Eigen::VectorXd &x, const double k_al) override;
		double compute_error(const Eigen::VectorXd &x) const override;

	private:
		StiffnessMatrix AtA;
		Eigen::VectorXd Atb;
	};
} // namespace polyfem::solver
