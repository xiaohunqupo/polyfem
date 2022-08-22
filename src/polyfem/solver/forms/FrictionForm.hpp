#pragma once

#include "Form.hpp"

#include <polyfem/State.hpp>

#include <polyfem/utils/Types.hpp>

#include <ipc/ipc.hpp>
#include <ipc/collision_mesh.hpp>
#include <ipc/friction/friction_constraint.hpp>

namespace polyfem::solver
{
	class ContactForm;

	/// @brief Form of the lagged friction disapative potential and forces
	class FrictionForm : public Form
	{
	public:
		/// @brief Construct a new Friction Form object
		/// @param state Reference to the simulation state
		/// @param epsv Smoothing factor between static and dynamic friction
		/// @param mu Global coefficient of friction
		/// @param dhat Barrier activation distance
		/// @param barrier_stiffness Barrier stiffness multiplier
		/// @param broad_phase_method Broad-phase method used for distance computation and collision detection
		/// @param dt Time step size
		/// @param contact_form Pointer to contact form; necessary to have the barrier stiffnes, maybe clean me
		FrictionForm(
			const State &state,
			const double epsv,
			const double mu,
			const double dhat,
			const ipc::BroadPhaseMethod broad_phase_method,
			const double dt,
			const ContactForm &contact_form);

		/// @brief Compute the value of the form
		/// @param x Current solution
		/// @return Computed value
		double value(const Eigen::VectorXd &x) const override;

		/// @brief Compute the first derivative of the value wrt x
		/// @param[in] x Current solution
		/// @param[out] gradv Output gradient of the value wrt x
		void first_derivative(const Eigen::VectorXd &x, Eigen::VectorXd &gradv) const override;

		/// @brief Compute the second derivative of the value wrt x
		/// @param[in] x Current solution
		/// @param[out] hessian Output Hessian of the value wrt x
		void second_derivative(const Eigen::VectorXd &x, StiffnessMatrix &hessian) override;

		/// @brief Initialize lagged fields
		/// @param x Current solution
		void init_lagging(const Eigen::VectorXd &x) override;

		/// @brief Update lagged fields
		/// @param x Current solution
		void update_lagging(const Eigen::VectorXd &x) override;

	private:
		const State &state_; ///< Reference to the simulation state

		const double epsv_;                              ///< Smoothing factor between static and dynamic friction
		const double mu_;                                ///< Global coefficient of friction
		const double dt_;                                ///< Time step size
		const double dhat_;                              ///< Barrier activation distance
		const ipc::BroadPhaseMethod broad_phase_method_; ///< Broad-phase method used for distance computation and collision detection

		ipc::FrictionConstraints friction_constraint_set_; ///< Lagged friction constraint set
		Eigen::MatrixXd displaced_surface_prev_;           ///< Displaced vertices at the start of the time-step.

		/// @brief Compute the displaced positions of the surface nodes
		Eigen::MatrixXd compute_displaced_surface(const Eigen::VectorXd &x) const;

		const ContactForm &contact_form_; ///> necessary to have the barrier stiffnes, maybe clean me
	};
} // namespace polyfem::solver
