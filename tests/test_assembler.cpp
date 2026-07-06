#include <polyfem/State.hpp>

#include <polyfem/assembler/NeoHookeanElasticity.hpp>
#include <polyfem/assembler/NeoHookeanElasticityAutodiff.hpp>
#include <polyfem/assembler/ModifiedNeoHookeanElasticity.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <iostream>

using namespace polyfem;
using namespace polyfem::assembler;
using namespace polyfem::basis;
using namespace polyfem::mesh;
using namespace polyfem::utils;

TEST_CASE("hessian_lin", "[assembler]")
{
	const std::string path = POLYFEM_DATA_DIR;
	json in_args = json({});
	in_args["geometry"] = {};
	in_args["geometry"]["mesh"] = path + "/plane_hole.obj";
	in_args["geometry"]["surface_selection"] = 7;
	// in_args["geometry"]["mesh"] = path + "/circle2.msh";
	// in_args["force_linear_geometry"] = true;

	in_args["preset_problem"] = {};
	in_args["preset_problem"]["type"] = "ElasticExact";

	in_args["materials"] = {};
	in_args["materials"]["type"] = "LinearElasticity";
	in_args["materials"]["E"] = 1e5;
	in_args["materials"]["nu"] = 0.3;

	State state;
	state.init_logger("", spdlog::level::err, spdlog::level::off, false);
	state.init(in_args, true);
	state.load_mesh();

	// state.compute_mesh_stats();
	state.build_basis();

	state.assemble_mass_mat();

	SparseMatrixCache mat_cache;
	StiffnessMatrix hessian, stiffness;
	Eigen::MatrixXd disp(state.n_bases * 2, 1);
	disp.setZero();

	state.build_stiffness_mat(stiffness);

	for (int rand = 0; rand < 10; ++rand)
	{
		state.assembler->assemble_hessian(false, state.n_bases, false,
										  state.bases, state.bases, state.ass_vals_cache, 0, 0, disp, Eigen::MatrixXd(), mat_cache, hessian);

		const StiffnessMatrix tmp = stiffness - hessian;
		const auto val = Catch::Approx(0).margin(1e-8);

		for (int k = 0; k < tmp.outerSize(); ++k)
		{
			for (StiffnessMatrix::InnerIterator it(tmp, k); it; ++it)
			{
				REQUIRE(it.value() == val);
			}
		}

		disp.setRandom();
	}
}

TEST_CASE("hessian_hooke", "[assembler]")
{
	const std::string path = POLYFEM_DATA_DIR;
	json in_args = json({});
	in_args["geometry"] = {};
	in_args["geometry"]["mesh"] = path + "/plane_hole.obj";
	in_args["geometry"]["surface_selection"] = 7;
	// in_args["geometry"]["mesh"] = path + "/circle2.msh";
	// in_args["force_linear_geometry"] = true;

	in_args["preset_problem"] = {};
	in_args["preset_problem"]["type"] = "ElasticExact";

	in_args["materials"] = {};
	in_args["materials"]["type"] = "HookeLinearElasticity";
	in_args["materials"]["E"] = 1e5;
	in_args["materials"]["nu"] = 0.3;

	State state;
	state.init_logger("", spdlog::level::err, spdlog::level::off, false);
	state.init(in_args, true);
	state.load_mesh();

	// state.compute_mesh_stats();
	state.build_basis();

	state.assemble_mass_mat();

	SparseMatrixCache mat_cache;
	StiffnessMatrix hessian, stiffness;
	Eigen::MatrixXd disp(state.n_bases * 2, 1);
	disp.setZero();

	state.build_stiffness_mat(stiffness);

	for (int rand = 0; rand < 10; ++rand)
	{
		state.assembler->assemble_hessian(false, state.n_bases, false,
										  state.bases, state.bases, state.ass_vals_cache, 0, 0, disp, Eigen::MatrixXd(), mat_cache, hessian);

		const StiffnessMatrix tmp = stiffness - hessian;
		const auto val = Catch::Approx(0).margin(1e-8);

		for (int k = 0; k < tmp.outerSize(); ++k)
		{
			for (StiffnessMatrix::InnerIterator it(tmp, k); it; ++it)
			{
				REQUIRE(it.value() == val);
			}
		}

		disp.setRandom();
	}
}

TEST_CASE("generic_elastic_assembler", "[assembler]")
{

	const std::string path = POLYFEM_DATA_DIR;
	json in_args = json({});
	in_args["geometry"] = {};
	in_args["geometry"]["mesh"] = path + "/plane_hole.obj";
	in_args["geometry"]["surface_selection"] = 7;
	// in_args["geometry"]["mesh"] = path + "/circle2.msh";
	// in_args["force_linear_geometry"] = true;

	in_args["preset_problem"] = {};
	in_args["preset_problem"]["type"] = "ElasticExact";

	in_args["materials"] = {};
	in_args["materials"]["type"] = "LinearElasticity";
	in_args["materials"]["E"] = 1e5;
	in_args["materials"]["nu"] = 0.3;

	State state;
	state.init_logger("", spdlog::level::err, spdlog::level::off, false);
	state.init(in_args, true);
	state.load_mesh();

	// state.compute_mesh_stats();
	state.build_basis();

	NeoHookeanAutodiff autodiff;
	NeoHookeanElasticity real;

	autodiff.set_size(2);
	real.set_size(2);

	autodiff.add_multimaterial(0, in_args["materials"], state.units, state.root_path());
	real.add_multimaterial(0, in_args["materials"], state.units, state.root_path());

	const int el_id = 0;
	const auto &bs = state.bases[el_id];
	Eigen::MatrixXd local_pts;
	Eigen::MatrixXi f;
	regular_2d_grid(10, true, local_pts, f);

	Eigen::MatrixXd displacement(state.n_bases, 1);

	ElementAssemblyValues vals;
	vals.compute(el_id, false, bs, bs);

	const auto &quadrature = vals.quadrature;
	const QuadratureVector da = vals.det.array() * quadrature.weights.array();

	for (int rand = 0; rand < 10; ++rand)
	{
		displacement.setRandom();

		// value
		{
			const NonLinearAssemblerData data(vals, 0, 0, displacement, displacement, da);

			const double ea = autodiff.compute_energy(data);
			const double e = real.compute_energy(data);

			if (std::isnan(e))
				REQUIRE(std::isnan(ea));
			else
				REQUIRE(ea == Catch::Approx(e).margin(1e-12));
		}

		// grad
		{
			const NonLinearAssemblerData data(vals, 0, 0, displacement, displacement, da);

			const Eigen::VectorXd grada = autodiff.assemble_gradient(data);
			const Eigen::VectorXd grad = real.assemble_gradient(data);

			for (int i = 0; i < grada.size(); ++i)
			{
				if (std::isnan(grad(i)))
					REQUIRE(std::isnan(grada(i)));
				else
					REQUIRE(grada(i) == Catch::Approx(grad(i)).margin(1e-12));
			}
		}

		// hessian
		{
			const NonLinearAssemblerData data(vals, 0, 0, displacement, displacement, da);

			const Eigen::MatrixXd hessiana = autodiff.assemble_hessian(data);
			const Eigen::MatrixXd hessian = real.assemble_hessian(data);

			for (int i = 0; i < hessiana.size(); ++i)
			{
				if (std::isnan(hessian(i)))
					REQUIRE(std::isnan(hessiana(i)));
				else
					REQUIRE(hessiana(i) == Catch::Approx(hessian(i)).margin(1e-12));
			}
		}

		// F stress
		{
			Eigen::MatrixXd stressa, stress;
			autodiff.compute_stress_tensor(OutputData(0, el_id, bs, bs, local_pts, displacement), ElasticityTensorType::F, stressa);
			real.compute_stress_tensor(OutputData(0, el_id, bs, bs, local_pts, displacement), ElasticityTensorType::F, stress);

			for (int i = 0; i < stressa.size(); ++i)
			{
				if (std::isnan(stress(i)))
					REQUIRE(std::isnan(stressa(i)));
				else
					REQUIRE(stressa(i) == Catch::Approx(stress(i)).margin(1e-12));
			}
		}

		// cauchy stress
		{
			Eigen::MatrixXd stressa, stress;
			autodiff.compute_stress_tensor(OutputData(0, el_id, bs, bs, local_pts, displacement), ElasticityTensorType::CAUCHY, stressa);
			real.compute_stress_tensor(OutputData(0, el_id, bs, bs, local_pts, displacement), ElasticityTensorType::CAUCHY, stress);

			for (int i = 0; i < stressa.size(); ++i)
			{
				if (std::isnan(stress(i)))
					REQUIRE(std::isnan(stressa(i)));
				else
					REQUIRE(stressa(i) == Catch::Approx(stress(i)).margin(1e-12));
			}
		}

		// pk1 stress
		{
			Eigen::MatrixXd stressa, stress;
			autodiff.compute_stress_tensor(OutputData(0, el_id, bs, bs, local_pts, displacement), ElasticityTensorType::PK1, stressa);
			real.compute_stress_tensor(OutputData(0, el_id, bs, bs, local_pts, displacement), ElasticityTensorType::PK1, stress);

			for (int i = 0; i < stressa.size(); ++i)
			{
				if (std::isnan(stress(i)))
					REQUIRE(std::isnan(stressa(i)));
				else
					REQUIRE(stressa(i) == Catch::Approx(stress(i)).margin(1e-12));
			}
		}

		// pk2 stress
		{
			Eigen::MatrixXd stressa, stress;
			autodiff.compute_stress_tensor(OutputData(0, el_id, bs, bs, local_pts, displacement), ElasticityTensorType::PK2, stressa);
			real.compute_stress_tensor(OutputData(0, el_id, bs, bs, local_pts, displacement), ElasticityTensorType::PK2, stress);

			for (int i = 0; i < stressa.size(); ++i)
			{
				if (std::isnan(stress(i)))
					REQUIRE(std::isnan(stressa(i)));
				else
					REQUIRE(stressa(i) == Catch::Approx(stress(i)).margin(1e-12));
			}
		}
	}
}

namespace
{
	// Build a State on plane_hole.obj (2D): mesh + P1 bases only; material type irrelevant.
	std::shared_ptr<State> make_state_2d()
	{
		const std::string path = POLYFEM_DATA_DIR;
		json in_args = json({});
		in_args["geometry"]["mesh"] = path + "/plane_hole.obj";
		in_args["geometry"]["surface_selection"] = 7;
		in_args["materials"]["type"] = "NeoHookean";
		in_args["materials"]["E"] = 1e5;
		in_args["materials"]["nu"] = 0.3;
		// Time block avoids the "Static problem needs Dirichlet nodes" check
		in_args["time"]["dt"] = 0.001;
		in_args["time"]["tend"] = 1.0;

		auto state = std::make_shared<State>();
		state->init_logger("", spdlog::level::err, spdlog::level::off, false);
		state->init(in_args, true);
		state->load_mesh();
		state->build_basis();
		return state;
	}

	// Set up a ModifiedNeoHookeanElasticity with E=1e5, nu=0.3.
	ModifiedNeoHookeanElasticity make_modified_assembler(const State &state)
	{
		ModifiedNeoHookeanElasticity a;
		a.set_size(2);
		json mat;
		mat["type"] = "ModifiedNeoHookean";
		mat["E"] = 1e5;
		mat["nu"] = 0.3;
		a.add_multimaterial(0, mat, state.units, state.root_path());
		return a;
	}
} // namespace

// Gradient consistency: compare analytic gradient against central-difference of energy.
TEST_CASE("modified-neohookean-gradient", "[assembler]")
{
	auto state = make_state_2d();
	auto assembler = make_modified_assembler(*state);

	const int dim = 2;
	const int el_id = 0;
	const auto &bs = state->bases[el_id];
	ElementAssemblyValues vals;
	vals.compute(el_id, false, bs, bs);
	const QuadratureVector da = vals.det.array() * vals.quadrature.weights.array();
	const int n_local = static_cast<int>(vals.basis_values.size());

	// Correct global DOF vector size
	Eigen::MatrixXd x0(state->n_bases * dim, 1);
	const double eps = 1e-7;

	for (int trial = 0; trial < 5; ++trial)
	{
		x0.setRandom();
		x0 /= 20.0; // keep J > 0 in all trials

		const Eigen::VectorXd grad = assembler.assemble_gradient(NonLinearAssemblerData(vals, 0, 0, x0, x0, da));

		for (int i = 0; i < n_local; ++i)
		{
			REQUIRE(vals.basis_values[i].global.size() == 1);
			const int g = vals.basis_values[i].global[0].index;
			for (int d = 0; d < dim; ++d)
			{
				Eigen::MatrixXd xp = x0, xm = x0;
				xp(g * dim + d) += eps;
				xm(g * dim + d) -= eps;
				const double ep = assembler.compute_energy(NonLinearAssemblerData(vals, 0, 0, xp, xp, da));
				const double em = assembler.compute_energy(NonLinearAssemblerData(vals, 0, 0, xm, xm, da));
				REQUIRE(grad(i * dim + d) == Catch::Approx((ep - em) / (2.0 * eps)).margin(1e-5));
			}
		}
	}
}

// Hessian consistency: compare analytic hessian against central-difference of gradient.
TEST_CASE("modified-neohookean-hessian", "[assembler]")
{
	auto state = make_state_2d();
	auto assembler = make_modified_assembler(*state);

	const int dim = 2;
	const int el_id = 0;
	const auto &bs = state->bases[el_id];
	ElementAssemblyValues vals;
	vals.compute(el_id, false, bs, bs);
	const QuadratureVector da = vals.det.array() * vals.quadrature.weights.array();
	const int n_local = static_cast<int>(vals.basis_values.size());

	Eigen::MatrixXd x0(state->n_bases * dim, 1);
	const double eps = 1e-5;

	for (int trial = 0; trial < 3; ++trial)
	{
		x0.setRandom();
		x0 /= 20.0;

		const Eigen::MatrixXd hess = assembler.assemble_hessian(NonLinearAssemblerData(vals, 0, 0, x0, x0, da));

		for (int i = 0; i < n_local; ++i)
		{
			const int g = vals.basis_values[i].global[0].index;
			for (int di = 0; di < dim; ++di)
			{
				Eigen::MatrixXd xp = x0, xm = x0;
				xp(g * dim + di) += eps;
				xm(g * dim + di) -= eps;
				const Eigen::VectorXd gp = assembler.assemble_gradient(NonLinearAssemblerData(vals, 0, 0, xp, xp, da));
				const Eigen::VectorXd gm = assembler.assemble_gradient(NonLinearAssemblerData(vals, 0, 0, xm, xm, da));
				const Eigen::VectorXd fd_col = (gp - gm) / (2.0 * eps);

				for (int j = 0; j < n_local; ++j)
				{
					for (int dj = 0; dj < dim; ++dj)
						REQUIRE(hess(j * dim + dj, i * dim + di) == Catch::Approx(fd_col(j * dim + dj)).margin(1e-4));
				}
			}
		}
	}
}

// Barrier is zero when J >= 0.5: at zero displacement, energy == NeoHookean energy.
TEST_CASE("modified-neohookean-barrier-zero", "[assembler]")
{
	auto state = make_state_2d();
	auto modified = make_modified_assembler(*state);

	NeoHookeanElasticity neo;
	neo.set_size(2);
	json mat;
	mat["type"] = "NeoHookean";
	mat["E"] = 1e5;
	mat["nu"] = 0.3;
	neo.add_multimaterial(0, mat, state->units, state->root_path());

	const int el_id = 0;
	const auto &bs = state->bases[el_id];
	ElementAssemblyValues vals;
	vals.compute(el_id, false, bs, bs);
	const QuadratureVector da = vals.det.array() * vals.quadrature.weights.array();

	// At zero displacement J = 1 >> 0.5 → barrier = 0 → energies match
	const Eigen::MatrixXd x_zero = Eigen::MatrixXd::Zero(state->n_bases * 2, 1);
	const NonLinearAssemblerData data(vals, 0, 0, x_zero, x_zero, da);

	REQUIRE(modified.compute_energy(data) == Catch::Approx(neo.compute_energy(data)).margin(1e-14));

	const Eigen::VectorXd g_mod = modified.assemble_gradient(data);
	const Eigen::VectorXd g_neo = neo.assemble_gradient(data);
	for (int i = 0; i < g_mod.size(); ++i)
		REQUIRE(g_mod(i) == Catch::Approx(g_neo(i)).margin(1e-14));
}

// Barrier is active when J < 0.5: compressed element has strictly higher energy than NeoHookean.
TEST_CASE("modified-neohookean-barrier-active", "[assembler]")
{
	auto state = make_state_2d();
	auto modified = make_modified_assembler(*state);

	NeoHookeanElasticity neo;
	neo.set_size(2);
	json mat;
	mat["type"] = "NeoHookean";
	mat["E"] = 1e5;
	mat["nu"] = 0.3;
	neo.add_multimaterial(0, mat, state->units, state->root_path());

	const int el_id = 0;
	const auto &bs = state->bases[el_id];
	ElementAssemblyValues vals;
	vals.compute(el_id, false, bs, bs);
	const QuadratureVector da = vals.det.array() * vals.quadrature.weights.array();
	const int n_local = static_cast<int>(vals.basis_values.size());

	// Scale entire mesh by 0.4 in x: u_x = -0.6 * X_x for every node.
	// For P1 this gives ∂u_x/∂X_x = -0.6 exactly → F = diag(0.4, 1) → J = 0.4 < 0.5.
	// mesh->point(v) is indexed by INPUT vertex; in_node_to_node[v] gives the DOF.
	Eigen::MatrixXd x_comp = Eigen::MatrixXd::Zero(state->n_bases * 2, 1);
	for (int v = 0; v < state->mesh->n_vertices(); ++v)
	{
		const int dof = state->in_node_to_node[v];
		if (dof < state->n_bases)
			x_comp(dof * 2) = -0.6 * state->mesh->point(v)(0); // u_x = -0.6 X_x
	}

	// Verify J at the first quadrature point so failures show the actual deformation gradient.
	{
		Eigen::MatrixXd local_disp(n_local, 2);
		local_disp.setZero();
		for (int i = 0; i < n_local; ++i)
			for (const auto &gv : vals.basis_values[i].global)
				for (int d = 0; d < 2; ++d)
					local_disp(i, d) += gv.val * x_comp(gv.index * 2 + d);

		Eigen::MatrixXd grad(n_local, 2);
		for (int i = 0; i < n_local; ++i)
			grad.row(i) = vals.basis_values[i].grad.row(0);

		const Eigen::MatrixXd def_grad =
			(local_disp.transpose() * grad) * vals.jac_it[0] + Eigen::MatrixXd::Identity(2, 2);
		const double J_check = def_grad.determinant();

		for (int i = 0; i < n_local; ++i)
		{
			const int g = vals.basis_values[i].global[0].index;
			const auto p = state->mesh->point(g);
			INFO("  node " << i << " (global " << g << "): ref=(" << p(0) << ", " << p(1)
				 << ")  u=(" << x_comp(g * 2) << ", " << x_comp(g * 2 + 1) << ")");
		}
		INFO("def_grad:\n" << def_grad);
		INFO("J at quad point 0 = " << J_check);
		REQUIRE(J_check < 0.5); // confirm the deformation actually inverts
	}

	const NonLinearAssemblerData data(vals, 0, 0, x_comp, x_comp, da);
	const double e_mod = modified.compute_energy(data);
	const double e_neo = neo.compute_energy(data);
	{
		INFO("e_neo = " << e_neo << "  e_mod = " << e_mod << "  diff = " << (e_mod - e_neo));
		REQUIRE(std::isfinite(e_mod));
		REQUIRE(e_mod > e_neo); // barrier term > 0 when J < 0.5
	}
}
