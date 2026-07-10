#include <polyfem/State.hpp>

#include <polyfem/Units.hpp>
#include <polyfem/assembler/Assembler.hpp>
#include <polyfem/assembler/AssemblyValsCache.hpp>
#include <polyfem/assembler/NeoHookeanElasticity.hpp>
#include <polyfem/assembler/NeoHookeanElasticityAutodiff.hpp>
#include <polyfem/assembler/InversionBarrier.hpp>
#include <polyfem/assembler/SumModel.hpp>
#include <polyfem/utils/RefElementSampler.hpp>
#include <polyfem/varforms/VarForm.hpp>

#include "VarFormTestAccess.hpp"

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
	test::VarFormTestAccess::prepare(*state.variational_formulation);

	SparseMatrixCache mat_cache;
	StiffnessMatrix hessian, stiffness;
	varform::VarForm &form = *state.variational_formulation;
	const test::VarFormDebugData debug = test::VarFormTestAccess::debug_data(form);
	REQUIRE(debug.assembler != nullptr);
	REQUIRE(debug.mesh != nullptr);
	REQUIRE(debug.bases != nullptr);
	REQUIRE(debug.geometry_bases != nullptr);
	AssemblyValsCache ass_vals_cache;
	ass_vals_cache.init_empty();
	Eigen::MatrixXd disp(debug.n_bases * debug.mesh->dimension(), 1);
	disp.setZero();

	REQUIRE(test::VarFormTestAccess::build_stiffness_mat(form, stiffness));

	for (int rand = 0; rand < 10; ++rand)
	{
		debug.assembler->assemble_hessian(
			debug.mesh->is_volume(), debug.n_bases, false,
			*debug.bases, *debug.geometry_bases, ass_vals_cache, 0, 0, disp, Eigen::MatrixXd(), mat_cache, hessian);

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
	test::VarFormTestAccess::prepare(*state.variational_formulation);

	SparseMatrixCache mat_cache;
	StiffnessMatrix hessian, stiffness;
	varform::VarForm &form = *state.variational_formulation;
	const test::VarFormDebugData debug = test::VarFormTestAccess::debug_data(form);
	REQUIRE(debug.assembler != nullptr);
	REQUIRE(debug.mesh != nullptr);
	REQUIRE(debug.bases != nullptr);
	REQUIRE(debug.geometry_bases != nullptr);
	AssemblyValsCache ass_vals_cache;
	ass_vals_cache.init_empty();
	Eigen::MatrixXd disp(debug.n_bases * debug.mesh->dimension(), 1);
	disp.setZero();

	REQUIRE(test::VarFormTestAccess::build_stiffness_mat(form, stiffness));

	for (int rand = 0; rand < 10; ++rand)
	{
		debug.assembler->assemble_hessian(
			debug.mesh->is_volume(), debug.n_bases, false,
			*debug.bases, *debug.geometry_bases, ass_vals_cache, 0, 0, disp, Eigen::MatrixXd(), mat_cache, hessian);

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
	test::VarFormTestAccess::prepare(*state.variational_formulation);

	NeoHookeanAutodiff autodiff;
	NeoHookeanElasticity real;

	autodiff.set_size(2);
	real.set_size(2);

	Units units;
	units.init(state.args["units"]);
	const test::VarFormDebugData debug = test::VarFormTestAccess::debug_data(*state.variational_formulation);
	REQUIRE(debug.mesh != nullptr);
	REQUIRE(debug.bases != nullptr);
	REQUIRE(debug.geometry_bases != nullptr);

	autodiff.add_multimaterial(0, in_args["materials"], units, debug.root_path);
	real.add_multimaterial(0, in_args["materials"], units, debug.root_path);

	const int el_id = 0;
	const auto &bs = (*debug.bases)[el_id];
	const auto &gbs = (*debug.geometry_bases)[el_id];
	Eigen::MatrixXd local_pts;
	Eigen::MatrixXi f;
	regular_2d_grid(10, true, local_pts, f);

	Eigen::MatrixXd displacement(debug.n_bases, 1);

	ElementAssemblyValues vals;
	vals.compute(el_id, debug.mesh->is_volume(), bs, gbs);

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
			autodiff.compute_stress_tensor(OutputData(0, el_id, bs, gbs, local_pts, displacement), ElasticityTensorType::F, stressa);
			real.compute_stress_tensor(OutputData(0, el_id, bs, gbs, local_pts, displacement), ElasticityTensorType::F, stress);

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
			autodiff.compute_stress_tensor(OutputData(0, el_id, bs, gbs, local_pts, displacement), ElasticityTensorType::CAUCHY, stressa);
			real.compute_stress_tensor(OutputData(0, el_id, bs, gbs, local_pts, displacement), ElasticityTensorType::CAUCHY, stress);

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
			autodiff.compute_stress_tensor(OutputData(0, el_id, bs, gbs, local_pts, displacement), ElasticityTensorType::PK1, stressa);
			real.compute_stress_tensor(OutputData(0, el_id, bs, gbs, local_pts, displacement), ElasticityTensorType::PK1, stress);

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
			autodiff.compute_stress_tensor(OutputData(0, el_id, bs, gbs, local_pts, displacement), ElasticityTensorType::PK2, stressa);
			real.compute_stress_tensor(OutputData(0, el_id, bs, gbs, local_pts, displacement), ElasticityTensorType::PK2, stress);

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
		test::VarFormTestAccess::prepare(*state->variational_formulation);
		return state;
	}

	// Set up a SumModel combining NeoHookean and InversionBarrier, both with E=1e5, nu=0.3.
	// This is the composable replacement for the old ModifiedNeoHookeanElasticity.
	std::shared_ptr<SumModel> make_modified_assembler(const Units &units, const std::string &root_path)
	{
		auto a = std::make_shared<SumModel>();
		a->set_size(2);
		json mat;
		mat["type"] = "MaterialSum";
		mat["models"] = json::array();
		json neo_mat;
		neo_mat["type"] = "NeoHookean";
		neo_mat["E"] = 1e5;
		neo_mat["nu"] = 0.3;
		mat["models"].push_back(neo_mat);
		json barrier_mat;
		barrier_mat["type"] = "InversionBarrier";
		barrier_mat["E"] = 1e5;
		barrier_mat["nu"] = 0.3;
		mat["models"].push_back(barrier_mat);
		a->add_multimaterial(0, mat, units, root_path);
		return a;
	}
} // namespace

// Gradient consistency: compare analytic gradient against central-difference of energy.
TEST_CASE("modified-neohookean-gradient", "[assembler]")
{
	auto state = make_state_2d();
	const test::VarFormDebugData debug = test::VarFormTestAccess::debug_data(*state->variational_formulation);
	Units units;
	units.init(state->args["units"]);
	auto assembler = make_modified_assembler(units, debug.root_path);

	const int dim = 2;
	const int el_id = 0;
	const auto &bs = (*debug.bases)[el_id];
	ElementAssemblyValues vals;
	vals.compute(el_id, false, bs, bs);
	const QuadratureVector da = vals.det.array() * vals.quadrature.weights.array();
	const int n_local = static_cast<int>(vals.basis_values.size());

	// Correct global DOF vector size
	Eigen::MatrixXd x0(debug.n_bases * dim, 1);
	const double eps = 1e-7;

	for (int trial = 0; trial < 5; ++trial)
	{
		x0.setRandom();
		x0 /= 20.0; // keep J > 0 in all trials

		const Eigen::VectorXd grad = assembler->assemble_gradient(NonLinearAssemblerData(vals, 0, 0, x0, x0, da));

		for (int i = 0; i < n_local; ++i)
		{
			REQUIRE(vals.basis_values[i].global.size() == 1);
			const int g = vals.basis_values[i].global[0].index;
			for (int d = 0; d < dim; ++d)
			{
				Eigen::MatrixXd xp = x0, xm = x0;
				xp(g * dim + d) += eps;
				xm(g * dim + d) -= eps;
				const double ep = assembler->compute_energy(NonLinearAssemblerData(vals, 0, 0, xp, xp, da));
				const double em = assembler->compute_energy(NonLinearAssemblerData(vals, 0, 0, xm, xm, da));
				REQUIRE(grad(i * dim + d) == Catch::Approx((ep - em) / (2.0 * eps)).margin(1e-5));
			}
		}
	}
}

// Hessian consistency: compare analytic hessian against central-difference of gradient.
TEST_CASE("modified-neohookean-hessian", "[assembler]")
{
	auto state = make_state_2d();
	const test::VarFormDebugData debug = test::VarFormTestAccess::debug_data(*state->variational_formulation);
	Units units;
	units.init(state->args["units"]);
	auto assembler = make_modified_assembler(units, debug.root_path);

	const int dim = 2;
	const int el_id = 0;
	const auto &bs = (*debug.bases)[el_id];
	ElementAssemblyValues vals;
	vals.compute(el_id, false, bs, bs);
	const QuadratureVector da = vals.det.array() * vals.quadrature.weights.array();
	const int n_local = static_cast<int>(vals.basis_values.size());

	Eigen::MatrixXd x0(debug.n_bases * dim, 1);
	const double eps = 1e-5;

	for (int trial = 0; trial < 3; ++trial)
	{
		x0.setRandom();
		x0 /= 20.0;

		const Eigen::MatrixXd hess = assembler->assemble_hessian(NonLinearAssemblerData(vals, 0, 0, x0, x0, da));

		for (int i = 0; i < n_local; ++i)
		{
			const int g = vals.basis_values[i].global[0].index;
			for (int di = 0; di < dim; ++di)
			{
				Eigen::MatrixXd xp = x0, xm = x0;
				xp(g * dim + di) += eps;
				xm(g * dim + di) -= eps;
				const Eigen::VectorXd gp = assembler->assemble_gradient(NonLinearAssemblerData(vals, 0, 0, xp, xp, da));
				const Eigen::VectorXd gm = assembler->assemble_gradient(NonLinearAssemblerData(vals, 0, 0, xm, xm, da));
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
	const test::VarFormDebugData debug = test::VarFormTestAccess::debug_data(*state->variational_formulation);
	Units units;
	units.init(state->args["units"]);
	auto modified = make_modified_assembler(units, debug.root_path);

	NeoHookeanElasticity neo;
	neo.set_size(2);
	json mat;
	mat["type"] = "NeoHookean";
	mat["E"] = 1e5;
	mat["nu"] = 0.3;
	neo.add_multimaterial(0, mat, units, debug.root_path);

	const int el_id = 0;
	const auto &bs = (*debug.bases)[el_id];
	ElementAssemblyValues vals;
	vals.compute(el_id, false, bs, bs);
	const QuadratureVector da = vals.det.array() * vals.quadrature.weights.array();

	// At zero displacement J = 1 >> 0.5 → barrier = 0 → energies match
	const Eigen::MatrixXd x_zero = Eigen::MatrixXd::Zero(debug.n_bases * 2, 1);
	const NonLinearAssemblerData data(vals, 0, 0, x_zero, x_zero, da);

	REQUIRE(modified->compute_energy(data) == Catch::Approx(neo.compute_energy(data)).margin(1e-14));

	const Eigen::VectorXd g_mod = modified->assemble_gradient(data);
	const Eigen::VectorXd g_neo = neo.assemble_gradient(data);
	for (int i = 0; i < g_mod.size(); ++i)
		REQUIRE(g_mod(i) == Catch::Approx(g_neo(i)).margin(1e-14));
}

// Barrier activates below J=0.5: at large compression e_mod > e_neo.
TEST_CASE("modified-neohookean-barrier-active", "[assembler]")
{
	// Single known triangle: vertices at (0,0),(1,0),(0,1). DOFs are 0,1,2.
	Eigen::MatrixXd V(3, 2);
	V << 0, 0,
		1, 0,
		0, 1;
	Eigen::MatrixXi F(1, 3);
	F << 0, 1, 2;

	json in_args;
	in_args["geometry"][0]["mesh"] = "";
	in_args["materials"]["type"] = "NeoHookean";
	in_args["materials"]["E"] = 1e5;
	in_args["materials"]["nu"] = 0.3;
	in_args["time"]["dt"] = 0.001;
	in_args["time"]["tend"] = 1.0;

	State state;
	state.init_logger("", spdlog::level::err, spdlog::level::off, false);
	state.init(in_args, false);
	state.load_mesh(V, F);
	test::VarFormTestAccess::prepare(*state.variational_formulation);

	const auto debug = test::VarFormTestAccess::debug_data(*state.variational_formulation);
	Units units;
	units.init(state.args["units"]);

	auto modified = make_modified_assembler(units, debug.root_path);
	NeoHookeanElasticity neo;
	neo.set_size(2);
	json mat;
	mat["type"] = "NeoHookean";
	mat["E"] = 1e5;
	mat["nu"] = 0.3;
	neo.add_multimaterial(0, mat, units, debug.root_path);

	const auto &bs = (*debug.bases)[0];
	ElementAssemblyValues vals;
	vals.compute(0, false, bs, bs);
	const QuadratureVector da = vals.det.array() * vals.quadrature.weights.array();

	// u = (s-1)*V, packed as [u0x,u0y, u1x,u1y, u2x,u2y].
	// DOFs 0,1,2 correspond directly to triangle vertices 0,1,2.
	auto make_x = [&](double s) {
		Eigen::MatrixXd x = Eigen::MatrixXd::Zero(6, 1);
		for (int i = 0; i < 3; ++i)
		{
			x(i * 2) = (s - 1.0) * V(i, 0);
			x(i * 2 + 1) = (s - 1.0) * V(i, 1);
		}
		return x;
	};

	// Zero displacement: barrier inactive, energies equal.
	{
		const auto x = make_x(1.0);
		const NonLinearAssemblerData data(vals, 0, 0, x, x, da);
		REQUIRE(modified->compute_energy(data) == Catch::Approx(neo.compute_energy(data)).margin(1e-14));
	}

	// s=0.9: J=0.81 > 0.5, barrier inactive, energies equal.
	{
		const auto x = make_x(0.9);
		const NonLinearAssemblerData data(vals, 0, 0, x, x, da);
		REQUIRE(modified->compute_energy(data) == Catch::Approx(neo.compute_energy(data)).margin(1e-10));
	}

	// s=0.6: J=0.36 < 0.5, barrier active, e_mod > e_neo.
	{
		const auto x = make_x(0.6);
		const NonLinearAssemblerData data(vals, 0, 0, x, x, da);
		REQUIRE(std::isfinite(modified->compute_energy(data)));
		REQUIRE(modified->compute_energy(data) > neo.compute_energy(data));
	}
}
