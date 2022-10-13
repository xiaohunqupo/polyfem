// #define EIGEN_STACK_ALLOCATION_LIMIT 0

#include "OgdenElasticity.hpp"
#include <polyfem/basis/Basis.hpp>
#include <polyfem/utils/MatrixUtils.hpp>

#include <polyfem/autogen/auto_eigs.hpp>
#include <polyfem/autogen/auto_elasticity_rhs.hpp>

#include <igl/Timer.h>

namespace polyfem::assembler
{
	namespace
	{
		void fill_mat_from_vect(const std::vector<double> vals, Eigen::VectorXd &vec)
		{
			vec.resize(vals.size());

			for (size_t i = 0; i < vals.size(); ++i)
				vec(i) = vals[i];
		}
	} // namespace

	OgdenElasticity::OgdenElasticity()
	{
		alphas_.resize(1);
		mus_.resize(1);
		Ds_.resize(1);

		alphas_.setOnes();
		mus_.setOnes();
		Ds_.setOnes();
	}

	void OgdenElasticity::add_multimaterial(const int index, const json &params)
	{
		assert(size_ == 2 || size_ == 3);

		if (params.count("alphas"))
		{
			const std::vector<double> tmp = params["alphas"];
			fill_mat_from_vect(tmp, alphas_);
		}
		if (params.count("mus"))
		{
			const std::vector<double> tmp = params["mus"];
			fill_mat_from_vect(tmp, mus_);
		}
		if (params.count("Ds"))
		{
			const std::vector<double> tmp = params["Ds"];
			fill_mat_from_vect(tmp, Ds_);
		}

		assert(alphas_.size() == mus_.size());
		assert(Ds_.size() == mus_.size());
	}

	void OgdenElasticity::set_size(const int size)
	{
		size_ = size;
	}

	Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 3, 1>
	OgdenElasticity::compute_rhs(const AutodiffHessianPt &pt) const
	{
		assert(pt.size() == size());
		Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 3, 1> res;

		// if(size() == 2)
		// 	autogen::saint_venant_2d_function(pt, elasticity_tensor_, res);
		// else if(size() == 3)
		// 	autogen::saint_venant_3d_function(pt, elasticity_tensor_, res);
		// else
		// 	assert(false);

		return res;
	}

	Eigen::VectorXd
	OgdenElasticity::assemble_grad(const NonLinearAssemblerData &data) const
	{
		// igl::Timer time; time.start();

		const int n_bases = vals.basis_values.size();

		return polyfem::gradient_from_energy(
			size(), n_bases, data,
			[&](const NonLinearAssemblerData &data) { return compute_energy_aux<DScalar1<double, Eigen::Matrix<double, 6, 1>>>(data); },
			[&](const NonLinearAssemblerData &data) { return compute_energy_aux<DScalar1<double, Eigen::Matrix<double, 8, 1>>>(data); },
			[&](const NonLinearAssemblerData &data) { return compute_energy_aux<DScalar1<double, Eigen::Matrix<double, 12, 1>>>(data); },
			[&](const NonLinearAssemblerData &data) { return compute_energy_aux<DScalar1<double, Eigen::Matrix<double, 18, 1>>>(data); },
			[&](const NonLinearAssemblerData &data) { return compute_energy_aux<DScalar1<double, Eigen::Matrix<double, 24, 1>>>(data); },
			[&](const NonLinearAssemblerData &data) { return compute_energy_aux<DScalar1<double, Eigen::Matrix<double, 30, 1>>>(data); },
			[&](const NonLinearAssemblerData &data) { return compute_energy_aux<DScalar1<double, Eigen::Matrix<double, 60, 1>>>(data); },
			[&](const NonLinearAssemblerData &data) { return compute_energy_aux<DScalar1<double, Eigen::Matrix<double, 81, 1>>>(data); },
			[&](const NonLinearAssemblerData &data) { return compute_energy_aux<DScalar1<double, Eigen::Matrix<double, Eigen::Dynamic, 1, 0, SMALL_N, 1>>>(data); },
			[&](const NonLinearAssemblerData &data) { return compute_energy_aux<DScalar1<double, Eigen::Matrix<double, Eigen::Dynamic, 1, 0, BIG_N, 1>>>(data); },
			[&](const NonLinearAssemblerData &data) { return compute_energy_aux<DScalar1<double, Eigen::VectorXd>>(data); });
	}

	Eigen::MatrixXd
	OgdenElasticity::assemble_hessian(const NonLinearAssemblerData &data) const
	{
		// igl::Timer time; time.start();

		const int n_bases = vals.basis_values.size();
		return polyfem::hessian_from_energy(
			size(), n_bases, data,
			[&](const NonLinearAssemblerData &data) { return compute_energy_aux<DScalar2<double, Eigen::Matrix<double, 6, 1>, Eigen::Matrix<double, 6, 6>>>(data); },
			[&](const NonLinearAssemblerData &data) { return compute_energy_aux<DScalar2<double, Eigen::Matrix<double, 8, 1>, Eigen::Matrix<double, 8, 8>>>(data); },
			[&](const NonLinearAssemblerData &data) { return compute_energy_aux<DScalar2<double, Eigen::Matrix<double, 12, 1>, Eigen::Matrix<double, 12, 12>>>(data); },
			[&](const NonLinearAssemblerData &data) { return compute_energy_aux<DScalar2<double, Eigen::Matrix<double, 18, 1>, Eigen::Matrix<double, 18, 18>>>(data); },
			[&](const NonLinearAssemblerData &data) { return compute_energy_aux<DScalar2<double, Eigen::Matrix<double, 24, 1>, Eigen::Matrix<double, 24, 24>>>(data); },
			[&](const NonLinearAssemblerData &data) { return compute_energy_aux<DScalar2<double, Eigen::Matrix<double, 30, 1>, Eigen::Matrix<double, 30, 30>>>(data); },
			[&](const NonLinearAssemblerData &data) { return compute_energy_aux<DScalar2<double, Eigen::Matrix<double, 60, 1>, Eigen::Matrix<double, 60, 60>>>(data); },
			[&](const NonLinearAssemblerData &data) { return compute_energy_aux<DScalar2<double, Eigen::Matrix<double, 81, 1>, Eigen::Matrix<double, 81, 81>>>(data); },
			[&](const NonLinearAssemblerData &data) { return compute_energy_aux<DScalar2<double, Eigen::Matrix<double, Eigen::Dynamic, 1, 0, SMALL_N, 1>, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, SMALL_N, SMALL_N>>>(data); },
			[&](const NonLinearAssemblerData &data) { return compute_energy_aux<DScalar2<double, Eigen::VectorXd, Eigen::MatrixXd>>(data); });
	}

	void OgdenElasticity::compute_stress_tensor(const int el_id, const ElementBases &bs, const ElementBases &gbs, const Eigen::MatrixXd &local_pts, const Eigen::MatrixXd &displacement, Eigen::MatrixXd &stresses) const
	{
		assign_stress_tensor(el_id, bs, gbs, local_pts, displacement, size() * size(), stresses, [&](const Eigen::MatrixXd &stress) {
			Eigen::MatrixXd tmp = stress;
			auto a = Eigen::Map<Eigen::MatrixXd>(tmp.data(), 1, size() * size());
			return Eigen::MatrixXd(a);
		});
	}

	void OgdenElasticity::compute_von_mises_stresses(const int el_id, const basis::ElementBases &bs, const basis::ElementBases &gbs, const Eigen::MatrixXd &local_pts, const Eigen::MatrixXd &displacement, Eigen::MatrixXd &stresses) const
	{
		assign_stress_tensor(el_id, bs, gbs, local_pts, displacement, 1, stresses, [&](const Eigen::MatrixXd &stress) {
			Eigen::Matrix<double, 1, 1> res;
			res.setConstant(von_mises_stress_for_stress_tensor(stress));
			return res;
		});
	}

	void OgdenElasticity::assign_stress_tensor(const int el_id, const ElementBases &bs, const ElementBases &gbs, const Eigen::MatrixXd &local_pts, const Eigen::MatrixXd &displacement, const int all_size, Eigen::MatrixXd &all, const std::function<Eigen::MatrixXd(const Eigen::MatrixXd &)> &fun) const
	{
		Eigen::MatrixXd displacement_grad(size(), size());

		assert(displacement.cols() == 1);
		ElementAssemblyValues vals;
		vals.compute(el_id, size() == 3, local_pts, bs, bs);

		all.resize(local_pts.rows(), all_size);

		// TODO implement stresses
		assert(false);

		// for(long p = 0; p < local_pts.rows(); ++p)
		// {
		// 	compute_diplacement_grad(size(), bs, vals, local_pts, p, displacement, displacement_grad);

		// 	Eigen::MatrixXd strain = strain_from_disp_grad(displacement_grad);
		// 	Eigen::MatrixXd stress_tensor(size(), size());

		// 	if(size() == 2)
		// 	{
		// 		std::array<double, 3> eps;
		// 		eps[0] = strain(0,0);
		// 		eps[1] = strain(1,1);
		// 		eps[2] = 2*strain(0,1);

		// 		stress_tensor <<
		// 		stress(eps, 0), stress(eps, 2),
		// 		stress(eps, 2), stress(eps, 1);
		// 	}
		// 	else
		// 	{
		// 		std::array<double, 6> eps;
		// 		eps[0] = strain(0,0);
		// 		eps[1] = strain(1,1);
		// 		eps[2] = strain(2,2);
		// 		eps[3] = 2*strain(1,2);
		// 		eps[4] = 2*strain(0,2);
		// 		eps[5] = 2*strain(0,1);

		// 		stress_tensor <<
		// 		stress(eps, 0), stress(eps, 5), stress(eps, 4),
		// 		stress(eps, 5), stress(eps, 1), stress(eps, 3),
		// 		stress(eps, 4), stress(eps, 3), stress(eps, 2);
		// 	}

		// 	stress_tensor = (Eigen::MatrixXd::Identity(size(), size()) + displacement_grad) * stress_tensor;

		// 	all.row(p) = fun(stress_tensor);
		// }
	}

	double OgdenElasticity::compute_energy(const NonLinearAssemblerData &data) const
	{
		return compute_energy_aux<double>(data);
	}

	// http://abaqus.software.polimi.it/v6.14/books/stm/default.htm?startat=ch04s06ath123.html Ogden form
	template <typename T>
	T OgdenElasticity::compute_energy_aux(const NonLinearAssemblerData &data) const
	{
		typedef Eigen::Matrix<T, Eigen::Dynamic, 1> AutoDiffVect;
		typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, 0, 3, 3> AutoDiffGradMat;

		AutoDiffVect local_disp;
		get_local_disp(data, size(), local_disp);

		AutoDiffGradMat def_grad(size(), size());

		T energy = T(0.0);

		const int n_pts = data.da.size();
		for (long p = 0; p < n_pts; ++p)
		{
			compute_disp_grad_at_quad(data, local_disp, p, size(), def_grad);
			// Id + grad d
			for (int d = 0; d < size; ++d)
				def_grad(d, d) += T(1);

			Eigen::Matrix<T, Eigen::Dynamic, 1, 0, 3, 1> eigs;

			if (size() == 2)
				autogen::eigs_2d<T>(def_grad, eigs);
			else if (size() == 3)
				autogen::eigs_3d<T>(def_grad, eigs);
			else
				assert(false);

			const T J = polyfem::determinant(def_grad);
			const T Jdenom = pow(J, -1. / size());

			for (long i = 0; i < eigs.size(); ++i)
				eigs(i) = eigs(i) * Jdenom;

			auto val = T(0);
			for (long N = 0; N < alphas_.size(); ++N)
			{
				auto tmp = T(-size());
				const double alpha = alphas_(N);
				const double mu = mus_(N);

				for (long i = 0; i < eigs.size(); ++i)
					tmp += pow(eigs(i), alpha);

				val += 2 * mu / (alpha * alpha) * tmp;
			}

			// std::cout<<val<<std::endl;

			for (long N = 0; N < Ds_.size(); ++N)
			{
				const double D = Ds_(N);

				val += 1. / D * pow(J - T(1), 2 * (N + 1));
			}

			energy += val * da(p);
		}
		// std::cout<<"\n\n------------\n\n\n"<<std::endl;
		return energy;
	}

} // namespace polyfem::assembler
