#pragma once

#define _USE_MATH_DEFINES
#include <cmath>  
#include "time.h"
#include <chrono>
#include "TriangleTile.h"
#include "Clusters.h"
#include "Particle_Dynamics.h"
#include "Wave_States.h"


class Compute
{
public:
	double Min_Distance;
	double Max_Distance;

	double ClassicalEnergy;

	Clusters * Hierarchical_Clusters;
	Eigen::VectorXd* Laplacian_Energy;
	Eigen::VectorXd* Laplacian_Energy_Derivative;
	Eigen::VectorXd* Laplacian_Energy_Derivative_smoothed;

	Eigen::VectorXd* Commutator_Energy;
	Eigen::VectorXd* Mass_Vector;
	Eigen::VectorXd* Energy_Vector;
	Eigen::VectorXd ClassicalHamiltonian_Energy;

	Eigen::MatrixXd* Laplacian_Orthonormal_Transformation;
	Eigen::MatrixXd* Energy_Orthonormal_Transformation;
	Eigen::MatrixXd* Commutator_Orthonormal_Transformation_Real;
	Eigen::MatrixXd* Commutator_Orthonormal_Transformation_Imag;
	Eigen::MatrixXd ClassicalHamiltonian_EigenStates;

	Multiplicity* Laplacian_Multiplicity;

	// Perturbation
	int Peturb_Order;

	void Perturb(double** dH, int num_points, long num_scale_bins, double* energy_vector, double** orthonormal_transformation);
	void Perturb(Eigen::MatrixXd dH, int num_points, long num_scale_bins, Eigen::VectorXd Energy_Vector_P, Eigen::MatrixXd Orthonormal_Transformation_P);
	// Smoothing
	void Smooth(Eigen::VectorXd* values, Eigen::MatrixXd* vectors, int scale_counter, int num_points);
	void Smooth(Eigen::VectorXd* values, Eigen::MatrixXcd* vectors, int scale_counter, int num_points);
	void Smooth(Eigen::VectorXd* values, Eigen::MatrixXd* vectors, Multiplicity* values_multiplicity, int num_points, int num_scale_bins);
	//

	void Run(double** positions, double** velocities, double* masses, int num_points, int dimension, double dt, bool eigenvectors_flag = false, bool perturb_flag = false);
	void Run(double** positions, double** velocities, double* masses, int num_points, int dimension, double dt, long num_scale_bins, bool eigenvectors_flag = false,  bool perturb_flag = false, bool smooth_flag = false);

	Compute(int num_points, long num_scale_bins, int perturb_order);
	Compute(int num_points, long num_scale_bins);
	~Compute();
private:
	long Number_Pairs;

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> Laplacian_Eigenstructure;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> MutualInteraction_Eigenstructure;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> Commutator_Eigenstructure;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ClassicalHamiltonian_Eigenstructure;

	Eigen::VectorXd Vac;
	Eigen::MatrixXd Correlation_Operator;
	Eigen::MatrixXd Mass_Operator;
	Eigen::MatrixXd ClassicalEnergy_Hamiltonian;
	Eigen::MatrixXd* Laplacian;

	Eigen::MatrixXcd* Commutator_Orthonormal_Transformation;

	//std::ofstream logFile;
};
