#ifndef REMESH_BOTSCH
#define REMESH_BOTSCH



#include <Eigen/Core>

// void adaptive_remesh_botsch(
// 	const Eigen::MatrixXd &Vattrs_in, const Eigen::MatrixXi &F_in, const Eigen::VectorXi &Vselection_in, 
// 	double selection_threshold, double epsilon, int iters, bool project, bool adaptive, int verbose, 
// 	Eigen::MatrixXd &Vattrs_out, Eigen::MatrixXi &F_out, Eigen::VectorXi &Vselection_out
// );


void remesh_botsch(
	const Eigen::MatrixXd &Vattrs_in, 
	const Eigen::MatrixXi &F_in,
	const Eigen::VectorXi &Vselection_in, 
	const Eigen::VectorXd &Vtargetlen_in,
	// ^ "user-specifiable sizing field"; if an item is negative, its negation is interpreted
	// as an epsilon for calc_sizing_field. if positive, it's used directly as the sizingField value
	double selection_threshold, 
	int iters, 
	bool project, 
	int verbose, 
	Eigen::MatrixXd &Vattrs_out, 
	Eigen::MatrixXi &F_out, 
	Eigen::VectorXi &Vselection_out);


// overload with constant targetlen
void remesh_botsch(
	const Eigen::MatrixXd &Vattrs_in, 
	const Eigen::MatrixXi &F_in,
	const Eigen::VectorXi &Vselection_in, 
	double targetlen,
	// ^ "user-specifiable sizing field"; if an item is negative, its negation is interpreted
	// as an epsilon for calc_sizing_field. if positive, it's used directly as the sizingField value
	double selection_threshold, 
	int iters, 
	bool project, 
	int verbose, 
	Eigen::MatrixXd &Vattrs_out, 
	Eigen::MatrixXi &F_out, 
	Eigen::VectorXi &Vselection_out);



void remesh_botsch_adaptive(
	const Eigen::MatrixXd &Vattrs_in, 
	const Eigen::MatrixXi &F_in,
	const Eigen::VectorXi &Vselection_in, 
	double epsilon,
	bool adaptive,
	double selection_threshold, 
	int iters, 
	bool project, 
	int verbose, 
	Eigen::MatrixXd &Vattrs_out, 
	Eigen::MatrixXi &F_out, 
	Eigen::VectorXi &Vselection_out);
	
#endif
