/*
s file is part of ACADO Toolkit.
 *
 *    ACADO Toolkit -- A Toolkit for Automatic Control and Dynamic Optimization.
 *    Copyright (C) 2008-2014 by Boris Houska, Hans Joachim Ferreau,
 *    Milan Vukov, Rien Quirynen, KU Leuven.
 *    Developed within the Optimization in Engineering Center (OPTEC)
 *    under supervision of Moritz Diehl. All rights reserved.
 *
 *    ACADO Toolkit is free s0ftware; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    ACADO Toolkit is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with ACADO Toolkit; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


 /**
  *    \file   examples/ocp/bioreactor.cpp
  *    \author Boris Houska, Filip Logist, Rien Quirynen
  *    \date   2014
  */
#include <acado_optimal_control.hpp>
#include <acado_gnuplot.hpp>
#include <acado_code_generation.hpp>

/* >>> start tutorial code >>> */
int main( ){

    USING_NAMESPACE_ACADO

    // INTRODUCE THE VARIABLES:
    // -------------------------
    DifferentialState     x,y,theta,s;
    Control               v,w;
    DifferentialEquation  f;

	OnlineData a_X1;
	OnlineData b_X1;
	OnlineData c_X1;
	OnlineData d_X1;
	OnlineData a_Y1;
	OnlineData b_Y1;
	OnlineData c_Y1;
	OnlineData d_Y1;

	OnlineData a_X2;
	OnlineData b_X2;
	OnlineData c_X2;
	OnlineData d_X2;
	OnlineData a_Y2;
	OnlineData b_Y2;
	OnlineData c_Y2;
	OnlineData d_Y2;

	OnlineData a_X3;
	OnlineData b_X3;
	OnlineData c_X3;
	OnlineData d_X3;
	OnlineData a_Y3;
	OnlineData b_Y3;
	OnlineData c_Y3;
	OnlineData d_Y3;

	OnlineData Wx;
	OnlineData Wy;
	OnlineData Wv;
	OnlineData Ww;

	OnlineData s01;
	OnlineData s02;
	OnlineData s03;

	OnlineData vref1;
	OnlineData vref2;
	OnlineData vref3;

	OnlineData delta1;
	OnlineData delta2;

	Expression lambda1 = 1/(1 + exp((s - delta1)/0.1));
	Expression lambda2 = 1/(1 + exp((s - delta2)/0.1));

	Expression x_path1 = (a_X1*(s-s01)*(s-s01)*(s-s01) + b_X1*(s-s01)*(s-s01) + c_X1*(s-s01) + d_X1) ;
	Expression y_path1 = (a_Y1*(s-s01)*(s-s01)*(s-s01) + b_Y1*(s-s01)*(s-s01) + c_Y1*(s-s01) + d_Y1) ;
	Expression dx_path1 = (3*a_X1*(s-s01)*(s-s01) + 2*b_X1*(s-s01) + c_X1) ;
	Expression dy_path1 = (3*a_Y1*(s-s01)*(s-s01) + 2*b_Y1*(s-s01) + c_Y1) ;

	Expression x_path2 = (a_X2*(s-s02)*(s-s02)*(s-s02) + b_X2*(s-s02)*(s-s02) + c_X2*(s-s02) + d_X2) ;
	Expression y_path2 = (a_Y2*(s-s02)*(s-s02)*(s-s02) + b_Y2*(s-s02)*(s-s02) + c_Y2*(s-s02) + d_Y2) ;
	Expression dx_path2 = (3*a_X2*(s-s02)*(s-s02) + 2*b_X2*(s-s02) + c_X2) ;
	Expression dy_path2 = (3*a_Y2*(s-s02)*(s-s02) + 2*b_Y2*(s-s02) + c_Y2) ;

	Expression x_path3 = (a_X3*(s-s03)*(s-s03)*(s-s03) + b_X3*(s-s03)*(s-s03) + c_X3*(s-s03) + d_X3) ;
	Expression y_path3 = (a_Y3*(s-s03)*(s-s03)*(s-s03) + b_Y3*(s-s03)*(s-s03) + c_Y3*(s-s03) + d_Y3) ;
	Expression dx_path3 = (3*a_X3*(s-s03)*(s-s03) + 2*b_X3*(s-s03) + c_X3) ;
	Expression dy_path3 = (3*a_Y3*(s-s03)*(s-s03) + 2*b_Y3*(s-s03) + c_Y3) ;

	Expression x_path = lambda1*lambda2*x_path1 + (1 - lambda1)*lambda2*x_path2 + (1 - lambda2)*x_path3;
	Expression y_path = lambda1*lambda2*y_path1 + (1 - lambda1)*lambda2*y_path2 + (1 - lambda2)*y_path3;
	Expression dx_path = lambda1*lambda2*dx_path1 + (1 - lambda1)*lambda2*dx_path2 + (1 - lambda2)*dx_path3;
	Expression dy_path = lambda1*lambda2*dy_path1 + (1 - lambda1)*lambda2*dy_path2 + (1 - lambda2)*dy_path3;

	Expression abs_grad = sqrt(dx_path.getPowInt(2) + dy_path.getPowInt(2));
	Expression dx_path_norm = dx_path/abs_grad;
	Expression dy_path_norm =  dy_path/abs_grad;

	Expression vref = lambda1*lambda2*vref1 + (1 - lambda1)*lambda2*vref2 + (1 - lambda2)*vref3;

    // DEFINE A DIFFERENTIAL EQUATION:
    // -------------------------------
    
    f << dot(x) == v*cos(theta);
    f << dot(y) == v*sin(theta);
    f << dot(theta) == w;
	f << dot(s) == v;

    // DEFINE AN OPTIMAL CONTROL PROBLEM:
    // ----------------------------------
    OCP ocp( 0.0, 2.5, 25.0 );

    // Need to set the number of online variables!
    ocp.setNOD(36);

	Expression error_contour   = dy_path_norm * (x - x_path) - dx_path_norm * (y - y_path);

	Expression error_lag       = -dx_path_norm * (x - x_path) - dy_path_norm * (y - y_path);

	ocp.minimizeLagrangeTerm( Wx*error_contour*error_contour + Wy*error_lag*error_lag + Ww*w*w +Wv*(v-vref)*(v-vref) ); // weight this with the physical cost!!!
	ocp.setModel(f);

    ocp.subjectTo( -2.0 <= v <= 2.0 );
    ocp.subjectTo( -1.0 <= w <= 1.0 );

    // DEFINE AN MPC EXPORT MODULE AND GENERATE THE CODE:
	// ----------------------------------------------------------
	OCPexport mpc( ocp );

	mpc.set( HESSIAN_APPROXIMATION,            EXACT_HESSIAN  		);
	mpc.set( DISCRETIZATION_TYPE,              MULTIPLE_SHOOTING 	);
	mpc.set( INTEGRATOR_TYPE,                  INT_RK4			    );
	mpc.set( NUM_INTEGRATOR_STEPS,             25            		);
	mpc.set( QP_SOLVER,                        QP_QPOASES    		);
	mpc.set( HOTSTART_QP,                      NO             		);
	mpc.set( GENERATE_TEST_FILE,               YES            		);
	mpc.set( GENERATE_MAKE_FILE,               YES            		);
	mpc.set( GENERATE_MATLAB_INTERFACE,        NO            		);
	mpc.set( SPARSE_QP_SOLUTION, 		       FULL_CONDENSING_N2	);
	mpc.set( DYNAMIC_SENSITIVITY, 		       SYMMETRIC			);
	mpc.set( CG_HARDCODE_CONSTRAINT_VALUES,    NO 					);
	mpc.set( CG_USE_VARIABLE_WEIGHTING_MATRIX, YES 					);

	mpc.exportCode( "generated_mpc" );
	mpc.printDimensionsQP( );
	// ----------------------------------------------------------

    return 0;
}
