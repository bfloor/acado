/*
s file is part of ACADO Toolkit.
 *
 *    ACADO Toolkit -- A Toolkit for Automatic Control and Dynamic Optimization.
 *    Copyright (C) 2008-2014 by Boris Houska, Hans Joachim Ferreau,
 *    Milan Vukov, Rien Quirynen, KU Leuven.
 *    Developed within the Optimization in Engineering Center (OPTEC)
 *    under supervision of Moritz Diehl. All rights reserved.
 *
 *    ACADO Toolkit is free software; you can redistribute it and/or
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
    Control               v,w   ;
    DifferentialEquation  f    ;

    //OnlineData goal_x;
    //OnlineData goal_y;
    //OnlineData goal_theta;


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

	OnlineData Wx;
	OnlineData Wy;
	OnlineData Wv;
	OnlineData Ww;

	OnlineData ws1;
	OnlineData ws2;

	OnlineData delta;

	Expression lambda = 1/(1 + exp((s - delta)/0.01));

	Expression x_path1 = (a_X1*(s-ws1)*(s-ws1)*(s-ws1) + b_X1*(s-ws1)*(s-ws1) + c_X1*(s-ws1) + d_X1) ;
	Expression y_path1 = (a_Y1*(s-ws1)*(s-ws1)*(s-ws1) + b_Y1*(s-ws1)*(s-ws1) + c_Y1*(s-ws1) + d_Y1) ;
	Expression dx_path1 = (3*a_X1*(s-ws1)*(s-ws1) + 2*b_X1*(s-ws1) + c_X1) ;
	Expression dy_path1 = (3*a_Y1*(s-ws1)*(s-ws1) + 2*b_Y1*(s-ws1) + c_Y1) ;

	Expression x_path2 = (a_X2*(s-ws2)*(s-ws2)*(s-ws2) + b_X2*(s-ws2)*(s-ws2) + c_X2*(s-ws2) + d_X2) ;
	Expression y_path2 = (a_Y2*(s-ws2)*(s-ws2)*(s-ws2) + b_Y2*(s-ws2)*(s-ws2) + c_Y2*(s-ws2) + d_Y2) ;
	Expression dx_path2 = (3*a_X2*(s-ws2)*(s-ws2) + 2*b_X2*(s-ws2) + c_X2) ;
	Expression dy_path2 = (3*a_Y2*(s-ws2)*(s-ws2) + 2*b_Y2*(s-ws2) + c_Y2) ;

	Expression x_path = lambda*x_path1 + (1 - lambda)*x_path2;
	Expression y_path = lambda*y_path1 + (1 - lambda)*y_path2;
	Expression dx_path = lambda*dx_path1 + (1 - lambda)*dx_path2;
	Expression dy_path = lambda*dy_path1 + (1 - lambda)*dy_path2;

	Expression abs_grad = sqrt(dx_path.getPowInt(2) + dy_path.getPowInt(2));
	Expression dx_path_norm = dx_path/abs_grad;
	Expression dy_path_norm =  dy_path/abs_grad;

    // DEFINE A DIFFERENTIAL EQUATION:
    // -------------------------------
    
    f << dot(x) == v*cos(theta);
    f << dot(y) == v*sin(theta);
    f << dot(theta) == w;
	f << dot(s) == v;

    // DEFINE AN OPTIMAL CONTROL PROBLEM:
    // ----------------------------------
    OCP ocp( 0.0, 5.0, 25.0 );

    // Need to set the number of online variables!
    ocp.setNOD(23);

	Expression error_contour   = dy_path_norm * (x - x_path) - dx_path_norm * (y - y_path);

	Expression error_lag       = -dx_path_norm * (x - x_path) - dy_path_norm * (y - y_path);


	ocp.minimizeLagrangeTerm(Wx*error_contour*error_contour + Wy*error_lag*error_lag + Ww*w*w +Wv*(v-0.5)*(v-0.5));// weight this with the physical cost!!!
    //ocp.subjectTo( f );
	ocp.setModel(f);

    //ocp.subjectTo( AT_END, s ==  2.5 );
    //ocp.subjectTo( AT_START, y == 0.0 );
    //ocp.subjectTo( AT_START, theta == 0.0 );
	//ocp.subjectTo( AT_START, s == 0.0 );

    ocp.subjectTo( -1 <= v <= 1.0 );
    ocp.subjectTo( -1.0 <= w <= 1.0 );


    // DEFINE A PLOT WINDOW:
    // ---------------------
    /*GnuplotWindow window;
        window.addSubplot( x ,"X"  );
        window.addSubplot( y ,"Y"  );
        window.addSubplot( theta ,"Theta"  );
        window.addSubplot( s,"V" );


    // DEFINE AN OPTIMIZATION ALGORITHM AND SOLVE THE OCP:
    // ---------------------------------------------------
	OptimizationAlgorithm algorithm(ocp);
	//RealTimeAlgorithm algorithm(ocp);
    algorithm.set( HESSIAN_APPROXIMATION, BLOCK_BFGS_UPDATE );
	algorithm.set(PRINTLEVEL, NONE);                       // default MEDIUM (NONE, MEDIUM, HIGH)
	algorithm.set(PRINT_SCP_METHOD_PROFILE, false);        // default false
	algorithm.set(PRINT_COPYRIGHT, false);                 // default true
	algorithm.set( DISCRETIZATION_TYPE, MULTIPLE_SHOOTING);
	Grid t(0,5.0,50);
	VariablesGrid s2(4,0,5.0,50),c2(2,0,5.0,50);

    algorithm.initializeDifferentialStates(s2);
    algorithm.initializeControls          (c2);
    
    algorithm.set( MAX_NUM_ITERATIONS, 100 );
    algorithm.set( KKT_TOLERANCE, 1e-8 );
    algorithm << window;
    //algorithm.solve(0.0,state_ini);
	algorithm.solve();
    VariablesGrid s3,c3;
    algorithm.getDifferentialStates(s3);
    algorithm.getControls          (c3);

	IntegratorRK45 integrator( f );

	integrator.set( INTEGRATOR_PRINTLEVEL, HIGH );
	integrator.set( INTEGRATOR_TOLERANCE, 1.0e-6 );

	// DEFINE INITIAL VALUES:
	// ----------------------

	double x_start[4] = { 0.0, 0.0 , 0.0, 0.0 };
	double u      [2] = {0.0,0.0};//c3.getFirstVector();
	double p      [1] = { 1.0      };

	double t_start    =  0.0        ;
	double t_end      =  2.5        ;


	// START THE INTEGRATION:
	// ----------------------

	//integrator.freezeAll();
	integrator.integrate( t_start, t_end, x_start, 0, p, u );


	// GET THE RESULTS
	// ---------------

	VariablesGrid differentialStates;
	integrator.getX( differentialStates );

	differentialStates.print( "x" );*/

	// DEFINE AN MPC EXPORT MODULE AND GENERATE THE CODE:
	// ----------------------------------------------------------
	OCPexport mpc( ocp );

	mpc.set( HESSIAN_APPROXIMATION,       EXACT_HESSIAN  		);
	mpc.set( DISCRETIZATION_TYPE,         MULTIPLE_SHOOTING 	);
	mpc.set( INTEGRATOR_TYPE,             INT_RK4			);
	mpc.set( NUM_INTEGRATOR_STEPS,        25            		);
	mpc.set( QP_SOLVER,                   QP_QPOASES    		);
	mpc.set( HOTSTART_QP,                 NO             		);
	mpc.set( GENERATE_TEST_FILE,          YES            		);
	mpc.set( GENERATE_MAKE_FILE,          YES            		);
	mpc.set( GENERATE_MATLAB_INTERFACE,   NO            		);
	mpc.set( SPARSE_QP_SOLUTION, 		  FULL_CONDENSING_N2	);
	mpc.set( DYNAMIC_SENSITIVITY, 		  SYMMETRIC				);
	mpc.set( CG_HARDCODE_CONSTRAINT_VALUES, NO 					);
	mpc.set( CG_USE_VARIABLE_WEIGHTING_MATRIX, YES 				);

	mpc.exportCode( "generated_mpc" );
	mpc.printDimensionsQP( );
	// ----------------------------------------------------------


    return 0;
}
