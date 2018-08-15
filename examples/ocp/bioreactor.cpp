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
#include <cmath>

/* >>> start tutorial code >>> */
int main( ){

    USING_NAMESPACE_ACADO

    // INTRODUCE THE VARIABLES:
    // -------------------------
    DifferentialState     x,y,theta, dummy, x_obst1, y_obst1, x_obst2, y_obst2, x_obst3, y_obst3, x_obst4, y_obst4;
    Control               v,w, sv;
    DifferentialEquation  f    ;

    OnlineData goal_x;
    OnlineData goal_y;
    OnlineData goal_theta;

    OnlineData wX;
    OnlineData wY;
    OnlineData wTheta;
    OnlineData wV;
    OnlineData wW;

    OnlineData wX_T;
    OnlineData wY_T;
    OnlineData wTheta_T;

    OnlineData ws;
    OnlineData wP;

    OnlineData r_disc;
    OnlineData disc_pos;

    OnlineData obst1_x;
    OnlineData obst1_y;
    OnlineData obst1_theta;
    OnlineData obst1_major;
    OnlineData obst1_minor;
    OnlineData obst1_vx;
    OnlineData obst1_vy;

    OnlineData obst2_x;
    OnlineData obst2_y;
    OnlineData obst2_theta;
    OnlineData obst2_major;
    OnlineData obst2_minor;
    OnlineData obst2_vx;
    OnlineData obst2_vy;

    OnlineData obst3_x;
    OnlineData obst3_y;
    OnlineData obst3_theta;
    OnlineData obst3_major;
    OnlineData obst3_minor;
    OnlineData obst3_vx;
    OnlineData obst3_vy;

    OnlineData obst4_x;
    OnlineData obst4_y;
    OnlineData obst4_theta;
    OnlineData obst4_major;
    OnlineData obst4_minor;
    OnlineData obst4_vx;
    OnlineData obst4_vy;

    // DEFINE A DIFFERENTIAL EQUATION:
    // -------------------------------

    f << dot(x) == v*cos(theta);
    f << dot(y) == v*sin(theta);
    f << dot(theta) == w;

    f << dot(dummy) == sv;

    f << dot(x_obst1) == obst1_vx;
    f << dot(y_obst1) == obst1_vy;
    f << dot(x_obst2) == obst2_vx;
    f << dot(y_obst2) == obst2_vy;
    f << dot(x_obst3) == obst3_vx;
    f << dot(y_obst3) == obst3_vy;
    f << dot(x_obst4) == obst4_vx;
    f << dot(y_obst4) == obst4_vy;

    // DEFINE AN OPTIMAL CONTROL PROBLEM:
    // ----------------------------------
    OCP ocp( 0.0, 10, 25.0 );

    ocp.setNOD(43);

    ocp.minimizeLagrangeTerm(wX*(x-goal_x)*(x-goal_x)+ wY*(y-goal_y)*(y-goal_y)+ wTheta*(theta-goal_theta)*(theta-goal_theta)+wV*v*v+wW*w*w + ws*sv*sv + wP*((1/((x-x_obst1)*(x-x_obst1)+(y-y_obst1)*(y-y_obst1)+0.0001)) + (1/((x-x_obst2)*(x-x_obst2)+(y-y_obst2)*(y-y_obst2)+0.0001)) + (1/((x-x_obst3)*(x-x_obst3)+(y-y_obst3)*(y-y_obst3)+0.0001)) + (1/((x-x_obst4)*(x-x_obst4)+(y-y_obst4)*(y-y_obst4)+0.0001))));  // weigh this with the physical cost!!!
//    ocp.minimizeLagrangeTerm(wX*(x-goal_x)*(x-goal_x)+ wY*(y-goal_y)*(y-goal_y)+ wTheta*(theta-goal_theta)*(theta-goal_theta)+wV*v*v+wW*w*w + ws*sv*sv + wP*((1/((x-obst1_x)*(x-obst1_x)+(y-obst1_y)*(y-obst1_y)+0.0001)) + (1/((x-obst2_x)*(x-obst2_x)+(y-obst2_y)*(y-obst2_y)+0.0001))));  // weigh this with the physical cost!!!

    ocp.minimizeMayerTerm(wX_T*(x-goal_x)*(x-goal_x)+ wY_T*(y-goal_y)*(y-goal_y)+ wTheta_T*(theta-goal_theta)*(theta-goal_theta));
    
    ocp.subjectTo( f );

    ocp.subjectTo( -0.5 <= v <= 0.5 );
    ocp.subjectTo( -0.5 <= w <= 0.5 );

    // DEFINE COLLISION CONSTRAINTS:
    // ---------------------------------------

    Expression ab_1(2,2);
    ab_1(0,0) = 1/((obst1_major + r_disc)*(obst1_major + r_disc));
    ab_1(0,1) = 0;
    ab_1(1,1) = 1/((obst1_minor + r_disc)*(obst1_minor + r_disc));
    ab_1(1,0) = 0;

    Expression ab_2(2,2);
    ab_2(0,0) = 1/((obst2_major + r_disc)*(obst2_major + r_disc));
    ab_2(0,1) = 0;
    ab_2(1,1) = 1/((obst2_minor + r_disc)*(obst2_minor + r_disc));
    ab_2(1,0) = 0;

    Expression ab_3(2,2);
    ab_3(0,0) = 1/((obst3_major + r_disc)*(obst3_major + r_disc));
    ab_3(0,1) = 0;
    ab_3(1,1) = 1/((obst3_minor + r_disc)*(obst3_minor + r_disc));
    ab_3(1,0) = 0;

    Expression ab_4(2,2);
    ab_4(0,0) = 1/((obst4_major + r_disc)*(obst4_major + r_disc));
    ab_4(0,1) = 0;
    ab_4(1,1) = 1/((obst4_minor + r_disc)*(obst4_minor + r_disc));
    ab_4(1,0) = 0;

    Expression R_obst_1(2,2);
    R_obst_1(0,0) = cos(obst1_theta);
    R_obst_1(0,1) = -sin(obst1_theta);
    R_obst_1(1,0) = sin(obst1_theta);
    R_obst_1(1,1) = cos(obst1_theta);

    Expression R_obst_2(2,2);
    R_obst_2(0,0) = cos(obst2_theta);
    R_obst_2(0,1) = -sin(obst2_theta);
    R_obst_2(1,0) = sin(obst2_theta);
    R_obst_2(1,1) = cos(obst2_theta);

    Expression R_obst_3(2,2);
    R_obst_3(0,0) = cos(obst3_theta);
    R_obst_3(0,1) = -sin(obst3_theta);
    R_obst_3(1,0) = sin(obst3_theta);
    R_obst_3(1,1) = cos(obst3_theta);

    Expression R_obst_4(2,2);
    R_obst_4(0,0) = cos(obst4_theta);
    R_obst_4(0,1) = -sin(obst4_theta);
    R_obst_4(1,0) = sin(obst4_theta);
    R_obst_4(1,1) = cos(obst4_theta);

    Expression deltaPos_disc_1(2,1);
    deltaPos_disc_1(0) =  x - x_obst1;
    deltaPos_disc_1(1) =  y - y_obst1;

    Expression deltaPos_disc_2(2,1);
    deltaPos_disc_2(0) =  x - x_obst2;
    deltaPos_disc_2(1) =  y - y_obst2;

    Expression deltaPos_disc_3(2,1);
    deltaPos_disc_3(0) =  x - x_obst3;
    deltaPos_disc_3(1) =  y - y_obst3;

    Expression deltaPos_disc_4(2,1);
    deltaPos_disc_4(0) =  x - x_obst4;
    deltaPos_disc_4(1) =  y - y_obst4;

    Expression c_obst_1, c_obst_2, c_obst_3, c_obst_4;
    c_obst_1 = deltaPos_disc_1.transpose() * R_obst_1.transpose() * ab_1 * R_obst_1 * deltaPos_disc_1;
    c_obst_2 = deltaPos_disc_2.transpose() * R_obst_2.transpose() * ab_2 * R_obst_2 * deltaPos_disc_2;
    c_obst_3 = deltaPos_disc_3.transpose() * R_obst_3.transpose() * ab_3 * R_obst_3 * deltaPos_disc_3;
    c_obst_4 = deltaPos_disc_4.transpose() * R_obst_4.transpose() * ab_4 * R_obst_4 * deltaPos_disc_4;

    ocp.subjectTo(c_obst_1 - sv >= 1);
    ocp.subjectTo(c_obst_2 - sv >= 1);
    ocp.subjectTo(c_obst_3 - sv >= 1);
    ocp.subjectTo(c_obst_4 - sv >= 1);

//    ocp.subjectTo(c_obst_1 >= 1);
//    ocp.subjectTo(c_obst_2 >= 1);
//    ocp.subjectTo(c_obst_3 >= 1);
//    ocp.subjectTo(c_obst_4 >= 1);

    // DEFINE AN MPC EXPORT MODULE AND GENERATE THE CODE:
	// ----------------------------------------------------------
	OCPexport mpc( ocp );

	mpc.set( HESSIAN_APPROXIMATION,       EXACT_HESSIAN  		);
    mpc.set( DISCRETIZATION_TYPE,         MULTIPLE_SHOOTING 	);
//	mpc.set( INTEGRATOR_TYPE,             INT_RK4			    );
    mpc.set( INTEGRATOR_TYPE,             INT_IRK_GL4			);
//    mpc.set( INTEGRATOR_TYPE,             INT_IRK_RIIA3			);
	mpc.set( NUM_INTEGRATOR_STEPS,        25            		);
	mpc.set( QP_SOLVER,                   QP_QPOASES    		);
	mpc.set( HOTSTART_QP,                 YES             		);
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
