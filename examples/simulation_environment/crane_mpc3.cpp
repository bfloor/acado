/*
 *    This file is part of ACADO Toolkit.
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
 *    \file crane_simulation.cpp
 *    \author Boris Houska, Hans Joachim Ferreau
 *    \date 2008
 */



#include <acado_toolkit.hpp>
#include <acado_gnuplot.hpp>


int main( ){

    USING_NAMESPACE_ACADO;

    // VARIABLES:
    // ------------------------
    DifferentialState        x,y,theta;   // Position of the trolley
	Control 				v,w;   // trolley accelaration

    double L = 1.0 ;              // length
	double m = 1.0 ;              // mass
	double g = 9.81;              // gravitational constant
	double b = 0.2 ;              // friction coefficient


    // DIFFERENTIAL EQUATION:
    // ------------------------
    DifferentialEquation     f, fSim;   // The model equations

    f << dot(x) ==  v*cos(theta);
	f << dot(y) ==  v*sin(theta);
    f << dot(theta) ==  w;

	L = 1.2;							// introduce model plant mismatch
	
	fSim << dot(x) ==  v*cos(theta);;
	fSim << dot(y) ==  v*sin(theta);
	fSim << dot(theta) ==  w;
	

    // DEFINE LEAST SQUARE FUNCTION:
    // -----------------------------
    Function h;

    h << x;
    h << y;
    h << theta;
	h << w;

    DMatrix Q(4,4); // LSQ coefficient matrix
    Q.setIdentity();
	Q(1,1)=10;
    DVector r(4); // Reference
	r(0)=1;
	r(1)=1;

    // DEFINE AN OPTIMAL CONTROL PROBLEM:
    // ----------------------------------
    const double t_start = 0.0;
    const double t_end   = 5.0;

    OCP ocp( t_start, t_end, 25 );

    ocp.minimizeLSQ( Q, h, r );
    ocp.subjectTo( f );
    ocp.subjectTo( -1.0 <= w <= 1.0 );
	ocp.subjectTo( 0.0 <= v <= 1.0 );

	//Collision avoidance constraints
	double obst1_major = 0.2;
	double obst1_minor = 0.2;
	double r_disc =0.5;
	double obst1_theta = 0;
	double obst1_x=0.5;
	double obst1_y =0;
	Expression ab(2,2);
	ab(0,0) = 1/((obst1_major + r_disc)*(obst1_major + r_disc));
	ab(0,1) = 0;
	ab(1,1) = 1/((obst1_minor + r_disc)*(obst1_minor + r_disc));
	ab(1,0) = 0;

	Expression R_obst(2,2);
	R_obst(0,0) = cos(obst1_theta);
	R_obst(0,1) = -sin(obst1_theta);
	R_obst(1,0) = sin(obst1_theta);
	R_obst(1,1) = cos(obst1_theta);

	Expression deltaPos_disc1(2,1);
	deltaPos_disc1(0) =  x - obst1_x;// - cos(obst1_theta)*disc_pos;
	deltaPos_disc1(1) =  y - obst1_y;// - sin(obst1_theta)*disc_pos;

//    Expression deltaPos_disc2(2,1);
//    deltaPos_disc2(0) = x - obst1_x + cos(obst1_theta)*disc_pos;
//    deltaPos_disc2(1) = y - obst1_y + sin(obst1_theta)*disc_pos;

	Expression c_obst1_1;
	c_obst1_1 = deltaPos_disc1.transpose() * R_obst.transpose() * ab * R_obst * deltaPos_disc1;
//    c_obst1_2 = deltaPos_disc2.transpose() * R_obst.transpose() * ab * R_obst * deltaPos_disc2;

//    c_obst1_1 = ((cos(obst1_theta)*(cos(obst1_theta)*(obst1_x - x + disc_pos*cos(obst1_theta)) - sin(obst1_theta)*(obst1_y - y + disc_pos*sin(obst1_theta))))/((obst1_major + r_disc)*(obst1_major + r_disc)) + (sin(obst1_theta)*(sin(obst1_theta)*(obst1_x - x + disc_pos*cos(obst1_theta)) + cos(obst1_theta)*(obst1_y - y + disc_pos*sin(obst1_theta))))/((obst1_minor + r_disc)*(obst1_minor + r_disc)))*(obst1_x - x + disc_pos*cos(obst1_theta)) + ((cos(obst1_theta)*(sin(obst1_theta)*(obst1_x - x + disc_pos*cos(obst1_theta)) + cos(obst1_theta)*(obst1_y - y + disc_pos*sin(obst1_theta))))/((obst1_minor + r_disc)*(obst1_minor + r_disc)) - (sin(obst1_theta)*(cos(obst1_theta)*(obst1_x - x + disc_pos*cos(obst1_theta)) - sin(obst1_theta)*(obst1_y - y + disc_pos*sin(obst1_theta))))/((obst1_major + r_disc)*(obst1_major + r_disc)))*(obst1_y - y + disc_pos*sin(obst1_theta));

	ocp.subjectTo(c_obst1_1 >= 1);


    // SETTING UP THE (SIMULATED) PROCESS:
    // -----------------------------------
	OutputFcn identity;
	DynamicSystem dynamicSystem( fSim,identity );

	Process process( dynamicSystem,INT_RK45 );

	//VariablesGrid disturbance; disturbance.read( "dist.txt" );
	//if (process.setProcessDisturbance( disturbance ) != SUCCESSFUL_RETURN)
	//	exit( EXIT_FAILURE );

    // SETTING UP THE MPC CONTROLLER:
    // ------------------------------

	RealTimeAlgorithm alg( ocp,0.1 );
//  	alg.set( USE_REALTIME_ITERATIONS,NO );
//  	alg.set( MAX_NUM_ITERATIONS,20 );
	VariablesGrid traj(4,t_start,t_end,25);

	for(int i=0;i<25;i++){
		traj(i,0)=2;
		traj(i,1)=2;
	}

	StaticReferenceTrajectory zeroReference(traj);

	Controller controller( alg,zeroReference );
	

    // SETTING UP THE SIMULATION ENVIRONMENT,  RUN THE EXAMPLE...
    // ----------------------------------------------------------
	SimulationEnvironment sim( 0.0,20.0,process,controller );

	DVector x0(3);
	x0.setZero();
	x0(0) = 0.0;

	if (sim.init( x0 ) != SUCCESSFUL_RETURN)
		exit( EXIT_FAILURE );
	if (sim.run( ) != SUCCESSFUL_RETURN)
		exit( EXIT_FAILURE );

    // ...AND PLOT THE RESULTS
    // ----------------------------------------------------------
	VariablesGrid diffStates;
	if (sim.getProcessDifferentialStates( diffStates ) != SUCCESSFUL_RETURN)
		exit( EXIT_FAILURE );

	VariablesGrid feedbackControl;
	if (sim.getFeedbackControl( feedbackControl ) != SUCCESSFUL_RETURN)
		exit( EXIT_FAILURE );

	GnuplotWindow window;
		window.addSubplot( diffStates(0),   "X-POSITION" );
		window.addSubplot( diffStates(1),   "Y-POSITION" );
		window.addSubplot( diffStates(2),   "THETA" );
		window.addSubplot( feedbackControl(0), "Veclocity [m/s]" );
		window.addSubplot( feedbackControl(1), "Veclocity [rad/s]" );
	window.plot();


    return EXIT_SUCCESS;
}

/* <<< end tutorial code <<< */
