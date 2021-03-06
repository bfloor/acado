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
 *    \author Hans Joachim Ferreau, Boris Houska, Milan Vukov
 *    \date   2011-2013
 */

#include <acado_code_generation.hpp>

int main( )
{
	USING_NAMESPACE_ACADO

	// Variables:
	DifferentialState   p    ;  // the trolley position
	DifferentialState   v    ;  // the trolley velocity 
	DifferentialState   phi  ;  // the excitation angle
	DifferentialState   omega;  // the angular velocity
	Control             a    ;  // the acc. of the trolley

	const double     g = 9.81;  // the gravitational constant 
	const double     b = 0.20;  // the friction coefficient

	// Model equations:
	DifferentialEquation f; 

	f << dot( p ) == v;
	f << dot( v ) == a;
	f << dot( phi ) == omega;
	f << dot( omega ) == -g * sin(phi) - a * cos(phi) - b * omega;

	// Reference functions and weighting matrices:
	Function h, hN;
	h << p << v << phi << omega << a;
	hN << p << v << phi << omega;

	BMatrix W = eye<bool>( h.getDim() );
	BMatrix WN = eye<bool>( hN.getDim() );

	//
	// Optimal Control Problem
	//
	OCP ocp(0.0, 3.0, 10);

	ocp.subjectTo( f );

	ocp.minimizeLSQ(W, h);
	ocp.minimizeLSQEndTerm(WN, hN);

	ocp.subjectTo( -1.0 <= a <= 1.0 );
	ocp.subjectTo( -0.5 <= v <= 1.5 );

	// Export the code:
	OCPexport mpc( ocp );

	mpc.set( HESSIAN_APPROXIMATION, GAUSS_NEWTON );
	mpc.set( DISCRETIZATION_TYPE, MULTIPLE_SHOOTING );
	mpc.set( INTEGRATOR_TYPE, INT_RK4 );
	mpc.set( NUM_INTEGRATOR_STEPS, 30 );

	mpc.set( SPARSE_QP_SOLUTION, CONDENSING );
	mpc.set( QP_SOLVER, QP_QPOASES );
	
// 	mpc.set( SPARSE_QP_SOLUTION, SPARSE_SOLVER );
// 	mpc.set( QP_SOLVER, QP_QPDUNES );
	
	mpc.set( HOTSTART_QP, YES );
	mpc.set( GENERATE_TEST_FILE, NO);
	mpc.set( GENERATE_MAKE_FILE, NO );
	mpc.set( GENERATE_MATLAB_INTERFACE, NO );
	mpc.set( GENERATE_SIMULINK_INTERFACE, YES );

	// Optionally set custom module name:
	mpc.set( CG_MODULE_NAME, "nmpc" );

	if (mpc.exportCode( "crane_cl_nmpc_export" ) != SUCCESSFUL_RETURN)
		exit( EXIT_FAILURE );

	mpc.printDimensionsQP( );

	return EXIT_SUCCESS;
}
