#!/bin/bash

cp -r ./ocp ~/predictive_workspace/src/lmpcc/src/

mv ~/predictive_workspace/src/lmpcc/src/ocp/acado_common.h ~/predictive_workspace/src/lmpcc/include/lmpcc/ocp/
mv ~/predictive_workspace/src/lmpcc/src/ocp/acado_auxiliary_functions.h ~/predictive_workspace/src/lmpcc/include/lmpcc/ocp/
mv ~/predictive_workspace/src/lmpcc/src/ocp/acado_qpoases_interface.hpp ~/predictive_workspace/src/lmpcc/include/lmpcc/ocp/

