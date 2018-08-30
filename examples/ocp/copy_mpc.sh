#!/bin/bash

cp -r ./generated_mpc ~/predictive_workspace/src/predictive_control/src/

mv ~/predictive_workspace/src/predictive_control/src/generated_mpc/acado_common.h ~/predictive_workspace/src/predictive_control/include/predictive_control/generated_mpc/
mv ~/predictive_workspace/src/predictive_control/src/generated_mpc/acado_auxiliary_functions.h ~/predictive_workspace/src/predictive_control/include/predictive_control/generated_mpc/
mv ~/predictive_workspace/src/predictive_control/src/generated_mpc/acado_qpoases_interface.hpp ~/predictive_workspace/src/predictive_control/include/predictive_control/generated_mpc/

