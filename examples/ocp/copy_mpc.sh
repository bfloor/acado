#!/bin/bash

cp -r ./generated_mpc ~/workspace/src/predictive_control/src/

mv ~/workspace/src/predictive_control/src/generated_mpc/acado_common.h ~/workspace/src/predictive_control/include/predictive_control/generated_mpc/
mv ~/workspace/src/predictive_control/src/generated_mpc/acado_auxiliary_functions.h ~/workspace/src/predictive_control/include/predictive_control/generated_mpc/
mv ~/workspace/src/predictive_control/src/generated_mpc/acado_qpoases_interface.hpp ~/workspace/src/predictive_control/include/predictive_control/generated_mpc/

