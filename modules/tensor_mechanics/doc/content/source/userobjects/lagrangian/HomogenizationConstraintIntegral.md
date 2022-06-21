# HomogenizationConstraintIntegral

!syntax description /UserObjects/HomogenizationConstraintIntegral

## Overview

This `UserObject` actually performs the volume integral over the domain to determine the
residual of the constraints imposed by the [Lagrangian kernel homogenization system](Homogenization.md).
A separate `UserObject` is required because `ScalarKernel` derived classes do not access each element
during residual and Jacobian assembly.
This user object then performs and stores the integrals for the kernel.
For more details refer to the homogenization system documentation.

The [TensorMechanics/MasterAction](/Modules/TensorMechanics/Master/index.md) can add this object
automatically, which is the recommended way to set up homogenization constraints.

## Example Input File Syntax

This example manually specifies the parameters required to initialize the object for a
large deformation, 3D example.
The `targets` parameters are the [Functions](syntax/Functions/index.md) specifying the constrained
cell-average values as a function of time.
The `constraint_types` parameters controls the type of constraint (deformation or stress) for each input.
The [homogenization system](Homogenization.md) documentation lists the order of these inputs
for each problem dimension/type.
The `large_kinematics` flag controls whether the constraints are on the deformation gradient or 
1st Piola Kirchhoff stress (`true`) or the engineering strain and stress (`false`).

!alert warning
This object should be run with `execute_on = 'initial linear'` to provide the 
updated integral values when required by the kernels.

!listing modules/tensor_mechanics/test/tests/lagrangian/total/homogenization/large-tests/3d-stress.i
         block=UserObjects

!syntax parameters /UserObjects/HomogenizationConstraintIntegral

!syntax inputs /UserObjects/HomogenizationConstraintIntegral

!syntax children /UserObjects/HomogenizationConstraintIntegral
