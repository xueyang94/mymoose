//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "SmoothMultiBoundingBoxBaseIC.h"
#include "InitialCondition.h"
#include "MooseRandom.h"

class CircleInBoxIC : public SmoothMultiBoundingBoxBaseIC
{
public:
  static InputParameters validParams();

  CircleInBoxIC(const InputParameters & parameters);

  Real value(const Point & p);
  // Real value(const Point & p, const Point & center, const Real & radius);

protected:
  // virtual Real computeCircleValue(const Point & p, const Point & center, const Real & radius);

  MooseMesh & _mesh;

  const Real _outside;
  Real _x1;
  Real _y1;
  Real _z1;
  const Real _radius;
  const Point _center;
  Real _inside_circle;

  // enum class ProfileType
  // {
  //   COS,
  //   TANH
  // } _profile;
};
