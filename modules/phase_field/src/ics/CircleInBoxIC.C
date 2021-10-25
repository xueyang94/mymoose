//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CircleInBoxIC.h"
#include "MooseMesh.h"

registerMooseObject("PhaseFieldApp", CircleInBoxIC);

InputParameters
CircleInBoxIC::validParams()
{
  InputParameters params = SmoothMultiBoundingBoxBaseIC::validParams();
  params.addClassDescription("Specify variable values inside a list of nested boxes shaped "
                             "axis-aligned regions defined by pairs of opposing corners,"
                             "and an inner circle defined by center coordinates and radius");
  params.addParam<Real>("outside", 0.0, "The value of the variable outside the largest boxes");
  params.addRequiredParam<Real>("x1", "The x coordinate of the circle center");
  params.addRequiredParam<Real>("y1", "The y coordinate of the circle center");
  params.addParam<Real>("z1", 0.0, "The z coordinate of the circle center");
  params.addRequiredParam<Real>("radius", "The radius of a circle");
  params.addRequiredParam<Real>("value_inside_circle", "The value inside the inner circle");
  return params;
}

CircleInBoxIC::CircleInBoxIC(const InputParameters & parameters)
  : SmoothMultiBoundingBoxBaseIC(parameters),
    _mesh(_fe_problem.mesh()),
    _outside(getParam<Real>("outside")),
    _x1(parameters.get<Real>("x1")),
    _y1(parameters.get<Real>("y1")),
    _z1(parameters.get<Real>("z1")),
    _radius(parameters.get<Real>("radius")),
    _center(_x1, _y1, _z1),
    _inside_circle(parameters.get<Real>("value_inside_circle"))
{
}

Real
CircleInBoxIC::value(const Point & p)
{
  Real value = SmoothMultiBoundingBoxBaseIC::value(p);

  // Compute the distance between the current point and the center
  Real dist = _mesh.minPeriodicDistance(_var.number(), p, _center);
  Real _outvalue = _inside[0]; // the value outside the circle inside the first nested box
  // Real _invalue = _inside_circle;

  if (dist <= _radius - _int_width / 2.0) // Inside circle
  {
    value = _inside_circle;
  }
  else if (dist < _radius + _int_width / 2.0) // Smooth interface
  {
    Real int_pos = (dist - _radius + _int_width / 2.0) / _int_width;
    value =
        _outvalue + (_inside_circle - _outvalue) * (1.0 + std::cos(int_pos * libMesh::pi)) / 2.0;
  }
  else
  {
    if (_int_width != 0.0)
    {
      for (unsigned int b = 0; b < _nbox; ++b)
      {
        for (unsigned int i = 0; i < _dim; ++i)
        {
          if (_c1[b](i) < _c2[b](i) && p(i) >= _c1[b](i) - _int_width &&
              p(i) <= _c2[b](i) + _int_width)
          {
            if (i != _dim - 1)
              continue;
            Real f_in = 1.0;
            for (unsigned int j = 0; j < _dim; ++j)
              f_in *= 0.5 * (std::tanh(2.0 * libMesh::pi * (p(j) - _c1[b](j)) / _int_width) -
                             std::tanh(2.0 * libMesh::pi * (p(j) - _c2[b](j)) / _int_width));
            if (b == _nbox - 1)
              value = _outside + (_inside[b] - _outside) * f_in;
            else
              value = _inside[b + 1] + (_inside[b] - _inside[b + 1]) * f_in;
            goto label;
          }
          else if (_c1[b](i) >= _c2[b](i))
            mooseError("The coordinates of the smaller_coordinate_corners are equal to or larger "
                       "than that of "
                       "the larger_coordinate_corners.");
          else
            break;
        }
      }
    }
    // label:
    //   return value;
  }
label:
  return value;
}

// Real
// CircleInBoxIC::computeCircleValue(const Point & p, const Point & center, const Real & radius)
// {
//   Point l_center = center;
//   Point l_p = p;
//   // Compute the distance between the current point and the center
//   Real dist = _mesh.minPeriodicDistance(_var.number(), l_p, l_center);
//   Real _outvalue = _inside[0]; // the value outside the circle inside the first nested box
//   Real _invalue = _inside_circle;
//
//   switch (_profile)
//   {
//     case ProfileType::COS:
//     {
//       // Return value
//       Real value = _inside[0]; // Outside circle
//
//       if (dist <= radius - _int_width / 2.0) // Inside circle
//         value = _inside_circle;
//       else if (dist < radius + _int_width / 2.0) // Smooth interface
//       {
//         Real int_pos = (dist - radius + _int_width / 2.0) / _int_width;
//         value = _outvalue +
//                 (_inside_circle - _outvalue) * (1.0 + std::cos(int_pos * libMesh::pi)) / 2.0;
//       }
//       return value;
//     }
//
//       // case ProfileType::TANH:
//       return (_invalue - _outvalue) * 0.5 * (std::tanh(2.0 * (radius - dist) / _int_width) + 1.0)
//       +
//              _outvalue;
//
//     default:
//       mooseError("Internal error.");
//   }
// }
