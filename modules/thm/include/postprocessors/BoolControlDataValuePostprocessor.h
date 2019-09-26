#pragma once

#include "GeneralPostprocessor.h"
#include "ControlData.h"

class BoolControlDataValuePostprocessor;
class THMProblem;

template <>
InputParameters validParams<BoolControlDataValuePostprocessor>();

/**
 * Reads a boolean control value data and prints it out
 */
class BoolControlDataValuePostprocessor : public GeneralPostprocessor
{
public:
  BoolControlDataValuePostprocessor(const InputParameters & parameters);

  virtual void initialize();
  virtual Real getValue();
  virtual void execute();

protected:
  THMProblem * _thm_problem;
  /// The name of the control data value
  const std::string & _control_data_name;
  /// The boolean value of the control data
  const ControlData<bool> * _control_data_value;
};
