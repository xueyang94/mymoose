#ifndef TIMEFUNCTIONCONTROL_H
#define TIMEFUNCTIONCONTROL_H

#include "THMControl.h"

class TimeFunctionControl;
class Function;

template <>
InputParameters validParams<TimeFunctionControl>();

class TimeFunctionControl : public THMControl
{
public:
  TimeFunctionControl(const InputParameters & parameters);

  virtual void execute();

protected:
  const std::string & _component_name;
  const std::string & _param_name;
  MooseObjectParameterName _ctrl_param_name;
  Function & _function;
};

#endif // TIMEFUNCTIONCONTROL_H
