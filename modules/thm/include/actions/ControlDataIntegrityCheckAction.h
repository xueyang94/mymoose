#ifndef CONTROLDATAINTEGRITYCHECKACTION_H
#define CONTROLDATAINTEGRITYCHECKACTION_H

#include "THMAction.h"

class ControlDataIntegrityCheckAction;

template <>
InputParameters validParams<ControlDataIntegrityCheckAction>();

/**
 * Action to trigger the check of control data integrity
 */
class ControlDataIntegrityCheckAction : public THMAction
{
public:
  ControlDataIntegrityCheckAction(InputParameters parameters);

  virtual void act();

protected:
};

#endif /* CONTROLDATAINTEGRITYCHECKACTION_H */
