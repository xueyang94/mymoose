//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "MooseTypes.h"
#include "RestartableData.h"

// Forward declarations
class PostprocessorData;
class SubProblem;
class InputParameters;
class MooseObject;
class MooseApp;
class MooseMesh;

/**
 * A class for creating restricted objects
 * \see BlockRestartable BoundaryRestartable
 */
class Restartable
{
public:
  /**
   * Class constructor
   *
   * @param moose_object The MooseObject that this interface is being implemented on.
   * @param system_name The name of the MOOSE system.  ie "Kernel", "BCs", etc.  Should roughly
   * correspond to the section in the input file so errors are easy to understand.
   *
   * This method will forward the thread id if it exists in the moose_object parameters. Delegates
   * to the "MooseApp &" constructor.
   */
  Restartable(const MooseObject * moose_object, const std::string & system_name);

  /**
   * Class constructor
   *
   * Similar to the other class constructor but also accepts an individual thread ID. If this
   * method is used, no thread ID in the parameters object is used. Delegates to the "MooseApp &"
   * constructor.
   */
  Restartable(const MooseObject * moose_object, const std::string & system_name, THREAD_ID tid);

  /**
   * This class constructor is used for non-Moose-based objects like interfaces. A name for the
   * storage as well as a system name must be passed in along with the thread ID explicitly.
   */
  Restartable(MooseApp & moose_app,
              const std::string & name,
              const std::string & system_name,
              THREAD_ID tid);

  /**
   * Emtpy destructor
   */
  virtual ~Restartable() = default;

protected:
  /**
   * Declare a piece of data as "restartable".
   * This means that in the event of a restart this piece of data
   * will be restored back to its previous value.
   *
   * NOTE: This returns a _reference_!  Make sure you store it in a _reference_!
   *
   * @param data_name The name of the data (usually just use the same name as the member variable)
   */
  template <typename T>
  T & declareRestartableData(const std::string & data_name);

  /**
   * Declare a piece of data as "restartable" and initialize it.
   * This means that in the event of a restart this piece of data
   * will be restored back to its previous value.
   *
   * NOTE: This returns a _reference_!  Make sure you store it in a _reference_!
   *
   * @param data_name The name of the data (usually just use the same name as the member variable)
   * @param init_value The initial value of the data
   */
  template <typename T>
  T & declareRestartableData(const std::string & data_name, const T & init_value);

  /**
   * Declare a piece of data as "restartable".
   * This means that in the event of a restart this piece of data
   * will be restored back to its previous value.
   *
   * NOTE: This returns a _reference_!  Make sure you store it in a _reference_!
   *
   * @param data_name The name of the data (usually just use the same name as the member variable)
   * @param context Context pointer that will be passed to the load and store functions
   */
  template <typename T>
  T & declareRestartableDataWithContext(const std::string & data_name, void * context);

  /**
   * Declare a piece of data as "restartable".
   * This means that in the event of a restart this piece of data
   * will be restored back to its previous value.
   *
   * NOTE: This returns a _reference_!  Make sure you store it in a _reference_!
   *
   * @param data_name The name of the data (usually just use the same name as the member variable)
   * @param prefix The prefix to prepend to the data_name, to retrieve data from another object.
   * @param context Context pointer that will be passed to the load and store functions
   */
  template <typename T>
  T & declareRestartableDataWithPrefixOverrideAndContext(const std::string & data_name,
                                                         const std::string & prefix,
                                                         void * context);

  /**
   * Declare a piece of data as "restartable" and initialize it.
   * This means that in the event of a restart this piece of data
   * will be restored back to its previous value.
   *
   * NOTE: This returns a _reference_!  Make sure you store it in a _reference_!
   *
   * @param data_name The name of the data (usually just use the same name as the member variable)
   * @param init_value The initial value of the data
   * @param context Context pointer that will be passed to the load and store functions
   */
  template <typename T>
  T & declareRestartableDataWithContext(const std::string & data_name,
                                        const T & init_value,
                                        void * context);

  /**
   * Declare a piece of data as "recoverable".
   * This means that in the event of a recovery this piece of data
   * will be restored back to its previous value.
   *
   * Note - this data will NOT be restored on _Restart_!
   *
   * NOTE: This returns a _reference_!  Make sure you store it in a _reference_!
   *
   * @param data_name The name of the data (usually just use the same name as the member variable)
   */
  template <typename T>
  T & declareRecoverableData(const std::string & data_name);

  /**
   * Declare a piece of data as "restartable" and initialize it.
   * This means that in the event of a restart this piece of data
   * will be restored back to its previous value.
   *
   * Note - this data will NOT be restored on _Restart_!
   *
   * NOTE: This returns a _reference_!  Make sure you store it in a _reference_!
   *
   * @param data_name The name of the data (usually just use the same name as the member variable)
   * @param init_value The initial value of the data
   */
  template <typename T>
  T & declareRecoverableData(const std::string & data_name, const T & init_value);

  /**
   * Declare a piece of data as "restartable".
   * This means that in the event of a restart this piece of data
   * will be restored back to its previous value.
   *
   * NOTE: This returns a _reference_!  Make sure you store it in a _reference_!
   *
   * @param data_name The name of the data (usually just use the same name as the member variable)
   * @param object_name A supplied name for the object that is declaring this data.
   */
  template <typename T>
  T & declareRestartableDataWithObjectName(const std::string & data_name,
                                           const std::string & object_name);

  /**
   * Declare a piece of data as "restartable".
   * This means that in the event of a restart this piece of data
   * will be restored back to its previous value.
   *
   * NOTE: This returns a _reference_!  Make sure you store it in a _reference_!
   *
   * @param data_name The name of the data (usually just use the same name as the member variable)
   * @param object_name A supplied name for the object that is declaring this data.
   * @param context Context pointer that will be passed to the load and store functions
   */
  template <typename T>
  T & declareRestartableDataWithObjectNameWithContext(const std::string & data_name,
                                                      const std::string & object_name,
                                                      void * context);

protected:
  /// Reference to the application
  MooseApp & _restartable_app;

  /// The system name this object is in
  const std::string _restartable_system_name;

  /// The thread ID for this object
  const THREAD_ID _restartable_tid;

  /// Flag for toggling read only status (see ReporterData)
  bool _restartable_read_only;

private:
  /// The name of the object
  std::string _restartable_name;

  /// Helper function for actually registering the restartable data.
  RestartableDataValue & registerRestartableDataOnApp(const std::string & name,
                                                      std::unique_ptr<RestartableDataValue> data,
                                                      THREAD_ID tid);

  /// Helper function for actually registering the restartable data.
  void registerRestartableNameWithFilterOnApp(const std::string & name,
                                              Moose::RESTARTABLE_FILTER filter);
};

template <typename T>
T &
Restartable::declareRestartableData(const std::string & data_name)
{
  return declareRestartableDataWithContext<T>(data_name, nullptr);
}

template <typename T>
T &
Restartable::declareRestartableData(const std::string & data_name, const T & init_value)
{
  return declareRestartableDataWithContext<T>(data_name, init_value, nullptr);
}

template <typename T>
T &
Restartable::declareRestartableDataWithContext(const std::string & data_name, void * context)
{
  std::string full_name = _restartable_system_name + "/" + _restartable_name + "/" + data_name;
  auto data_ptr = std::make_unique<RestartableData<T>>(full_name, context);

  // See comment in overloaded version of this function with "init_value"
  auto & restartable_data_ref = static_cast<RestartableData<T> &>(
      registerRestartableDataOnApp(full_name, std::move(data_ptr), _restartable_tid));

  return restartable_data_ref.set();
}

template <typename T>
T &
Restartable::declareRestartableDataWithContext(const std::string & data_name,
                                               const T & init_value,
                                               void * context)
{
  std::string full_name = _restartable_system_name + "/" + _restartable_name + "/" + data_name;

  // Here we will create the RestartableData even though we may not use this instance.
  // If it's already in use, the App will return a reference to the existing instance and we'll
  // return that one instead. We might refactor this to have the app create the RestartableData
  // at a later date.
  auto data_ptr = std::make_unique<RestartableData<T>>(full_name, context);
  auto & restartable_data_ref = static_cast<RestartableData<T> &>(
      registerRestartableDataOnApp(full_name, std::move(data_ptr), _restartable_tid));

  restartable_data_ref.set() = init_value;
  return restartable_data_ref.set();
}

template <typename T>
T &
Restartable::declareRestartableDataWithObjectName(const std::string & data_name,
                                                  const std::string & object_name)
{
  return declareRestartableDataWithObjectNameWithContext<T>(data_name, object_name, nullptr);
}

template <typename T>
T &
Restartable::declareRestartableDataWithObjectNameWithContext(const std::string & data_name,
                                                             const std::string & object_name,
                                                             void * context)
{
  std::string old_name = _restartable_name;

  _restartable_name = object_name;

  T & value = declareRestartableDataWithContext<T>(data_name, context);

  _restartable_name = old_name;

  return value;
}

template <typename T>
T &
Restartable::declareRecoverableData(const std::string & data_name)
{
  std::string full_name = _restartable_system_name + "/" + _restartable_name + "/" + data_name;

  registerRestartableNameWithFilterOnApp(full_name, Moose::RESTARTABLE_FILTER::RECOVERABLE);

  return declareRestartableDataWithContext<T>(data_name, nullptr);
}

template <typename T>
T &
Restartable::declareRecoverableData(const std::string & data_name, const T & init_value)
{
  std::string full_name = _restartable_system_name + "/" + _restartable_name + "/" + data_name;

  registerRestartableNameWithFilterOnApp(full_name, Moose::RESTARTABLE_FILTER::RECOVERABLE);

  return declareRestartableDataWithContext<T>(data_name, init_value, nullptr);
}
