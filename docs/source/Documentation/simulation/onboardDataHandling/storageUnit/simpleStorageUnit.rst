.. _simpleStorageUnit:

Module: simpleStorageUnit
=========================

Executive Summary
-----------------

The SimpleStorageUnit class is a model of storage unit functionality that keeps only one "partition" or buffer that all data is stored in. This module considers:

1. Integrated net input data of the attached dataNodes in one partition.
2. The storage unit's maximum storage capacity as defined by the ``storageCapacity`` attribute.

Integration of the net input data is performed with a simple Euler method for each partition.

    :math:`Data_{stored} = (baudRate) (t_{current} - t_{previous})`


Module Assumptions and Limitations
----------------------------------
See :ref:`DataStorageUnitBase` class for inherited assumption and limitations.

Message Connection Descriptions
-------------------------------
This module only uses the input and output messages of the :ref:`DataStorageUnitBase` base class.

User Guide
----------

To set up this module users must create a SimpleStorageUnit instance::

   storageUnit = simpleStorageUnit.SimpleStorageUnit()
   storageUnit.ModelTag = "storageUnit"

In addition to the variables that must be set for the :ref:`DataStorageUnitBase` base class, this module requires the ``storageCapacity`` attribute to be specified.  The total data stored in the storageUnit will be limited to not exceed this capacity value::

   storageUnit.storageCapacity = 1E5 # Given in bits

The next step is to attach one or more :ref:`DataNodeUsageMsgPayload` instances to it using the ``addDataNodeToModel()`` method::

   storageUnit.addDataNodeToModel(dataMsg)

For more information on how to set up and use this module, see the simple data system example :ref:`scenarioDataDemo`.


----

.. autodoxygenfile:: simpleStorageUnit.h
   :project: storageUnit

