.. _facetDragDynamicEffector:

Module: facetDragDynamicEffector
================================


Executive Summary
-----------------

Drag dynamics class used to compute drag effects on spacecraft bodies

This class is used to implement drag dynamic effects on spacecraft using a variety of simple or complex models, which will include
cannonball (attitude-independent) drag, single flat-plate drag, faceted drag models, and an interface to full-CAD GPU-accellerated
drag models.
For more information see the
:download:`PDF Description </../../src/simulation/dynamics/facetDragEffector/_Documentation/Basilisk-facet_drag-20190515.pdf>`.



Message Connection Descriptions
-------------------------------
The following table lists all the module input and output messages.  The module msg variable name is set by the
user from python.  The msg type contains a link to the message structure definition, while the description
provides information on what this message is used for.

.. list-table:: Module I/O Messages
    :widths: 25 25 50
    :header-rows: 1

    * - Msg Variable Name
      - Msg Type
      - Description
    * - atmoDensInMsg
      - :ref:`AtmoPropsMsgPayload`
      - input message for atmospheric density information



----

.. autodoxygenfile:: facetDragDynamicEffector.h
   :project: facetDragEffector

