/*
 ISC License

 Copyright (c) 2016, Autonomous Vehicle Systems Lab, University of Colorado at Boulder

 Permission to use, copy, modify, and/or distribute this software for any
 purpose with or without fee is hereby granted, provided that the above
 copyright notice and this permission notice appear in all copies.

 THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

 */

%module simpleStorageUnit
%{
#include "simpleStorageUnit.h"
%}


%pythoncode %{
from Basilisk.architecture.swig_common_model import *
%}
%include "std_string.i"
%include "swig_conly_data.i"
%include "swig_eigen.i"
%include "std_vector.i"
%include "sys_model.i"
%include "stdint.i"

//When using scientific notation in Python (1E9), it is interpreted as float
// giving a type error when assigning storageCapacity or using setDataBuffer.
// This maps that float to long int in C++ in this module.
%typemap(in) long long int {
    $1 = static_cast<long long int>(PyFloat_AsDouble($input));
}

%include "simulation/onboardDataHandling/_GeneralModuleFiles/dataStorageUnitBase.h"
%include "simpleStorageUnit.h"
%include "architecture/msgPayloadDefC/DataNodeUsageMsgPayload.h"
struct DataNodeUsageMsg_C;
%include "architecture/msgPayloadDefCpp/DataStorageStatusMsgPayload.h"

%pythoncode %{
import sys
protectAllClasses(sys.modules[__name__])
%}
