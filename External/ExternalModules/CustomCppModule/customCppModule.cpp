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
#include "customCppModule.h"
#include <iostream>
#include <cstring>
#include "architecture/utilities/avsEigenSupport.h"
#include "architecture/utilities/linearAlgebra.h"

/*! This is the constructor for the module class.  It sets default variable
    values and initializes the various parts of the model */
CustomCppModule::CustomCppModule()
{
}

/*! Module Destructor.  */
CustomCppModule::~CustomCppModule()
{
    return;
}


/*! This method is used to reset the module.
 @return void
 */
void CustomCppModule::Reset(uint64_t CurrentSimNanos)
{
    /*! - reset any required variables */
    this->dummy = 0.0;
    bskLogger.bskLog(BSK_INFORMATION, "Variable dummy set to %f in reset.",this->dummy);
}


void CustomCppModule::UpdateState(uint64_t CurrentSimNanos)
{
    double Lr[3];                                   /*!< [unit] variable description */
    CustomModuleCppMsgPayload outMsgBuffer;       /*!< local output message copy */
    CustomModuleMsgPayload inMsgBuffer;        /*!< local copy of input message */
    double  inputVector[3];

    // always zero the output buffer first
    outMsgBuffer = this->dataOutMsg.zeroMsgPayload;
    v3SetZero(inputVector);

    /*! - Read the input messages */
    if (this->dataInMsg.isLinked()) {
        inMsgBuffer = this->dataInMsg();
        v3Copy(inMsgBuffer.dataVector, inputVector);
    }

    /*! - Add the module specific code */
    v3Copy(inputVector, Lr);
    this->dummy += 1.0;
    Lr[0] += this->dummy;

    /*! - store the output message */
    v3Copy(Lr, outMsgBuffer.dataVector);

    /*! - write the module output message */
    this->dataOutMsg.write(&outMsgBuffer, moduleID, CurrentSimNanos);

    /* this logging statement is not typically required.  It is done here to see in the
     quick-start guide which module is being executed */
    bskLogger.bskLog(BSK_INFORMATION, "C++ Module ID %lld ran Update at %fs", moduleID, (double) CurrentSimNanos/(1e9));

}
