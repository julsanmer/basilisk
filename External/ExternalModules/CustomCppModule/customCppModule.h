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

#ifndef _CUSTOM_CPP_MODULE_H
#define _CUSTOM_CPP_MODULE_H

#include "architecture/_GeneralModuleFiles/sys_model.h"
#include "../../msgPayloadDefCpp/CustomModuleCppMsgPayload.h"
#include "../../msgPayloadDefC/CustomModuleMsgPayload.h"
#include "architecture/utilities/bskLogging.h"
#include "architecture/messaging/messaging.h"

/*! @brief basic Basilisk C++ module class */
class CustomCppModule: public SysModel {
public:
    CustomCppModule();
    ~CustomCppModule();

    void Reset(uint64_t CurrentSimNanos);
    void UpdateState(uint64_t CurrentSimNanos);

public:

    double dummy;                                   //!< [units] sample module variable declaration
    double dumVector[3];                            //!< [units] sample vector variable

    Message<CustomModuleCppMsgPayload> dataOutMsg;     //!< attitude navigation output msg
    ReadFunctor<CustomModuleMsgPayload> dataInMsg;  //!< translation navigation output msg

    BSKLogger bskLogger;              //!< -- BSK Logging

};


#endif
