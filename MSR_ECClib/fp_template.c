/**************************************************************************
* MSR ECClib, an efficient and secure elliptic curve cryptographic library
*
*    Copyright (c) Microsoft Corporation. All rights reserved.
*    Licensed under the Apache License, Version 2.0 (the "License"); 
*    you may not use these files except in compliance with the License. 
*    You may obtain a copy of the License at
*                http://www.apache.org/licenses/LICENSE-2.0
*    Unless required by applicable law or agreed to in writing, software 
*    distributed under the License is distributed on an "AS IS" BASIS, 
*    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or 
*    implied. See the License for the specific language governing 
*    permissions and limitations under the License.
*
*
* Abstract: template for field functions
*
* This software is based on the article by Joppe Bos, Craig Costello, 
* Patrick Longa and Michael Naehrig, "Selecting elliptic curves for
* cryptography: an efficiency and security analysis", preprint available
* at http://eprint.iacr.org/2014/130.
***************************************************************************/  

#include "msr_ecclib.h"
#include "msr_ecclib_priv.h"


void FP_MUL(FELM a, FELM b, FELM c)
{ // Modular multiplication, c=a*b mod p

    FP_MUL_LOW(a, b, c);
    return;
}


void FP_SQR(FELM a, FELM c)
{ // Modular squaring, c=a^2 mod p

    FP_SQR_LOW(a, c);
    return;
}


void FP_ADD(FELM a, FELM b, FELM c)
{ // Modular addition, c=a+b mod p

    FP_ADD_LOW(a, b, c);
    return;
}


void FP_SUB(FELM a, FELM b, FELM c)
{ // Modular subtraction, c=a-b mod p

    FP_SUB_LOW(a, b, c);
    return;
}


void FP_DIV2(FELM a, FELM c)
{ // Modular division by two, c=a/2 mod p

    FP_DIV2_LOW(a, c);
    return;
}


BOOL FP_NEG(FELM modulus, FELM a)
{ // Subtraction, a=modulus-a
  // If modulus=p then it performs a modular negation a=-a mod p
  // If a <= modulus then eval = 1 (TRUE)
    BOOL eval = FALSE;

    eval = FP_NEG_LOW(modulus, a);
    return eval;
}


void FP_MULC(FELM a, FELM b, FELM c)
{ // Modular multiplication by a single-word multiplier, c=a*b mod p, where b is a single-word constant

    FP_MULC_LOW(a, b, c);
    return;
}


void FP_ZERO(FELM a)
{ // Zero a field element, a=0

    FP_ZERO_LOW(a);
    return;
}


void FP_INV(FELM a)
{ // Modular inversion, a=a^-1 mod p

    FP_INV_LOW(a);
    return;
}