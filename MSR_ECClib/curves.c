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
* Abstract: supported elliptic curves
*
* This software is based on the article by Joppe Bos, Craig Costello, 
* Patrick Longa and Michael Naehrig, "Selecting elliptic curves for
* cryptography: an efficiency and security analysis", preprint available
* at http://eprint.iacr.org/2014/130.
***************************************************************************/  

#include "msr_ecclib.h"


//
// "Jac256": Weierstrass curve a=-3, E: y^2 = x^3 - 3x + 152961, p = 2^256-189
//
CurveStruct curve_Jac256 = {
     // Curve ID, 2 x targeted security level, order bitlength, prime bitlength
     Jac256, 256, 256, 256,
     // Prime p = 2^256-189
     {0xFFFFFFFFFFFFFF43, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF},
     // Parameter "a"
     {0xFFFFFFFFFFFFFF40, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF},
     // Parameter "b"
     {0x25581},
     // Order of the group
     {0x20AB20294751A825, 0xE43C8275EA265C60, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF},
     // x(generator)
     {0x01},
     // y(generator)
     {0x0BF46306C2B56C77, 0xD02C2F9375894EC1, 0xFC82C96CCEEEDD6B, 0x696F1853C1E466D7},
     // co-factor
     1 };


//
// "Ted256": twisted Edwards curve a=-1, E: -x^2 + y^2 = 1 + 15342x^2y^2, p = 2^256-189
//
CurveStruct curve_Ted256 = {
     // Curve ID, 2 x targeted security level, order bitlength, prime bitlength
     Ted256, 256, 254, 256,
     // Prime p = 2^256-189
     {0xFFFFFFFFFFFFFF43, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF},
     // Parameter "a"
     {0xFFFFFFFFFFFFFF42, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF},
     // Parameter "d"
     {0x3BEE},
     // Order of the subgroup
     {0xE5B84E6F1122B4AD, 0xBE6AA55AD0A6BC64, 0xFFFFFFFFFFFFFFFF, 0x3FFFFFFFFFFFFFFF},
     // x(generator)
     {0x0D},
     // y(generator)
     {0x707F6FB5331CADBA, 0xBE2A6D63824D303F, 0xA3D330B39FA046BF, 0x7D0AB41E2A1276DB},
     // co-factor
     4 };


//
// "Jac384": Weierstrass curve a=-3, E: y^2 = x^3 - 3x - 34568, p = 2^384-317
//
CurveStruct curve_Jac384 = {
     // Curve ID, 2 x targeted security level, order bitlength, prime bitlength
     Jac384, 384, 384, 384,
     // Prime p = 2^384-317
     {0xFFFFFFFFFFFFFEC3, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF},
     // Parameter "a"
     {0xFFFFFFFFFFFFFEC0, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF},
     // Parameter "b"
     {0xFFFFFFFFFFFF77BB, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF},    
     // Order of the group
     {0x604D81F67B0E61B9, 0xBEDA9D3D4C37E27A, 0xD61EAF1EEB5D6881, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF},
     // x(generator)
     {0x02},     
     // y(generator)
     {0x5BCBBF503EA66F43, 0xB4A3D0BAD5D30F86, 0x73330DC58D15FFA2, 0x8034ED422F04F826, 0x71E763E0663E5DBD, 0x3C9F82CB4B87B4DC}, 
     // co-factor
     1 };


//
// "Ted384": twisted Edwards curve a=-1, E: -x^2 + y^2 = 1 + 333194x^2y^2, p = 2^384-317
//
CurveStruct curve_Ted384 = {
     // Curve ID, 2 x targeted security level, order bitlength, prime bitlength
     Ted384, 384, 382, 384,
     // Prime p = 2^2^384-317
     {0xFFFFFFFFFFFFFEC3, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF},
     // Parameter "a"
     {0xFFFFFFFFFFFFFEC2, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF},
     // Parameter "d"
     {0x5158A},     
     // Order of the subgroup
     {0x51D6D71F70426E25, 0x5A13A0458E39F4E4, 0xECD7D11ED5A259A2, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0x3FFFFFFFFFFFFFFF}, 
     // x(generator)
     {0x08},     
     // y(generator)
     {0x34A4ACF5FA7F5BEE, 0x2F003860FEBABAD5, 0xBD798A8AE753C6D7, 0xAA5C7B4C930BFF8E, 0x5BD4471794AA619D, 0x749CDABA136CE9B6}, 
     // co-factor
     4 };


//
// "Jac512": Weierstrass curve a=-3, E: y^2 = x^3 - 3x + 121243, p = 2^512-569
//
CurveStruct curve_Jac512 = {
     // Curve ID, 2 x targeted security level, order bitlength, prime bitlength
     Jac512, 512, 512, 512,
     // Prime p = 2^512-569
     {0xFFFFFFFFFFFFFDC7, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF}, 
     // Parameter "a"
     {0xFFFFFFFFFFFFFDC4, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF}, 
     // Parameter "b"
     {0x1D99B},
     // Order of the group
     {0xCE153F390433555D, 0x3B568B36607CD243, 0x4FC258ED97D0BDC6, 0x5B3CA4FB94E7831B, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF},
     // x(generator)
     {0x02},              
     // y(generator)
     {0x205CC30D7F83CF28, 0x00BB3D08DC83755B, 0x724DDE537C2B0ADB, 0xCC5A1C1F0C716FDC, 0xE0B9328433AFBDD8, 0xFCC13031CF6DD336, 0x1952C250EA61AD53, 0x1C282EB23327F971},              
     // co-factor
     1 };


//
// "Ted512": twisted Edwards curve a=-1, E: -x^2 + y^2 = 1 + 637608x^2y^2, p = 2^512-569
//
CurveStruct curve_Ted512 = {
     // Curve ID, , 2 x targeted security level, order bitlength, prime bitlength
     Ted512, 512, 510, 512,
     // Prime p = 2^512-569
     {0xFFFFFFFFFFFFFDC7, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF}, 
     // Parameter "a"
     {0xFFFFFFFFFFFFFDC6, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF}, 
     // Parameter "d"
     {0x9BAA8}, 
     // Order of the subgroup
     {0x4E78D1CB0B5F0189, 0xDCEA5FF0CB800F89, 0xA624784F449545F0, 0xA7E50809EFDABBB9, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0x3FFFFFFFFFFFFFFF}, 
     // x(generator)
     {0x20},              
     // y(generator)
     {0x0CB91027543B1C5E, 0x9261134638750F4F, 0xD96B7A897C1D7279, 0x62AE2ECE5057B5DA, 0x9FF048779E1D614E, 0x9CEB124BF726973F, 0x605091D80869212F, 0x7D67E841DC4C467B},              
     // co-factor
     4 };