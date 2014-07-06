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
* Abstract: internal declarations for MSR ECClib
*
* This software is based on the article by Joppe Bos, Craig Costello, 
* Patrick Longa and Michael Naehrig, "Selecting elliptic curves for
* cryptography: an efficiency and security analysis", preprint available
* at http://eprint.iacr.org/2014/130.
***************************************************************************/  

#ifndef __MSR_ECCLIB_PRIV_H__
#define __MSR_ECCLIB_PRIV_H__

// For C++
#ifdef __cplusplus
extern "C" {
#endif


#include <stdint.h>
#include "msr_ecclib.h"
    

#define UNREFERENCED_PARAMETER(PAR)    (PAR)

    
// Additional point representations used internally:

// Types for point representations for Weierstrass a=-3 curve "Jac256" 
typedef struct { dig256 X; dig256 Y; dig256 Z; dig256 Z2; dig256 Z3; } point_chu_precomp_Jac256t;  // Point representation in Chudnovsky coordinates (X:Y:Z:Z^2:Z^3) (used for precomputed points).
typedef point_chu_precomp_Jac256t point_chu_precomp_Jac256[1];                      

// Types for point representations for twisted Edwards a=-1 curve "Ted256"                  
typedef struct { dig256 XY; dig256 YX; dig256 Z2; dig256 T2; } point_extproj_precomp_Ted256t;      // Point representation in extended homogeneous coordinates (X+Y:Y-Z:2Z:2T) (used for precomputed points).
typedef point_extproj_precomp_Ted256t point_extproj_precomp_Ted256[1]; 


// Types for point representations for Weierstrass a=-3 curve "Jac384"  
typedef struct { dig384 X; dig384 Y; dig384 Z; dig384 Z2; dig384 Z3; } point_chu_precomp_Jac384t;  // Point representation in Chudnovsky coordinates (X:Y:Z:Z^2:Z^3) (used for precomputed points).
typedef point_chu_precomp_Jac384t point_chu_precomp_Jac384[1];                       

// Types for point representations for twisted Edwards a=-1 curve "Ted384"                  
typedef struct { dig384 XY; dig384 YX; dig384 Z2; dig384 T2; } point_extproj_precomp_Ted384t;      // Point representation in extended homogeneous coordinates (X+Y:Y-Z:2Z:2T) (used for precomputed points).
typedef point_extproj_precomp_Ted384t point_extproj_precomp_Ted384[1];  


// Types for point representations for Weierstrass a=-3 curve "Jac512"   
typedef struct { dig512 X; dig512 Y; dig512 Z; dig512 Z2; dig512 Z3; } point_chu_precomp_Jac512t;  // Point representation in Chudnovsky coordinates (X:Y:Z:Z^2:Z^3) (used for precomputed points).
typedef point_chu_precomp_Jac512t point_chu_precomp_Jac512[1];                        

// Types for point representations for twisted Edwards a=-1 curve "Ted512"                  
typedef struct { dig512 XY; dig512 YX; dig512 Z2; dig512 T2; } point_extproj_precomp_Ted512t;      // Point representation in extended homogeneous coordinates (X+Y:Y-Z:2Z:2T) (used for precomputed points).
typedef point_extproj_precomp_Ted512t point_extproj_precomp_Ted512[1];  



/*************************** Function prototypes *****************************/

/********** Field functions ***********/

// Low-level field multiplication c=a*b mod p in x64 assembly
void fpmul256_a(dig256 a, dig256 b, dig256 c);
void fpmul384_a(dig384 a, dig384 b, dig384 c);
void fpmul512_a(dig512 a, dig512 b, dig512 c);

// Low-level field squaring c=a^2 mod p in x64 assembly
void fpsqr256_a(dig256 a, dig256 c);
void fpsqr384_a(dig384 a, dig384 c);
void fpsqr512_a(dig512 a, dig512 c);

// Low-level subtraction a = modulus-a, or field negation, a = -a (mod p) if modulus=p, in x64 assembly
BOOL fpneg256_a(dig256 modulus, dig256 a);
BOOL fpneg384_a(dig384 modulus, dig384 a);
BOOL fpneg512_a(dig512 modulus, dig512 a);

// Low-level field addition c = a+b mod p in x64 assembly
void fpadd256_a(dig256 a, dig256 b, dig256 c);
void fpadd384_a(dig384 a, dig384 b, dig384 c);
void fpadd512_a(dig512 a, dig512 b, dig512 c);

// Low-level field subtraction c = a-b mod p in x64 assembly
void fpsub256_a(dig256 a, dig256 b, dig256 c);
void fpsub384_a(dig384 a, dig384 b, dig384 c);
void fpsub512_a(dig512 a, dig512 b, dig512 c);

// Low-level field division by two c = a/2 mod p in x64 assembly
void fpdiv2_256_a(dig256 a, dig256 c);
void fpdiv2_384_a(dig384 a, dig384 c);
void fpdiv2_512_a(dig512 a, dig512 c);

// Low-level field zeroing in x64 assembly
void fpzero256_a(dig256 a);
void fpzero384_a(dig384 a);
void fpzero512_a(dig512 a);

// Low-level field inversion, a = a^-1 mod p (= a^(p-2) mod p)
void fpinv256_fixedchain(dig256 a);
void fpinv384_fixedchain(dig384 a);
void fpinv512_fixedchain(dig512 a);


/********** Curve functions ***********/

// SECURITY NOTE: the following functions are used for internal computations only and their correctness depends on the context 
//                in which they are used.   

// Weierstrass a=-3 curves (Jac curves)   

// "Complete" addition P = P+Q using mixed Jacobian-affine coordinates, Weierstrass a=-3 curve
void eccadd_mixed_jac_Jac256(point_Jac256 Q, point_jac_Jac256 P, point_jac_Jac256 *table, PCurveStruct JacCurve); 
void eccadd_mixed_jac_Jac384(point_Jac384 Q, point_jac_Jac384 P, point_jac_Jac384 *table, PCurveStruct JacCurve); 
void eccadd_mixed_jac_Jac512(point_Jac512 Q, point_jac_Jac512 P, point_jac_Jac512 *table, PCurveStruct JacCurve);  

// Point mixed addition P = P+Q using conditionals (if-statements), Weierstrass a=-3 curve
static void eccadd_mixed_jac_conditionals_Jac256(point_Jac256 Q, point_jac_Jac256 P, PCurveStruct JacCurve); 
static void eccadd_mixed_jac_conditionals_Jac384(point_Jac384 Q, point_jac_Jac384 P, PCurveStruct JacCurve); 
static void eccadd_mixed_jac_conditionals_Jac512(point_Jac512 Q, point_jac_Jac512 P, PCurveStruct JacCurve); 

// "Complete" addition P = P+Q, Weierstrass a=-3 curve. Table initialization is not included.
static void eccadd_jac_no_init_Jac256(point_jac_Jac256 Q, point_jac_Jac256 P, point_jac_Jac256 *table, PCurveStruct JacCurve); 
static void eccadd_jac_no_init_Jac384(point_jac_Jac384 Q, point_jac_Jac384 P, point_jac_Jac384 *table, PCurveStruct JacCurve); 
static void eccadd_jac_no_init_Jac512(point_jac_Jac512 Q, point_jac_Jac512 P, point_jac_Jac512 *table, PCurveStruct JacCurve); 

// Point doubling-addition P = 2P+Q using conditionals (if-statements), Weierstrass a=-3 curve
static void eccdoubleadd_jac_conditionals_Jac256(point_chu_precomp_Jac256 Q, point_jac_Jac256 P, PCurveStruct JacCurve);
static void eccdoubleadd_jac_conditionals_Jac384(point_chu_precomp_Jac384 Q, point_jac_Jac384 P, PCurveStruct JacCurve);
static void eccdoubleadd_jac_conditionals_Jac512(point_chu_precomp_Jac512 Q, point_jac_Jac512 P, PCurveStruct JacCurve); 

// Point doubling-addition P = 2P+Q using Jacobian coordinates, Weierstrass a=-3 curve
void eccdoubleadd_jac_Jac256(point_chu_precomp_Jac256 Q, point_jac_Jac256 P, PCurveStruct JacCurve);
void eccdoubleadd_jac_Jac384(point_chu_precomp_Jac384 Q, point_jac_Jac384 P, PCurveStruct JacCurve);
void eccdoubleadd_jac_Jac512(point_chu_precomp_Jac512 Q, point_jac_Jac512 P, PCurveStruct JacCurve);

// Special point addition R = P+Q with identical Z-coordinate for the precomputation, Weierstrass a=-3 curve
static void eccadd_jac_precomp_Jac256(point_jac_Jac256 P, point_chu_precomp_Jac256 Q, point_chu_precomp_Jac256 R);
static void eccadd_jac_precomp_Jac384(point_jac_Jac384 P, point_chu_precomp_Jac384 Q, point_chu_precomp_Jac384 R);
static void eccadd_jac_precomp_Jac512(point_jac_Jac512 P, point_chu_precomp_Jac512 Q, point_chu_precomp_Jac512 R);

// Precomputation scheme using Jacobian coordinates, Weierstrass a=-3 curve
void ecc_precomp_jac_Jac256(point_Jac256 P, point_chu_precomp_Jac256 *T, unsigned int npoints, PCurveStruct JacCurve);
void ecc_precomp_jac_Jac384(point_Jac384 P, point_chu_precomp_Jac384 *T, unsigned int npoints, PCurveStruct JacCurve);
void ecc_precomp_jac_Jac512(point_Jac512 P, point_chu_precomp_Jac512 *T, unsigned int npoints, PCurveStruct JacCurve);

// Constant-time table lookup to extract a Chudnovsky point from the precomputed table, Weierstrass a=-3 curve
void lut_chu_Jac256(point_chu_precomp_Jac256* table, point_chu_precomp_Jac256 P, int digit, unsigned int npoints, PCurveStruct JacCurve);
void lut_chu_Jac384(point_chu_precomp_Jac384* table, point_chu_precomp_Jac384 P, int digit, unsigned int npoints, PCurveStruct JacCurve);
void lut_chu_Jac512(point_chu_precomp_Jac512* table, point_chu_precomp_Jac512 P, int digit, unsigned int npoints, PCurveStruct JacCurve);

// Evaluation for the complete addition: determines the index for table lookup and the mask for element selections using complete_select_Jacxxx
unsigned int complete_eval_Jac256(dig256 val1, dig256 val2, dig256 val3, sdig *mask);
unsigned int complete_eval_Jac384(dig384 val1, dig384 val2, dig384 val3, sdig *mask);
unsigned int complete_eval_Jac512(dig512 val1, dig512 val2, dig512 val3, sdig *mask);

// Field element selection for the complete addition using mask
void complete_select_Jac256(dig256 in1, dig256 in2, dig256 out, sdig mask);
void complete_select_Jac384(dig384 in1, dig384 in2, dig384 out, sdig mask);
void complete_select_Jac512(dig512 in1, dig512 in2, dig512 out, sdig mask);

// Point extraction from 4-LUT for the complete mixed addition
void complete_lut4_Jac256(point_jac_Jac256 *table, unsigned int index, point_jac_Jac256 P);
void complete_lut4_Jac384(point_jac_Jac384 *table, unsigned int index, point_jac_Jac384 P);
void complete_lut4_Jac512(point_jac_Jac512 *table, unsigned int index, point_jac_Jac512 P);

// Point extraction from 5-LUT for the complete addition
void complete_lut5_Jac256(point_jac_Jac256 *table, unsigned int index, point_jac_Jac256 P);
void complete_lut5_Jac384(point_jac_Jac384 *table, unsigned int index, point_jac_Jac384 P);
void complete_lut5_Jac512(point_jac_Jac512 *table, unsigned int index, point_jac_Jac512 P);

// (Internal) fixed-base scalar multiplication Q = k.P, where P = P_table, using the Modified LSB-set method, Weierstrass a=-3 curve
BOOL ecc_scalar_mul_fixed_internal_Jac256(point_Jac256 *P_table, dig *k, point_Jac256 Q, unsigned int w, unsigned int v, PCurveStruct JacCurve);
BOOL ecc_scalar_mul_fixed_internal_Jac384(point_Jac384 *P_table, dig *k, point_Jac384 Q, unsigned int w, unsigned int v, PCurveStruct JacCurve);
BOOL ecc_scalar_mul_fixed_internal_Jac512(point_Jac512 *P_table, dig *k, point_Jac512 Q, unsigned int w, unsigned int v, PCurveStruct JacCurve);

// (Internal) precomputation function using affine coordinates for fixed-base scalar multiplication, Weierstrass a=-3 curve
point_Jac256* ecc_precomp_fixed_internal_Jac256(point_Jac256 P, unsigned int w, unsigned int v, unsigned int d, unsigned int e, PCurveStruct JacCurve);
point_Jac384* ecc_precomp_fixed_internal_Jac384(point_Jac384 P, unsigned int w, unsigned int v, unsigned int d, unsigned int e, PCurveStruct JacCurve);
point_Jac512* ecc_precomp_fixed_internal_Jac512(point_Jac512 P, unsigned int w, unsigned int v, unsigned int d, unsigned int e, PCurveStruct JacCurve);

// Constant-time table lookup to extract an affine point from the fixed-base precomputed table, Weierstrass a=-3 curve
void lut_aff_Jac256(point_Jac256* table, point_Jac256 P, int digit, int sign, unsigned int npoints, PCurveStruct JacCurve);
void lut_aff_Jac384(point_Jac384* table, point_Jac384 P, int digit, int sign, unsigned int npoints, PCurveStruct JacCurve);
void lut_aff_Jac512(point_Jac512* table, point_Jac512 P, int digit, int sign, unsigned int npoints, PCurveStruct JacCurve);

// (Internal) double-scalar multiplication R = k.P+l.Q, where P = P_table, using wNAF with Interleaving, Weierstrass a=-3 curve
BOOL ecc_double_scalar_mul_internal_Jac256(point_Jac256 *P_table, dig *k, point_Jac256 Q, dig *l, point_Jac256 R, unsigned int w_P, PCurveStruct JacCurve);
BOOL ecc_double_scalar_mul_internal_Jac384(point_Jac384 *P_table, dig *k, point_Jac384 Q, dig *l, point_Jac384 R, unsigned int w_P, PCurveStruct JacCurve);
BOOL ecc_double_scalar_mul_internal_Jac512(point_Jac512 *P_table, dig *k, point_Jac512 Q, dig *l, point_Jac512 R, unsigned int w_P, PCurveStruct JacCurve);

// (Internal)  precomputation function using affine coordinates for the fixed-base of double-scalar multiplication, Weierstrass a=-3 curve
point_Jac256* ecc_precomp_dblmul_internal_Jac256(point_Jac256 P, unsigned int w, PCurveStruct JacCurve);
point_Jac384* ecc_precomp_dblmul_internal_Jac384(point_Jac384 P, unsigned int w, PCurveStruct JacCurve);
point_Jac512* ecc_precomp_dblmul_internal_Jac512(point_Jac512 P, unsigned int w, PCurveStruct JacCurve);

// Additional functions based on macros

// Zero Chudnovsky point (X:Y:Z:Z^2:Z^3): P = (0:0:0:0:0)
#define ecczero_chu_Jac256(P);    fpzero256(P->X);  \
                                  fpzero256(P->Y);  \
                                  fpzero256(P->Z);  \
                                  fpzero256(P->Z2); \
                                  fpzero256(P->Z3); \

#define ecczero_chu_Jac384(P);    fpzero384(P->X);  \
                                  fpzero384(P->Y);  \
                                  fpzero384(P->Z);  \
                                  fpzero384(P->Z2); \
                                  fpzero384(P->Z3); \

#define ecczero_chu_Jac512(P);    fpzero512(P->X);  \
                                  fpzero512(P->Y);  \
                                  fpzero512(P->Z);  \
                                  fpzero512(P->Z2); \
                                  fpzero512(P->Z3); \


// Twisted Edwards a=-1 curves (Ted curves) 

// (Internal) complete point addition P = P+Q or P = P+P using extended projective coordinates (X:Y:Z:Ta:Tb) for P and (X+Y:Y-Z:2*Z:2*d*T) for Q, twisted Edwards a=-1 curve
void eccadd_extproj_internal_Ted256(point_extproj_precomp_Ted256 Q, point_extproj_Ted256 P, PCurveStruct TedCurve);
void eccadd_extproj_internal_Ted384(point_extproj_precomp_Ted384 Q, point_extproj_Ted384 P, PCurveStruct TedCurve);
void eccadd_extproj_internal_Ted512(point_extproj_precomp_Ted512 Q, point_extproj_Ted512 P, PCurveStruct TedCurve);

// Complete mixed addition P = P+Q or P = P+P using mixed extended projective-affine coordinates, twisted Edwards a=-1 curve
void eccadd_mixed_extproj_Ted256(point_extaff_precomp_Ted256 Q, point_extproj_Ted256 P, PCurveStruct TedCurve);
void eccadd_mixed_extproj_Ted384(point_extaff_precomp_Ted384 Q, point_extproj_Ted384 P, PCurveStruct TedCurve);
void eccadd_mixed_extproj_Ted512(point_extaff_precomp_Ted512 Q, point_extproj_Ted512 P, PCurveStruct TedCurve);

// Special point addition R = P+Q for the precomputation, twisted Edwards a=-1 curve
static void eccadd_extproj_precomp_Ted256(point_extproj_precomp_Ted256 P, point_extproj_precomp_Ted256 Q, point_extproj_precomp_Ted256 R, PCurveStruct TedCurve);
static void eccadd_extproj_precomp_Ted384(point_extproj_precomp_Ted384 P, point_extproj_precomp_Ted384 Q, point_extproj_precomp_Ted384 R, PCurveStruct TedCurve);
static void eccadd_extproj_precomp_Ted512(point_extproj_precomp_Ted512 P, point_extproj_precomp_Ted512 Q, point_extproj_precomp_Ted512 R, PCurveStruct TedCurve);

// Precomputation scheme using extended projective coordinates, twisted Edwards a=-1 curve
void ecc_precomp_extproj_Ted256(point_extproj_Ted256 P, point_extproj_precomp_Ted256 *T, unsigned int npoints, PCurveStruct TedCurve);
void ecc_precomp_extproj_Ted384(point_extproj_Ted384 P, point_extproj_precomp_Ted384 *T, unsigned int npoints, PCurveStruct TedCurve);
void ecc_precomp_extproj_Ted512(point_extproj_Ted512 P, point_extproj_precomp_Ted512 *T, unsigned int npoints, PCurveStruct TedCurve);

// Constant-time table lookup to extract an extended projective point from the precomputed table, twisted Edwards a=-1 curve (it does not use coordinates Ta, Tb)
void lut_extproj_Ted256(point_extproj_precomp_Ted256* table, point_extproj_precomp_Ted256 P, int digit, unsigned int npoints, PCurveStruct TedCurve);
void lut_extproj_Ted384(point_extproj_precomp_Ted384* table, point_extproj_precomp_Ted384 P, int digit, unsigned int npoints, PCurveStruct TedCurve);
void lut_extproj_Ted512(point_extproj_precomp_Ted512* table, point_extproj_precomp_Ted512 P, int digit, unsigned int npoints, PCurveStruct TedCurve);

// (Internal) fixed-base scalar multiplication Q = k.P, where P = P_table, using affine coordinates and the Modified LSB-set method, twisted Edwards a=-1 curve
BOOL ecc_scalar_mul_fixed_internal_Ted256(point_extaff_precomp_Ted256 *P_table, dig *k, point_Ted256 Q, unsigned int w, unsigned int v, PCurveStruct TedCurve);
BOOL ecc_scalar_mul_fixed_internal_Ted384(point_extaff_precomp_Ted384 *P_table, dig *k, point_Ted384 Q, unsigned int w, unsigned int v, PCurveStruct TedCurve);
BOOL ecc_scalar_mul_fixed_internal_Ted512(point_extaff_precomp_Ted512 *P_table, dig *k, point_Ted512 Q, unsigned int w, unsigned int v, PCurveStruct TedCurve);

// (Internal) precomputation function using extended affine coordinates for fixed-base scalar multiplication, twisted Edwards a=-1 curve
point_extaff_precomp_Ted256* ecc_precomp_fixed_internal_Ted256(point_Ted256 P, unsigned int w, unsigned int v, unsigned int d, unsigned int e, PCurveStruct TedCurve);
point_extaff_precomp_Ted384* ecc_precomp_fixed_internal_Ted384(point_Ted384 P, unsigned int w, unsigned int v, unsigned int d, unsigned int e, PCurveStruct TedCurve);
point_extaff_precomp_Ted512* ecc_precomp_fixed_internal_Ted512(point_Ted512 P, unsigned int w, unsigned int v, unsigned int d, unsigned int e, PCurveStruct TedCurve);

// Constant-time table lookup to extract an extended affine point from the fixed-base precomputed table, twisted Edwards a=-1 curve
void lut_extaff_Ted256(point_extaff_precomp_Ted256* table, point_extaff_precomp_Ted256 P, int digit, int sign, unsigned int npoints, PCurveStruct TedCurve);
void lut_extaff_Ted384(point_extaff_precomp_Ted384* table, point_extaff_precomp_Ted384 P, int digit, int sign, unsigned int npoints, PCurveStruct TedCurve);
void lut_extaff_Ted512(point_extaff_precomp_Ted512* table, point_extaff_precomp_Ted512 P, int digit, int sign, unsigned int npoints, PCurveStruct TedCurve);

// (Internal) double-scalar multiplication R = k.P+l.Q, where P = P_table, using wNAF with Interleaving, twisted Edwards a=-1 curve
BOOL ecc_double_scalar_mul_internal_Ted256(point_extaff_precomp_Ted256 *P_table, dig *k, point_Ted256 Q, dig *l, point_Ted256 R, unsigned int w_P, PCurveStruct TedCurve);
BOOL ecc_double_scalar_mul_internal_Ted384(point_extaff_precomp_Ted384 *P_table, dig *k, point_Ted384 Q, dig *l, point_Ted384 R, unsigned int w_P, PCurveStruct TedCurve);
BOOL ecc_double_scalar_mul_internal_Ted512(point_extaff_precomp_Ted512 *P_table, dig *k, point_Ted512 Q, dig *l, point_Ted512 R, unsigned int w_P, PCurveStruct TedCurve);

// (Internal) precomputation function using extended affine coordinates for the fixed-base of double-scalar multiplication, twisted Edwards a=-1 curve
point_extaff_precomp_Ted256* ecc_precomp_dblmul_internal_Ted256(point_Ted256 P, unsigned int w_P, PCurveStruct TedCurve);
point_extaff_precomp_Ted384* ecc_precomp_dblmul_internal_Ted384(point_Ted384 P, unsigned int w_P, PCurveStruct TedCurve);
point_extaff_precomp_Ted512* ecc_precomp_dblmul_internal_Ted512(point_Ted512 P, unsigned int w_P, PCurveStruct TedCurve);

// Additional functions based on macros

// Zero extended projective point (X+Y:Y-X:2*Z:2*d*T): P = (0:0:0:0)
#define ecczero_extproj_precomp_Ted256(P); fpzero256(P->XY);  \
                                           fpzero256(P->YX);  \
                                           fpzero256(P->Z2);  \
                                           fpzero256(P->T2);  \

#define ecczero_extproj_precomp_Ted384(P); fpzero384(P->XY);  \
                                           fpzero384(P->YX);  \
                                           fpzero384(P->Z2);  \
                                           fpzero384(P->T2);  \

#define ecczero_extproj_precomp_Ted512(P); fpzero512(P->XY);  \
                                           fpzero512(P->YX);  \
                                           fpzero512(P->Z2);  \
                                           fpzero512(P->T2);  \
                                           
// Zero extended affine point (x+y:y-x:2*d*t): P = (0:0:0)
#define ecczero_extaff_precomp_Ted256(P); fpzero256(P->xy);  \
                                          fpzero256(P->yx);  \
                                          fpzero256(P->t2);  \

#define ecczero_extaff_precomp_Ted384(P); fpzero384(P->xy);  \
                                          fpzero384(P->yx);  \
                                          fpzero384(P->t2);  \

#define ecczero_extaff_precomp_Ted512(P); fpzero512(P->xy);  \
                                          fpzero512(P->yx);  \
                                          fpzero512(P->t2);  \


// Recoding functions

void fixed_window_recode(dig *scalar, unsigned int nbit, unsigned int w, int *digits);
void mLSB_set_recode(dig *scalar, unsigned int nbit, unsigned int l, unsigned int d, int *digits);
void wNAF_recode(dig *scalar, unsigned int nbits, unsigned int w, int *digits);


#ifdef __cplusplus
}
#endif

#endif