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
* Abstract: main header file
*
* This software is based on the article by Joppe Bos, Craig Costello, 
* Patrick Longa and Michael Naehrig, "Selecting elliptic curves for
* cryptography: an efficiency and security analysis", preprint available
* at http://eprint.iacr.org/2014/130.
***************************************************************************/  

#ifndef __MSR_ECCLIB_H__
#define __MSR_ECCLIB_H__

// For C++
#ifdef __cplusplus
extern "C" {
#endif


#include <stdint.h>


// Definition of compiler

#define COMPILER_VC      1
#define COMPILER_GCC     2

#if defined(_MSC_VER)      // Microsoft Visual C compiler
    #define COMPILER COMPILER_VC
#elif defined(__GNUC__)    // GNU GCC compiler
    #define COMPILER COMPILER_GCC
    #error -- "Unsupported COMPILER"    
#else
    #error -- "Unknown COMPILER"
#endif


// Definition of the targeted architecture and basic data types. "dig" and "sdig" 
// are data types for representing unsigned and signed computer words, respectively.
    
#define TARGET_AMD64     1

#if defined(_AMD64_)
    #define TARGET TARGET_AMD64
    #define RADIX  64
    typedef uint64_t dig;       // unsigned 64-bit digit
    typedef int64_t  sdig;      // signed 64-bit digit
#else
    #error -- "Unsupported ARCHITECTURE"
#endif
    

// Definition of BOOL type

typedef int BOOL; 
#ifndef TRUE
    #define TRUE 1
#endif
#ifndef FALSE
    #define FALSE 0
#endif


// Select supported bitlengths (security levels)

#define ECCURVES_256
#define ECCURVES_384
#define ECCURVES_512


// Definition of memory footprint type for precomputed tables in the scalar multiplication

typedef enum {
    MEM_DEFAULT,
    MEM_COMPACT,
    MEM_LARGE,
    MEM_END_OF_LIST    // List size indicator
} MemType;

#define MemTypeSize (MEM_END_OF_LIST)


// Window width for variable-base precomputed table in variable-base scalar multiplication and double-scalar multiplication

#define W_VARBASE         6    // Memory requirement: 2.5KB for Jac256, 3.75KB for Jac384 and 5KB for Jac512 
                               //                     2KB for Ted256, 3KB for Ted384 and 4KB for Ted512 


// Window width for fixed-base precomputed table in fixed-base scalar multiplication

#define W_MEM_LARGE       6    // Memory requirement: 6KB for Jac256, 9KB for Jac384 and 12KB for Jac512 
#define V_MEM_LARGE       3    //                     9KB for Ted256, 13.5KB for Ted384 and 18KB for Ted512 

#define W_MEM_COMPACT_W   3    // Memory requirement: 768 bytes for Jac256, 1.125KB for Jac384 and 1.5KB for Jac512 
#define V_MEM_COMPACT_W   3    //                       

#define W_MEM_COMPACT_TE  2    // Memory requirement: 384 bytes for Ted256, 576 bytes for Ted384 and 768 bytes for Ted512 
#define V_MEM_COMPACT_TE  2    //                      


// Window width for fixed-base precomputed table in double-scalar multiplication

#define W_P_MEM_LARGE     7    // Memory requirement: 2KB for Jac256, 3KB for Jac384 and 4KB for Jac512 
                               //                     3KB for Ted256, 4.5KB for Ted384 and 6KB for Ted512 

#define W_P_MEM_COMPACT   2    // Memory requirement: 64 bytes for Jac256, 96 bytes for Jac384 and 128 bytes for Jac512 
                               //                     96 bytes for Ted256, 128 bytes for Ted384 and 192 bytes for Ted512 


// Basic parameters for supported curves

#define MAXBITLENGTH    512                                 // Max. field bitlength supported
#define ML_WORD         (sizeof(dig)*8)                     // Number of bits in a computer word
#define MAXWORDS_FIELD  ((MAXBITLENGTH+ML_WORD-1)/ML_WORD)  // Max. number of words needed to represent supported fields
#define MAXPOINTS       64                                  // Max. number of precomputed points for variable-base scalar multiplication
#define WMAX            8                                   // Max. window size (number of rows) for fixed-base scalar multiplication
#define VMAX            8                                   // Max. table size (number of tables) for fixed-base scalar multiplication

#if RADIX >= 32
    //#define ML_COUNT                                        // Defined when operation counting is required (used for debugging purposes)
#endif
#ifdef ML_COUNT
    // Variables storing operation counts (used for debugging purposes):
    // number of inversions, multiplications, adds/subs/div2/negations, table lookups, two-field element selections (in the complete addition), 4-point table lookups (in the complete addition)
    int ninv, nmul, nsqr, nadd, nlut, ncsel, nclut;       
#endif


// Definition of data types for multi-precision field elements

#define ML_WORDS256 ((256+ML_WORD-1)/ML_WORD)               // Number of words to represent 256-bit field elements
#define ML_WORDS384 ((384+ML_WORD-1)/ML_WORD)               // Number of words to represent 384-bit field elements
#define ML_WORDS512 ((512+ML_WORD-1)/ML_WORD)               // Number of words to represent 512-bit field elements

typedef dig dig256[ML_WORDS256];                            // Multiprecision type to represent 256-bit field elements or elements in Z_r
typedef dig dig384[ML_WORDS384];                            // Multiprecision type to represent 384-bit field elements or elements in Z_r
typedef dig dig512[ML_WORDS512];                            // Multiprecision type to represent 512-bit field elements or elements in Z_r


// Types for point representations for Weierstrass a=-3 curve "Jac256"  
typedef struct { dig256 X; dig256 Y; dig256 Z; } point_jac_Jac256t;                // Point representation in Jacobian coordinates (X:Y:Z) such that x = X/Z^2, y = Y/Z^3.
typedef point_jac_Jac256t point_jac_Jac256[1];                                                              
typedef struct { dig256 x; dig256 y; } point_Jac256t;                              // Point representation in affine coordinates (x,y).
typedef point_Jac256t point_Jac256[1];                                               

// Types for point representations for twisted Edwards a=-1 curve "Ted256"
typedef struct { dig256 X; dig256 Y; dig256 Z; dig256 Ta; dig256 Tb; } point_extproj_Ted256t; // Point representation in extended homogeneous coordinates (X:Y:Z:Ta:Tb) such that
typedef point_extproj_Ted256t point_extproj_Ted256[1];                                        // x = X/Z, y = Y/Z, T = Ta*Tb = X*Y/Z.
typedef struct { dig256 xy; dig256 yx; dig256 t2; } point_extaff_precomp_Ted256t;             // Point representation in extended affine coordinates (x+y,y-x,2dt) (used for precomputed points).
typedef point_extaff_precomp_Ted256t point_extaff_precomp_Ted256[1];  
typedef struct { dig256 x; dig256 y; } point_Ted256t;                                         // Point representation in affine coordinates (x,y)
typedef point_Ted256t point_Ted256[1]; 


// Types for point representations for Weierstrass a=-3 curve "Jac384" 
typedef struct { dig384 X; dig384 Y; dig384 Z; } point_jac_Jac384t;                // Point representation in Jacobian coordinates (X:Y:Z) such that x = X/Z^2, y = Y/Z^3.
typedef point_jac_Jac384t point_jac_Jac384[1];                                                                 
typedef struct { dig384 x; dig384 y; } point_Jac384t;                              // Point representation in affine coordinates (x,y).
typedef point_Jac384t point_Jac384[1];                                               

// Types for point representations for twisted Edwards a=-1 curve "Ted384"
typedef struct { dig384 X; dig384 Y; dig384 Z; dig384 Ta; dig384 Tb; } point_extproj_Ted384t; // Point representation in extended homogeneous coordinates (X:Y:Z:Ta:Tb) such that
typedef point_extproj_Ted384t point_extproj_Ted384[1];                                        // x = X/Z, y = Y/Z, T = Ta*Tb = X*Y/Z.
typedef struct { dig384 xy; dig384 yx; dig384 t2; } point_extaff_precomp_Ted384t;             // Point representation in extended affine coordinates (x+y,y-x,2dt) (used for precomputed points).
typedef point_extaff_precomp_Ted384t point_extaff_precomp_Ted384[1];  
typedef struct { dig384 x; dig384 y; } point_Ted384t;                                         // Point representation in affine coordinates (x,y)
typedef point_Ted384t point_Ted384[1]; 


// Types for point representations for Weierstrass a=-3 curve "Jac512"    
typedef struct { dig512 X; dig512 Y; dig512 Z; } point_jac_Jac512t;                // Point representation in Jacobian coordinates (X:Y:Z) such that x = X/Z^2, y = Y/Z^3.
typedef point_jac_Jac512t point_jac_Jac512[1];                                                                    
typedef struct { dig512 x; dig512 y; } point_Jac512t;                              // Point representation in affine coordinates (x,y).
typedef point_Jac512t point_Jac512[1];                                               

// Types for point representations for twisted Edwards a=-1 curve "Ted512"
typedef struct { dig512 X; dig512 Y; dig512 Z; dig512 Ta; dig512 Tb; } point_extproj_Ted512t; // Point representation in extended homogeneous coordinates (X:Y:Z:Ta:Tb) such that
typedef point_extproj_Ted512t point_extproj_Ted512[1];                                        // x = X/Z, y = Y/Z, T = Ta*Tb = X*Y/Z.
typedef struct { dig512 xy; dig512 yx; dig512 t2; } point_extaff_precomp_Ted512t;             // Point representation in extended affine coordinates (x+y,y-x,2dt) (used for precomputed points).
typedef point_extaff_precomp_Ted512t point_extaff_precomp_Ted512[1];  
typedef struct { dig512 x; dig512 y; } point_Ted512t;                                         // Point representation in affine coordinates (x,y)
typedef point_Ted512t point_Ted512[1]; 


// Elliptic curve structure

typedef struct
{    
    unsigned int     Curve;                             // Curve ID, curve defined over GF(prime)  
    unsigned int     nbits;                             // Two times the targeted security level 
    unsigned int     rbits;                             // Bitlength of the order of the curve (sub)group 
    unsigned int     pbits;                             // Bitlength of the prime 
    dig              prime[MAXWORDS_FIELD];             // Prime
    dig              parameter1[MAXWORDS_FIELD];        // Curve parameter ("a" for Weierstrass, "a" for twisted Edwards)
    dig              parameter2[MAXWORDS_FIELD];        // Curve parameter ("b" for Weierstrass, "d" for twisted Edwards)
    dig              order[MAXWORDS_FIELD];             // Prime order of the curve (sub)group 
    dig              generator_x[MAXWORDS_FIELD];       // x-coordinate of generator
    dig              generator_y[MAXWORDS_FIELD];       // y-coordinate of generator
    unsigned int     cofactor;                          // Co-factor of the curve group
} CurveStruct, *PCurveStruct;


// Curve IDs of supported curves
typedef enum 
{
    Jac256,        // Weierstrass curves defined over a pseudo-Mersenne prime
    Jac384,
    Jac512,
    Ted256,        // Twisted Edwards curves defined over a pseudo-Mersenne prime
    Ted384,
    Ted512
} Curve_ID;


// Supported curves:

// "Jac256": Weierstrass curve a=-3, E: y^2 = x^3 - 3x + 152961, p = 2^256-189
extern CurveStruct curve_Jac256; 

// "Ted256": twisted Edwards curve a=-1, E: -x^2 + y^2 = 1 + 15342x^2y^2, p = 2^256-189
extern CurveStruct curve_Ted256;

// "Jac384": Weierstrass curve a=-3, E: y^2 = x^3 - 3x - 34568, p = 2^384-317
extern CurveStruct curve_Jac384;

// "Ted384": twisted Edwards curve a=-1, E: -x^2 + y^2 = 1 + 333194x^2y^2, p = 2^384-317
extern CurveStruct curve_Ted384;

// "Jac512": Weierstrass curve a=-3, E: y^2 = x^3 - 3x + 121243, p = 2^512-569
extern CurveStruct curve_Jac512;

// "Ted512": twisted Edwards curve a=-1, E: -x^2 + y^2 = 1 + 637608x^2y^2, p = 2^512-569
extern CurveStruct curve_Ted512;



/*************************** Function prototypes *****************************/

/********** Field functions ***********/

// Copy of a field element, c = a 
void fpcopy256(dig256 a, dig256 c);
void fpcopy384(dig384 a, dig384 c);
void fpcopy512(dig512 a, dig512 c);

// Zero a field element, a=0 
void fpzero256(dig256 a);
void fpzero384(dig384 a);
void fpzero512(dig512 a);

// Is field element zero, a=0?  
BOOL fp_iszero256(dig256 a);
BOOL fp_iszero384(dig384 a);
BOOL fp_iszero512(dig512 a);

// Field multiplication c=a*b mod p
void fpmul256(dig256 a, dig256 b, dig256 c);
void fpmul384(dig384 a, dig384 b, dig384 c);
void fpmul512(dig512 a, dig512 b, dig512 c);

// Field squaring c=a^2 mod p
void fpsqr256(dig256 a, dig256 c);
void fpsqr384(dig384 a, dig384 c);
void fpsqr512(dig512 a, dig512 c);

// Field inversion, a = a^-1 mod p (= a^(p-2) mod p)
void fpinv256(dig256 a);
void fpinv384(dig384 a);
void fpinv512(dig512 a);

// Subtraction a = modulus-a, or field negation, a = -a (mod p) if modulus=p
BOOL fpneg256(dig256 modulus, dig256 a);
BOOL fpneg384(dig384 modulus, dig384 a);
BOOL fpneg512(dig512 modulus, dig512 a);

// Evaluate if an element is in [0, modulus-1]
BOOL mod_eval256(dig256 a, dig256 modulus);
BOOL mod_eval384(dig384 a, dig384 modulus);
BOOL mod_eval512(dig512 a, dig512 modulus);

// Field addition, c = a+b mod p
void fpadd256(dig256 a, dig256 b, dig256 c);
void fpadd384(dig384 a, dig384 b, dig384 c);
void fpadd512(dig512 a, dig512 b, dig512 c);

// Field subtraction, c = a-b mod p
void fpsub256(dig256 a, dig256 b, dig256 c);
void fpsub384(dig384 a, dig384 b, dig384 c);
void fpsub512(dig512 a, dig512 b, dig512 c);

// Field division by 2, c = a/2 mod p
void fpdiv2_256(dig256 a, dig256 c);
void fpdiv2_384(dig384 a, dig384 c);
void fpdiv2_512(dig512 a, dig512 c);



/********** Main curve functions ***********/

// NOTE: (Most of) the following functions accept input/output points in standard affine representation (x,y).
//       These are the functions that are commonly required by most ECC protocols. 

// Weierstrass a=-3 curves ("Jac" curves)

// Set generator P = (x,y) on Weierstrass a=-3 curve
void eccset_Jac256(point_Jac256 P, PCurveStruct JacCurve);
void eccset_Jac384(point_Jac384 P, PCurveStruct JacCurve);
void eccset_Jac512(point_Jac512 P, PCurveStruct JacCurve);

// Check if point P = (x,y) on Weierstrass a=-3 curve is the point at infinity (0,0)
BOOL ecc_is_infinity_Jac256(point_Jac256 P, PCurveStruct JacCurve);
BOOL ecc_is_infinity_Jac384(point_Jac384 P, PCurveStruct JacCurve);
BOOL ecc_is_infinity_Jac512(point_Jac512 P, PCurveStruct JacCurve);

// Variable-base scalar multiplication Q = k.P using fixed-window method, Weierstrass a=-3 curve
// P is the input point and k is the scalar
BOOL ecc_scalar_mul_Jac256(point_Jac256 P, dig *k, point_Jac256 Q, PCurveStruct JacCurve);
BOOL ecc_scalar_mul_Jac384(point_Jac384 P, dig *k, point_Jac384 Q, PCurveStruct JacCurve);
BOOL ecc_scalar_mul_Jac512(point_Jac512 P, dig *k, point_Jac512 Q, PCurveStruct JacCurve);

// Fixed-base scalar multiplication Q = k.P, where P = P_table, using the Modified LSB-set method, Weierstrass a=-3 curve
// P_table is the input point table and k is the scalar. P_table is precalculated by calling ecc_precomp_fixed_Jacxxx
// Parameter "memory_use" (see possible values in the definition of "MemType") determines the size of P_table
BOOL ecc_scalar_mul_fixed_Jac256(point_Jac256 *P_table, dig *k, point_Jac256 Q, MemType memory_use, PCurveStruct JacCurve);
BOOL ecc_scalar_mul_fixed_Jac384(point_Jac384 *P_table, dig *k, point_Jac384 Q, MemType memory_use, PCurveStruct JacCurve);
BOOL ecc_scalar_mul_fixed_Jac512(point_Jac512 *P_table, dig *k, point_Jac512 Q, MemType memory_use, PCurveStruct JacCurve);

// Function that outputs the precomputed table "P_table" to be used by fixed-base scalar multiplications ecc_scalar_mul_fixed_Jacxxx
// P is the input point used to create the table
// Parameter "memory_use" should have the same value used by the corresponding invocation to ecc_scalar_mul_fixed_Jacxxx
point_Jac256* ecc_precomp_fixed_Jac256(point_Jac256 P, MemType memory_use, PCurveStruct JacCurve);
point_Jac384* ecc_precomp_fixed_Jac384(point_Jac384 P, MemType memory_use, PCurveStruct JacCurve);
point_Jac512* ecc_precomp_fixed_Jac512(point_Jac512 P, MemType memory_use, PCurveStruct JacCurve);

// Double-scalar multiplication R = k.P+l.Q, where P = P_table, using wNAF with Interleaving, Weierstrass a=-3 curve
// P_table is the fixed-base input point table, with corresponding scalar k, and Q is the variable-base input point, with corresponding scalar l
// P_table is precalculated by calling ecc_precomp_dblmul_Jacxxx
// Parameter "memory_use" (see possible values in the definition of "MemType") determines the size of P_table
BOOL ecc_double_scalar_mul_Jac256(point_Jac256 *P_table, dig *k, point_Jac256 Q, dig *l, point_Jac256 R, MemType memory_use, PCurveStruct JacCurve);
BOOL ecc_double_scalar_mul_Jac384(point_Jac384 *P_table, dig *k, point_Jac384 Q, dig *l, point_Jac384 R, MemType memory_use, PCurveStruct JacCurve);
BOOL ecc_double_scalar_mul_Jac512(point_Jac512 *P_table, dig *k, point_Jac512 Q, dig *l, point_Jac512 R, MemType memory_use, PCurveStruct JacCurve);

// Function that outputs the precomputed table "P_table" to be used by double-scalar multiplications ecc_double_scalar_mul_Jacxxx
// P is the input point used to create the table
// Parameter "memory_use" should have the same value used by the corresponding invocation to ecc_double_scalar_mul_Jacxxx
point_Jac256* ecc_precomp_dblmul_Jac256(point_Jac256 P, MemType memory_use, PCurveStruct JacCurve);
point_Jac384* ecc_precomp_dblmul_Jac384(point_Jac384 P, MemType memory_use, PCurveStruct JacCurve);
point_Jac512* ecc_precomp_dblmul_Jac512(point_Jac512 P, MemType memory_use, PCurveStruct JacCurve);

// Frees memory occupied by precomputation tables used during fixed-base and double-scalar multiplications.
// This function must be called once done using a table generated by ecc_precomp_fixed_Jacxxx or ecc_precomp_dblmul_Jacxxx. 
BOOL ecc_destroy_precomp_Jac256(point_Jac256* T_fixed);
BOOL ecc_destroy_precomp_Jac384(point_Jac384* T_fixed);
BOOL ecc_destroy_precomp_Jac512(point_Jac512* T_fixed);

// Main functions based on macros

// Copy point Q = (x,y) to P: P = Q
#define ecccopy_Jac256(Q, P); fpcopy256(Q->x, P->x);  \
                              fpcopy256(Q->y, P->y);  \

#define ecccopy_Jac384(Q, P); fpcopy384(Q->x, P->x);  \
                              fpcopy384(Q->y, P->y);  \

#define ecccopy_Jac512(Q, P); fpcopy512(Q->x, P->x);  \
                              fpcopy512(Q->y, P->y);  \

// Zero point (x,y): P = (0,0)
#define ecczero_Jac256(P);    fpzero256(P->x);  \
                              fpzero256(P->y);  \

#define ecczero_Jac384(P);    fpzero384(P->x);  \
                              fpzero384(P->y);  \

#define ecczero_Jac512(P);    fpzero512(P->x);  \
                              fpzero512(P->y);  \


// Twisted Edwards a=-1 curves ("Ted" Curves)

// Set generator P = (x,y) on twisted Edwards a=-1 curve
void eccset_Ted256(point_Ted256 P, PCurveStruct TedCurve);
void eccset_Ted384(point_Ted384 P, PCurveStruct TedCurve);
void eccset_Ted512(point_Ted512 P, PCurveStruct TedCurve);

// Check if point P = (x,y) on twisted Edwards a=-1 curve is the neutral point (0,1) 
BOOL ecc_is_neutral_Ted256(point_Ted256 P, PCurveStruct TedCurve);
BOOL ecc_is_neutral_Ted384(point_Ted384 P, PCurveStruct TedCurve);
BOOL ecc_is_neutral_Ted512(point_Ted512 P, PCurveStruct TedCurve);

// Variable-base scalar multiplication Q = k.P using fixed-window method, twisted Edwards a=-1 curve
// P is the input point and k is the scalar
BOOL ecc_scalar_mul_Ted256(point_Ted256 P, dig *k, point_Ted256 Q, PCurveStruct TedCurve);
BOOL ecc_scalar_mul_Ted384(point_Ted384 P, dig *k, point_Ted384 Q, PCurveStruct TedCurve);
BOOL ecc_scalar_mul_Ted512(point_Ted512 P, dig *k, point_Ted512 Q, PCurveStruct TedCurve);

// Fixed-base scalar multiplication Q = k.P, where P = P_table, using the Modified LSB-set method, twisted Edwards a=-1 curve
// P_table is the input point table and k is the scalar. P_table is precalculated by calling ecc_precomp_fixed_Tedxxx
// Parameter "memory_use" (see possible values in the definition of "MemType") determines the size of P_table
BOOL ecc_scalar_mul_fixed_Ted256(point_extaff_precomp_Ted256 *P_table, dig *k, point_Ted256 Q, MemType memory_use, PCurveStruct TedCurve);
BOOL ecc_scalar_mul_fixed_Ted384(point_extaff_precomp_Ted384 *P_table, dig *k, point_Ted384 Q, MemType memory_use, PCurveStruct TedCurve);
BOOL ecc_scalar_mul_fixed_Ted512(point_extaff_precomp_Ted512 *P_table, dig *k, point_Ted512 Q, MemType memory_use, PCurveStruct TedCurve);

// Function that outputs the precomputed table "P_table" to be used by fixed-base scalar multiplications ecc_scalar_mul_fixed_Tedxxx
// P is the input point used to create the table
// Parameter "memory_use" should have the same value used by the corresponding invocation to ecc_scalar_mul_fixed_Tedxxx
point_extaff_precomp_Ted256* ecc_precomp_fixed_Ted256(point_Ted256 P, MemType memory_use, PCurveStruct TedCurve);
point_extaff_precomp_Ted384* ecc_precomp_fixed_Ted384(point_Ted384 P, MemType memory_use, PCurveStruct TedCurve);
point_extaff_precomp_Ted512* ecc_precomp_fixed_Ted512(point_Ted512 P, MemType memory_use, PCurveStruct TedCurve);

// Double-scalar multiplication R = k.P+l.Q, where P = P_table, using wNAF with Interleaving, twisted Edwards a=-1 curve
// P_table is the fixed-base input point table, with corresponding scalar k, and Q is the variable-base input point, with corresponding scalar l
// P_table is precalculated by calling ecc_precomp_dblmul_Tedxxx
// Parameter "memory_use" (see possible values in the definition of "MemType") determines the size of P_table
BOOL ecc_double_scalar_mul_Ted256(point_extaff_precomp_Ted256 *P_table, dig *k, point_Ted256 Q, dig *l, point_Ted256 R, MemType memory_use, PCurveStruct TedCurve);
BOOL ecc_double_scalar_mul_Ted384(point_extaff_precomp_Ted384 *P_table, dig *k, point_Ted384 Q, dig *l, point_Ted384 R, MemType memory_use, PCurveStruct TedCurve);
BOOL ecc_double_scalar_mul_Ted512(point_extaff_precomp_Ted512 *P_table, dig *k, point_Ted512 Q, dig *l, point_Ted512 R, MemType memory_use, PCurveStruct TedCurve);

// Function that outputs the precomputed table "P_table" to be used by double-scalar multiplications ecc_double_scalar_mul_Tedxxx
// P is the input point used to create the table
// Parameter "memory_use" should have the same value used by the corresponding invocation to ecc_double_scalar_mul_Tedxxx
point_extaff_precomp_Ted256* ecc_precomp_dblmul_Ted256(point_Ted256 P, MemType memory_use, PCurveStruct TedCurve);
point_extaff_precomp_Ted384* ecc_precomp_dblmul_Ted384(point_Ted384 P, MemType memory_use, PCurveStruct TedCurve);
point_extaff_precomp_Ted512* ecc_precomp_dblmul_Ted512(point_Ted512 P, MemType memory_use, PCurveStruct TedCurve);

// Frees memory occupied by precomputation tables used during fixed-base and double-scalar multiplications.
// This function must be called once done using a table generated by ecc_precomp_fixed_Tedxxx or ecc_precomp_dblmul_Tedxxx. 
BOOL ecc_destroy_precomp_Ted256(point_extaff_precomp_Ted256* T_fixed);
BOOL ecc_destroy_precomp_Ted384(point_extaff_precomp_Ted384* T_fixed);
BOOL ecc_destroy_precomp_Ted512(point_extaff_precomp_Ted512* T_fixed);

// Main functions based on macros

// Copy point Q = (x,y) to P: P = Q
#define ecccopy_Ted256(Q, P)  ecccopy_Jac256(Q, P) 
#define ecccopy_Ted384(Q, P)  ecccopy_Jac384(Q, P) 
#define ecccopy_Ted512(Q, P)  ecccopy_Jac512(Q, P)

// Zero point (x,y): P = (0,0)
#define ecczero_Ted256(P)     ecczero_Jac256(P) 
#define ecczero_Ted384(P)     ecczero_Jac384(P) 
#define ecczero_Ted512(P)     ecczero_Jac512(P) 



/********** Additional curve functions ***********/

// NOTE: the following functions also accept input/output points in some projective coordinate system (X:Y:Z) or variants.
//       These functions can be needed in some special cases in which access to projective coordinates is required.   

// Weierstrass a=-3 curves ("Jac" curves) 

// Copy Jacobian point (X:Y:Z) on Weierstrass a=-3 curve, P = Q
void ecccopy_jac_Jac256(point_jac_Jac256 Q, point_jac_Jac256 P, PCurveStruct JacCurve);
void ecccopy_jac_Jac384(point_jac_Jac384 Q, point_jac_Jac384 P, PCurveStruct JacCurve);
void ecccopy_jac_Jac512(point_jac_Jac512 Q, point_jac_Jac512 P, PCurveStruct JacCurve);

// Check if Jacobian point P = (X:Y:Z) on Weierstrass a=-3 curve is the point at infinity (0:Y:0)
BOOL ecc_is_infinity_jac_Jac256(point_jac_Jac256 P, PCurveStruct JacCurve);
BOOL ecc_is_infinity_jac_Jac384(point_jac_Jac384 P, PCurveStruct JacCurve);
BOOL ecc_is_infinity_jac_Jac512(point_jac_Jac512 P, PCurveStruct JacCurve);

// Normalize a Jacobian point Q = (X:Y:Z) -> P = (x,y)
BOOL eccnorm_Jac256(point_jac_Jac256 Q, point_Jac256 P, PCurveStruct JacCurve);
BOOL eccnorm_Jac384(point_jac_Jac384 Q, point_Jac384 P, PCurveStruct JacCurve);
BOOL eccnorm_Jac512(point_jac_Jac512 Q, point_Jac512 P, PCurveStruct JacCurve);

// Point doubling P = 2P using Jacobian coordinates (X:Y:Z), Weierstrass a=-3 curve
BOOL eccdouble_jac_Jac256(point_jac_Jac256 P, PCurveStruct JacCurve);
BOOL eccdouble_jac_Jac384(point_jac_Jac384 P, PCurveStruct JacCurve);
BOOL eccdouble_jac_Jac512(point_jac_Jac512 P, PCurveStruct JacCurve);

// "Complete" addition P = P+Q using Jacobian coordinates (X:Y:Z), Weierstrass a=-3 curve
void eccadd_jac_Jac256(point_jac_Jac256 Q, point_jac_Jac256 P, PCurveStruct JacCurve); 
void eccadd_jac_Jac384(point_jac_Jac384 Q, point_jac_Jac384 P, PCurveStruct JacCurve); 
void eccadd_jac_Jac512(point_jac_Jac512 Q, point_jac_Jac512 P, PCurveStruct JacCurve);

// Additional functions based on macros

// Copy Jacobian point Q = (X:Y:Z) to P: P = Q
#define ecccopy_jac_Jac256(Q, P); fpcopy256(Q->X, P->X);  \
                                  fpcopy256(Q->Y, P->Y);  \
                                  fpcopy256(Q->Z, P->Z);  \

#define ecccopy_jac_Jac384(Q, P); fpcopy384(Q->X, P->X);  \
                                  fpcopy384(Q->Y, P->Y);  \
                                  fpcopy384(Q->Z, P->Z);  \

#define ecccopy_jac_Jac512(Q, P); fpcopy512(Q->X, P->X);  \
                                  fpcopy512(Q->Y, P->Y);  \
                                  fpcopy512(Q->Z, P->Z);  \

// Zero Jacobian point (X:Y:Z): P = (0:0:0)
#define ecczero_jac_Jac256(P);    fpzero256(P->X);  \
                                  fpzero256(P->Y);  \
                                  fpzero256(P->Z);  \

#define ecczero_jac_Jac384(P);    fpzero384(P->X);  \
                                  fpzero384(P->Y);  \
                                  fpzero384(P->Z);  \

#define ecczero_jac_Jac512(P);    fpzero512(P->X);  \
                                  fpzero512(P->Y);  \
                                  fpzero512(P->Z);  \

// Convert affine point Q = (x,y) to Jacobian P = (X:Y:1), where X=x, Y=y
#define eccconvert_aff_to_jac_Jac256(Q, P); fpcopy256(Q->x, P->X);        \
                                            fpcopy256(Q->y, P->Y);        \
                                            fpzero256(P->Z); P->Z[0] = 1; \
                              
#define eccconvert_aff_to_jac_Jac384(Q, P); fpcopy384(Q->x, P->X);        \
                                            fpcopy384(Q->y, P->Y);        \
                                            fpzero384(P->Z); P->Z[0] = 1; \

#define eccconvert_aff_to_jac_Jac512(Q, P); fpcopy512(Q->x, P->X);        \
                                            fpcopy512(Q->y, P->Y);        \
                                            fpzero512(P->Z); P->Z[0] = 1; \

// Twisted Edwards a=-1 curves (Ted Curves)

// Copy extended projective point on twisted Edwards a=-1 curve using projective coordinates (X:Y:Z:Ta:Tb), P = Q
void ecccopy_extproj_Ted256(point_extproj_Ted256 Q, point_extproj_Ted256 P, PCurveStruct TedCurve);
void ecccopy_extproj_Ted384(point_extproj_Ted384 Q, point_extproj_Ted384 P, PCurveStruct TedCurve);
void ecccopy_extproj_Ted512(point_extproj_Ted512 Q, point_extproj_Ted512 P, PCurveStruct TedCurve);

// Check if extended projective point P = (X:Y:Z:Ta:Tb) on twisted Edwards a=-1 curve is the neutral point (0:1:1) 
BOOL ecc_is_neutral_extproj_Ted256(point_extproj_Ted256 P, PCurveStruct TedCurve);
BOOL ecc_is_neutral_extproj_Ted384(point_extproj_Ted384 P, PCurveStruct TedCurve);
BOOL ecc_is_neutral_extproj_Ted512(point_extproj_Ted512 P, PCurveStruct TedCurve);

// Normalize a twisted Edwards point Q = (X:Y:Z) -> P = (x,y)
void eccnorm_Ted256(point_extproj_Ted256 Q, point_Ted256 P, PCurveStruct TedCurve);
void eccnorm_Ted384(point_extproj_Ted384 Q, point_Ted384 P, PCurveStruct TedCurve);
void eccnorm_Ted512(point_extproj_Ted512 Q, point_Ted512 P, PCurveStruct TedCurve);

// Point doubling 2P using extended projective coordinates (X:Y:Z:Ta:Tb), twisted Edwards a=-1 curve
void eccdouble_extproj_Ted256(point_extproj_Ted256 P, PCurveStruct TedCurve);
void eccdouble_extproj_Ted384(point_extproj_Ted384 P, PCurveStruct TedCurve);
void eccdouble_extproj_Ted512(point_extproj_Ted512 P, PCurveStruct TedCurve);

// Complete point addition P = P+Q or P = P+P using extended projective coordinates (X:Y:Z:Ta:Tb), twisted Edwards a=-1 curve
void eccadd_extproj_Ted256(point_extproj_Ted256 Q, point_extproj_Ted256 P, PCurveStruct TedCurve);
void eccadd_extproj_Ted384(point_extproj_Ted384 Q, point_extproj_Ted384 P, PCurveStruct TedCurve);
void eccadd_extproj_Ted512(point_extproj_Ted512 Q, point_extproj_Ted512 P, PCurveStruct TedCurve);

// Additional functions based on macros

// Copy extended projective point Q = (X:Y:Z:Ta:Tb) to P: P = Q
#define ecccopy_extproj_Ted256(Q, P); fpcopy256(Q->X, P->X);  \
                                      fpcopy256(Q->Y, P->Y);  \
                                      fpcopy256(Q->Z, P->Z);  \
                                      fpcopy256(Q->Ta, P->Ta);\
                                      fpcopy256(Q->Tb, P->Tb);\

#define ecccopy_extproj_Ted384(Q, P); fpcopy384(Q->X, P->X);  \
                                      fpcopy384(Q->Y, P->Y);  \
                                      fpcopy384(Q->Z, P->Z);  \
                                      fpcopy384(Q->Ta, P->Ta);\
                                      fpcopy384(Q->Tb, P->Tb);\

#define ecccopy_extproj_Ted512(Q, P); fpcopy512(Q->X, P->X);  \
                                      fpcopy512(Q->Y, P->Y);  \
                                      fpcopy512(Q->Z, P->Z);  \
                                      fpcopy512(Q->Ta, P->Ta);\
                                      fpcopy512(Q->Tb, P->Tb);\

// Zero extended projective point (X:Y:Z:Ta:Tb): P = (0:0:0:0:0)
#define ecczero_extproj_Ted256(P); fpzero256(P->X);  \
                                   fpzero256(P->Y);  \
                                   fpzero256(P->Z);  \
                                   fpzero256(P->Ta); \
                                   fpzero256(P->Tb); \

#define ecczero_extproj_Ted384(P); fpzero384(P->X);  \
                                   fpzero384(P->Y);  \
                                   fpzero384(P->Z);  \
                                   fpzero384(P->Ta); \
                                   fpzero384(P->Tb); \

#define ecczero_extproj_Ted512(P); fpzero512(P->X);  \
                                   fpzero512(P->Y);  \
                                   fpzero512(P->Z);  \
                                   fpzero512(P->Ta); \
                                   fpzero512(P->Tb); \

// Convert affine point Q = (x,y) to extended projective P = (X:Y:1:Ta:Tb), where X=x, Y=y, Ta=x, Ty=y
#define eccconvert_aff_to_extproj_Ted256(Q, P); fpcopy256(Q->x, P->X);        \
                                                fpcopy256(Q->y, P->Y);        \
                                                fpzero256(P->Z); P->Z[0] = 1; \
                                                fpcopy256(Q->x, P->Ta);       \
                                                fpcopy256(Q->y, P->Tb);       \
                              
#define eccconvert_aff_to_extproj_Ted384(Q, P); fpcopy384(Q->x, P->X);        \
                                                fpcopy384(Q->y, P->Y);        \
                                                fpzero384(P->Z); P->Z[0] = 1; \
                                                fpcopy384(Q->x, P->Ta);       \
                                                fpcopy384(Q->y, P->Tb);       \

#define eccconvert_aff_to_extproj_Ted512(Q, P); fpcopy512(Q->x, P->X);        \
                                                fpcopy512(Q->y, P->Y);        \
                                                fpzero512(P->Z); P->Z[0] = 1; \
                                                fpcopy512(Q->x, P->Ta);       \
                                                fpcopy512(Q->y, P->Tb);       \


#ifdef __cplusplus
}
#endif

#endif