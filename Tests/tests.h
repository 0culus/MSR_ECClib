/**************************************************************************
* Benchmarking/testing suite for MSR ECClib
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
* Abstract: header file for benchmarking/testing suite
*
* This software is based on the article by Joppe Bos, Craig Costello, 
* Patrick Longa and Michael Naehrig, "Selecting elliptic curves for
* cryptography: an efficiency and security analysis", preprint available
* at http://eprint.iacr.org/2014/130.
***************************************************************************/    

#include "msr_ecclib.h"


// Benchmark and test parameters  
#define ML_BENCH_LOOPS       100000     // Number of iterations per bench
#define ML_SHORT_BENCH_LOOPS 10000      // Number of iterations per bench (for expensive operations)
#define ML_TEST_LOOPS        100        // Number of iterations per test
#define ML_SHORT_TEST_LOOPS  1          // Number of iterations per test (for expensive operations)
#define MAIN_TESTS_ONLY                 // Select if tests need to be performed on "main functions" only

// This instruction queries the processor for information about the supported features and CPU type
#define cpuid(cpuinfo, infotype)    __cpuidex(cpuinfo, infotype, 0)


// Access timestamp counter for benchmarking
__int64 cpucycles(void);

// Detecting support for AVX instructions
BOOL is_avx_supported(void);


/********** Prototypes of utility functions ***********/

// Compare two field elements, a = b? (1) a>b, (0) a=b, (-1) a<b 
// SECURITY NOTE: this function does not have constant-time execution. TO BE USED FOR TESTING ONLY.
int fpcompare256(dig256 a, dig256 b);
int fpcompare384(dig384 a, dig384 b);
int fpcompare512(dig512 a, dig512 b);

// Get prime value p
void fp_prime256(dig256 p, PCurveStruct PCurve);
void fp_prime384(dig384 p, PCurveStruct PCurve);
void fp_prime512(dig512 p, PCurveStruct PCurve);

// Generate a pseudo-random element in [0, modulus-1]
// SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.
void random256(dig256 a, dig256 modulus);
void random384(dig384 a, dig384 modulus);
void random512(dig512 a, dig512 modulus);

// Subtraction without borrow, c = a-b where a>b
// SECURITY NOTE: this function does not have constant-time execution. TO BE USED FOR TESTING ONLY.
void sub256_test(dig256 a, dig256 b, dig256 c);
void sub384_test(dig384 a, dig384 b, dig384 c);
void sub512_test(dig512 a, dig512 b, dig512 c);

// Point doubling P = 2P using affine coordinates, Weierstrass curve, generic "a"
// SECURITY NOTE: this function does not have constant-time execution. TO BE USED FOR TESTING ONLY.
void eccdouble_waff_256(point_Jac256 P, PCurveStruct PCurve);
void eccdouble_waff_384(point_Jac384 P, PCurveStruct PCurve);
void eccdouble_waff_512(point_Jac512 P, PCurveStruct PCurve);

// Point addition P = P+Q using affine coordinates, Weierstrass a=-3 curve
// SECURITY NOTE: this function does not have constant-time execution. TO BE USED FOR TESTING ONLY.
void eccadd_waff_256(point_Jac256 Q, point_Jac256 P, PCurveStruct PCurve);
void eccadd_waff_384(point_Jac384 Q, point_Jac384 P, PCurveStruct PCurve);
void eccadd_waff_512(point_Jac512 Q, point_Jac512 P, PCurveStruct PCurve);

// Variable-base scalar multiplication Q = k.P using affine coordinates and the binary representation, Weierstrass a=-3 curve
// SECURITY NOTE: this function does not have constant-time execution. TO BE USED FOR TESTING ONLY.
void ecc_mul_waff_256(point_Jac256 P, dig *k, point_Jac256 Q, PCurveStruct PCurve);
void ecc_mul_waff_384(point_Jac384 P, dig *k, point_Jac384 Q, PCurveStruct PCurve);
void ecc_mul_waff_512(point_Jac512 P, dig *k, point_Jac512 Q, PCurveStruct PCurve);

// Convert point on one of the Ted curves to its corresponding isomorphic Weierstrass curve
void ecc_Ted256_to_weierstrass(point_Ted256 Q, point_Jac256 P, PCurveStruct TedCurve);
void ecc_Ted384_to_weierstrass(point_Ted384 Q, point_Jac384 P, PCurveStruct TedCurve);
void ecc_Ted512_to_weierstrass(point_Ted512 Q, point_Jac512 P, PCurveStruct TedCurve);

// Convert point on isomorphic Weierstrass curve to its corresponding Ted curve
void ecc_weierstrass_to_Ted256(point_Jac256 Q, point_Ted256 P, PCurveStruct TedCurve);
void ecc_weierstrass_to_Ted384(point_Jac384 Q, point_Ted384 P, PCurveStruct TedCurve);
void ecc_weierstrass_to_Ted512(point_Jac512 Q, point_Ted512 P, PCurveStruct TedCurve);


// Other testing functions used for verifying correctness of recoding algorithms

BOOL verify_mLSB_recoding(dig *scalar, int *digits, unsigned int nbits, unsigned int l, unsigned int d);
BOOL verify_recoding(dig *scalar, int *digits, unsigned int nbits);