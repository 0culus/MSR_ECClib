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
* Abstract: utility functions 
*
* This software is based on the article by Joppe Bos, Craig Costello, 
* Patrick Longa and Michael Naehrig, "Selecting elliptic curves for
* cryptography: an efficiency and security analysis", preprint available
* at http://eprint.iacr.org/2014/130.
***************************************************************************/  

#include <windows.h>
#include <intrin.h>
#include "tests.h"


__int64 cpucycles(void)
{ // Access timestamp counter for benchmarking
    return __rdtsc();
}


BOOL is_avx_supported()
{ // Determining if AVX instructions are supported
    int info[4];
    int nIds;
    BOOL AVX_SUPPORT = FALSE, AVX_available, XSAVE_XRSTORE_available;
    dig XCR_mask;
        
// Check if Visual Studio 2010 SP1 or later
#if (_MSC_FULL_VER >= 160040219)

    cpuid(info, 0);
    nIds = info[0];
    
    if (nIds >= 1){
        cpuid(info, 1);
        XSAVE_XRSTORE_available = ((info[2] & ((int)1 << 28)) != 0);  // XSAVE/XRSTORE used by OS?
        AVX_available = ((info[2] & ((int)1 << 28)) != 0);            // AVX supported by processor?
    } else {
        return FALSE;
    }

    if (AVX_available && XSAVE_XRSTORE_available) {
        XCR_mask = _xgetbv(_XCR_XFEATURE_ENABLED_MASK);                
        AVX_SUPPORT = ((XCR_mask & 0x06) != 0);                       // AVX registers restored at context switch?
    } else {
        return FALSE;
    }

#endif
        
    return AVX_SUPPORT;
}


#ifdef ECCURVES_256
//
// Utility functions for 256-bit prime fields   
//

int fpcompare256(dig256 a, dig256 b)
{ // Comparing two 256-bit field elements: (1) a>b, (0) a=b, (-1) a<b
  // SECURITY NOTE: this function does not have constant-time execution. It is for TESTING ONLY.
    int i;

    for (i = ML_WORDS256-1; i >= 0; i--)
    {
        if (a[i]>b[i]) return 1;
        else if (a[i]<b[i]) return -1;
    }
    return 0; 
}


void fp_prime256(dig256 p, PCurveStruct PCurve)
{ // Extracting 256-bit prime  
    unsigned int i;

    for (i = 0; i < ML_WORDS256; i++)
    {
        p[i] = PCurve->prime[i]; 
    }
    return;
}


void random256(dig256 a, dig256 modulus)
{ // Generating a pseudo-random 256-bit element in [0, modulus-1] 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.
    int i;
    dig *ptr0, *ptr1, *ptr2, *ptr3;
    unsigned char cadena[32];

    ptr0 = (dig*) &cadena[0];
    ptr1 = (dig*) &cadena[8];
    ptr2 = (dig*) &cadena[16];
    ptr3 = (dig*) &cadena[24];

    for (i=0; i<32; i++) cadena[i]=(unsigned char)rand();   // Obtain 256-bit number 
    a[0] = *ptr0; a[1] = *ptr1; a[2] = *ptr2; a[3] = *ptr3;
    while (fpcompare256(modulus, a)<1) {                    // Force it to [0, modulus-1]
        sub256_test(a, modulus, a);
    }
}


void sub256_test(dig256 a, dig256 b, dig256 c)
{ // 256-bit subtraction without borrow, c = a-b where a>b
  // SECURITY NOTE: this function does not have constant-time execution. It is for TESTING ONLY.
     
    unsigned int i;
    dig res, carry, borrow = 0;
  
    for (i = 0; i < ML_WORDS256; i++)
    {
        res = a[i] - b[i];
        carry = (a[i] < b[i]);
        c[i] = res - borrow;
        borrow = carry || (res < borrow);
    } 

    return;
}
#endif


#ifdef ECCURVES_384
//
// Utility functions for 384-bit prime fields
//

int fpcompare384(dig384 a, dig384 b)
{// Comparing two 384-bit field elements: (1) a>b, (0) a=b, (-1) a<b
  // SECURITY NOTE: this function does not have constant-time execution. It is for TESTING ONLY.
    int i;

    for (i = ML_WORDS384-1; i >= 0; i--)
    {
        if (a[i]>b[i]) return 1;
        else if (a[i]<b[i]) return -1;
    }
    return 0; 
}


void fp_prime384(dig384 p, PCurveStruct PCurve)
{// Extracting 384-bit prime  
    unsigned int i;

    for (i = 0; i < ML_WORDS384; i++)
    {
        p[i] = PCurve->prime[i]; 
    }
    return;
}


void random384(dig384 a, dig384 mopulus)
{ // Generating a random 384-bit element in [0, modulus-1]
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.
    int i;
    dig *ptr0, *ptr1, *ptr2, *ptr3, *ptr4, *ptr5;
    unsigned char cadena[48];

    ptr0 = (dig*) &cadena[0];
    ptr1 = (dig*) &cadena[8];
    ptr2 = (dig*) &cadena[16];
    ptr3 = (dig*) &cadena[24];
    ptr4 = (dig*) &cadena[32];
    ptr5 = (dig*) &cadena[40];

    for (i=0; i<48; i++) cadena[i]=(unsigned char)rand();   // Obtain 384-bit number
    a[0] = *ptr0; a[1] = *ptr1; a[2] = *ptr2; a[3] = *ptr3; a[4] = *ptr4; a[5] = *ptr5;
    while (fpcompare384(mopulus, a)<1) {                    // Force it to [0, modulus-1]
        sub384_test(a, mopulus, a);
    }
}


void sub384_test(dig384 a, dig384 b, dig384 c)
{ // 384-bit subtraction without borrow, c = a-b where a>b
  // SECURITY NOTE: this function does not have constant-time execution. It is for TESTING ONLY.
    unsigned int i;
    dig res, carry, borrow = 0;
  
    for (i = 0; i < ML_WORDS384; i++)
    {
        res = a[i] - b[i];
        carry = (a[i] < b[i]);
        c[i] = res - borrow;
        borrow = carry || (res < borrow);
    } 

    return;
}
#endif


#ifdef ECCURVES_512
//
// Utility functions for 512-bit prime fields
//

int fpcompare512(dig512 a, dig512 b)
{ // Comparing two 512-bit field elements: (1) a>b, (0) a=b, (-1) a<b
  // SECURITY NOTE: this function does not have constant-time execution. It is for TESTING ONLY.
    int i;

    for (i = ML_WORDS512-1; i >= 0; i--)
    {
        if (a[i]>b[i]) return 1;
        else if (a[i]<b[i]) return -1;
    }
    return 0; 
}


void fp_prime512(dig512 p, PCurveStruct PCurve)
{// Extracting 512-bit prime  
    unsigned int i;

    for (i = 0; i < ML_WORDS512; i++)
    {
        p[i] = PCurve->prime[i]; 
    }
    return;
}


void random512(dig512 a, dig512 mopulus)
{ // Generating a random 512-bit element in [0, modulus-1]
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.
    int i;
    dig *ptr0, *ptr1, *ptr2, *ptr3, *ptr4, *ptr5, *ptr6, *ptr7;
    unsigned char cadena[64];

    ptr0 = (dig*) &cadena[0];
    ptr1 = (dig*) &cadena[8];
    ptr2 = (dig*) &cadena[16];
    ptr3 = (dig*) &cadena[24];
    ptr4 = (dig*) &cadena[32];
    ptr5 = (dig*) &cadena[40];
    ptr6 = (dig*) &cadena[48];
    ptr7 = (dig*) &cadena[56];

    for (i=0; i<64; i++) cadena[i]=(unsigned char)rand();   // Obtain 512-bit number  
    a[0] = *ptr0; a[1] = *ptr1; a[2] = *ptr2; a[3] = *ptr3; a[4] = *ptr4; a[5] = *ptr5, a[6] = *ptr6; a[7] = *ptr7;
    while (fpcompare512(mopulus, a)<1) {                    // Force it to [0, modulus-1]
        sub512_test(a, mopulus, a);
    }
}


void sub512_test(dig512 a, dig512 b, dig512 c)
{ // 512-bit subtraction without borrow, c = a-b where a>b
  // SECURITY NOTE: this function does not have constant-time execution. It is for TESTING ONLY.
    unsigned int i;
    dig res, carry, borrow = 0;
  
    for (i = 0; i < ML_WORDS512; i++)
    {
        res = a[i] - b[i];
        carry = (a[i] < b[i]);
        c[i] = res - borrow;
        borrow = carry || (res < borrow);
    } 

    return;
}
#endif


BOOL verify_mLSB_recoding(dig *scalar, int *digits, unsigned int nbits, unsigned int l, unsigned int d)
{ // Verification of mLSB-set's recoding algorithm used in fixed-base scalar multiplication 
    unsigned int j, cwords = (l+ML_WORD-1)/ML_WORD;    // Number of computer words to represent scalar
    int i, digit;
    dig *generated_scalar, temp, temp2, carry, borrow;
    BOOL passed = TRUE;

    generated_scalar = (dig*)calloc(1, cwords*sizeof(*generated_scalar));
    if (generated_scalar == NULL)
        return FALSE;

    for (i = (l-1); i >= 0; i--)
    {
        // Shift generated scalar to the left by 1 (multiply by 2)
        temp = ((generated_scalar[0] >> (ML_WORD-1)) & 1) ;
        generated_scalar[0] = generated_scalar[0] << 1;

        for (j = 1; j < cwords; j++) {
            temp2 = ((generated_scalar[j] >> (ML_WORD-1)) & 1) ;
            generated_scalar[j] = (generated_scalar[j] << 1) | temp;
            temp = temp2;
        }
     
        // generated scalar + digit_i
        if (i < (int)d) {
            digit = digits[i] | 1;
            if (digit >= 0) {
                generated_scalar[0] = generated_scalar[0] + digit;
                carry = (generated_scalar[0] < digit);
                for (j = 1; j < cwords; j++)
                {
                    generated_scalar[j] = generated_scalar[j] + carry;    
                    carry = (generated_scalar[j] < carry);
                }
            } else {
                borrow = 0;
                temp = (dig)(-digit);
                for (j = 0; j < cwords; j++)
                {
                    temp2 = generated_scalar[j] - temp;
                    carry = (generated_scalar[j] < temp);
                    generated_scalar[j] = temp2 - borrow;
                    borrow = carry || (temp2 < borrow);
                    temp = 0;
                }
            } 
        } else {
            digit = digits[i]*(digits[i-(i/d)*d] | 1);
            if (digit >= 0) {
                generated_scalar[0] = generated_scalar[0] + digit;
                carry = (generated_scalar[0] < digit);
                for (j = 1; j < cwords; j++)
                {
                    generated_scalar[j] = generated_scalar[j] + carry;    
                    carry = (generated_scalar[j] < carry);
                }
            } else {
                borrow = 0;
                temp = (dig)(-digit);
                for (j = 0; j < cwords; j++)
                {
                    temp2 = generated_scalar[j] - temp;
                    carry = (generated_scalar[j] < temp);
                    generated_scalar[j] = temp2 - borrow;
                    borrow = carry || (temp2 < borrow);
                    temp = 0;
                }
            } 
        }
    }

    for (j = 0; j < (nbits+ML_WORD-1)/(ML_WORD); j++)
    {
        if (scalar[j] != generated_scalar[j]) 
            passed = FALSE;
    }

// cleanup
    if (generated_scalar != NULL)
        free(generated_scalar);

    return passed;
}


BOOL verify_recoding(dig *scalar, int *digits, unsigned int nbits)
{ // Verification of wNAF's recoding algorithm used in double-scalar multiplication 
    unsigned int j, cwords = (nbits+ML_WORD-1)/ML_WORD;    // Number of computer words to represent scalar
    int i, digit;
    dig *generated_scalar, temp, temp2, carry, borrow;
    BOOL passed = TRUE;

    generated_scalar = (dig *)calloc(1, cwords*sizeof(*generated_scalar));
    if (generated_scalar == NULL)
        return FALSE;

    for (i = (nbits-1); i >= 0; i--)
    {
        // Shift generated scalar to the left by 1 (multiply by 2)
        temp = ((generated_scalar[0] >> (ML_WORD-1)) & 1) ;
        generated_scalar[0] = generated_scalar[0] << 1;

        for (j = 1; j < cwords; j++) {
            temp2 = ((generated_scalar[j] >> (ML_WORD-1)) & 1) ;
            generated_scalar[j] = (generated_scalar[j] << 1) | temp;
            temp = temp2;
        }
     
        // generated scalar + digit_i
        digit = digits[i];
        if (digit >= 0) {
            generated_scalar[0] = generated_scalar[0] + digit;
            carry = (generated_scalar[0] < digit);
            for (j = 1; j < cwords; j++)
            {
                generated_scalar[j] = generated_scalar[j] + carry;    
                carry = (generated_scalar[j] < carry);
            }
        } else {
            borrow = 0;
            temp = (dig)(-digit);
            for (j = 0; j < cwords; j++)
            {
                temp2 = generated_scalar[j] - temp;
                carry = (generated_scalar[j] < temp);
                generated_scalar[j] = temp2 - borrow;
                borrow = carry || (temp2 < borrow);
                temp = 0;
            }
        } 
    }

    for (j = 0; j < (nbits+ML_WORD-1)/(ML_WORD); j++)
    {
        if (scalar[j]!=generated_scalar[j]) 
            passed = FALSE;
    }

// cleanup
    if (generated_scalar != NULL)
        free(generated_scalar);

    return passed;
}


#ifdef ECCURVES_256
//
// Utility functions for elliptic curves   
//

void eccdouble_waff_256(point_Jac256 P, PCurveStruct PCurve)      
{ // Point doubling in affine coordinates, P = 2P
  // Weierstrass curve, generic "a"
  // Input:  P = (x,y) in affine coordinates
  // Output: 2P = (x,y) in affine coordinates
  // SECURITY NOTE: this function does not have constant-time execution. It is for TESTING ONLY.
    dig i, num_words = (PCurve->pbits+ML_WORD-1)/ML_WORD;
    dig256 t1, t2, t3, a;

    for (i = 0; i < num_words; i++) 
    {
         a[i] = PCurve->parameter1[i]; 
    }

    // Check if P is the point at infinity (0,0)
    if (ecc_is_infinity_Jac256(P, PCurve) == TRUE) return;
    
    fpsqr256(P->x, t1);               // t1 = x^2
    fpadd256(t1, t1, t2);             // t2 = 2x^2
    fpadd256(t1, t2, t1);             // t1 = 3x^2
    fpadd256(t1, a, t1);              // t1 = 3x^2+a
    fpadd256(P->y, P->y, t3);         // t3 = 2y
    fpinv256(t3);                     // t3 = 1/2y
    fpmul256(t1, t3, t2);             // t2 = alpha = (3x^2+a)/2y
    fpsqr256(t2, t3);                 // t3 = alpha^2
    fpsub256(t3, P->x, t3);           // t3 = alpha^2 - x
    fpsub256(t3, P->x, t3);           // t3 = alpha^2 - 2x
    fpsub256(P->x, t3, t1);           // t1 = x-xfinal
    fpcopy256(t3, P->x);               // xfinal = alpha^2 - 2x
    fpmul256(t1, t2, t3);             // t3 = alpha.(x-xfinal)
    fpsub256(t3, P->y, P->y);         // yfinal = alpha.(x-xfinal) - y

    return;
}


void eccadd_waff_256(point_Jac256 Q, point_Jac256 P, PCurveStruct PCurve)      
{ // Point addition in affine coordinates P = P+Q
  // Weierstrass a=-3 curve
  // Input:  P = (x1,y1) and Q = (x2,y2) in affine coordinates
  // Output: P+Q = (x1,y1) in affine coordinates
  // SECURITY NOTE: this function does not have constant-time execution. It is for TESTING ONLY.
    dig256 t1, t2, t3, t4; 
    
    // Check if P is the point at infinity (0,0)
    if (ecc_is_infinity_Jac256(P, PCurve) == TRUE) {     
        fpcopy256(Q->x, P->x); fpcopy256(Q->y, P->y);
        return;   
    }    
    // Check if Q is the point at infinity (0,0)
    if (ecc_is_infinity_Jac256(Q, PCurve) == TRUE) return;    

    fpsub256(Q->x, P->x, t1);           // t1 = x2-x1
    fpsub256(Q->y, P->y, t2);           // t2 = y2-y1
    
    if (fp_iszero256(t1) == TRUE) {
        if (fp_iszero256(t2) == TRUE) {
            eccdouble_waff_256(P, PCurve);
            return;
        } else {
            fpzero256(P->x); fpzero256(P->y);                                 
            return;
        }
    }
    
    fpinv256(t1);                     // t2 = 1/(x2-x1)
    fpmul256(t1, t2, t3);             // t3 = alpha = (y2-y1)/(x2-x1)
    fpsqr256(t3, t2);                 // t2 = alpha^2
    fpsub256(t2, P->x, t2);           // t2 = alpha^2 - x1
    fpsub256(t2, Q->x, t2);           // t2 = alpha^2 - x1 - x2
    fpsub256(P->x, t2, t4);           // t4 = x1-xfinal
    fpcopy256(t2, P->x);              // xfinal = alpha^2 - x1 - x2
    fpmul256(t3, t4, t1);             // t1 = alpha.(x1-xfinal)
    fpsub256(t1, P->y, P->y);         // yfinal = alpha.(x1-xfinal) - y1

    return;
}


void ecc_mul_waff_256(point_Jac256 P, dig *k, point_Jac256 Q, PCurveStruct PCurve)                                
{ // Variable-base scalar multiplication Q = k.P using affine coordinates and the binary method 
  // Weierstrass a=-3 curve
  // SECURITY NOTE: this function does not have regular execution. It should be used for TESTING ONLY.
    unsigned int i, j=0, scalar_digit, scalar_bit=0, nbit=PCurve->nbits;
    unsigned int* pscalar_digit, nwords;

    nwords = (nbit+8*sizeof(dig)-1)/(8*sizeof(dig));
    nbit = 8*sizeof(dig)*nwords;
    pscalar_digit = (unsigned int*)k + (nbit+8*sizeof(unsigned int)-1)/(8*sizeof(unsigned int)) - 1;
    scalar_digit = *pscalar_digit; 

    while (scalar_bit==0 && nbit>0) {
        if (j==sizeof(unsigned int)*8) {
            j = 0;
            pscalar_digit--;
            scalar_digit = *pscalar_digit; 
        }
        scalar_bit = (scalar_digit & 0xF0000000) >> 31;
        scalar_digit = scalar_digit << 1;
        j++;
        nbit--;
    }
    
    fpcopy256(P->x, Q->x);
    fpcopy256(P->y, Q->y);
    for (i = 0; i < nbit; i++)
    {
        if (j==sizeof(unsigned int)*8) {
            j = 0;
            pscalar_digit--;
            scalar_digit = *pscalar_digit; 
        }
        scalar_bit = (scalar_digit & 0xF0000000) >> 31;
        scalar_digit = scalar_digit << 1;
        
        eccdouble_waff_256(Q, PCurve); 
        if (scalar_bit==1) {
            eccadd_waff_256(P, Q, PCurve); 
        }
        j++;
    }
        
    return;
}

/******** Point conversion functions for "Ted256" Curves ********/

void ecc_Ted256_to_weierstrass(point_Ted256 Q, point_Jac256 P, PCurveStruct TedCurve)
{ // Convert point on one of the Ted curves to its corresponding isomorphic Weierstrass curve
	// Input:  Q = (xTE,yTE) on twisted Edwards curve
	//         Twisted Edwards curve struct TedCurve    
	// Output: P = (xW,yW) on Weierstrass curve
	dig i, num_words = (TedCurve->pbits+ML_WORD-1)/ML_WORD;
	dig256 t1, t2, t3, t4, a, d;

	for (i = 0; i < num_words; i++) 
	{
		a[i] = TedCurve->parameter1[i]; 
		d[i] = TedCurve->parameter2[i]; 
	}

	fpzero256(t1); t1[0] = 5;             // t1 = 5
	fpmul256(a, t1, t2);                  // t2 = 5a
	fpsub256(t2, d, t2);                  // t2 = 5a-d
	fpmul256(d, t1, t3);                  // t3 = 5d
	fpsub256(a, t3, t1);                  // t1 = a-5d
	fpmul256(Q->y, t1, t3);               // t3 = yTE*(a-5d)
	fpadd256(t3, t2, t2);                 // t2 = (5a-d) + yTE*(a-5d)
	fpzero256(t1); t1[0] = 1;             // t1 = 1
	fpsub256(t1, Q->y, t3);               // t3 = 1-yTE
	fpzero256(t1); t1[0] = 12;            // t1 = 12
	fpmul256(t1, t3, t4);                 // t4 = 12(1-yTE)
	fpinv256(t4);                         // t4 = 1/12(1-yTE)
	fpmul256(Q->x, t3, t1);               // t1 = xTE*(1-yTE)
	fpadd256(t1, t1, t3);                 // t3 = 2xTE*(1-yTE)
	fpadd256(t3, t3, t3);                 // t3 = 4xTE*(1-yTE)
	fpinv256(t3);                         // t3 = 1/4xTE*(1-yTE)
	fpmul256(t4, t2, P->x);               // Xfinal = ((5a-d) + yTE*(a-5d))/12(1-yTE)
	fpzero256(t1); t1[0] = 1;             // t1 = 1
	fpadd256(Q->y, t1, t1);               // t1 = yTE+1
	fpsub256(a, d, t2);                   // t2 = a-d
	fpmul256(t1, t2, t4);                 // t4 = (a-d)*(yTE+1)
	fpmul256(t4, t3, P->y);               // Yfinal = ((a-d)*(yTE+1))/4xTE*(1-yTE)

	// cleanup
	fpzero256(t1);
	fpzero256(t2);
	fpzero256(t3);
	fpzero256(t4);

	return;
}


void ecc_weierstrass_to_Ted256(point_Jac256 Q, point_Ted256 P, PCurveStruct TedCurve)
{ // Convert point on isomorphic Weierstrass curve to its corresponding Ted curve
	// Input:  Q = (xW,yW) on Weierstrass curve
	//         Twisted Edwards curve struct TedCurve    
	// Output: P = (xTE,yTE) on twisted Edwards curve
	dig i, num_words = (TedCurve->pbits+ML_WORD-1)/ML_WORD;
	dig256 t1, t2, t3, a, d;

	for (i = 0; i < num_words; i++) 
	{
		a[i] = TedCurve->parameter1[i]; 
		d[i] = TedCurve->parameter2[i]; 
	}

	fpadd256(Q->x, Q->x, t1);
	fpadd256(Q->x, t1, t1);
	fpadd256(t1, t1, t1);             // t1 = 6xW
	fpsub256(t1, a, t2);              // t2 = 6xW - a
	fpsub256(t2, d, t2);              // t2 = 6xW - a - d
	fpadd256(Q->y, Q->y, t3);
	fpadd256(Q->y, t3, t3);
	fpadd256(t3, t3, t3);             // t3 = 6yW
	fpinv256(t3);                     // t3 = 1/6yW
	fpmul256(t2, t3, P->x);           // Xfinal = (6xW - a - d)/6yW
	fpadd256(t1, t1, t1);             // t1 = 12xW
	fpadd256(t1, d, t2);              // t2 = 12xW + d
	fpadd256(t1, a, t1);              // t1 = 12xW + a
	fpadd256(a, a, t3);  
	fpsub256(t2, t3, t2);             // t2 = 12xW + d - 2a   
	fpsub256(t2, t3, t2);             // t2 = 12xW + d - 4a 
	fpsub256(t2, a, t2);              // t2 = 12xW + d - 5a  
	fpadd256(d, d, t3);  
	fpsub256(t1, t3, t1);             // t1 = 12xW + a - 2d   
	fpsub256(t1, t3, t1);             // t1 = 12xW + a - 4d 
	fpsub256(t1, d, t1);              // t1 = 12xW + a - 5d          
	fpinv256(t1);                     // t1 = 1/(12xW + a - 5d)
	fpmul256(t1, t2, P->y);           // Yfinal = (12xW + d - 5a)/(12xW + a - 5d)

	// cleanup
	fpzero256(t1);
	fpzero256(t2);
	fpzero256(t3);

	return;
}
#endif


#ifdef ECCURVES_384

void eccdouble_waff_384(point_Jac384 P, PCurveStruct PCurve)      
{ // Point doubling in affine coordinates, P = 2P
  // Weierstrass curve, generic "a"
  // Input:  P = (x,y) in affine coordinates
  // Output: 2P = (x,y) in affine coordinates
  // SECURITY NOTE: this function does not have constant-time execution. It is for TESTING ONLY.
    dig i, num_words = (PCurve->pbits+ML_WORD-1)/ML_WORD;
    dig384 t1, t2, t3, a;

    for (i = 0; i < num_words; i++) 
    {
         a[i] = PCurve->parameter1[i]; 
    }

    // Check if P is the point at infinity (0,0)
    if (ecc_is_infinity_Jac384(P, PCurve) == TRUE) return;
    
    fpsqr384(P->x, t1);               // t1 = x^2
    fpadd384(t1, t1, t2);             // t2 = 2x^2
    fpadd384(t1, t2, t1);             // t1 = 3x^2
    fpadd384(t1, a, t1);              // t1 = 3x^2+a
    fpadd384(P->y, P->y, t3);         // t3 = 2y
    fpinv384(t3);                     // t3 = 1/2y
    fpmul384(t1, t3, t2);             // t2 = alpha = (3x^2+a)/2y
    fpsqr384(t2, t3);                 // t3 = alpha^2
    fpsub384(t3, P->x, t3);           // t3 = alpha^2 - x
    fpsub384(t3, P->x, t3);           // t3 = alpha^2 - 2x
    fpsub384(P->x, t3, t1);           // t1 = x-xfinal
    fpcopy384(t3, P->x);               // xfinal = alpha^2 - 2x
    fpmul384(t1, t2, t3);             // t3 = alpha.(x-xfinal)
    fpsub384(t3, P->y, P->y);         // yfinal = alpha.(x-xfinal) - y

    return;
}


void eccadd_waff_384(point_Jac384 Q, point_Jac384 P, PCurveStruct PCurve)      
{ // Point addition in affine coordinates P = P+Q
  // Weierstrass a=-3 curve
  // Input:  P = (x1,y1) and Q = (x2,y2) in affine coordinates
  // Output: P+Q = (x1,y1) in affine coordinates
  // SECURITY NOTE: this function does not have constant-time execution. It is for TESTING ONLY.
    dig384 t1, t2, t3, t4; 
    
    // Check if P is the point at infinity (0,0)
    if (ecc_is_infinity_Jac384(P, PCurve) == TRUE) {     
        fpcopy384(Q->x, P->x); fpcopy384(Q->y, P->y);
        return;   
    }    
    // Check if Q is the point at infinity (0,0)
    if (ecc_is_infinity_Jac384(Q, PCurve) == TRUE) return;    

    fpsub384(Q->x, P->x, t1);           // t1 = x2-x1
    fpsub384(Q->y, P->y, t2);           // t2 = y2-y1
    
    if (fp_iszero384(t1) == TRUE) {
        if (fp_iszero384(t2) == TRUE) {
            eccdouble_waff_384(P, PCurve);
            return;
        } else {
            fpzero384(P->x); fpzero384(P->y);                                 
            return;
        }
    }
    
    fpinv384(t1);                     // t2 = 1/(x2-x1)
    fpmul384(t1, t2, t3);             // t3 = alpha = (y2-y1)/(x2-x1)
    fpsqr384(t3, t2);                 // t2 = alpha^2
    fpsub384(t2, P->x, t2);           // t2 = alpha^2 - x1
    fpsub384(t2, Q->x, t2);           // t2 = alpha^2 - x1 - x2
    fpsub384(P->x, t2, t4);           // t4 = x1-xfinal
    fpcopy384(t2, P->x);              // xfinal = alpha^2 - x1 - x2
    fpmul384(t3, t4, t1);             // t1 = alpha.(x1-xfinal)
    fpsub384(t1, P->y, P->y);         // yfinal = alpha.(x1-xfinal) - y1

    return;
}


void ecc_mul_waff_384(point_Jac384 P, dig *k, point_Jac384 Q, PCurveStruct PCurve)                                
{ // Variable-base scalar multiplication Q = k.P using affine coordinates and the binary method 
  // Weierstrass a=-3 curve
  // SECURITY NOTE: this function does not have regular execution. It should be used for TESTING ONLY.
    unsigned int i, j=0, scalar_digit, scalar_bit=0, nbit=PCurve->nbits;
    unsigned int* pscalar_digit, nwords;

    nwords = (nbit+8*sizeof(dig)-1)/(8*sizeof(dig));
    nbit = 8*sizeof(dig)*nwords;
    pscalar_digit = (unsigned int*)k + (nbit+8*sizeof(unsigned int)-1)/(8*sizeof(unsigned int)) - 1;
    scalar_digit = *pscalar_digit; 

    while (scalar_bit==0 && nbit>0) {
        if (j==sizeof(unsigned int)*8) {
            j = 0;
            pscalar_digit--;
            scalar_digit = *pscalar_digit; 
        }
        scalar_bit = (scalar_digit & 0xF0000000) >> 31;
        scalar_digit = scalar_digit << 1;
        j++;
        nbit--;
    }
    
    fpcopy384(P->x, Q->x);
    fpcopy384(P->y, Q->y);
    for (i = 0; i < nbit; i++)
    {
        if (j==sizeof(unsigned int)*8) {
            j = 0;
            pscalar_digit--;
            scalar_digit = *pscalar_digit; 
        }
        scalar_bit = (scalar_digit & 0xF0000000) >> 31;
        scalar_digit = scalar_digit << 1;
        
        eccdouble_waff_384(Q, PCurve); 
        if (scalar_bit==1) {
            eccadd_waff_384(P, Q, PCurve); 
        }
        j++;
    }
        
    return;
}

/******** Point conversion functions for "Ted384" Curves ********/

void ecc_Ted384_to_weierstrass(point_Ted384 Q, point_Jac384 P, PCurveStruct TedCurve)
{ // Convert point on one of the Ted curves to its corresponding isomorphic Weierstrass curve
	// Input:  Q = (xTE,yTE) on twisted Edwards curve
	//         Twisted Edwards curve struct TedCurve    
	// Output: P = (xW,yW) on Weierstrass curve
	dig i, num_words = (TedCurve->pbits+ML_WORD-1)/ML_WORD;
	dig384 t1, t2, t3, t4, a, d;

	for (i = 0; i < num_words; i++) 
	{
		a[i] = TedCurve->parameter1[i]; 
		d[i] = TedCurve->parameter2[i]; 
	}

	fpzero384(t1); t1[0] = 5;             // t1 = 5
	fpmul384(a, t1, t2);                  // t2 = 5a
	fpsub384(t2, d, t2);                  // t2 = 5a-d
	fpmul384(d, t1, t3);                  // t3 = 5d
	fpsub384(a, t3, t1);                  // t1 = a-5d
	fpmul384(Q->y, t1, t3);               // t3 = yTE*(a-5d)
	fpadd384(t3, t2, t2);                 // t2 = (5a-d) + yTE*(a-5d)
	fpzero384(t1); t1[0] = 1;             // t1 = 1
	fpsub384(t1, Q->y, t3);               // t3 = 1-yTE
	fpzero384(t1); t1[0] = 12;            // t1 = 12
	fpmul384(t1, t3, t4);                 // t4 = 12(1-yTE)
	fpinv384(t4);                         // t4 = 1/12(1-yTE)
	fpmul384(Q->x, t3, t1);               // t1 = xTE*(1-yTE)
	fpadd384(t1, t1, t3);                 // t3 = 2xTE*(1-yTE)
	fpadd384(t3, t3, t3);                 // t3 = 4xTE*(1-yTE)
	fpinv384(t3);                         // t3 = 1/4xTE*(1-yTE)
	fpmul384(t4, t2, P->x);               // Xfinal = ((5a-d) + yTE*(a-5d))/12(1-yTE)
	fpzero384(t1); t1[0] = 1;             // t1 = 1
	fpadd384(Q->y, t1, t1);               // t1 = yTE+1
	fpsub384(a, d, t2);                   // t2 = a-d
	fpmul384(t1, t2, t4);                 // t4 = (a-d)*(yTE+1)
	fpmul384(t4, t3, P->y);               // Yfinal = ((a-d)*(yTE+1))/4xTE*(1-yTE)

	// cleanup
	fpzero384(t1);
	fpzero384(t2);
	fpzero384(t3);
	fpzero384(t4);

	return;
}


void ecc_weierstrass_to_Ted384(point_Jac384 Q, point_Ted384 P, PCurveStruct TedCurve)
{ // Convert point on isomorphic Weierstrass curve to its corresponding Ted curve
	// Input:  Q = (xW,yW) on Weierstrass curve
	//         Twisted Edwards curve struct TedCurve    
	// Output: P = (xTE,yTE) on twisted Edwards curve
	dig i, num_words = (TedCurve->pbits+ML_WORD-1)/ML_WORD;
	dig384 t1, t2, t3, a, d;

	for (i = 0; i < num_words; i++) 
	{
		a[i] = TedCurve->parameter1[i]; 
		d[i] = TedCurve->parameter2[i]; 
	}

	fpadd384(Q->x, Q->x, t1);
	fpadd384(Q->x, t1, t1);
	fpadd384(t1, t1, t1);             // t1 = 6xW
	fpsub384(t1, a, t2);              // t2 = 6xW - a
	fpsub384(t2, d, t2);              // t2 = 6xW - a - d
	fpadd384(Q->y, Q->y, t3);
	fpadd384(Q->y, t3, t3);
	fpadd384(t3, t3, t3);             // t3 = 6yW
	fpinv384(t3);                     // t3 = 1/6yW
	fpmul384(t2, t3, P->x);           // Xfinal = (6xW - a - d)/6yW
	fpadd384(t1, t1, t1);             // t1 = 12xW
	fpadd384(t1, d, t2);              // t2 = 12xW + d
	fpadd384(t1, a, t1);              // t1 = 12xW + a
	fpadd384(a, a, t3);  
	fpsub384(t2, t3, t2);             // t2 = 12xW + d - 2a   
	fpsub384(t2, t3, t2);             // t2 = 12xW + d - 4a 
	fpsub384(t2, a, t2);              // t2 = 12xW + d - 5a  
	fpadd384(d, d, t3);  
	fpsub384(t1, t3, t1);             // t1 = 12xW + a - 2d   
	fpsub384(t1, t3, t1);             // t1 = 12xW + a - 4d 
	fpsub384(t1, d, t1);              // t1 = 12xW + a - 5d          
	fpinv384(t1);                     // t1 = 1/(12xW + a - 5d)
	fpmul384(t1, t2, P->y);           // Yfinal = (12xW + d - 5a)/(12xW + a - 5d)

	// cleanup
	fpzero384(t1);
	fpzero384(t2);
	fpzero384(t3);

	return;
}

#endif


#ifdef ECCURVES_512

void eccdouble_waff_512(point_Jac512 P, PCurveStruct PCurve)      
{ // Point doubling in affine coordinates, P = 2P
  // Weierstrass curve, generic "a"
  // Input:  P = (x,y) in affine coordinates
  // Output: 2P = (x,y) in affine coordinates
  // SECURITY NOTE: this function does not have constant-time execution. It is for TESTING ONLY.
    dig i, num_words = (PCurve->pbits+ML_WORD-1)/ML_WORD;
    dig512 t1, t2, t3, a;

    for (i = 0; i < num_words; i++) 
    {
         a[i] = PCurve->parameter1[i]; 
    }

    // Check if P is the point at infinity (0,0)
    if (ecc_is_infinity_Jac512(P, PCurve) == TRUE) return;
    
    fpsqr512(P->x, t1);               // t1 = x^2
    fpadd512(t1, t1, t2);             // t2 = 2x^2
    fpadd512(t1, t2, t1);             // t1 = 3x^2
    fpadd512(t1, a, t1);              // t1 = 3x^2+a
    fpadd512(P->y, P->y, t3);         // t3 = 2y
    fpinv512(t3);                     // t3 = 1/2y
    fpmul512(t1, t3, t2);             // t2 = alpha = (3x^2+a)/2y
    fpsqr512(t2, t3);                 // t3 = alpha^2
    fpsub512(t3, P->x, t3);           // t3 = alpha^2 - x
    fpsub512(t3, P->x, t3);           // t3 = alpha^2 - 2x
    fpsub512(P->x, t3, t1);           // t1 = x-xfinal
    fpcopy512(t3, P->x);               // xfinal = alpha^2 - 2x
    fpmul512(t1, t2, t3);             // t3 = alpha.(x-xfinal)
    fpsub512(t3, P->y, P->y);         // yfinal = alpha.(x-xfinal) - y

    return;
}


void eccadd_waff_512(point_Jac512 Q, point_Jac512 P, PCurveStruct PCurve)      
{ // Point addition in affine coordinates P = P+Q
  // Weierstrass a=-3 curve
  // Input:  P = (x1,y1) and Q = (x2,y2) in affine coordinates
  // Output: P+Q = (x1,y1) in affine coordinates
  // SECURITY NOTE: this function does not have constant-time execution. It is for TESTING ONLY.
    dig512 t1, t2, t3, t4; 
    
    // Check if P is the point at infinity (0,0)
    if (ecc_is_infinity_Jac512(P, PCurve) == TRUE) {     
        fpcopy512(Q->x, P->x); fpcopy512(Q->y, P->y);
        return;   
    }    
    // Check if Q is the point at infinity (0,0)
    if (ecc_is_infinity_Jac512(Q, PCurve) == TRUE) return;    

    fpsub512(Q->x, P->x, t1);           // t1 = x2-x1
    fpsub512(Q->y, P->y, t2);           // t2 = y2-y1
    
    if (fp_iszero512(t1) == TRUE) {
        if (fp_iszero512(t2) == TRUE) {
            eccdouble_waff_512(P, PCurve);
            return;
        } else {
            fpzero512(P->x); fpzero512(P->y);                                 
            return;
        }
    }
    
    fpinv512(t1);                     // t2 = 1/(x2-x1)
    fpmul512(t1, t2, t3);             // t3 = alpha = (y2-y1)/(x2-x1)
    fpsqr512(t3, t2);                 // t2 = alpha^2
    fpsub512(t2, P->x, t2);           // t2 = alpha^2 - x1
    fpsub512(t2, Q->x, t2);           // t2 = alpha^2 - x1 - x2
    fpsub512(P->x, t2, t4);           // t4 = x1-xfinal
    fpcopy512(t2, P->x);              // xfinal = alpha^2 - x1 - x2
    fpmul512(t3, t4, t1);             // t1 = alpha.(x1-xfinal)
    fpsub512(t1, P->y, P->y);         // yfinal = alpha.(x1-xfinal) - y1

    return;
}


void ecc_mul_waff_512(point_Jac512 P, dig *k, point_Jac512 Q, PCurveStruct PCurve)                                
{ // Variable-base scalar multiplication Q = k.P using affine coordinates and the binary method 
  // Weierstrass a=-3 curve
  // SECURITY NOTE: this function does not have regular execution. It should be used for TESTING ONLY.
    unsigned int i, j=0, scalar_digit, scalar_bit=0, nbit=PCurve->nbits;
    unsigned int* pscalar_digit, nwords;

    nwords = (nbit+8*sizeof(dig)-1)/(8*sizeof(dig));
    nbit = 8*sizeof(dig)*nwords;
    pscalar_digit = (unsigned int*)k + (nbit+8*sizeof(unsigned int)-1)/(8*sizeof(unsigned int)) - 1;
    scalar_digit = *pscalar_digit; 

    while (scalar_bit==0 && nbit>0) {
        if (j==sizeof(unsigned int)*8) {
            j = 0;
            pscalar_digit--;
            scalar_digit = *pscalar_digit; 
        }
        scalar_bit = (scalar_digit & 0xF0000000) >> 31;
        scalar_digit = scalar_digit << 1;
        j++;
        nbit--;
    }
    
    fpcopy512(P->x, Q->x);
    fpcopy512(P->y, Q->y);
    for (i = 0; i < nbit; i++)
    {
        if (j==sizeof(unsigned int)*8) {
            j = 0;
            pscalar_digit--;
            scalar_digit = *pscalar_digit; 
        }
        scalar_bit = (scalar_digit & 0xF0000000) >> 31;
        scalar_digit = scalar_digit << 1;
        
        eccdouble_waff_512(Q, PCurve); 
        if (scalar_bit==1) {
            eccadd_waff_512(P, Q, PCurve); 
        }
        j++;
    }
        
    return;
}

/******** Point conversion functions for "Ted512" Curves ********/

void ecc_Ted512_to_weierstrass(point_Ted512 Q, point_Jac512 P, PCurveStruct TedCurve)
{ // Convert point on one of the Ted curves to its corresponding isomorphic Weierstrass curve
  // Input:  Q = (xTE,yTE) on twisted Edwards curve
  //         Twisted Edwards curve struct TedCurve    
  // Output: P = (xW,yW) on Weierstrass curve
    dig i, num_words = (TedCurve->pbits+ML_WORD-1)/ML_WORD;
    dig512 t1, t2, t3, t4, a, d;

    for (i = 0; i < num_words; i++) 
    {
        a[i] = TedCurve->parameter1[i]; 
        d[i] = TedCurve->parameter2[i]; 
    }
    
    fpzero512(t1); t1[0] = 5;             // t1 = 5
    fpmul512(a, t1, t2);                  // t2 = 5a
    fpsub512(t2, d, t2);                  // t2 = 5a-d
    fpmul512(d, t1, t3);                  // t3 = 5d
    fpsub512(a, t3, t1);                  // t1 = a-5d
    fpmul512(Q->y, t1, t3);               // t3 = yTE*(a-5d)
    fpadd512(t3, t2, t2);                 // t2 = (5a-d) + yTE*(a-5d)
    fpzero512(t1); t1[0] = 1;             // t1 = 1
    fpsub512(t1, Q->y, t3);               // t3 = 1-yTE
    fpzero512(t1); t1[0] = 12;            // t1 = 12
    fpmul512(t1, t3, t4);                 // t4 = 12(1-yTE)
    fpinv512(t4);                         // t4 = 1/12(1-yTE)
    fpmul512(Q->x, t3, t1);               // t1 = xTE*(1-yTE)
    fpadd512(t1, t1, t3);                 // t3 = 2xTE*(1-yTE)
    fpadd512(t3, t3, t3);                 // t3 = 4xTE*(1-yTE)
    fpinv512(t3);                         // t3 = 1/4xTE*(1-yTE)
    fpmul512(t4, t2, P->x);               // Xfinal = ((5a-d) + yTE*(a-5d))/12(1-yTE)
    fpzero512(t1); t1[0] = 1;             // t1 = 1
    fpadd512(Q->y, t1, t1);               // t1 = yTE+1
    fpsub512(a, d, t2);                   // t2 = a-d
    fpmul512(t1, t2, t4);                 // t4 = (a-d)*(yTE+1)
    fpmul512(t4, t3, P->y);               // Yfinal = ((a-d)*(yTE+1))/4xTE*(1-yTE)
    
// cleanup
    fpzero512(t1);
    fpzero512(t2);
    fpzero512(t3);
    fpzero512(t4);
    
    return;
}


void ecc_weierstrass_to_Ted512(point_Jac512 Q, point_Ted512 P, PCurveStruct TedCurve)
{ // Convert point on isomorphic Weierstrass curve to its corresponding Ted curve
  // Input:  Q = (xW,yW) on Weierstrass curve
  //         Twisted Edwards curve struct TedCurve    
  // Output: P = (xTE,yTE) on twisted Edwards curve
    dig i, num_words = (TedCurve->pbits+ML_WORD-1)/ML_WORD;
    dig512 t1, t2, t3, a, d;

    for (i = 0; i < num_words; i++) 
    {
        a[i] = TedCurve->parameter1[i]; 
        d[i] = TedCurve->parameter2[i]; 
    }

    fpadd512(Q->x, Q->x, t1);
    fpadd512(Q->x, t1, t1);
    fpadd512(t1, t1, t1);             // t1 = 6xW
    fpsub512(t1, a, t2);              // t2 = 6xW - a
    fpsub512(t2, d, t2);              // t2 = 6xW - a - d
    fpadd512(Q->y, Q->y, t3);
    fpadd512(Q->y, t3, t3);
    fpadd512(t3, t3, t3);             // t3 = 6yW
    fpinv512(t3);                     // t3 = 1/6yW
    fpmul512(t2, t3, P->x);           // Xfinal = (6xW - a - d)/6yW
    fpadd512(t1, t1, t1);             // t1 = 12xW
    fpadd512(t1, d, t2);              // t2 = 12xW + d
    fpadd512(t1, a, t1);              // t1 = 12xW + a
    fpadd512(a, a, t3);  
    fpsub512(t2, t3, t2);             // t2 = 12xW + d - 2a   
    fpsub512(t2, t3, t2);             // t2 = 12xW + d - 4a 
    fpsub512(t2, a, t2);              // t2 = 12xW + d - 5a  
    fpadd512(d, d, t3);  
    fpsub512(t1, t3, t1);             // t1 = 12xW + a - 2d   
    fpsub512(t1, t3, t1);             // t1 = 12xW + a - 4d 
    fpsub512(t1, d, t1);              // t1 = 12xW + a - 5d          
    fpinv512(t1);                     // t1 = 1/(12xW + a - 5d)
    fpmul512(t1, t2, P->y);           // Yfinal = (12xW + d - 5a)/(12xW + a - 5d)
    
// cleanup
    fpzero512(t1);
    fpzero512(t2);
    fpzero512(t3);
    
    return;
}
#endif
