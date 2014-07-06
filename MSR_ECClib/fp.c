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
* Abstract: additional field and recoding functions
*
* This software is based on the article by Joppe Bos, Craig Costello, 
* Patrick Longa and Michael Naehrig, "Selecting elliptic curves for
* cryptography: an efficiency and security analysis", preprint available
* at http://eprint.iacr.org/2014/130.
***************************************************************************/  

#include "msr_ecclib.h"
#include "msr_ecclib_priv.h"
#include "immintrin.h"


#ifdef ECCURVES_256
//
// Specialized 256-bit field operations for curves "Jac256" and "Ted256"
//

void fpcopy256(dig256 a, dig256 c)
{ // Copy of a 256-bit field element, c = a 
    unsigned int i;

    for (i = 0; i < ML_WORDS256; i++)
    {
        c[i] = a[i]; 
    }
    return;
}


BOOL fp_iszero256(dig256 a)
{ // Is 256-bit field element zero, a=0?  
    unsigned int i;
    dig c;
    BOOL answer = TRUE;

    c = a[0];
    for (i = 1; i < ML_WORDS256; i++)
    {
        c = c | a[i]; 
    }
    answer = (BOOL)((~(0 - c) >> (ML_WORD-1)) & (~c >> (ML_WORD-1)));

    return answer;
}


BOOL mod_eval256(dig256 a, dig256 modulus)
{ // Evaluate if 256-bit element is in [0, modulus-1] 
  // eval = TRUE if 0 <= a < modulus, else eval = FALSE 
    BOOL eval = FALSE;
    dig256 t1;
    
    fpcopy256(a, t1);
    eval = fpneg256(modulus, t1);            // eval = TRUE if a <= modulus
    eval = (eval & (~fp_iszero256(t1) & 1)); // eval = TRUE if a < modulus

// cleanup
    fpzero256(t1);

    return eval;
}


void fpinv256_fixedchain(dig256 a)
{ // Inverse of field element, af = a^-1 = a^(p-2) mod p
  // Hardwired for p = 2^256-189
    int i, j;
    dig256 t1, t2, t3, t4, t5;
     
    fpsqr256(a, t1);                   // t1 = a^2                  
    fpmul256(a, t1, t2);               // t2 = a^3   
    fpsqr256(t2, t3);                  // t3 = a^6   
    fpsqr256(t3, t4);                  // t4 = a^12                  
    fpmul256(t2, t4, t5);              // t5 = a^15                   
    fpsqr256(t1, t2);                  // t2 = a^4      
    fpsqr256(t2, t1);                                      
    fpsqr256(t1, t2);                                      
    fpsqr256(t2, t1);                                      
    fpsqr256(t1, t2);                  // t2 = a^64                 
    fpmul256(a, t2, t3);               // t3 = a^65                    
    fpsqr256(t5, t2);                  // t2 = a^30                     
    fpsqr256(t2, t1);                   
    fpsqr256(t1, t2);                   
    fpsqr256(t2, t4);                  // t4 = a^240                  
    fpmul256(t5, t4, t1);              // t1 = a^255                    
    fpsqr256(t1, t2);                  // t2 = a^510                      
    fpsqr256(t2, t4);                                  
    fpsqr256(t4, t2);                                          
    fpsqr256(t2, t4);                                          
    fpsqr256(t4, t2);                                          
    fpsqr256(t2, t4);                                          
    fpsqr256(t4, t2);                                          
    fpsqr256(t2, t4);                  // t4 = a^65280                 
    fpmul256(t3, t4, t2);              // t2 = a^65345                 
    fpmul256(t1, t4, t3);              // t3 = a^65535      
    fpcopy256(t3, a);                  // af = a^65535

    for (i=0; i<14; i++) {
        for (j=0; j<16; j++) { 
            fpsqr256(a, t1); 
            fpcopy256(t1, a); 
        }                              // af = af^65536
        fpmul256(t3, t1, a);           // af = af * a^65535    
    }
    for (i=0; i<16; i++) { 
        fpsqr256(a, t1); 
        fpcopy256(t1, a); 
    }                                  // af = af^65536
    fpmul256(t2, t1, a);               // af = af * a^65345 = a^(2^256-191)
    
// cleanup
    fpzero256(t1);
    fpzero256(t2);
    fpzero256(t3);
    fpzero256(t4);
    fpzero256(t5);

    return;
}


void lut_chu_Jac256(point_chu_precomp_Jac256* table, point_chu_precomp_Jac256 P, int digit, unsigned int npoints, PCurveStruct PCurve)
{ // Constant-time table lookup to extract a Chudnovsky point (X:Y:Z:Z^2:Z^3) from the precomputed table
  // Weierstrass a=-3 curve over p = 2^256-189
  // Operation: P = sign * table[(|digit|-1)/2], where sign=1 if digit>0 and sign=-1 if digit<0
    unsigned int i;
    int sign=0, pos, mask;
    __m256d point[5], temp_point[5], full_mask;
    
    sign = ((unsigned int)digit >> (sizeof(unsigned int)*8-1)) - 1;    // if digit<0 then sign = 0x00...0 else sign = 0xFF...F
    pos = ((sign & (digit ^ -digit)) ^ -digit) >> 1;                   // position = (|digit|-1)/2  
    point[0] = _mm256_loadu_pd ((double const *) table[0]->X);         // point = table[0] 
    point[1] = _mm256_loadu_pd ((double const *) table[0]->Y);    
    point[2] = _mm256_loadu_pd ((double const *) table[0]->Z);    
    point[3] = _mm256_loadu_pd ((double const *) table[0]->Z2);   
    point[4] = _mm256_loadu_pd ((double const *) table[0]->Z3);    

    for (i=1; i<npoints; i++) 
    { 
        pos--;
        // If match then mask = 0xFF...F else sign = 0x00...0
        mask = (((unsigned int)-pos >> (sizeof(unsigned int)*8-1)) || ((unsigned int)pos >> (sizeof(unsigned int)*8-1))) - 1; 
        full_mask = _mm256_set1_pd ((double)mask);
        temp_point[0] = _mm256_loadu_pd ((double const *) table[i]->X);    // temp_point = table[i+1]
        temp_point[1] = _mm256_loadu_pd ((double const *) table[i]->Y);
        temp_point[2] = _mm256_loadu_pd ((double const *) table[i]->Z);
        temp_point[3] = _mm256_loadu_pd ((double const *) table[i]->Z2);
        temp_point[4] = _mm256_loadu_pd ((double const *) table[i]->Z3);
        // If mask = 0x00...0 then point = point, else if mask = 0xFF...F then point = temp_point
        point[0] = _mm256_blendv_pd (point[0], temp_point[0], full_mask);     
        point[1] = _mm256_blendv_pd (point[1], temp_point[1], full_mask);  
        point[2] = _mm256_blendv_pd (point[2], temp_point[2], full_mask);  
        point[3] = _mm256_blendv_pd (point[3], temp_point[3], full_mask);  
        point[4] = _mm256_blendv_pd (point[4], temp_point[4], full_mask); 
    }
    
    _mm256_storeu_pd ((double *)P->X, point[0]);    
    _mm256_storeu_pd ((double *)P->Y, point[1]);  
    _mm256_storeu_pd ((double *)P->Z, point[2]);    
    _mm256_storeu_pd ((double *)P->Z2, point[3]);    
    _mm256_storeu_pd ((double *)P->Z3, point[4]);    
    fpneg256(PCurve->prime, P->Y);                                    // point[1]: y coordinate  
    temp_point[1] = _mm256_loadu_pd ((double const *)P->Y);           // temp_point[1]: -y coordinate
    full_mask = _mm256_set1_pd ((double)sign);
    point[1] = _mm256_blendv_pd (temp_point[1], point[1], full_mask); // if mask = 0x00...0 then choose negative of the point
    _mm256_storeu_pd ((double *)P->Y, point[1]); 

    return;
}


void lut_aff_Jac256(point_Jac256* table, point_Jac256 P, int digit, int sign, unsigned int npoints, PCurveStruct PCurve)
{ // Constant-time table lookup to extract an affine point from the precomputed table
  // If (sign = 0x00...0) then final digit is positive, else if (sign = 0xFF...F) then final digit is negative
  // Weierstrass a=-3 curve over p = 2^256-189
  // Operation: if sign=0 then P = table[digit], else if (sign=-1) then P = -table[digit]
    unsigned int i;
    int pos, mask;
    __m256d point[2], temp_point[2], full_mask;
    
    pos = digit;                                                       // Load digit position.  
    point[0] = _mm256_loadu_pd ((double const *) table[0]->x);         // point = table[0] 
    point[1] = _mm256_loadu_pd ((double const *) table[0]->y);  

    for (i=1; i<npoints; i++) 
    { 
        pos--;
        // If match then mask = 0xFF...F else sign = 0x00...0
        mask = (((unsigned int)-pos >> (sizeof(unsigned int)*8-1)) || ((unsigned int)pos >> (sizeof(unsigned int)*8-1))) - 1; 
        full_mask = _mm256_set1_pd ((double)mask);
        temp_point[0] = _mm256_loadu_pd ((double const *) table[i]->x);    // temp_point = table[i+1]
        temp_point[1] = _mm256_loadu_pd ((double const *) table[i]->y);
        // If mask = 0x00...0 then point = point, else if mask = 0xFF...F then point = temp_point
        point[0] = _mm256_blendv_pd (point[0], temp_point[0], full_mask);     
        point[1] = _mm256_blendv_pd (point[1], temp_point[1], full_mask); 
    }
    
    _mm256_storeu_pd ((double *)P->x, point[0]);    
    _mm256_storeu_pd ((double *)P->y, point[1]); 
    fpneg256(PCurve->prime, P->y);                                    // point[1]: y coordinate  
    temp_point[1] = _mm256_loadu_pd ((double const *)P->y);           // temp_point[1]: -y coordinate
    full_mask = _mm256_set1_pd ((double)sign);
    point[1] = _mm256_blendv_pd (point[1], temp_point[1], full_mask); // if mask = 0xFF...F then choose negative of the point
    _mm256_storeu_pd ((double *)P->y, point[1]); 
    
    return;
}


void lut_extproj_Ted256(point_extproj_precomp_Ted256* table, point_extproj_precomp_Ted256 P, int digit, unsigned int npoints, PCurveStruct PCurve)
{ // Constant-time table lookup to extract an extended twisted Edwards point (X+Y:Y-X:2Z:2T) from the precomputed table
  // Twisted Edwards a=-1 curve over p = 2^256-189
  // Operation: P = sign * table[(|digit|-1)/2], where sign=1 if digit>0 and sign=-1 if digit<0
    unsigned int i;
    int sign=0, pos, mask;
    __m256d point[4], temp_point[4], full_mask;
    dig256 t1;

    sign = ((unsigned int)digit >> (sizeof(unsigned int)*8-1)) - 1; // if digit<0 then sign = 0x00...0 else sign = 0xFF...F
    pos = ((sign & (digit ^ -digit)) ^ -digit) >> 1;                // position = (|digit|-1)/2  
    point[0] = _mm256_loadu_pd ((double const *) table[0]->XY);     // point = table[0] 
    point[1] = _mm256_loadu_pd ((double const *) table[0]->YX);    
    point[2] = _mm256_loadu_pd ((double const *) table[0]->Z2);    
    point[3] = _mm256_loadu_pd ((double const *) table[0]->T2);    

    for (i=1; i<npoints; i++) 
    { 
        pos--;
        // If match then mask = 0xFF...F else sign = 0x00...0
        mask = (((unsigned int)-pos >> (sizeof(unsigned int)*8-1)) || ((unsigned int)pos >> (sizeof(unsigned int)*8-1))) - 1; 
        full_mask = _mm256_set1_pd ((double)mask);
        temp_point[0] = _mm256_loadu_pd ((double const *) table[i]->XY);    // temp_point = table[i+1]
        temp_point[1] = _mm256_loadu_pd ((double const *) table[i]->YX);
        temp_point[2] = _mm256_loadu_pd ((double const *) table[i]->Z2);
        temp_point[3] = _mm256_loadu_pd ((double const *) table[i]->T2);
        // If mask = 0x00...0 then point = point, else if mask = 0xFF...F then point = temp_point
        point[0] = _mm256_blendv_pd (point[0], temp_point[0], full_mask);     
        point[1] = _mm256_blendv_pd (point[1], temp_point[1], full_mask);  
        point[2] = _mm256_blendv_pd (point[2], temp_point[2], full_mask);  
        point[3] = _mm256_blendv_pd (point[3], temp_point[3], full_mask);  
    }
    
    _mm256_storeu_pd ((double *)P->XY, point[0]);    
    _mm256_storeu_pd ((double *)P->YX, point[1]);  
    _mm256_storeu_pd ((double *)P->Z2, point[2]);    
    _mm256_storeu_pd ((double *)P->T2, point[3]);    
    fpcopy256(P->XY, t1);
    fpcopy256(P->YX, P->XY);
    fpcopy256(t1, P->YX);
    fpneg256(PCurve->prime, P->T2);                                     // point: x, t coordinate  
    temp_point[0] = _mm256_loadu_pd ((double const *)P->XY);            // temp_point: -x, -t coordinate
    temp_point[1] = _mm256_loadu_pd ((double const *)P->YX);        
    temp_point[3] = _mm256_loadu_pd ((double const *)P->T2);  
    full_mask = _mm256_set1_pd ((double)sign);
    point[0] = _mm256_blendv_pd (temp_point[0], point[0], full_mask);   // if mask = 0x00...0 then choose negative of the point
    point[1] = _mm256_blendv_pd (temp_point[1], point[1], full_mask);
    point[3] = _mm256_blendv_pd (temp_point[3], point[3], full_mask);
    _mm256_storeu_pd ((double *)P->XY, point[0]); 
    _mm256_storeu_pd ((double *)P->YX, point[1]); 
    _mm256_storeu_pd ((double *)P->T2, point[3]); 

    return;
}


void lut_extaff_Ted256(point_extaff_precomp_Ted256* table, point_extaff_precomp_Ted256 P, int digit, int sign, unsigned int npoints, PCurveStruct PCurve)
{ // Constant-time table lookup to extract an extended affine point (x+y,y-x,2t) from the precomputed table
  // Twisted Edwards a=-1 curve over p = 2^256-189
  // Operation: if sign=0 then P = table[digit], else if (sign=-1) then P = -table[digit]
    unsigned int i;
    int pos, mask;
    __m256d point[3], temp_point[3], full_mask;

    pos = digit;                                                       // Load digit position.  
    point[0] = _mm256_loadu_pd ((double const *) table[0]->xy);        // point = table[0] 
    point[1] = _mm256_loadu_pd ((double const *) table[0]->yx);  
    point[2] = _mm256_loadu_pd ((double const *) table[0]->t2);  

    for (i=1; i<npoints; i++) 
    { 
        pos--;
        // If match then mask = 0xFF...F else sign = 0x00...0
        mask = (((unsigned int)-pos >> (sizeof(unsigned int)*8-1)) || ((unsigned int)pos >> (sizeof(unsigned int)*8-1))) - 1; 
        full_mask = _mm256_set1_pd ((double)mask);
        temp_point[0] = _mm256_loadu_pd ((double const *) table[i]->xy);    // temp_point = table[i+1]
        temp_point[1] = _mm256_loadu_pd ((double const *) table[i]->yx);
        temp_point[2] = _mm256_loadu_pd ((double const *) table[i]->t2);
        // If mask = 0x00...0 then point = point, else if mask = 0xFF...F then point = temp_point
        point[0] = _mm256_blendv_pd (point[0], temp_point[0], full_mask);     
        point[1] = _mm256_blendv_pd (point[1], temp_point[1], full_mask);    
        point[2] = _mm256_blendv_pd (point[2], temp_point[2], full_mask); 
    }
    
    _mm256_storeu_pd ((double *)P->xy, point[1]);    
    _mm256_storeu_pd ((double *)P->yx, point[0]);   
    _mm256_storeu_pd ((double *)P->t2, point[2]); 
    fpneg256(PCurve->prime, P->t2);                                   // point negated: -x, -t coordinate 
    temp_point[0] = _mm256_loadu_pd ((double const *)P->xy);          // temp_point is point negated
    temp_point[1] = _mm256_loadu_pd ((double const *)P->yx);       
    temp_point[2] = _mm256_loadu_pd ((double const *)P->t2);       
    full_mask = _mm256_set1_pd ((double)sign);
    point[0] = _mm256_blendv_pd (point[0], temp_point[0], full_mask); // if mask = 0xFF...F then choose negative of the point
    point[1] = _mm256_blendv_pd (point[1], temp_point[1], full_mask);
    point[2] = _mm256_blendv_pd (point[2], temp_point[2], full_mask);
    _mm256_storeu_pd ((double *)P->xy, point[0]); 
    _mm256_storeu_pd ((double *)P->yx, point[1]); 
    _mm256_storeu_pd ((double *)P->t2, point[2]); 

    return;
}
#endif


#ifdef ECCURVES_384
//
// Specialized 384-bit field operations for curves "Jac384" and "Ted384"
//

void fpcopy384(dig384 a, dig384 c)
{ // Copy of a 384-bit field element, c = a 
    unsigned int i;

    for (i = 0; i < ML_WORDS384; i++)
    {
        c[i] = a[i]; 
    }
    return;
}


BOOL fp_iszero384(dig384 a)
{ // Is 384-bit field element zero, a=0?  
    unsigned int i;
    dig c;
    BOOL answer = TRUE;

    c = a[0];
    for (i = 1; i < ML_WORDS384; i++)
    {
        c = c | a[i]; 
    }
    answer = (BOOL)((~(0 - c) >> (ML_WORD-1)) & (~c >> (ML_WORD-1)));

    return answer;
}


BOOL mod_eval384(dig384 a, dig384 modulus)
{ // Evaluate if 384-bit element is in [0, modulus-1] 
  // eval = TRUE if 0 <= a < modulus, else eval = FALSE 
    BOOL eval = FALSE;
    dig384 t1;
    
    fpcopy384(a, t1);
    eval = fpneg384(modulus, t1);            // eval = TRUE if a <= modulus
    eval = (eval & (~fp_iszero384(t1) & 1)); // eval = TRUE if a < modulus

// cleanup
    fpzero384(t1);

    return eval;
}


void fpinv384_fixedchain(dig384 a)
{ // Inverse of field element, af = a^-1 = a^(p-2) mod p
  // Hardwired for p = 2^384-319
    int j;
    dig384 t3, t12, tF, T, o10, o40, aux;

    fpsqr384(a, aux);           // a = t1
    fpmul384(a, aux, t3);       // t3
    fpsqr384(t3, aux); 
    fpsqr384(aux, t12);         // t12
    fpmul384(t3, t12, tF);      // tF
    fpsqr384(tF, T); 
    fpsqr384(T, aux);  
    fpsqr384(aux, T);  
    fpsqr384(T, aux); 
    fpmul384(tF, aux, T);       // T
    fpsqr384(T, aux);  
    fpsqr384(aux, T); 
    fpmul384(t3, T, o10);       // o10
    fpcopy384(o10, aux);
    for (j=0; j<10; j++) { 
        fpsqr384(aux, T); 
        fpcopy384(T, aux); 
    }
    fpmul384(o10, aux, T);
    fpcopy384(T, aux);
    for (j=0; j<20; j++) { 
        fpsqr384(aux, o40); 
        fpcopy384(o40, aux); 
    }
    fpmul384(T, aux, o40);
    fpcopy384(o40, aux);
    for (j=0; j<40; j++) { 
        fpsqr384(aux, T); 
        fpcopy384(T, aux); 
    }
    fpmul384(o40, aux, T);
    fpcopy384(T, aux);
    for (j=0; j<80; j++) { 
        fpsqr384(aux, t12); 
        fpcopy384(t12, aux); 
    }
    fpmul384(T, t12, aux);
    fpcopy384(aux, T);
    for (j=0; j<160; j++) { 
        fpsqr384(aux, t12); 
        fpcopy384(t12, aux); 
    }
    fpmul384(T, t12, aux);
    for (j=0; j<10; j++) { 
        fpsqr384(aux, T); 
        fpcopy384(T, aux); 
    }
    fpmul384(o10, T, aux);
    for (j=0; j<40; j++) { 
        fpsqr384(aux, T); 
        fpcopy384(T, aux); 
    }
    fpmul384(o40, aux, T);
    fpsqr384(T, aux); 
    fpsqr384(aux, T);  
    fpsqr384(T, aux);  
    fpsqr384(aux, T); 
    fpmul384(tF, T, aux);
    fpsqr384(aux, T); 
    fpmul384(a, T, aux); 
    fpsqr384(aux, T);  
    fpsqr384(T, aux);  
    fpsqr384(aux, T); 
    fpmul384(t3, T, aux);
    for (j=0; j<6; j++) { 
        fpsqr384(aux, T); 
        fpcopy384(T, aux); 
    }
    fpcopy384(a, aux);
    fpmul384(T, aux, a);
    
// cleanup
    fpzero384(t3);
    fpzero384(t12);
    fpzero384(tF);
    fpzero384(T);
    fpzero384(o10);
    fpzero384(o40);
    fpzero384(aux);

    return;
}


void lut_chu_Jac384(point_chu_precomp_Jac384* table, point_chu_precomp_Jac384 P, int digit, unsigned int npoints, PCurveStruct PCurve)
{ // Constant-time table lookup to extract a Chudnovsky point (X:Y:Z:Z^2:Z^3) from the precomputed table
  // Weierstrass a=-3 curve over p = 2^384-319
  // Operation: P = sign * table[(|digit|-1)/2], where sign=1 if digit>0 and sign=-1 if digit<0
    unsigned int i;
    int sign=0, pos, mask;
    __m256d point[5], temp_point[5], full_mask;
    __m128d pointb[5], temp_pointb[5], full_maskb;
    
    sign = ((unsigned int)digit >> (sizeof(unsigned int)*8-1)) - 1;    // if digit<0 then sign = 0x00...0 else sign = 0xFF...F
    pos = ((sign & (digit ^ -digit)) ^ -digit) >> 1;                   // position = (|digit|-1)/2  
    point[0] = _mm256_loadu_pd ((double const *) table[0]->X);         // point = table[0] 
    pointb[0] = _mm_loadu_pd ((double const *) table[0]->X+4);          
    point[1] = _mm256_loadu_pd ((double const *) table[0]->Y);  
    pointb[1] = _mm_loadu_pd ((double const *) table[0]->Y+4);            
    point[2] = _mm256_loadu_pd ((double const *) table[0]->Z);  
    pointb[2] = _mm_loadu_pd ((double const *) table[0]->Z+4);            
    point[3] = _mm256_loadu_pd ((double const *) table[0]->Z2);  
    pointb[3] = _mm_loadu_pd ((double const *) table[0]->Z2+4);           
    point[4] = _mm256_loadu_pd ((double const *) table[0]->Z3);  
    pointb[4] = _mm_loadu_pd ((double const *) table[0]->Z3+4);            

    for (i=1; i<npoints; i++) 
    { 
        pos--;
        // If match then mask = 0xFF...F else sign = 0x00...0
        mask = (((unsigned int)-pos >> (sizeof(unsigned int)*8-1)) || ((unsigned int)pos >> (sizeof(unsigned int)*8-1))) - 1; 
        full_mask = _mm256_set1_pd ((double)mask);
        full_maskb = _mm_set1_pd ((double)mask);
        temp_point[0] = _mm256_loadu_pd ((double const *) table[i]->X);    // temp_point = table[i+1]
        temp_pointb[0] = _mm_loadu_pd ((double const *) table[i]->X+4); 
        temp_point[1] = _mm256_loadu_pd ((double const *) table[i]->Y);
        temp_pointb[1] = _mm_loadu_pd ((double const *) table[i]->Y+4); 
        temp_point[2] = _mm256_loadu_pd ((double const *) table[i]->Z);
        temp_pointb[2] = _mm_loadu_pd ((double const *) table[i]->Z+4); 
        temp_point[3] = _mm256_loadu_pd ((double const *) table[i]->Z2);
        temp_pointb[3] = _mm_loadu_pd ((double const *) table[i]->Z2+4); 
        temp_point[4] = _mm256_loadu_pd ((double const *) table[i]->Z3);
        temp_pointb[4] = _mm_loadu_pd ((double const *) table[i]->Z3+4); 
        // If mask = 0x00...0 then point = point, else if mask = 0xFF...F then point = temp_point
        point[0] = _mm256_blendv_pd (point[0], temp_point[0], full_mask);    
        pointb[0] = _mm_blendv_pd (pointb[0], temp_pointb[0], full_maskb);      
        point[1] = _mm256_blendv_pd (point[1], temp_point[1], full_mask);   
        pointb[1] = _mm_blendv_pd (pointb[1], temp_pointb[1], full_maskb);   
        point[2] = _mm256_blendv_pd (point[2], temp_point[2], full_mask);    
        pointb[2] = _mm_blendv_pd (pointb[2], temp_pointb[2], full_maskb);  
        point[3] = _mm256_blendv_pd (point[3], temp_point[3], full_mask);    
        pointb[3] = _mm_blendv_pd (pointb[3], temp_pointb[3], full_maskb);  
        point[4] = _mm256_blendv_pd (point[4], temp_point[4], full_mask);   
        pointb[4] = _mm_blendv_pd (pointb[4], temp_pointb[4], full_maskb);  
    }
    
    _mm256_storeu_pd ((double *)P->X, point[0]);  
    _mm_storeu_pd ((double *)P->X+4, pointb[0]);   
    _mm256_storeu_pd ((double *)P->Y, point[1]);   
    _mm_storeu_pd ((double *)P->Y+4, pointb[1]);   
    _mm256_storeu_pd ((double *)P->Z, point[2]);  
    _mm_storeu_pd ((double *)P->Z+4, pointb[2]);      
    _mm256_storeu_pd ((double *)P->Z2, point[3]);  
    _mm_storeu_pd ((double *)P->Z2+4, pointb[3]);      
    _mm256_storeu_pd ((double *)P->Z3, point[4]);   
    _mm_storeu_pd ((double *)P->Z3+4, pointb[4]);     
    fpneg384(PCurve->prime, P->Y);                                     // point[1]: y coordinate  
    temp_point[1] = _mm256_loadu_pd ((double const *)P->Y);            // temp_point[1]: -y coordinate
    temp_pointb[1] = _mm_loadu_pd ((double const *)P->Y+4);     
    full_mask = _mm256_set1_pd ((double)sign);    
    full_maskb = _mm_set1_pd ((double)sign);
    point[1] = _mm256_blendv_pd (temp_point[1], point[1], full_mask);  // if mask = 0x00...0 then choose negative of the point
    pointb[1] = _mm_blendv_pd (temp_pointb[1], pointb[1], full_maskb);
    _mm256_storeu_pd ((double *)P->Y, point[1]); 
    _mm_storeu_pd ((double *)P->Y+4, pointb[1]);  

    return;
}


void lut_aff_Jac384(point_Jac384* table, point_Jac384 P, int digit, int sign, unsigned int npoints, PCurveStruct PCurve)
{ // Constant-time table lookup to extract an affine point from the precomputed table
  // If (sign = 0x00...0) then final digit is positive, else if (sign = 0xFF...F) then final digit is negative
  // Weierstrass a=-3 curve over p = 2^384-319
  // Operation: if sign=0 then P = table[digit], else if (sign=-1) then P = -table[digit]
    unsigned int i;
    int pos, mask;
    __m256d point[2], temp_point[2], full_mask;
    __m128d pointb[2], temp_pointb[2], full_maskb;
            
    pos = digit;                                                       // Load digit position.  
    point[0] = _mm256_loadu_pd ((double const *) table[0]->x);         // point = table[0] 
    pointb[0] = _mm_loadu_pd ((double const *) table[0]->x+4);          
    point[1] = _mm256_loadu_pd ((double const *) table[0]->y);  
    pointb[1] = _mm_loadu_pd ((double const *) table[0]->y+4);            

    for (i=1; i<npoints; i++) 
    { 
        pos--;
        // If match then mask = 0xFF...F else sign = 0x00...0
        mask = (((unsigned int)-pos >> (sizeof(unsigned int)*8-1)) || ((unsigned int)pos >> (sizeof(unsigned int)*8-1))) - 1; 
        full_mask = _mm256_set1_pd ((double)mask);
        full_maskb = _mm_set1_pd ((double)mask);
        temp_point[0] = _mm256_loadu_pd ((double const *) table[i]->x);    // temp_point = table[i+1]
        temp_pointb[0] = _mm_loadu_pd ((double const *) table[i]->x+4); 
        temp_point[1] = _mm256_loadu_pd ((double const *) table[i]->y);
        temp_pointb[1] = _mm_loadu_pd ((double const *) table[i]->y+4);
        // If mask = 0x00...0 then point = point, else if mask = 0xFF...F then point = temp_point
        point[0] = _mm256_blendv_pd (point[0], temp_point[0], full_mask);    
        pointb[0] = _mm_blendv_pd (pointb[0], temp_pointb[0], full_maskb);      
        point[1] = _mm256_blendv_pd (point[1], temp_point[1], full_mask);   
        pointb[1] = _mm_blendv_pd (pointb[1], temp_pointb[1], full_maskb); 
    }
    
    _mm256_storeu_pd ((double *)P->x, point[0]);  
    _mm_storeu_pd ((double *)P->x+4, pointb[0]);   
    _mm256_storeu_pd ((double *)P->y, point[1]);   
    _mm_storeu_pd ((double *)P->y+4, pointb[1]); 
    fpneg384(PCurve->prime, P->y);                                    // point[1]: y coordinate  
    temp_point[1] = _mm256_loadu_pd ((double const *)P->y);           // temp_point[1]: -y coordinate
    temp_pointb[1] = _mm_loadu_pd ((double const *)P->y+4);     
    full_mask = _mm256_set1_pd ((double)sign);    
    full_maskb = _mm_set1_pd ((double)sign);
    point[1] = _mm256_blendv_pd (point[1], temp_point[1], full_mask); // if mask = 0xFF...F then choose negative of the point
    pointb[1] = _mm_blendv_pd (pointb[1], temp_pointb[1], full_maskb);
    _mm256_storeu_pd ((double *)P->y, point[1]); 
    _mm_storeu_pd ((double *)P->y+4, pointb[1]); 

    return;
}


void lut_extproj_Ted384(point_extproj_precomp_Ted384* table, point_extproj_precomp_Ted384 P, int digit, unsigned int npoints, PCurveStruct PCurve)
{ // Constant-time table lookup to extract an extended twisted Edwards point (X+Y:Y-X:2Z:2T) from the precomputed table
  // Twisted Edwards a=-1 curve over p = 2^384-319
  // Operation: P = sign * table[(|digit|-1)/2], where sign=1 if digit>0 and sign=-1 if digit<0
    unsigned int i;
    int sign=0, pos, mask;
    __m256d point[6], temp_point[6], full_mask;
    __m128d pointb[3], temp_pointb[3], full_maskb;
    dig384 t1;

    sign = ((unsigned int)digit >> (sizeof(unsigned int)*8-1)) - 1; // if digit<0 then sign = 0x00...0 else sign = 0xFF...F
    pos = ((sign & (digit ^ -digit)) ^ -digit) >> 1;                // position = (|digit|-1)/2  
    point[0] = _mm256_loadu_pd ((double const *) table[0]->XY);     // point = table[0]  
    point[1] = _mm256_loadu_pd ((double const *) table[0]->XY+4);  
    point[2] = _mm256_loadu_pd ((double const *) table[0]->XY+2*4);  
    point[3] = _mm256_loadu_pd ((double const *) table[0]->XY+3*4);  
    point[4] = _mm256_loadu_pd ((double const *) table[0]->XY+4*4);  
    point[5] = _mm256_loadu_pd ((double const *) table[0]->XY+5*4);    // 384*4 coord = 1536 and 1536/256 = 6   

    for (i=1; i<npoints; i++) 
    { 
        pos--;
        // If match then mask = 0xFF...F else sign = 0x00...0
        mask = (((unsigned int)-pos >> (sizeof(unsigned int)*8-1)) || ((unsigned int)pos >> (sizeof(unsigned int)*8-1))) - 1; 
        full_mask = _mm256_set1_pd ((double)mask);
        temp_point[0] = _mm256_loadu_pd ((double const *) table[i]->XY);    // temp_point = table[i+1]
        temp_point[1] = _mm256_loadu_pd ((double const *) table[i]->XY+4);
        temp_point[2] = _mm256_loadu_pd ((double const *) table[i]->XY+2*4);
        temp_point[3] = _mm256_loadu_pd ((double const *) table[i]->XY+3*4);
        temp_point[4] = _mm256_loadu_pd ((double const *) table[i]->XY+4*4);
        temp_point[5] = _mm256_loadu_pd ((double const *) table[i]->XY+5*4);
        // If mask = 0x00...0 then point = point, else if mask = 0xFF...F then point = temp_point
        point[0] = _mm256_blendv_pd (point[0], temp_point[0], full_mask);     
        point[1] = _mm256_blendv_pd (point[1], temp_point[1], full_mask);  
        point[2] = _mm256_blendv_pd (point[2], temp_point[2], full_mask);  
        point[3] = _mm256_blendv_pd (point[3], temp_point[3], full_mask);   
        point[4] = _mm256_blendv_pd (point[4], temp_point[4], full_mask);  
        point[5] = _mm256_blendv_pd (point[5], temp_point[5], full_mask); 
    }
    
    _mm256_storeu_pd ((double *)P->XY, point[0]);    
    _mm256_storeu_pd ((double *)P->XY+4, point[1]); 
    _mm256_storeu_pd ((double *)P->XY+2*4, point[2]);    
    _mm256_storeu_pd ((double *)P->XY+3*4, point[3]);      
    _mm256_storeu_pd ((double *)P->XY+4*4, point[4]);    
    _mm256_storeu_pd ((double *)P->XY+5*4, point[5]);  
    pointb[0] = _mm_loadu_pd ((double const *)P->XY+4); 
    point[1] = _mm256_loadu_pd ((double const *)P->YX);   
    pointb[1] = _mm_loadu_pd ((double const *)P->YX+4);  
    point[2] = _mm256_loadu_pd ((double const *)P->T2); 
    pointb[2] = _mm_loadu_pd ((double const *)P->T2+4);  
    fpcopy384(P->XY, t1);
    fpcopy384(P->YX, P->XY);
    fpcopy384(t1, P->YX);
    fpneg384(PCurve->prime, P->T2);                                     // point: x, t coordinate  
    temp_point[0] = _mm256_loadu_pd ((double const *)P->XY);            // temp_point: -x, -t coordinate
    temp_pointb[0] = _mm_loadu_pd ((double const *)P->XY+4);       
    temp_point[1] = _mm256_loadu_pd ((double const *)P->YX);        
    temp_pointb[1] = _mm_loadu_pd ((double const *)P->YX+4);       
    temp_point[2] = _mm256_loadu_pd ((double const *)P->T2);  
    temp_pointb[2] = _mm_loadu_pd ((double const *)P->T2+4);       
    full_mask = _mm256_set1_pd ((double)sign);      
    full_maskb = _mm_set1_pd ((double)sign);
    point[0] = _mm256_blendv_pd (temp_point[0], point[0], full_mask);   // if mask = 0x00...0 then choose negative of the point
    pointb[0] = _mm_blendv_pd (temp_pointb[0], pointb[0], full_maskb); 
    point[1] = _mm256_blendv_pd (temp_point[1], point[1], full_mask);
    pointb[1] = _mm_blendv_pd (temp_pointb[1], pointb[1], full_maskb); 
    point[2] = _mm256_blendv_pd (temp_point[2], point[2], full_mask);
    pointb[2] = _mm_blendv_pd (temp_pointb[2], pointb[2], full_maskb); 
    _mm256_storeu_pd ((double *)P->XY, point[0]); 
    _mm_storeu_pd ((double *)P->XY+4, pointb[0]); 
    _mm256_storeu_pd ((double *)P->YX, point[1]); 
    _mm_storeu_pd ((double *)P->YX+4, pointb[1]);
    _mm256_storeu_pd ((double *)P->T2, point[2]); 
    _mm_storeu_pd ((double *)P->T2+4, pointb[2]);

    return;
}


void lut_extaff_Ted384(point_extaff_precomp_Ted384* table, point_extaff_precomp_Ted384 P, int digit, int sign, unsigned int npoints, PCurveStruct PCurve)
{ // Constant-time table lookup to extract an extended affine point (x+y,y-x,2t) from the precomputed table
  // Twisted Edwards a=-1 curve over p = 2^384-319
  // Operation: if sign=0 then P = table[digit], else if (sign=-1) then P = -table[digit]
    unsigned int i;
    int pos, mask;
    __m256d point[3], temp_point[3], full_mask;
    __m128d pointb[3], temp_pointb[3], full_maskb;
            
    pos = digit;                                                       // Load digit position.  
    point[0] = _mm256_loadu_pd ((double const *) table[0]->xy);        // point = table[0] 
    pointb[0] = _mm_loadu_pd ((double const *) table[0]->xy+4);          
    point[1] = _mm256_loadu_pd ((double const *) table[0]->yx);  
    pointb[1] = _mm_loadu_pd ((double const *) table[0]->yx+4);          
    point[2] = _mm256_loadu_pd ((double const *) table[0]->t2);  
    pointb[2] = _mm_loadu_pd ((double const *) table[0]->t2+4);            

    for (i=1; i<npoints; i++) 
    { 
        pos--;
        // If match then mask = 0xFF...F else sign = 0x00...0
        mask = (((unsigned int)-pos >> (sizeof(unsigned int)*8-1)) || ((unsigned int)pos >> (sizeof(unsigned int)*8-1))) - 1; 
        full_mask = _mm256_set1_pd ((double)mask);
        full_maskb = _mm_set1_pd ((double)mask);
        temp_point[0] = _mm256_loadu_pd ((double const *) table[i]->xy);    // temp_point = table[i+1]
        temp_pointb[0] = _mm_loadu_pd ((double const *) table[i]->xy+4); 
        temp_point[1] = _mm256_loadu_pd ((double const *) table[i]->yx);
        temp_pointb[1] = _mm_loadu_pd ((double const *) table[i]->yx+4);
        temp_point[2] = _mm256_loadu_pd ((double const *) table[i]->t2);
        temp_pointb[2] = _mm_loadu_pd ((double const *) table[i]->t2+4);
        // If mask = 0x00...0 then point = point, else if mask = 0xFF...F then point = temp_point
        point[0] = _mm256_blendv_pd (point[0], temp_point[0], full_mask);    
        pointb[0] = _mm_blendv_pd (pointb[0], temp_pointb[0], full_maskb);      
        point[1] = _mm256_blendv_pd (point[1], temp_point[1], full_mask);   
        pointb[1] = _mm_blendv_pd (pointb[1], temp_pointb[1], full_maskb);     
        point[2] = _mm256_blendv_pd (point[2], temp_point[2], full_mask);   
        pointb[2] = _mm_blendv_pd (pointb[2], temp_pointb[2], full_maskb); 
    }
    
    _mm256_storeu_pd ((double *)P->xy, point[1]);  
    _mm_storeu_pd ((double *)P->xy+4, pointb[1]);   
    _mm256_storeu_pd ((double *)P->yx, point[0]);   
    _mm_storeu_pd ((double *)P->yx+4, pointb[0]);   
    _mm256_storeu_pd ((double *)P->t2, point[2]);   
    _mm_storeu_pd ((double *)P->t2+4, pointb[2]); 
    fpneg384(PCurve->prime, P->t2);                                    // point negated: -x, -t coordinate 
    temp_point[0] = _mm256_loadu_pd ((double const *)P->xy);           // temp_point is point negated
    temp_pointb[0] = _mm_loadu_pd ((double const *)P->xy+4);     
    temp_point[1] = _mm256_loadu_pd ((double const *)P->yx);           
    temp_pointb[1] = _mm_loadu_pd ((double const *)P->yx+4);    
    temp_point[2] = _mm256_loadu_pd ((double const *)P->t2);           
    temp_pointb[2] = _mm_loadu_pd ((double const *)P->t2+4);    
    full_mask = _mm256_set1_pd ((double)sign);    
    full_maskb = _mm_set1_pd ((double)sign);
    point[0] = _mm256_blendv_pd (point[0], temp_point[0], full_mask);  // if mask = 0xFF...F then choose negative of the point
    pointb[0] = _mm_blendv_pd (pointb[0], temp_pointb[0], full_maskb);
    point[1] = _mm256_blendv_pd (point[1], temp_point[1], full_mask); 
    pointb[1] = _mm_blendv_pd (pointb[1], temp_pointb[1], full_maskb);
    point[2] = _mm256_blendv_pd (point[2], temp_point[2], full_mask); 
    pointb[2] = _mm_blendv_pd (pointb[2], temp_pointb[2], full_maskb);
    _mm256_storeu_pd ((double *)P->xy, point[0]); 
    _mm_storeu_pd ((double *)P->xy+4, pointb[0]); 
    _mm256_storeu_pd ((double *)P->yx, point[1]); 
    _mm_storeu_pd ((double *)P->yx+4, pointb[1]); 
    _mm256_storeu_pd ((double *)P->t2, point[2]); 
    _mm_storeu_pd ((double *)P->t2+4, pointb[2]); 

    return;
}
#endif


#ifdef ECCURVES_512
//
// Specialized 512-bit field operations for curves "Jac512" and "Ted512"
//

void fpcopy512(dig512 a, dig512 c)
{ // Copy of a 512-bit field element, c = a 
    unsigned int i;

    for (i = 0; i < ML_WORDS512; i++)
    {
        c[i] = a[i]; 
    }
    return;
}


BOOL fp_iszero512(dig512 a)
{ // Is 512-bit field element zero, a=0?  
    unsigned int i;
    dig c;
    BOOL answer = TRUE;

    c = a[0];
    for (i = 1; i < ML_WORDS512; i++)
    {
        c = c | a[i]; 
    }
    answer = (BOOL)((~(0 - c) >> (ML_WORD-1)) & (~c >> (ML_WORD-1)));

    return answer;
}


BOOL mod_eval512(dig512 a, dig512 modulus)
{ // Evaluate if 512-bit element is in [0, modulus-1] 
  // eval = TRUE if 0 <= a < modulus, else eval = FALSE 
    BOOL eval = FALSE;
    dig512 t1;
    
    fpcopy512(a, t1);
    eval = fpneg512(modulus, t1);            // eval = TRUE if a <= modulus
    eval = (eval & (~fp_iszero512(t1) & 1)); // eval = TRUE if a < modulus

// cleanup
    fpzero512(t1);

    return eval;
}


void fpinv512_fixedchain(dig512 a)
{ // Inverse of field element, af = a^-1 = a^(p-2) mod p
  // Hardwired for p = 2^512-319
    int j;
    dig512 t2, T, t5, t7, t10, t80, aux, aux2;

    fpsqr512(a, t2);
    fpsqr512(t2, aux);
    fpmul512(a, aux, t5);
    fpmul512(t2, t5, t7);
    fpmul512(t2, a, T);
    fpcopy512(T, t80);
    fpsqr512(T, aux2);
    fpsqr512(aux2, aux);
    fpmul512(t80, aux, T);
    fpsqr512(T, aux);
    fpmul512(a, aux, T);
    fpcopy512(T, aux);
    for (j=0; j<5; j++) { 
        fpsqr512(aux, t10); 
        fpcopy512(t10, aux); 
    }
    fpmul512(T, aux, t10);
    fpcopy512(t10, aux);
    for (j=0; j<10; j++) { 
        fpsqr512(aux, T); 
        fpcopy512(T, aux); 
    }
    fpmul512(t10, aux, T);
    fpcopy512(T, aux);
    fpcopy512(T, aux2);
    for (j=0; j<20; j++) { 
        fpsqr512(aux, T); 
        fpcopy512(T, aux); 
    }
    fpmul512(aux2, aux, T);
    fpcopy512(T, aux);
    for (j=0; j<40; j++) { 
        fpsqr512(aux, t80); 
        fpcopy512(t80, aux); 
    }
    fpmul512(T, aux, t80);
    fpcopy512(t80, aux);
    for (j=0; j<80; j++) { 
        fpsqr512(aux, aux2); 
        fpcopy512(aux2, aux); 
    }
    fpmul512(t80, aux, T);
    fpcopy512(T, aux);
    for (j=0; j<80; j++) { 
        fpsqr512(aux, aux2); 
        fpcopy512(aux2, aux); 
    }
    fpmul512(t80, aux, T);
    fpcopy512(T, aux);
    for (j=0; j<10; j++) { 
        fpsqr512(aux, aux2); 
        fpcopy512(aux2, aux); 
    }
    fpmul512(t10, aux, T);
    fpsqr512(T, aux);
    fpmul512(a, aux, T);
    fpcopy512(T, aux);
    fpcopy512(T, aux2);
    for (j=0; j<251; j++) { 
        fpsqr512(aux, t80); 
        fpcopy512(t80, aux); 
    }
    fpmul512(aux2, aux, T);
    for (j=0; j<4; j++) { 
        fpsqr512(T, aux); 
        fpcopy512(aux, T); 
    }
    fpmul512(t7, aux, T);
    for (j=0; j<6; j++) { 
        fpsqr512(T, aux); 
        fpcopy512(aux, T); 
    }
    fpmul512(t5, aux, a);
    
// cleanup
    fpzero512(t2);
    fpzero512(T);
    fpzero512(t5);
    fpzero512(t7);
    fpzero512(t10);
    fpzero512(t80);
    fpzero512(aux);
    fpzero512(aux2);

    return;
}


void lut_chu_Jac512(point_chu_precomp_Jac512* table, point_chu_precomp_Jac512 P, int digit, unsigned int npoints, PCurveStruct PCurve)
{ // Constant-time table lookup to extract a Chudnovsky point (X:Y:Z:Z^2:Z^3) from the precomputed table
  // Weierstrass a=-3 curve over p = 2^512-319
  // Operation: P = sign * table[(|digit|-1)/2], where sign=1 if digit>0 and sign=-1 if digit<0
    unsigned int i;
    int sign=0, pos, mask;
    __m256d point[5], temp_point[5], point2[5], temp_point2[5], full_mask;
    
    sign = ((unsigned int)digit >> (sizeof(unsigned int)*8-1)) - 1;    // if digit<0 then sign = 0x00...0 else sign = 0xFF...F
    pos = ((sign & (digit ^ -digit)) ^ -digit) >> 1;                   // position = (|digit|-1)/2  
    point[0] = _mm256_loadu_pd ((double const *) table[0]->X);         // point = table[0] 
    point2[0] = _mm256_loadu_pd ((double const *) table[0]->X+4);          
    point[1] = _mm256_loadu_pd ((double const *) table[0]->Y);  
    point2[1] = _mm256_loadu_pd ((double const *) table[0]->Y+4);            
    point[2] = _mm256_loadu_pd ((double const *) table[0]->Z);  
    point2[2] = _mm256_loadu_pd ((double const *) table[0]->Z+4);            
    point[3] = _mm256_loadu_pd ((double const *) table[0]->Z2);  
    point2[3] = _mm256_loadu_pd ((double const *) table[0]->Z2+4);           
    point[4] = _mm256_loadu_pd ((double const *) table[0]->Z3);  
    point2[4] = _mm256_loadu_pd ((double const *) table[0]->Z3+4);            

    for (i=1; i<npoints; i++) 
    { 
        pos--;
        // If match then mask = 0xFF...F else sign = 0x00...0
        mask = (((unsigned int)-pos >> (sizeof(unsigned int)*8-1)) || ((unsigned int)pos >> (sizeof(unsigned int)*8-1))) - 1; 
        full_mask = _mm256_set1_pd ((double)mask);
        temp_point[0] = _mm256_loadu_pd ((double const *) table[i]->X);    // temp_point = table[i+1]
        temp_point2[0] = _mm256_loadu_pd ((double const *) table[i]->X+4); 
        temp_point[1] = _mm256_loadu_pd ((double const *) table[i]->Y);
        temp_point2[1] = _mm256_loadu_pd ((double const *) table[i]->Y+4); 
        temp_point[2] = _mm256_loadu_pd ((double const *) table[i]->Z);
        temp_point2[2] = _mm256_loadu_pd ((double const *) table[i]->Z+4); 
        temp_point[3] = _mm256_loadu_pd ((double const *) table[i]->Z2);
        temp_point2[3] = _mm256_loadu_pd ((double const *) table[i]->Z2+4); 
        temp_point[4] = _mm256_loadu_pd ((double const *) table[i]->Z3);
        temp_point2[4] = _mm256_loadu_pd ((double const *) table[i]->Z3+4); 
        // If mask = 0x00...0 then point = point, else if mask = 0xFF...F then point = temp_point
        point[0] = _mm256_blendv_pd (point[0], temp_point[0], full_mask);    
        point2[0] = _mm256_blendv_pd (point2[0], temp_point2[0], full_mask);      
        point[1] = _mm256_blendv_pd (point[1], temp_point[1], full_mask);   
        point2[1] = _mm256_blendv_pd (point2[1], temp_point2[1], full_mask);   
        point[2] = _mm256_blendv_pd (point[2], temp_point[2], full_mask);    
        point2[2] = _mm256_blendv_pd (point2[2], temp_point2[2], full_mask);  
        point[3] = _mm256_blendv_pd (point[3], temp_point[3], full_mask);    
        point2[3] = _mm256_blendv_pd (point2[3], temp_point2[3], full_mask);  
        point[4] = _mm256_blendv_pd (point[4], temp_point[4], full_mask);   
        point2[4] = _mm256_blendv_pd (point2[4], temp_point2[4], full_mask);  
    }
    
    _mm256_storeu_pd ((double *)P->X, point[0]);  
    _mm256_storeu_pd ((double *)P->X+4, point2[0]);   
    _mm256_storeu_pd ((double *)P->Y, point[1]);   
    _mm256_storeu_pd ((double *)P->Y+4, point2[1]);   
    _mm256_storeu_pd ((double *)P->Z, point[2]);  
    _mm256_storeu_pd ((double *)P->Z+4, point2[2]);      
    _mm256_storeu_pd ((double *)P->Z2, point[3]);  
    _mm256_storeu_pd ((double *)P->Z2+4, point2[3]);      
    _mm256_storeu_pd ((double *)P->Z3, point[4]);   
    _mm256_storeu_pd ((double *)P->Z3+4, point2[4]);     
    fpneg512(PCurve->prime, P->Y);                                    // point[1]: y coordinate  
    temp_point[1] = _mm256_loadu_pd ((double const *)P->Y);           // temp_point[1]: -y coordinate
    temp_point2[1] = _mm256_loadu_pd ((double const *)P->Y+4);     
    full_mask = _mm256_set1_pd ((double)sign); 
    point[1] = _mm256_blendv_pd (temp_point[1], point[1], full_mask); // if mask = 0x00...0 then choose negative of the point
    point2[1] = _mm256_blendv_pd (temp_point2[1], point2[1], full_mask);
    _mm256_storeu_pd ((double *)P->Y, point[1]); 
    _mm256_storeu_pd ((double *)P->Y+4, point2[1]); 

    return;
}


void lut_aff_Jac512(point_Jac512* table, point_Jac512 P, int digit, int sign, unsigned int npoints, PCurveStruct PCurve)
{ // Constant-time table lookup to extract an affine point from the precomputed table
  // If (sign = 0x00...0) then final digit is positive, else if (sign = 0xFF...F) then final digit is negative
  // Weierstrass a=-3 curve over p = 2^512-319
  // Operation: if sign=0 then P = table[digit], else if (sign=-1) then P = -table[digit]
    unsigned int i;
    int pos, mask;
    __m256d point[25], temp_point[2], point2[2], temp_point2[2], full_mask;
                           
    pos = digit;                                                       // Load digit position.  
    point[0] = _mm256_loadu_pd ((double const *) table[0]->x);         // point = table[0] 
    point2[0] = _mm256_loadu_pd ((double const *) table[0]->x+4);          
    point[1] = _mm256_loadu_pd ((double const *) table[0]->y);  
    point2[1] = _mm256_loadu_pd ((double const *) table[0]->y+4);            

    for (i=1; i<npoints; i++) 
    { 
        pos--;
        // If match then mask = 0xFF...F else sign = 0x00...0
        mask = (((unsigned int)-pos >> (sizeof(unsigned int)*8-1)) || ((unsigned int)pos >> (sizeof(unsigned int)*8-1))) - 1; 
        full_mask = _mm256_set1_pd ((double)mask);
        temp_point[0] = _mm256_loadu_pd ((double const *) table[i]->x);    // temp_point = table[i+1]
        temp_point2[0] = _mm256_loadu_pd ((double const *) table[i]->x+4); 
        temp_point[1] = _mm256_loadu_pd ((double const *) table[i]->y);
        temp_point2[1] = _mm256_loadu_pd ((double const *) table[i]->y+4);
        // If mask = 0x00...0 then point = point, else if mask = 0xFF...F then point = temp_point
        point[0] = _mm256_blendv_pd (point[0], temp_point[0], full_mask);    
        point2[0] = _mm256_blendv_pd (point2[0], temp_point2[0], full_mask);      
        point[1] = _mm256_blendv_pd (point[1], temp_point[1], full_mask);   
        point2[1] = _mm256_blendv_pd (point2[1], temp_point2[1], full_mask); 
    }
    
    _mm256_storeu_pd ((double *)P->x, point[0]);  
    _mm256_storeu_pd ((double *)P->x+4, point2[0]);   
    _mm256_storeu_pd ((double *)P->y, point[1]);   
    _mm256_storeu_pd ((double *)P->y+4, point2[1]); 
    fpneg512(PCurve->prime, P->y);                                    // point[1]: y coordinate  
    temp_point[1] = _mm256_loadu_pd ((double const *)P->y);           // temp_point[1]: -y coordinate
    temp_point2[1] = _mm256_loadu_pd ((double const *)P->y+4);     
    full_mask = _mm256_set1_pd ((double)sign);  
    point[1] = _mm256_blendv_pd (point[1], temp_point[1], full_mask); // if mask = 0xFF...F then choose negative of the point
    point2[1] = _mm256_blendv_pd (point2[1], temp_point2[1], full_mask);
    _mm256_storeu_pd ((double *)P->y, point[1]); 
    _mm256_storeu_pd ((double *)P->y+4, point2[1]);

    return;
}


void lut_extproj_Ted512(point_extproj_precomp_Ted512* table, point_extproj_precomp_Ted512 P, int digit, unsigned int npoints, PCurveStruct PCurve)
{ // Constant-time table lookup to extract an extended twisted Edwards point (X+Y:Y-X:2Z:2T) from the precomputed table
  // Twisted Edwards a=-1 curve over p = 2^512-319
  // Operation: P = sign * table[(|digit|-1)/2], where sign=1 if digit>0 and sign=-1 if digit<0
    unsigned int i;
    int sign=0, pos, mask;
    __m256d point[8], temp_point[8], point2[3], temp_point2[3], full_mask;
    dig512 t1;

    sign = ((unsigned int)digit >> (sizeof(unsigned int)*8-1)) - 1; // if digit<0 then sign = 0x00...0 else sign = 0xFF...F
    pos = ((sign & (digit ^ -digit)) ^ -digit) >> 1;                // position = (|digit|-1)/2  
    point[0] = _mm256_loadu_pd ((double const *) table[0]->XY);     // point = table[0]  
    point[1] = _mm256_loadu_pd ((double const *) table[0]->XY+4);  
    point[2] = _mm256_loadu_pd ((double const *) table[0]->XY+2*4);  
    point[3] = _mm256_loadu_pd ((double const *) table[0]->XY+3*4);  
    point[4] = _mm256_loadu_pd ((double const *) table[0]->XY+4*4);  
    point[5] = _mm256_loadu_pd ((double const *) table[0]->XY+5*4);       
    point[6] = _mm256_loadu_pd ((double const *) table[0]->XY+6*4);     
    point[7] = _mm256_loadu_pd ((double const *) table[0]->XY+7*4);    // 512*4 coord = 2048 and 2048/256 = 8  

    for (i=1; i<npoints; i++) 
    { 
        pos--;
        // If match then mask = 0xFF...F else sign = 0x00...0
        mask = (((unsigned int)-pos >> (sizeof(unsigned int)*8-1)) || ((unsigned int)pos >> (sizeof(unsigned int)*8-1))) - 1; 
        full_mask = _mm256_set1_pd ((double)mask);
        temp_point[0] = _mm256_loadu_pd ((double const *) table[i]->XY);    // temp_point = table[i+1]
        temp_point[1] = _mm256_loadu_pd ((double const *) table[i]->XY+4);
        temp_point[2] = _mm256_loadu_pd ((double const *) table[i]->XY+2*4);
        temp_point[3] = _mm256_loadu_pd ((double const *) table[i]->XY+3*4);
        temp_point[4] = _mm256_loadu_pd ((double const *) table[i]->XY+4*4);
        temp_point[5] = _mm256_loadu_pd ((double const *) table[i]->XY+5*4);
        temp_point[6] = _mm256_loadu_pd ((double const *) table[i]->XY+6*4);
        temp_point[7] = _mm256_loadu_pd ((double const *) table[i]->XY+7*4);
        // If mask = 0x00...0 then point = point, else if mask = 0xFF...F then point = temp_point
        point[0] = _mm256_blendv_pd (point[0], temp_point[0], full_mask);     
        point[1] = _mm256_blendv_pd (point[1], temp_point[1], full_mask);  
        point[2] = _mm256_blendv_pd (point[2], temp_point[2], full_mask);  
        point[3] = _mm256_blendv_pd (point[3], temp_point[3], full_mask);   
        point[4] = _mm256_blendv_pd (point[4], temp_point[4], full_mask);  
        point[5] = _mm256_blendv_pd (point[5], temp_point[5], full_mask);  
        point[6] = _mm256_blendv_pd (point[6], temp_point[6], full_mask);  
        point[7] = _mm256_blendv_pd (point[7], temp_point[7], full_mask); 
    }
    
    _mm256_storeu_pd ((double *)P->XY, point[0]);    
    _mm256_storeu_pd ((double *)P->XY+4, point[1]); 
    _mm256_storeu_pd ((double *)P->XY+2*4, point[2]);    
    _mm256_storeu_pd ((double *)P->XY+3*4, point[3]);      
    _mm256_storeu_pd ((double *)P->XY+4*4, point[4]);    
    _mm256_storeu_pd ((double *)P->XY+5*4, point[5]);      
    _mm256_storeu_pd ((double *)P->XY+6*4, point[6]);     
    _mm256_storeu_pd ((double *)P->XY+7*4, point[7]); 
    point2[0] = _mm256_loadu_pd ((double const *)P->XY+4); 
    point[1] = _mm256_loadu_pd ((double const *)P->YX);   
    point2[1] = _mm256_loadu_pd ((double const *)P->YX+4);  
    point[2] = _mm256_loadu_pd ((double const *)P->T2); 
    point2[2] = _mm256_loadu_pd ((double const *)P->T2+4);  
    fpcopy512(P->XY, t1);
    fpcopy512(P->YX, P->XY);
    fpcopy512(t1, P->YX);
    fpneg512(PCurve->prime, P->T2);                                     // point: x, t coordinate  
    temp_point[0] = _mm256_loadu_pd ((double const *)P->XY);            // temp_point: -x, -t coordinate
    temp_point2[0] = _mm256_loadu_pd ((double const *)P->XY+4);       
    temp_point[1] = _mm256_loadu_pd ((double const *)P->YX);        
    temp_point2[1] = _mm256_loadu_pd ((double const *)P->YX+4);       
    temp_point[2] = _mm256_loadu_pd ((double const *)P->T2);  
    temp_point2[2] = _mm256_loadu_pd ((double const *)P->T2+4);       
    full_mask = _mm256_set1_pd ((double)sign); 
    point[0] = _mm256_blendv_pd (temp_point[0], point[0], full_mask);   // if mask = 0x00...0 then choose negative of the point
    point2[0] = _mm256_blendv_pd (temp_point2[0], point2[0], full_mask); 
    point[1] = _mm256_blendv_pd (temp_point[1], point[1], full_mask);
    point2[1] = _mm256_blendv_pd (temp_point2[1], point2[1], full_mask); 
    point[2] = _mm256_blendv_pd (temp_point[2], point[2], full_mask);
    point2[2] = _mm256_blendv_pd (temp_point2[2], point2[2], full_mask); 
    _mm256_storeu_pd ((double *)P->XY, point[0]); 
    _mm256_storeu_pd ((double *)P->XY+4, point2[0]); 
    _mm256_storeu_pd ((double *)P->YX, point[1]); 
    _mm256_storeu_pd ((double *)P->YX+4, point2[1]);
    _mm256_storeu_pd ((double *)P->T2, point[2]); 
    _mm256_storeu_pd ((double *)P->T2+4, point2[2]);

    return;
}


void lut_extaff_Ted512(point_extaff_precomp_Ted512* table, point_extaff_precomp_Ted512 P, int digit, int sign, unsigned int npoints, PCurveStruct PCurve)
{ // Constant-time table lookup to extract an extended affine point (x+y,y-x,2t) from the precomputed table
  // Twisted Edwards a=-1 curve over p = 2^512-319
  // Operation: if sign=0 then P = table[digit], else if (sign=-1) then P = -table[digit]
    unsigned int i;
    int pos, mask;
    __m256d point[3], temp_point[3], point2[3], temp_point2[3], full_mask;
            
    pos = digit;                                                       // Load digit position.  
    point[0] = _mm256_loadu_pd ((double const *) table[0]->xy);        // point = table[0] 
    point2[0] = _mm256_loadu_pd ((double const *) table[0]->xy+4);          
    point[1] = _mm256_loadu_pd ((double const *) table[0]->yx);  
    point2[1] = _mm256_loadu_pd ((double const *) table[0]->yx+4);          
    point[2] = _mm256_loadu_pd ((double const *) table[0]->t2);  
    point2[2] = _mm256_loadu_pd ((double const *) table[0]->t2+4);            

    for (i=1; i<npoints; i++) 
    { 
        pos--;
        // If match then mask = 0xFF...F else sign = 0x00...0
        mask = (((unsigned int)-pos >> (sizeof(unsigned int)*8-1)) || ((unsigned int)pos >> (sizeof(unsigned int)*8-1))) - 1; 
        full_mask = _mm256_set1_pd ((double)mask);
        temp_point[0] = _mm256_loadu_pd ((double const *) table[i]->xy);    // temp_point = table[i+1]
        temp_point2[0] = _mm256_loadu_pd ((double const *) table[i]->xy+4); 
        temp_point[1] = _mm256_loadu_pd ((double const *) table[i]->yx);
        temp_point2[1] = _mm256_loadu_pd ((double const *) table[i]->yx+4);
        temp_point[2] = _mm256_loadu_pd ((double const *) table[i]->t2);
        temp_point2[2] = _mm256_loadu_pd ((double const *) table[i]->t2+4);
        // If mask = 0x00...0 then point = point, else if mask = 0xFF...F then point = temp_point
        point[0] = _mm256_blendv_pd (point[0], temp_point[0], full_mask);    
        point2[0] = _mm256_blendv_pd (point2[0], temp_point2[0], full_mask);      
        point[1] = _mm256_blendv_pd (point[1], temp_point[1], full_mask);   
        point2[1] = _mm256_blendv_pd (point2[1], temp_point2[1], full_mask);     
        point[2] = _mm256_blendv_pd (point[2], temp_point[2], full_mask);   
        point2[2] = _mm256_blendv_pd (point2[2], temp_point2[2], full_mask); 
    }
    
    _mm256_storeu_pd ((double *)P->xy, point[1]);  
    _mm256_storeu_pd ((double *)P->xy+4, point2[1]);   
    _mm256_storeu_pd ((double *)P->yx, point[0]);   
    _mm256_storeu_pd ((double *)P->yx+4, point2[0]);   
    _mm256_storeu_pd ((double *)P->t2, point[2]);   
    _mm256_storeu_pd ((double *)P->t2+4, point2[2]); 
    fpneg512(PCurve->prime, P->t2);                                    // point negated: -x, -t coordinate 
    temp_point[0] = _mm256_loadu_pd ((double const *)P->xy);           // temp_point is point negated
    temp_point2[0] = _mm256_loadu_pd ((double const *)P->xy+4);     
    temp_point[1] = _mm256_loadu_pd ((double const *)P->yx);           
    temp_point2[1] = _mm256_loadu_pd ((double const *)P->yx+4);    
    temp_point[2] = _mm256_loadu_pd ((double const *)P->t2);           
    temp_point2[2] = _mm256_loadu_pd ((double const *)P->t2+4);    
    full_mask = _mm256_set1_pd ((double)sign); 
    point[0] = _mm256_blendv_pd (point[0], temp_point[0], full_mask);  // if mask = 0xFF...F then choose negative of the point
    point2[0] = _mm256_blendv_pd (point2[0], temp_point2[0], full_mask);
    point[1] = _mm256_blendv_pd (point[1], temp_point[1], full_mask); 
    point2[1] = _mm256_blendv_pd (point2[1], temp_point2[1], full_mask);
    point[2] = _mm256_blendv_pd (point[2], temp_point[2], full_mask); 
    point2[2] = _mm256_blendv_pd (point2[2], temp_point2[2], full_mask);
    _mm256_storeu_pd ((double *)P->xy, point[0]); 
    _mm256_storeu_pd ((double *)P->xy+4, point2[0]); 
    _mm256_storeu_pd ((double *)P->yx, point[1]); 
    _mm256_storeu_pd ((double *)P->yx+4, point2[1]); 
    _mm256_storeu_pd ((double *)P->t2, point[2]); 
    _mm256_storeu_pd ((double *)P->t2+4, point2[2]); 

    return;
}
#endif



/******** Recoding function for variable-base scalar multiplication ********/

void fixed_window_recode(dig *scalar, unsigned int nbit, unsigned int w, int *digits)
{ //  Computes the fixed window representation of scalar, where nonzero digits are in set {+-1,+-3,...,+-(2^(w-1)-1)}
    unsigned int val1, val2, t, cwords, i, j;
    int k;
    sdig temp, temp2;
    dig res, borrow;

    cwords = (nbit+ML_WORD-1)/ML_WORD;    // Number of computer words to represent scalar
    t = (nbit+(w-2))/(w-1);               // Fixed length of the fixed window representation
    val1 = (1 << w) - 1;
    val2 = (1 << (w-1));

    for (i = 0; i <= (t-1); i++)
    {
        temp = (scalar[0] & val1) - val2; // ki = (k mod 2^w)/2^(w-1)
        *digits = (int)temp;
        digits++;
                 
        res = scalar[0] - temp;           // k = (k-ki)/2^(w-1)
        borrow = (temp > 0) & (scalar[0] < (unsigned int)temp);
        scalar[0] = res;
  
        for (j = 1; j < cwords; j++)
        {
            res = scalar[j];
            scalar[j] = res - borrow;
            borrow = (res < borrow); 
        }    
  
        temp = (scalar[cwords-1] & (val2-1)) << (ML_WORD-(w-1));
        scalar[cwords-1] = scalar[cwords-1] >> (w-1);

        for (k = cwords-2; k >= 0; k--) {
            temp2 = (scalar[k] & (val2-1)) << (ML_WORD-(w-1));
            scalar[k] = (scalar[k] >> (w-1)) | temp;
            temp = temp2;
        }
    } 
    *digits = (int)scalar[0];             // kt = k  (t+1 digits)

    return;
}


/******** Recoding function for fixed-base scalar multiplication ********/

void mLSB_set_recode(dig *scalar, unsigned int nbit, unsigned int l, unsigned int d, int *digits)
{ //  Computes the modified LSB-set representation of scalar
    unsigned int cwords, i, j;
    int k;
    dig temp, temp2, carry;
    
    cwords = (nbit+ML_WORD-1)/ML_WORD;                    // Number of computer words to represent scalar
    digits[d-1] = 0;

    // Shift scalar to the right by 1
    temp = (scalar[cwords-1] & 1) << (ML_WORD-1);
    scalar[cwords-1] = scalar[cwords-1] >> 1;

    for (k = cwords-2; k >= 0; k--) {
        temp2 = (scalar[k] & 1) << (ML_WORD-1);
        scalar[k] = (scalar[k] >> 1) | temp;
        temp = temp2;
    }

    for (i = 0; i <= (d-2); i++)
    {
        digits[i] = (int)((scalar[0] & 1) - 1);           // Convention for the "sign" row: 
                                                          // if k_(i+1) = 0 then digit_i = -1 (negative), else if k_(i+1) = 1 then digit_i = 0 (positive)
        // Shift scalar to the right by 1
        temp = (scalar[cwords-1] & 1) << (ML_WORD-1);
        scalar[cwords-1] = scalar[cwords-1] >> 1;

        for (k = cwords-2; k >= 0; k--) {
            temp2 = (scalar[k] & 1) << (ML_WORD-1);
            scalar[k] = (scalar[k] >> 1) | temp;
            temp = temp2;
        }  
    } 

    for (i = d; i <= (l-1); i++)
    {
        digits[i] = scalar[0] & 1;                        // digits_i = k mod 2. Sign is determined by the "sign" row

        // Shift scalar to the right by 1
        temp = (scalar[cwords-1] & 1) << (ML_WORD-1);
        scalar[cwords-1] = scalar[cwords-1] >> 1;

        for (k = cwords-2; k >= 0; k--) {
            temp2 = (scalar[k] & 1) << (ML_WORD-1);
            scalar[k] = (scalar[k] >> 1) | temp;
            temp = temp2;
        }
                        
        temp = -digits[i-(i/d)*d] & digits[i];            // if (digits_i=0 \/ 1) then temp = 0, else if (digits_i=-1) then temp = 1 
            
        // floor(scalar/2) + temp
        scalar[0] = scalar[0] + temp;
        carry = (temp & (~(0-scalar[0]) >> (ML_WORD-1)) & (~scalar[0] >> (ML_WORD-1)));      //carry = (scalar[0] < temp);
        for (j = 1; j < cwords; j++)
        {
            scalar[j] = scalar[j] + carry; 
            carry = (temp & (~(0-scalar[j]) >> (ML_WORD-1)) & (~scalar[j] >> (ML_WORD-1)));  //carry = (scalar[j] < temp);
        }
    } 
    return;              
}


/******** Non-constant time recoding function for double-scalar multiplication ********/

void wNAF_recode(dig *scalar, unsigned int nbits, unsigned int w, int *digits)
{ //  Computes wNAF of scalar, where digits are in set {0,+-1,+-3,...,+-(2^(w-1)-1)}
    unsigned int j, cwords, mask;
    int digit, index=0, k, val1, val2;
    dig temp=0, temp2, carry;
    
    for (j = 0; j<=nbits; j++) digits[j] = 0;              // Initialize digit output to zero
    cwords = (nbits+ML_WORD-1)/ML_WORD;                    // Number of computer words to represent scalar
    val1 = (int)(1 << (w-1)) - 1;                          // 2^(w-1) - 1
    val2 = (int)(1 << w);                                  // 2^w
    mask = (unsigned int)val2 - 1;                         // 2^w - 1
    
    for (j = 0; j < cwords; j++) {
        temp = temp | scalar[j];
    }

    while (temp != 0)
    {
        digit = (scalar[0] & 1); 

        if (digit==0) 
        {                                                
            // Shift scalar to the right by 1
            temp = (scalar[cwords-1] & 1) << (ML_WORD-1);
            scalar[cwords-1] = scalar[cwords-1] >> 1;

            for (k = cwords-2; k >= 0; k--) {
                temp2 = (scalar[k] & 1) << (ML_WORD-1);
                scalar[k] = (scalar[k] >> 1) | temp;
                temp = temp2;
            }  
            digits[index] = 0;
        }
        else
        {
            digit = (scalar[0] & mask); 

            // Shift scalar to the right by w
            temp = (scalar[cwords-1] & mask) << (ML_WORD-w);
            scalar[cwords-1] = scalar[cwords-1] >> w;

            for (k = cwords-2; k >= 0; k--) {
                temp2 = (scalar[k] & mask) << (ML_WORD-w);
                scalar[k] = (scalar[k] >> w) | temp;
                temp = temp2;
            }

            if (digit > val1) {
                digit = digit - val2; 
            }
            if (digit < 0) 
            {
                // scalar + 1
                scalar[0] = scalar[0] + 1;
                carry = (scalar[0] < 1);
                for (j = 1; j < cwords; j++)
                {
                    scalar[j] = scalar[j] + carry;                
                    carry = (scalar[j] < carry);
                }
            }
            digits[index] = digit; 
            
            temp = 0;
            for (j = 0; j < cwords; j++) {
                temp = temp | scalar[j];
            }
            if (temp != 0)              // Check if scalar != 0
            {
                for (j = 0; j < (w-1); j++) 
                {     
                    index++; 
                    digits[index] = 0;
                }
            }
        }
        
        // To check if scalar != 0 
        temp = 0;
        for (j = 0; j < cwords; j++) {
            temp = temp | scalar[j];
        }
        index++;
    } 
    return;
}