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
* Abstract: template for elliptic curve and scalar arithmetic functions
*
* This software is based on the article by Joppe Bos, Craig Costello, 
* Patrick Longa and Michael Naehrig, "Selecting elliptic curves for
* cryptography: an efficiency and security analysis", preprint available
* at http://eprint.iacr.org/2014/130.
***************************************************************************/  

#include <malloc.h>
#include <tchar.h>
#include "msr_ecclib.h"
#include "immintrin.h"


/***************************************************************************************/
/********** POINT/SCALAR FUNCTIONS FOR WEIERSTRASS a=-3 CURVES ("JAC" CURVES) **********/

void ECCSET_W(POINT_WAFF P, PCurveStruct JacCurve)
{ // Set generator P = (x,y) on Weierstrass a=-3 curve
    dig i, num_words = (JacCurve->pbits+ML_WORD-1)/ML_WORD; 

    for (i = 0; i < num_words; i++) 
    {
        P->x[i] = JacCurve->generator_x[i]; 
        P->y[i] = JacCurve->generator_y[i]; 
    }
    return;
}


BOOL ECC_IS_INFINITY_WAFF(POINT_WAFF P, PCurveStruct JacCurve)
{ // Check if point P is the point at infinity (0,0) 
    dig i, c, num_words = (JacCurve->pbits+ML_WORD-1)/ML_WORD;
    BOOL answer = TRUE;

    c = P->x[0] | P->y[0];
    for (i = 1; i < num_words; i++)
    {
        c = c | P->x[i] | P->y[i]; 
    }
    answer = (BOOL)((~(0 - c) >> (ML_WORD-1)) & (~c >> (ML_WORD-1)));

    return answer;
}


BOOL ECC_IS_INFINITY_WJAC(POINT_WJAC P, PCurveStruct JacCurve)
{ // Check if Jacobian point P is the point at infinity (0:Y:0) 
    dig i, c, num_words = (JacCurve->pbits+ML_WORD-1)/ML_WORD;
    BOOL answer = TRUE;

    c = P->X[0] | P->Z[0];
    for (i = 1; i < num_words; i++)
    {
        c = c | P->X[i] | P->Z[i]; 
    }
    answer = (BOOL)((~(0 - c) >> (ML_WORD-1)) & (~c >> (ML_WORD-1)));

    return answer;
}


BOOL ECCNORM_W(POINT_WJAC Q, POINT_WAFF P, PCurveStruct JacCurve)
{ // Normalize a Jacobian point Q = (X:Y:Z) -> P = (x,y)
    FELM t1, t2, t3;  
#ifdef ML_COUNT
    ninv++; nmul+=3; nsqr++;
#endif
    
    // Check if Q is the point at infinity (0:Y:0)
    // SECURITY NOTE: this if-statement evaluates over public information when the function is called from constant-time scalar multiplications, i.e.,
    //                Q is never the point at infinity when the call is from ecc_scalar_mul_Jacxxx() or ecc_scalar_mul_fixed_internal_Jacxxx(). For more information,
    //                refer to "Selecting elliptic curves for cryptography: an efficiency and security analysis", http://eprint.iacr.org/2014/130
    if (ECC_IS_INFINITY_WJAC(Q, JacCurve) == TRUE) {
        FP_ZERO(P->x); FP_ZERO(P->y);    // Output the point at infinity P = (0,0)   
        return TRUE;
    }
    
    FP_COPY(Q->Z, t1);  
    FP_INV(t1);                      // t1 = Z^-1
    FP_SQR(t1, t2);                  // t2 = Z^-2
    FP_MUL(Q->X, t2, t3);            // t3 = X/Z^2
    FP_COPY(t3, P->x);               // x = X/Z^2
    FP_MUL(t1, t2, t3);              // t3 = Z^-3
    FP_MUL(Q->Y, t3, t1);            // t1 = Y/Z^3 
    FP_COPY(t1, P->y);               // y = Y/Z^3

// cleanup
    FP_ZERO(t1);
    FP_ZERO(t2);
    FP_ZERO(t3);

    return TRUE;
}


BOOL ECCDOUBLE_WJAC(POINT_WJAC P, PCurveStruct JacCurve)      
{ // Point doubling P = 2P
  // Weierstrass a=-3 curve
  // Input:  P = (X,Y,Z) in Jacobian coordinates
  // Output: 2P = (X,Y,Z) in Jacobian coordinates
    FELM t1, t2, t3, t4;  
    BOOL OK = TRUE;
#ifdef ML_COUNT
    nmul+=4; nsqr+=4; nadd+=8;
#endif  
    UNREFERENCED_PARAMETER(JacCurve);
     
    // SECURITY NOTE: this function does not produce exceptions on prime-order Weierstrass curves. For more information, refer to
    //                "Selecting elliptic curves for cryptography: an efficiency and security analysis", http://eprint.iacr.org/2014/130

    FP_SQR(P->Z, t1);           // t1 = z^2
    FP_MUL(P->Z, P->Y, t4);     // t4 = zy
    FP_ADD(P->X, t1, t2);       // t2 = x + z^2
    FP_SUB(P->X, t1, t1);       // t1 = x - z^2
    FP_COPY(t4, P->Z);          // Zfinal = zy
    FP_MUL(t1, t2, t3);         // t3 = (x + z^2)(x - z^2)
    FP_DIV2(t3, t2);            // t2 = (x + u.z^2)(x - u.z^2)/2
    FP_ADD(t3, t2, t1);         // t1 = alpha = 3(x + u.z^2)(x - u.z^2)/2               
    FP_SQR(P->Y, t2);           // t2 = y^2
    FP_SQR(t1, t4);             // t4 = alpha^2
    FP_MUL(P->X, t2, t3);       // t3 = beta = xy^2
    FP_SUB(t4, t3, t4);         // t4 = alpha^2-beta
    FP_SUB(t4, t3, P->X);       // Xfinal = alpha^2-2beta
    FP_SUB(t3, P->X, t4);       // t4 = beta-Xfinal
    FP_SQR(t2, t3);             // t3 = y^4
    FP_MUL(t1, t4, t2);         // t2 = alpha.(beta-Xfinal)
    FP_SUB(t2, t3, P->Y);       // Yfinal = alpha.(beta-Xfinal)-y^4
    
// cleanup
    FP_ZERO(t1);
    FP_ZERO(t2);
    FP_ZERO(t3);
    FP_ZERO(t4);

    return OK;
}


void ECCDOUBLEADD_WJAC(POINT_PRECOMP_WCHU Q, POINT_WJAC P, PCurveStruct JacCurve)      
{ // Point addition P = 2P+Q
  // Weierstrass a=-3 curve
  // Inputs: P = (X1,Y1,Z1) in Jacobian coordinates
  //         Q = (X2,Y2,Z2,Z2^2,Z2^3) in Chudnovsky coordinates 
  // Output: P = (X1,Y1,Z1) in Jacobian coordinates
    FELM t1, t2, t3, t4, t5, t6, t7; 
#ifdef ML_COUNT
    nmul+=16; nsqr+=5; nadd+=13;
#endif
    UNREFERENCED_PARAMETER(JacCurve);
     
    // SECURITY NOTE: this function does not produce exceptions when P!=inf, Q!=inf, P!=Q, P!=-Q or Q!=-2P. In particular, it works when called from scalar multiplications 
    //                ecc_scalar_mul_Jacxxx() and ecc_scalar_mul_fixed_Jacxxx(). For more information, refer to "Selecting elliptic curves for cryptography: an efficiency 
    //                and security analysis", http://eprint.iacr.org/2014/130
    
    FP_SQR(P->Z, t2);               // t2 = z1^2
    FP_MUL(Q->Z3, P->Y, t3);        // t3 = z2^3*y1
    FP_MUL(P->Z, t2, t4);           // t4 = z1^3
    FP_MUL(t2, Q->X, t1);           // t1 = z1^2*x2                      
    FP_MUL(Q->Y, t4, t2);           // t2 = z1^3*y2
    FP_MUL(Q->Z2, P->X, t6);        // t6 = z2^2*x1
    FP_SUB(t2, t3, t2);             // t2 = alpha = z1^3*y2-z2^3*y1
    FP_SUB(t1, t6, t1);             // t1 = beta = z1^2*x2-z2^2*x1
    FP_SQR(t2, t4);                 // t4 = alpha^2
    FP_SQR(t1, t5);                 // t5 = beta^2
    FP_MUL(P->Z, Q->Z, t7);         // t5 = z1*z2
    FP_MUL(t6, t5, P->X);           // x1 = x1' = z2^2*x1*beta^2
    FP_MUL(t1, t5, t6);             // t6 = beta^3
    FP_SUB(t4, t6, t4);             // t4 = alpha^2 - beta^3
    FP_SUB(t4, P->X, t4);           // t4 = alpha^2 - beta^3 - x1'
    FP_SUB(t4, P->X, t4);           // t4 = alpha^2 - beta^3 - 2*x1'
    FP_SUB(t4, P->X, t4);           // t4 = omega = alpha^2 - beta^3 - 3*x1'
    FP_MUL(t6, t3, P->Y);           // y1 = y1' = z2^3*y1*beta^3 
    FP_MUL(t1, t7, t3);             // t3 = z1' = z1*z2*beta
    FP_MUL(t2, t4, t1);             // t1 = alpha.omega
    FP_SQR(t4, t2);                 // t2 = omega^2                      
    FP_ADD(t1, P->Y, t1);           // t1 = alpha.omega + y1'
    FP_ADD(t1, P->Y, t1);           // t1 = theta = alpha.omega + 2y1'    
    FP_MUL(t3, t4, P->Z);           // Zfinal = z1'*omega
    FP_MUL(t2, t4, t5);             // t5 = omega^3
    FP_MUL(t2, P->X, t4);           // t4 = x1'*omega^2
    FP_SQR(t1, t3);                 // t3 = theta^2
    FP_SUB(t3, t5, t3);             // t3 = theta^2 - omega^3
    FP_SUB(t3, t4, t3);             // t3 = theta^2 - omega^3 - x1'*omega^2
    FP_SUB(t3, t4, P->X);           // Xfinal = theta^2 - omega^3 - 2*x1'*omega^2
    FP_SUB(P->X, t4, t3);           // t3 = Xfinal-x1'*omega^2
    FP_MUL(P->Y, t5, t2);           // t2 = y1'*omega^3
    FP_MUL(t3, t1, t5);             // t5 = theta.(Xfinal-x1'*omega^2)
    FP_SUB(t5, t2, P->Y);           // Yfinal = theta.(Xfinal-x1'*omega^2) - y1'*omega^3
    
// cleanup
    FP_ZERO(t1);
    FP_ZERO(t2);
    FP_ZERO(t3);
    FP_ZERO(t4);
    FP_ZERO(t5);
    FP_ZERO(t6);
    FP_ZERO(t7);

    return;
}


// Functions for the precomputation (Weierstrass a=-3 curves)

static void ECCADD_PRECOMP_WJAC(POINT_WJAC P, POINT_PRECOMP_WCHU Q, POINT_PRECOMP_WCHU R)      
{ // Special point addition R = P+Q with identical Z-coordinate for the precomputation
  // Weierstrass a=-3 curve
  // Inputs:  P = (X1,Y1,Z) in Jacobian coordinates with the same Z-coordinate
  //          Q = (X2,Y2,Z,Z^2,Z^3) in Chudnovsky coordinates with the same Z-coordinate
  //          Values (X1',Y1')
  // Outputs: R = (X3,Y3,Z3,Z3^2,Z3^2) in Chudnovsky coordinates 
  //          new representation P  = (X1',Y1',Z1') = (X1.(X2-X1)^2, X1.(X2-X1)^3, Z.(X2-X1)) in Jacobian coordinates
    FELM t1, t2, t3, t4; 
#ifdef ML_COUNT
    nmul+=6; nsqr+=3; nadd+=7;
#endif
     
    // SECURITY NOTE: this function does not produce exceptions in the context of variable-base precomputation. For more information,
    //                refer to "Selecting elliptic curves for cryptography: an efficiency and security analysis", http://eprint.iacr.org/2014/130
    
    FP_SUB(Q->X, P->X, t1);         // t1 = x2-x1
    FP_MUL(P->Z, t1, R->Z);         // Zfinal = z.(x2-x1)
    FP_COPY(R->Z, P->Z);            // Z1' = z.(x2-x1)
    FP_SQR(t1, t2);                 // t2 = (x2-x1)^2
    FP_SQR(R->Z, R->Z2);            // Z2final = Zfinal^2
    FP_MUL(t1, t2, t3);             // t3 = (x2-x1)^3
    FP_MUL(P->X, t2, t4);           // t4 = X1' = x1.(x2-x1)^2
    FP_COPY(t4, P->X);              // X1'
    FP_SUB(Q->Y, P->Y, t1);         // t1 = y2-y1
    FP_SQR(t1, R->X);               // X3 = (y2-y1)^2
    FP_MUL(R->Z, R->Z2, R->Z3);     // Z3final = Zfinal^3
    FP_SUB(R->X, t3, R->X);         // X3 = (y2-y1)^2 - (x2-x1)^3
    FP_SUB(R->X, t4, R->X);         // X3 = (y2-y1)^2 - (x2-x1)^3 - x1.(x2-x1)^2
    FP_SUB(R->X, t4, R->X);         // X3final = (y2-y1)^2 - (x2-x1)^3 - 2*x1.(x2-x1)^2
    FP_SUB(t4, R->X, t2);           // t2 = x1.(x2-x1)^2-X3
    FP_MUL(t1, t2, t4);             // t4 = (y2-y1)[x1.(x2-x1)^2-X3]
    FP_MUL(P->Y, t3, t2);           // t2 = Y1' = y1*(x2-x1)^3
    FP_COPY(t2, P->Y);              // Y1'
    FP_SUB(t4, t2, R->Y);           // Yfinal = (y2-y1)[x1.(x2-x1)^2-X3] - y1*(x2-x1)^3
    
// cleanup
    FP_ZERO(t1);
    FP_ZERO(t2);
    FP_ZERO(t3);
    FP_ZERO(t4);

    return;
}


void ECC_PRECOMP_WJAC(POINT_WAFF P, POINT_PRECOMP_WCHU *T, unsigned int npoints, PCurveStruct JacCurve)
{ // Precomputation scheme using Jacobian coordinates
  // Weierstrass a=-3 curve
  // Input:   P = (x,y)
  // Outputs: T[0] = P, T[1] = 3*P, ... , T[npoints-1] = (2*npoints-1)*P in coordinates (X:Y:Z:Z^2:Z^3)
    POINT_WJAC P2;
    FELM t1, t2, t3; 
    unsigned int i;
#ifdef ML_COUNT
    nmul+=3; nsqr+=4; nadd+=7;
#endif
    UNREFERENCED_PARAMETER(JacCurve);
    
    // SECURITY NOTE: this function does not produce exceptions in the context of variable-base scalar multiplication and double-scalar multiplication.
    //                For more information, refer to "Selecting elliptic curves for cryptography: an efficiency and security analysis", http://eprint.iacr.org/2014/130

    // Generating 2P = 2(x,y) = (X2,Y2,Z2) and P = (x,y) = (X1',Y1',Z1',Z1^2',Z1^3') = (x*y^2, y*y^3, y, y^2, y^3)
    FP_ZERO(t2); t2[0] = 1;             // t2 = 1
    FP_SQR(P->x, t1);                   // t1 = x^2
    FP_SUB(t1, t2, t1);                 // t1 = x^2-1
    FP_DIV2(t1, t2);                    // t2 = (x^2-1)/2
    FP_ADD(t1, t2, t1);                 // t1 = alpha = 3(x^2-1)/2
    FP_SQR(P->y, T[0]->Z2);             // Z1^2' = y^2
    FP_MUL(T[0]->Z2, P->x, T[0]->X);    // X1' = beta = xy^2
    FP_MUL(T[0]->Z2, P->y, T[0]->Z3);   // Z1^3' = y^3
    FP_SQR(t1, t2);                     // t2 = alpha^2
    FP_SUB(t2, T[0]->X, t2);            // t2 = alpha^2-beta
    FP_SUB(t2, T[0]->X, P2->X);         // X2final = alpha^2-2beta
    FP_COPY(P->y, P2->Z);               // Z2final = y
    FP_COPY(P->y, T[0]->Z);             // Z1' = y
    FP_SQR(T[0]->Z2, T[0]->Y);          // Y1' = y^4
    FP_SUB(T[0]->X, P2->X, t2);         // t2 = beta-Xfinal
    FP_MUL(t1, t2, t3);                 // t3 = alpha.(beta-Xfinal)
    FP_SUB(t3, T[0]->Y, P2->Y);         // Y2final = alpha.(beta-Xfinal)-y^4
    
    for (i=1; i<npoints; i++) {
        // T[i] = 2P'+T[i-1] = (2*i+1)P = (X_(2*i+1),Y_(2*i+1),Z_(2*i+1),Z_(2*i+1)^2,Z_(2*i+1)^3) and new 2P' s.t. Z(2P')=Z_(2*i+1)
        ECCADD_PRECOMP_WJAC(P2, T[i-1], T[i]);
    }
    
// cleanup
    ECCZERO_WJAC(P2);
    FP_ZERO(t1);
    FP_ZERO(t2);
    FP_ZERO(t3);

    return;
}


static void ECCUADD_NO_INIT_WJAC(POINT_WJAC Q, POINT_WJAC P, POINT_WJAC *table, PCurveStruct JacCurve)      
{ // Complete point addition: if P=-Q then P=0, else if P=0 then P=Q, else if P=Q then P=2P, else P=P+Q 
  // Constant-time extraction over 5-LUT: table[0] = inf, table[1] = Q, table[2] = 2P, table[3] = P+Q, table[4] = P. First two entries and last one are pre-loaded. 
  // Weierstrass a=-3 curve
  // Inputs: P = (X1,Y1,Z1) in Jacobian coordinates
  //         Q = (X2,Y2,Z2) in Jacobian coordinates 
  // Output: P = P+Q = (X1,Y1,Z1) + (X2,Y2,Z2) in Jacobian coordinates
    FELM t1, t2, t3, t4, t5, t6, t7, t8; 
    unsigned int index = 0;
    sdig mask = 0;
#ifdef ML_COUNT
    nmul+=12; nsqr+=4; nadd+=11; ncsel+=8; nclut++; 
#endif
    UNREFERENCED_PARAMETER(JacCurve);
        
    // SECURITY NOTE: this constant-time addition function is complete (i.e., it works for any possible inputs, including the cases P!=Q, P=Q, P=-Q and P=inf) on prime-order Weierstrass curves.
    //                For more information, refer to "Selecting elliptic curves for cryptography: an efficiency and security analysis", http://eprint.iacr.org/2014/130
    
    FP_SQR(P->Z, t2);                         // t2 = z1^2
    FP_MUL(P->Z, t2, t3);                     // t3 = z1^3
    FP_MUL(t2, Q->X, t1);                     // t1 = z1^2*x2
    FP_MUL(t3, Q->Y, t4);                     // t4 = z1^3*y2
    FP_SQR(Q->Z, t3);                         // t3 = z2^2
    FP_MUL(Q->Z, t3, t5);                     // t5 = z2^3
    FP_MUL(t3, P->X, t7);                     // t7 = z2^2*x1
    FP_MUL(t5, P->Y, t8);                     // t8 = z2^3*y1
    FP_SUB(t1, t7, t1);                       // t1 = beta2 = z1^2*x2-z2^2*x1               
    FP_SUB(t4, t8, t4);                       // t4 = alpha2 = z1^3*y2-z2^3*y1
    index = COMPLETE_EVAL(t1, P->Z, t4, &mask);    // if t1=0 (P=-Q) then index=0, if Z1=0 (P inf) then index=1, if t4=0 (P=Q) then index=2, else index=3
                                                   // if index=3 then mask = 0xff...ff, else mask = 0
    mask = ~(-FP_ISZERO(Q->Z)) ;                   // if Z2=0 (Q inf) then mask = 0, else mask = 0xff...ff
    index = (mask & (index ^ 4)) ^ 4;         // if mask = 0 then index=4, else if mask = 0xff...ff then keep previous index  
    FP_ADD(P->X, t2, t3);                     // t3 = x1+z1^2
    FP_SUB(P->X, t2, t6);                     // t6 = x1-z1^2
    COMPLETE_SELECT(P->Y, t1, t2, mask);      // If mask=0 (DBL) then t2=y1, else if mask=-1 (ADD) then t2=beta2 
    FP_SQR(t2, t5);                           // t5 = y1^2 (DBL) or beta2^2 (ADD)
    COMPLETE_SELECT(P->X, t7, t7, mask);      // If mask=0 (DBL) then t7=x1, else if mask=-1 (ADD) then t7=z2^2*x1 
    FP_MUL(t5, t7, t1);                       // t1 = x1y1^2 = beta1 (DBL) or z2^2*x1*beta2^2 (ADD)
    FP_MUL(P->Z, t2, table[2]->Z);            // Z2Pfinal = z1y1
    FP_MUL(Q->Z, table[2]->Z, table[3]->Z);   // ZPQfinal = z1*z2*beta2
    COMPLETE_SELECT(t3, t2, t3, mask);        // If mask=0 (DBL) then t3=x1+z1^2, else if mask=-1 (ADD) then t3=beta2 
    COMPLETE_SELECT(t6, t5, t6, mask);        // If mask=0 (DBL) then t6=x1-z1^2, else if mask=-1 (ADD) then t6=beta2^2
    FP_MUL(t3, t6, t2);                       // t2 = (x1+z1^2)(x1-z1^2) (DBL) or beta2^3 (ADD)
    FP_DIV2(t2, t3);                          // t3 = (x1+z1^2)(x1-z1^2)/2
    FP_ADD(t2, t3, t3);                       // t3 = alpha1 = 3(x1+z1^2)(x1-z1^2)/2
    COMPLETE_SELECT(t3, t4, t3, mask);        // If mask=0 (DBL) then t3=alpha1, else if mask=-1 (ADD) then t3=alpha2
    FP_SQR(t3, t4);                           // t4 = alpha1^2 (DBL) or alpha2^2 (ADD)
    FP_SUB(t4, t1, t4);                       // t4 = alpha1^2-beta1 (DBL) or alpha2^2-z2^2*x1*beta2^2
    FP_SUB(t4, t1, table[2]->X);              // X2Pfinal = alpha1^2-2beta1 (DBL) or alpha2^2-2z2^2*x1*beta2^2 (ADD)
    FP_SUB(table[2]->X, t2, table[3]->X);     // XPQfinal = alpha^2-beta2^3-2z2^2*x1*beta2^2
    COMPLETE_SELECT(table[2]->X, table[3]->X, t4, mask);   // If mask=0 (DBL) then t4=X2Pfinal, else if mask=-1 (ADD) then t4=XPQfinal
    FP_SUB(t1, t4, t1);                       // t1 = beta1-X2Pfinal (DBL) or (ADD) z2^2*x1*beta2^2-XPQfinal
    FP_MUL(t3, t1, t4);                       // t4 = alpha1.(beta1-X2Pfinal) or alpha2.(z2^2*x1*beta2^2-XPQfinal)
    COMPLETE_SELECT(t5, t8, t1, mask);        // If mask=0 (DBL) then t1=y1^2, else if mask=-1 (ADD) then t1=z2^3*y1
    COMPLETE_SELECT(t5, t2, t2, mask);        // If mask=0 (DBL) then t2=y1^2, else if mask=-1 (ADD) then t2=beta2^3
    FP_MUL(t1, t2, t3);                       // t3 = y1^4 (DBL) or z2^3*y1*beta2^3 (ADD)
    FP_SUB(t4, t3, table[2]->Y);              // Y2Pfinal = alpha1.(beta1-X2Pfinal)-y1^4 (DBL) or alpha2.(z2^2*x1*beta2^2-XPQfinal)-z2^3*y1*beta2^3 (ADD)
    FP_COPY(table[2]->Y, table[3]->Y);        // YPQfinal = alpha2.(z2^2*x1*beta2^2-XPQfinal)-z2^3*y1*beta2^3
    COMPLETE_LUT5(table, index, P);           // P = table[index]
    
// cleanup
    FP_ZERO(t1);
    FP_ZERO(t2);
    FP_ZERO(t3);
    FP_ZERO(t4);
    FP_ZERO(t5);
    FP_ZERO(t6);
    FP_ZERO(t7);
    FP_ZERO(t8);
        
    return;
}


void ECCUADD_WJAC(POINT_WJAC Q, POINT_WJAC P, PCurveStruct JacCurve)      
{ // Complete point addition: if P=-Q then P=0, else if P=0 then P=Q, else if P=Q then P=2P, else P=P+Q 
  // Constant-time extraction over 5-LUT: table[0] = inf, table[1] = Q, table[2] = 2P, table[3] = P+Q, table[4] = P. 
  // Weierstrass a=-3 curve
  // Inputs: P = (X1,Y1,Z1) in Jacobian coordinates
  //         Q = (X2,Y2,Z2) in Jacobian coordinates 
  // Output: P = P+Q = (X1,Y1,Z1) + (X2,Y2,Z2) in Jacobian coordinates
    POINT_WJAC table[5] = {0};
    dig j;
         
    table[0]->Y[0] = 1;                       // Initialize table[0] with the point at infinity (0:1:0)
    ECCCOPY_WJAC(Q, table[1]);                // Initialize table[1] with Q
    ECCCOPY_WJAC(P, table[4]);                // Initialize table[4] with P
    ECCUADD_NO_INIT_WJAC(Q, P, table, JacCurve);
    
// cleanup
    for (j = 0; j < 5; j++) {
        ECCZERO_WJAC(table[j]);
    }
    
    return;
}


BOOL ECC_MUL_W(POINT_WAFF P, dig *k, POINT_WAFF Q, PCurveStruct JacCurve)
{ // Variable-base scalar multiplication Q = k.P using fixed-window method 
  // Weierstrass a=-3 curve
    unsigned int j, t = (JacCurve->rbits+(W_VARBASE-2))/(W_VARBASE-1);           // Fixed length of the fixed window representation  
    unsigned int npoints = 1 << (W_VARBASE-2);
    unsigned int num_words = (JacCurve->nbits+ML_WORD-1)/(ML_WORD);              // Number of words to represent field elements and elements in Z_r
    dig k_temp[((sizeof(FELM)*8)+ML_WORD-1)/ML_WORD] = {0};
    int digits[((sizeof(FELM)*8)+W_VARBASE-2)/(W_VARBASE-1) + 1]={0};
    sdig i, odd = 0;
    POINT_WJAC T, TT;
    POINT_PRECOMP_WCHU table[1 << (W_VARBASE-2)], R;
    FELM t1, t2, t3;
    BOOL OK = FALSE;
#ifdef ML_COUNT
    nmul++; nsqr+=2; nadd+=7; nlut+=(t+1);
#endif
        
    // SECURITY NOTE: the crypto sensitive part of this function is protected against timing attacks and runs in constant-time on prime-order Weierstrass curves. 
    //                Conditional if-statements evaluate public data only and the number of iterations for all loops is public. Refer to "Selecting elliptic curves 
    //                for cryptography: an efficiency and security analysis", http://eprint.iacr.org/2014/130, for the full proof demonstrating exception-less, 
    //                constant-time execution of scalar multiplication.
    // DISCLAIMER:    the protocol designer is responsible for guarantying that early termination produced after detecting errors during input validation
    //                (of scalar k or base point P) does not leak any secret information. 

    
    /***** Input validation: *****/
    // Check if P is the point at infinity (0,0)              
    if (ECC_IS_INFINITY_WAFF(P, JacCurve) == TRUE) {     
        return FALSE;   
    } 
    // Is scalar k in [1,r-1]?               
    if ((FP_ISZERO(k) == TRUE) || (MOD_EVAL(k, JacCurve->prime) == FALSE)) {
        return FALSE;
    }
    // Are (x,y) in [0,p-1]?
    if (MOD_EVAL(P->x, JacCurve->prime) == FALSE || MOD_EVAL(P->y, JacCurve->prime) == FALSE) {  
        return FALSE;
    }
    // Does P lie on the curve?
    FP_SQR(P->y, t1);           // y^2
    FP_SQR(P->x, t2);
    FP_MUL(P->x, t2, t3); 
    FP_ADD(t3, JacCurve->parameter2, t2);
    FP_ADD(P->x, P->x, t3); 
    FP_ADD(P->x, t3, t3); 
    FP_SUB(t2, t3, t2);         // x^3 - 3x + b
    FP_SUB(t1, t2, t1); 
    if (FP_ISZERO(t1) == FALSE) {
        OK = FALSE;
        goto cleanup;
    }
    /*****************************/

    ECC_PRECOMP_WJAC(P, table, npoints, JacCurve);                // Precomputation of points T[0],...,T[npoints-1]

    odd = -((sdig)k[0] & 1);     
    FP_SUB(JacCurve->order, k, k_temp);                           // Converting scalar to odd (r-k if even)
    for (j = 0; j < num_words; j++)                               // If (odd) then k = k_temp else k = k 
    {
        k_temp[j] = (odd & (k[j] ^ k_temp[j])) ^ k_temp[j];
    }
    
    fixed_window_recode(k_temp, JacCurve->rbits, W_VARBASE, digits); 

    LUT_WCHU(table, R, digits[t], npoints, JacCurve); 
    FP_COPY(R->X, T->X);                                          // Initialize T = (X_T:Y_T:Z_T) with a point from the precomputed table
    FP_COPY(R->Y, T->Y); 
    FP_COPY(R->Z, T->Z); 

    for (i = (t-1); i >= 1; i--)
    {
        for (j = 0; j < (W_VARBASE-2); j++)
        {
            ECCDOUBLE_WJAC(T, JacCurve);                         // Double (X_T:Y_T:Z_T) = 2(X_T:Y_T:Z_T) 
        }
        LUT_WCHU(table, R, digits[i], npoints, JacCurve);        // Load R = (X_R:Y_R:Z_R:Z_R^2:Z_R^3) with a point from the precomputed table
        ECCDOUBLEADD_WJAC(R, T, JacCurve);                       // Double-add (X_T:Y_T:Z_T) = 2(X_T:Y_T:Z_T) + (X_R:Y_R:Z_R:Z_R^2:Z_R^3) 
    }
    
    // Perform last iteration
    for (j = 0; j < (W_VARBASE-1); j++)
    {
        ECCDOUBLE_WJAC(T, JacCurve);                             // Double (X_T:Y_T:Z_T) = 2(X_T:Y_T:Z_T)
    }
    LUT_WCHU(table, R, digits[0], npoints, JacCurve);            // Load R = (X_R:Y_R:Z_R:Z_R^2:Z_R^3) with a point from the precomputed table
    FP_COPY(R->X, TT->X);                                        // TT = R = (X_R:Y_R:Z_R)
    FP_COPY(R->Y, TT->Y);
    FP_COPY(R->Z, TT->Z);
    ECCUADD_WJAC(TT, T, JacCurve);                               // Complete addition (X_T:Y_T:Z_T) = (X_T:Y_T:Z_T) + (X_R:Y_R:Z_R)

    FP_COPY(T->Y, k_temp);
    FP_NEG(JacCurve->prime, k_temp);                             // Correcting scalar (-Ty if even)
    for (j = 0; j < num_words; j++)                              // If (even) then Ty = -Ty 
    {
        T->Y[j] = (odd & (T->Y[j] ^ k_temp[j])) ^ k_temp[j];
    }
    OK = ECCNORM_W(T, Q, JacCurve);                              // Output Q = (x,y)
    
cleanup:
    for (j = 0; j < ((sizeof(FELM)*8)+ML_WORD-1)/ML_WORD; j++) {
        k_temp[j] = 0;
    }
    for (j = 0; j < ((sizeof(FELM)*8)+W_VARBASE-2)/(W_VARBASE-1) + 1; j++) {
        digits[j] = 0;
    }
    ECCZERO_WJAC(T);
    ECCZERO_WJAC(TT);
    ECCZERO_WCHU(R);
    for (j = 0; j < (1 << (W_VARBASE-2)); j++) {
        ECCZERO_WCHU(table[j]);
    }
    FP_ZERO(t1);
    FP_ZERO(t2);
    FP_ZERO(t3);

    return OK;
}


void ECCUMADD_WJAC(POINT_WAFF Q, POINT_WJAC P, POINT_WJAC *table, PCurveStruct JacCurve)      
{ // Complete mixed point addition: if P=-Q then P=0, else if P=0 then P=Q, else if P=Q then P=2P, else P=P+Q 
  // Constant-time extraction over 4-LUT: table[0] = inf, table[1] = Q, table[2] = 2P, table[3] = P+Q. First two entries are pre-loaded. 
  // Weierstrass a=-3 curve
  // Inputs: P = (X1,Y1,Z1) in Jacobian coordinates
  //         Q = (x,y) in affine coordinates 
  // Output: P = P+Q = (X1,Y1,Z1) + (x,y) in Jacobian coordinates
    FELM t1, t2, t3, t4, t5, t6; 
    unsigned int index = 0;
    sdig mask = 0;
#ifdef ML_COUNT
    nmul+=8; nsqr+=3; nadd+=11; ncsel+=7; nclut++; 
#endif
    UNREFERENCED_PARAMETER(JacCurve);
        
    // SECURITY NOTE: this constant-time addition function is complete (i.e., it works for any possible inputs, including the cases P!=Q, P=Q, P=-Q and P=inf) on prime-order Weierstrass curves.
    //                For more information, refer to "Selecting elliptic curves for cryptography: an efficiency and security analysis", http://eprint.iacr.org/2014/130
    
    FP_SQR(P->Z, t2);                         // t2 = z1^2
    FP_MUL(P->Z, t2, t3);                     // t3 = z1^3
    FP_MUL(t2, Q->x, t1);                     // t1 = z1^2*x2
    FP_MUL(t3, Q->y, t4);                     // t4 = z1^3*y2
    FP_SUB(t1, P->X, t1);                     // t1 = beta2 = z1^2*x2-x1                
    FP_SUB(t4, P->Y, t4);                     // t4 = alpha2 = z1^3*y2-y1
    index = COMPLETE_EVAL(t1, P->Z, t4, &mask);
    FP_ADD(P->X, t2, t3);                     // t3 = x1+z1^2
    FP_SUB(P->X, t2, t6);                     // t6 = x1-z1^2
    COMPLETE_SELECT(P->Y, t1, t2, mask);      // If mask=0 (DBL) then t2=y1, else if mask=-1 (ADD) then t2=beta2 
    FP_SQR(t2, t5);                           // t5 = y1^2 (DBL) or beta2^2 (ADD)
    FP_MUL(P->X, t5, t1);                     // t1 = x1y1^2 = beta1 (DBL) or x1*beta2^2 (ADD)
    FP_MUL(P->Z, t2, table[2]->Z);            // Z2Pfinal = z1y1
    FP_COPY(table[2]->Z, table[3]->Z);        // ZPQfinal = z1*beta2
    COMPLETE_SELECT(t3, t2, t3, mask);        // If mask=0 (DBL) then t3=x1+z1^2, else if mask=-1 (ADD) then t3=beta2 
    COMPLETE_SELECT(t6, t5, t6, mask);        // If mask=0 (DBL) then t6=x1-z1^2, else if mask=-1 (ADD) then t6=beta2^2
    FP_MUL(t3, t6, t2);                       // t2 = (x1+z1^2)(x1-z1^2) (DBL) or beta2^3 (ADD)
    FP_DIV2(t2, t3);                          // t3 = (x1+z1^2)(x1-z1^2)/2
    FP_ADD(t2, t3, t3);                       // t3 = alpha1 = 3(x1+z1^2)(x1-z1^2)/2
    COMPLETE_SELECT(t3, t4, t3, mask);        // If mask=0 (DBL) then t3=alpha1, else if mask=-1 (ADD) then t3=alpha2
    FP_SQR(t3, t4);                           // t4 = alpha1^2 (DBL) or alpha2^2 (ADD)
    FP_SUB(t4, t1, t4);                       // t4 = alpha1^2-beta1 (DBL) or alpha2^2-x1*beta2^2
    FP_SUB(t4, t1, table[2]->X);              // X2Pfinal = alpha1^2-2beta1 (DBL) or alpha2^2-2x1*beta2^2 (ADD)
    FP_SUB(table[2]->X, t2, table[3]->X);     // XPQfinal = alpha^2-beta2^3-2x1*beta2^2
    COMPLETE_SELECT(table[2]->X, table[3]->X, t4, mask);   // If mask=0 (DBL) then t4=X2Pfinal, else if mask=-1 (ADD) then t4=XPQfinal
    FP_SUB(t1, t4, t1);                       // t1 = beta1-X2Pfinal (DBL) or (ADD) x1*beta2^2-XPQfinal
    FP_MUL(t3, t1, t4);                       // t4 = alpha1.(beta1-X2Pfinal) or alpha2.(x1*beta2^2-XPQfinal)
    COMPLETE_SELECT(t5, P->Y, t1, mask);      // If mask=0 (DBL) then t1=y1^2, else if mask=-1 (ADD) then t1=y1
    COMPLETE_SELECT(t5, t2, t2, mask);        // If mask=0 (DBL) then t2=y1^2, else if mask=-1 (ADD) then t2=beta2^3
    FP_MUL(t1, t2, t3);                       // t3 = y1^4 (DBL) or y1*beta2^3 (ADD)
    FP_SUB(t4, t3, table[2]->Y);              // Y2Pfinal = alpha1.(beta1-X2Pfinal)-y1^4 (DBL) or alpha2.(x1*beta2^2-XPQfinal)-y1*beta2^3 (ADD)
    FP_COPY(table[2]->Y, table[3]->Y);        // YPQfinal = alpha2.(x1*beta2^2-XPQfinal)-y1*beta2^3
    COMPLETE_LUT4(table, index, P);           // P = table[index]
    
// cleanup
    FP_ZERO(t1);
    FP_ZERO(t2);
    FP_ZERO(t3);
    FP_ZERO(t4);
    FP_ZERO(t5);
    FP_ZERO(t6);
        
    return;
}


BOOL ECC_MUL_FIXED_W(POINT_WAFF *P_table, dig *k, POINT_WAFF Q, MemType memory_use, PCurveStruct JacCurve)
{ // Wrapper for fixed-base scalar multiplication Q = k.P, where P = P_table 
  // Weierstrass a=-3 curve
    unsigned int w, v;
    
    if (memory_use < 0 || memory_use >= MemTypeSize) {
        return FALSE;
    }

    if (memory_use == MEM_LARGE) {
        w = W_MEM_LARGE;
        v = V_MEM_LARGE;
    } else if (memory_use == MEM_COMPACT) {
        w = W_MEM_COMPACT_W;
        v = V_MEM_COMPACT_W;
    } else {
        w = W_MEM_LARGE;
        v = V_MEM_LARGE;
    }

    return ECC_MUL_FIXED_INTERNAL_W(P_table, k, Q, w, v, JacCurve);
}


BOOL ECC_MUL_FIXED_INTERNAL_W(POINT_WAFF *P_table, dig *k, POINT_WAFF Q, unsigned int w, unsigned int v, PCurveStruct JacCurve)
{ // Fixed-base scalar multiplication Q = k.P, where P = P_table, using the Modified LSB-set Comb method 
  // Weierstrass a=-3 curve
    unsigned int j, npoints, d, e, l, num_words = (JacCurve->nbits+ML_WORD-1)/ML_WORD;  // Number of words to represent field elements and elements in Z_r
    int i, ii, digits[2*MAXBITLENGTH] = {0};                                            /////
    dig k_temp[((sizeof(FELM)*8)+ML_WORD-1)/ML_WORD] = {0};
    POINT_WJAC T, complete_table[4] = {0};                // Table to store {inf, Q, 2P, P+Q}. This is used in the "complete" addition.
    POINT_WAFF R;
    sdig odd = 0;
    signed long digit = 0;
    BOOL OK = FALSE;
#ifdef ML_COUNT
    nadd+=2; nlut++;
#endif
        
    // SECURITY NOTE: the crypto sensitive part of this function is protected against timing attacks and runs in constant-time on prime-order Weierstrass curves. 
    //                Conditional if-statements evaluate public data only and the number of iterations for all loops is public. Refer to "Selecting elliptic curves 
    //                for cryptography: an efficiency and security analysis", http://eprint.iacr.org/2014/130, for the full proof demonstrating exception-less, 
    //                constant-time execution of scalar multiplication.
    // DISCLAIMER:    the protocol designer is responsible for guarantying that early termination produced after detecting errors during input validation
    //                (of scalar k) does not leak any secret information. 

    
    /***** Input validation: *****/          
    // Is scalar k in [1,r-1]?                
    if ((FP_ISZERO(k) == TRUE) || (MOD_EVAL(k, JacCurve->prime) == FALSE)) {
        return FALSE;
    }
    if (P_table == NULL) {    // Full point validation is done during offline precomputation
        return FALSE;
    }
    /*****************************/ 

    e = (JacCurve->rbits+w*v-1)/(w*v);     
    if (JacCurve->rbits-e*w*v == 0) {
        return FALSE;                                  // This parameter selection is not allowed
    }
    d = e*v;     
    l = d*w;                                           // Fixed length of the mLSB-set representation
    npoints = v*(1 << (w-1));

    odd = -((sdig)k[0] & 1);     
    FP_SUB(JacCurve->order, k, k_temp);                // Converting scalar to odd (r-k if even)
    for (j = 0; j < num_words; j++)                    // If (odd) then k = k_temp else k = k 
    {
        k_temp[j] = (odd & (k[j] ^ k_temp[j])) ^ k_temp[j];
    }
    mLSB_set_recode(k_temp, JacCurve->rbits, l, d, digits); 

    // Extracting initial digit 
    digit = digits[w*d-1];
    for (i = (int)((w-1)*d-1); i >= (int)(2*d-1); i = i-d)           
    {
        digit = 2*digit + digits[i];
    }
    // Initialize T = (X_T:Y_T:1) with a point from the precomputed table
    LUT_WAFF(P_table+(v-1)*(1 << (w-1)), R, digit, digits[d-1], 1 << (w-1), JacCurve);
    ECCCONVERT_AFF_TO_JAC_W(R, T);

    // Initialize complete_table[0] with the point at infinity (0:1:0)
    complete_table[0]->Y[0] = 1;     

    for (j = 0; j < (v-1); j++)
    {
        digit = digits[w*d-(j+1)*e-1];
        for (i=(int)((w-1)*d-(j+1)*e-1); i>=(int)(2*d-(j+1)*e-1); i=i-d)           
        {
            digit = 2*digit + digits[i];
        }
        LUT_WAFF(P_table+(v-j-2)*(1 << (w-1)), R, digit, digits[d-(j+1)*e-1], 1 << (w-1), JacCurve);  // Load R = (X_R:Y_R) with point from precomputed table
#ifdef ML_COUNT
        nlut++;
#endif
        ECCCONVERT_AFF_TO_JAC_W(R, complete_table[1]);                 // Load complete_table[1] with (X_R:Y_R:1)
        ECCUMADD_WJAC(R, T, complete_table, JacCurve);                 // "Complete" mixed addition (X_T:Y_T:Z_T) = (X_T:Y_T:Z_T) + (X_R:Y_R)
    }

    for (ii = (e-2); ii >= 0; ii--)
    {
        ECCDOUBLE_WJAC(T, JacCurve);                                   // Double (X_T:Y_T:Z_T) = 2(X_T:Y_T:Z_T)
        for (j = 0; j <= (v-1); j++)
        {
            digit = digits[w*d-j*e+ii-e];

            for (i = (int)((w-1)*d-j*e+ii-e); i >= (int)(2*d-j*e+ii-e); i = i-d)           
            {
                digit = 2*digit + digits[i];
            }
            LUT_WAFF(P_table+(v-j-1)*(1 << (w-1)), R, digit, digits[d-j*e+ii-e], 1 << (w-1), JacCurve);  // Load R = (X_R:Y_R) with point from precomputed table
#ifdef ML_COUNT
            nlut++;
#endif
            ECCCONVERT_AFF_TO_JAC_W(R, complete_table[1]);            // Load complete_table[1] with (X_R:Y_R:1)
            ECCUMADD_WJAC(R, T, complete_table, JacCurve);            // "Complete" mixed addition (X_T:Y_T:Z_T) = (X_T:Y_T:Z_T) + (X_R:Y_R)
        }        
    } 

    FP_COPY(T->Y, k_temp);
    FP_NEG(JacCurve->prime, k_temp);                                 // Correcting scalar (-Ty if even)
    for (j = 0; j < num_words; j++)                                  // If (even) then Ty = -Ty 
    {
        T->Y[j] = (odd & (T->Y[j] ^ k_temp[j])) ^ k_temp[j];
    }
    OK = ECCNORM_W(T, Q, JacCurve);                                  // Output Q = (x,y)
    
// cleanup
    for (j = 0; j < ((sizeof(FELM)*8)+ML_WORD-1)/ML_WORD; j++) {
        k_temp[j] = 0;
    }
    for (j = 0; j < (2*MAXBITLENGTH); j++) {
        digits[j] = 0;
    }
    ECCZERO_WJAC(T);
    for (j = 0; j < 4; j++) {
        ECCZERO_WJAC(complete_table[j]);
    }
    ECCZERO_WAFF(R);
    
    return OK;
}


POINT_WAFF* ECC_PRECOMP_FIXED_W(POINT_WAFF P, MemType memory_use, PCurveStruct JacCurve)
{ // Precomputation scheme using affine coordinates for fixed-base scalar multiplication
  // Weierstrass a=-3 curve
    unsigned int w, v, d, e;
    
    if (memory_use < 0 || memory_use >= MemTypeSize) {
        return FALSE;
    }

    if (memory_use == MEM_LARGE) {  
        w = W_MEM_LARGE;
        v = V_MEM_LARGE;
    } else if (memory_use == MEM_COMPACT) {
        w = W_MEM_COMPACT_W;
        v = V_MEM_COMPACT_W;
    } else {
        w = W_MEM_LARGE;
        v = V_MEM_LARGE;
    }
    e = (JacCurve->rbits+w*v-1)/(w*v);    
    if (JacCurve->rbits-e*w*v == 0) {    // This parameter selection is not allowed
        return FALSE;      
    }
    d = e*v;                             

    return ECC_PRECOMP_FIXED_INTERNAL_W(P, w, v, d, e, JacCurve);
}


POINT_WAFF* ECC_PRECOMP_FIXED_INTERNAL_W(POINT_WAFF P, unsigned int w, unsigned int v, unsigned int d, unsigned int e, PCurveStruct JacCurve)
{ // Precomputation scheme using affine coordinates for fixed-base scalar multiplication. Function returns NULL on error.
  // Weierstrass a=-3 curve
	POINT_WJAC R, base[WMAX], complete_table[5] = {0};                // Table to store {inf, Q, 2P, P+Q, P}. This is used in the "complete" addition.
    POINT_WAFF *T = NULL;
    unsigned int i, j, k, npoints, index;
    unsigned long index_group;  
    FELM t1, t2, t3;   
                
    // SECURITY NOTE: precomputation for fixed-base scalar multiplication uses public inputs. 

    /***** Input validation: *****/
    // Check if P is the point at infinity (0,0)
    if (ECC_IS_INFINITY_WAFF(P, JacCurve) == TRUE) { 
        return NULL;   
    } 
    // Are (x,y) in [0,p-1]?
    if (MOD_EVAL(P->x, JacCurve->prime) == FALSE || MOD_EVAL(P->y, JacCurve->prime) == FALSE) {  
        return NULL;
    }
    // Does P lie on the curve?
    FP_SQR(P->y, t1);           // y^2
    FP_SQR(P->x, t2);
    FP_MUL(P->x, t2, t3); 
    FP_ADD(t3, JacCurve->parameter2, t2);
    FP_ADD(P->x, P->x, t3); 
    FP_ADD(P->x, t3, t3); 
    FP_SUB(t2, t3, t2);         // x^3 - 3x + b
    FP_SUB(t1, t2, t1); 
    if (FP_ISZERO(t1)==FALSE) {
        T = NULL;
        goto cleanup;
    }
    /*****************************/
        
    if (JacCurve->rbits-e*w*v == 0) return FALSE;                     // This parameter selection is not allowed

    npoints = v*(1 << (w-1));
    T = (POINT_WAFF*)calloc(npoints, sizeof(POINT_WAFF));             // Allocating space for table
    if (T == NULL) {
        goto cleanup;
    }

    // Initialize complete_table[0] with the point at infinity (0:1:0) 
    complete_table[0]->Y[0] = 1;  
    
    ECCCONVERT_AFF_TO_JAC_W(P, base[0]);                             // base[0] = P in coordinates (X:Y:1)

    // Compute base point for each w (or row)
	for (i = 0; i < (w-1); i++) {
        ECCCOPY_WJAC(base[i], R);
		for (j = 0; j < d; j++) ECCDOUBLE_WJAC(R, JacCurve);         // base[i+1] = 2^d base[i] in coordinates (X:Y:Z)
        ECCCOPY_WJAC(R, base[i+1]);
	}
    ECCCOPY_W(P, T[0]);                                              // T[0] = P
    
    // Compute precomputed points for the first table
    index = 0;
    index_group = 1;
    for (i = 0; i < (w-1); i++)                                      // T[index] = (1 + u_0.2^d + ... + u_{w-2}.2^((w-1)d)) P
    {
        for (j = 0; j < index_group; j++)
        {
            ECCCONVERT_AFF_TO_JAC_W(T[j], R);
            ECCCOPY_WJAC(base[i+1], complete_table[1]);                      // Load complete_table[1] with base[i+1]
            ECCCOPY_WJAC(R, complete_table[4]);                              // Load complete_table[4] with (X_R:Y_R:Z_R)
            ECCUADD_NO_INIT_WJAC(base[i+1], R, complete_table, JacCurve);    // Complete addition R = R + base[i+1]  
            index++;
            ECCNORM_W(R, T[index], JacCurve);
        }
        index_group = 2*index_group;
    }
        
    // Compute precomputed points for the remaining tables
    index++;
    for (i = 0; i < (v-1); i++)                                      // T[index] = 2^(ev) (1 + u_0.2^d + ... + u_{w-2}.2^((w-1)d)) P
    {
        for (j = 0; j < index; j++)
        {
            ECCCONVERT_AFF_TO_JAC_W(T[i*index + j], R);
		    for (k = 0; k < e; k++) ECCDOUBLE_WJAC(R, JacCurve);     // 2^(ev) * X * P
            ECCNORM_W(R, T[(i+1)*index + j], JacCurve);
        }
    }                                 
    
cleanup:
    ECCZERO_WJAC(R);
    for (j = 0; j < WMAX; j++) {
        ECCZERO_WJAC(base[j]);
    }
    for (j = 0; j < 5; j++) {
        ECCZERO_WJAC(complete_table[j]);
    }
    FP_ZERO(t1);
    FP_ZERO(t2);
    FP_ZERO(t3);

    return T;
}


static void ECCDOUBLEADD_CONDITIONALS_WJAC(POINT_PRECOMP_WCHU Q, POINT_WJAC P, PCurveStruct JacCurve)      
{ // Point addition P = 2P+Q containing conditionals (if-statements)
  // Weierstrass a=-3 curve
  // Inputs: P = (X1,Y1,Z1) in Jacobian coordinates
  //         Q = (X2,Y2,Z2,Z2^2,Z2^3) in Chudnovsky coordinates 
  // Output: P = (X1,Y1,Z1) in Jacobian coordinates
    FELM t1, t2, t3, t4, t5, t6, t7; 
    POINT_WJAC T;
#ifdef ML_COUNT
    nmul+=16; nsqr+=5; nadd+=13;
#endif
    
    // SECURITY NOTE: this function should only be called from functions not requiring constant-time execution such as double-scalar multiplications 
    //                ecc_double_scalar_mul_internal_Jacxxx() used in signature verification.

    // Check if P is the point at infinity (0:Y:0)
    if (ECC_IS_INFINITY_WJAC(P, JacCurve) == TRUE) { 
        FP_COPY(Q->X, P->X); FP_COPY(Q->Y, P->Y); FP_COPY(Q->Z, P->Z);  // Output P = Q
        return;   
    }    
    // Check if Q is the point at infinity (0:Y:0)
    FP_COPY(Q->X, T->X); FP_COPY(Q->Y, T->Y); FP_COPY(Q->Z, T->Z);
    if (ECC_IS_INFINITY_WJAC(T, JacCurve) == TRUE) {        
        ECCDOUBLE_WJAC(P, JacCurve);                                    // Output P = 2*P
        return;
    }
    
    FP_SQR(P->Z, t2);               // t2 = z1^2
    FP_MUL(Q->Z3, P->Y, t3);        // t3 = z2^3*y1
    FP_MUL(P->Z, t2, t4);           // t4 = z1^3
    FP_MUL(t2, Q->X, t1);           // t1 = z1^2*x2                      
    FP_MUL(Q->Y, t4, t2);           // t2 = z1^3*y2
    FP_MUL(Q->Z2, P->X, t6);        // t6 = z2^2*x1
    FP_SUB(t2, t3, t2);             // t2 = alpha = z1^3*y2-z2^3*y1
    FP_SUB(t1, t6, t1);             // t1 = beta = z1^2*x2-z2^2*x1
    
    if ((FP_ISZERO(t1) & FP_ISZERO(t2)) == TRUE) {
        FP_COPY(P->X, T->X); FP_COPY(P->Y, T->Y); FP_COPY(P->Z, T->Z);
        FP_NEG(JacCurve->prime, T->Y);         // T = -P  
        ECCDOUBLE_WJAC(P, JacCurve);
        ECCDOUBLEADD_WJAC(Q, P, JacCurve);     // Output P = 2*(2P)-P = 3*P
        goto cleanup;
    } 
    if (FP_ISZERO(t1) == TRUE) goto cleanup;   // Output P

    FP_SQR(t2, t4);                 // t4 = alpha^2
    FP_SQR(t1, t5);                 // t5 = beta^2
    FP_MUL(P->Z, Q->Z, t7);         // t5 = z1*z2
    FP_MUL(t6, t5, P->X);           // x1 = x1' = z2^2*x1*beta^2
    FP_MUL(t1, t5, t6);             // t6 = beta^3
    FP_SUB(t4, t6, t4);             // t4 = alpha^2 - beta^3
    FP_SUB(t4, P->X, t4);           // t4 = alpha^2 - beta^3 - x1'
    FP_SUB(t4, P->X, t4);           // t4 = alpha^2 - beta^3 - 2*x1'
    FP_SUB(t4, P->X, t4);           // t4 = omega = alpha^2 - beta^3 - 3*x1'
    
    if (FP_ISZERO(t4) == TRUE) {
        FP_ZERO(P->X); FP_ZERO(P->Z);   // Output point at infinity (0:Y:0)  
        goto cleanup;                                    
    }

    FP_MUL(t6, t3, P->Y);           // y1 = y1' = z2^3*y1*beta^3 
    FP_MUL(t1, t7, t3);             // t3 = z1' = z1*z2*beta
    FP_MUL(t2, t4, t1);             // t1 = alpha.omega
    FP_SQR(t4, t2);                 // t2 = omega^2                      
    FP_ADD(t1, P->Y, t1);           // t1 = alpha.omega + y1'
    FP_ADD(t1, P->Y, t1);           // t1 = theta = alpha.omega + 2y1'    
    FP_MUL(t3, t4, P->Z);           // Zfinal = z1'*omega
    FP_MUL(t2, t4, t5);             // t5 = omega^3
    FP_MUL(t2, P->X, t4);           // t4 = x1'*omega^2
    FP_SQR(t1, t3);                 // t3 = theta^2
    FP_SUB(t3, t5, t3);             // t3 = theta^2 - omega^3
    FP_SUB(t3, t4, t3);             // t3 = theta^2 - omega^3 - x1'*omega^2
    FP_SUB(t3, t4, P->X);           // Xfinal = theta^2 - omega^3 - 2*x1'*omega^2
    FP_SUB(P->X, t4, t3);           // t3 = Xfinal-x1'*omega^2
    FP_MUL(P->Y, t5, t2);           // t2 = y1'*omega^3
    FP_MUL(t3, t1, t5);             // t5 = theta.(Xfinal-x1'*omega^2)
    FP_SUB(t5, t2, P->Y);           // Yfinal = theta.(Xfinal-x1'*omega^2) - y1'*omega^3
    
cleanup:
    ECCZERO_WJAC(T);
    FP_ZERO(t1);
    FP_ZERO(t2);
    FP_ZERO(t3);
    FP_ZERO(t4);
    FP_ZERO(t5);
    FP_ZERO(t6);
    FP_ZERO(t7);

    return;
}


static void ECCMADD_CONDITIONALS_WJAC(POINT_WAFF Q, POINT_WJAC P, PCurveStruct JacCurve)      
{ // Point addition P = P+Q
  // Weierstrass a=-3 curve
  // Inputs: P = (X1,Y1,Z1) in Jacobian coordinates
  //         Q = (x,y) in affine coordinates 
  // Output: P = (X1,Y1,Z1) in Jacobian coordinates
    FELM t1, t2, t3, t4; 
#ifdef ML_COUNT
    nmul+=8; nsqr+=3; nadd+=7; 
#endif
    
    // SECURITY NOTE: this function should only be called from functions not requiring constant-time execution such as double-scalar multiplications 
    //                ecc_double_scalar_mul_internal_Jacxxx() used in signature verification.

    // Check if P is the point at infinity (0:Y:0)
    if (ECC_IS_INFINITY_WJAC(P, JacCurve) == TRUE) {   
        ECCCONVERT_AFF_TO_JAC_W(Q, P);                       // Output P = Q = (X:Y:1)
        return;   
    }       
    // Check if Q is the point at infinity (0,0)
    if (ECC_IS_INFINITY_WAFF(Q, JacCurve) == TRUE) {
        return;                                             // Output P
    }

    FP_SQR(P->Z, t2);               // t2 = z1^2
    FP_MUL(P->Z, t2, t3);           // t3 = z1^3
    FP_MUL(t2, Q->x, t1);           // t1 = z1^2*x2
    FP_MUL(Q->y, t3, t2);           // t2 = z1^3*y2
    FP_SUB(t2, P->Y, t2);           // t2 = alpha = z1^3*y2-y1
    FP_SUB(t1, P->X, t1);           // t1 = beta = z1^2*x2-x1

    if ((FP_ISZERO(t1) & FP_ISZERO(t2)) == TRUE) {
        ECCDOUBLE_WJAC(P, JacCurve);     // Output P = 2*P since P=Q   
        goto cleanup;
    } 
    if (FP_ISZERO(t1) == TRUE) {
        FP_ZERO(P->X); FP_ZERO(P->Z);    // Output the point at infinity P = (0:Y:0) since P=-Q  
        goto cleanup;
    }

    FP_COPY(P->Z, t4);              // t4 = z1
    FP_MUL(t1, t4, P->Z);           // Zfinal = z1*beta
    FP_SQR(t1, t4);                 // t4 = beta^2
    FP_MUL(t1, t4, t3);             // t3 = beta^3
    FP_MUL(P->X, t4, t1);           // t1 = x1*beta^2
    FP_SQR(t2, t4);                 // t4 = alpha^2
    FP_SUB(t4, t3, t4);             // t4 = alpha^2 - beta^3
    FP_SUB(t4, t1, t4);             // t4 = alpha^2 - beta^3 - z2^2*x1*beta^2
    FP_SUB(t4, t1, P->X);           // Xfinal = alpha^2 - beta^3 - 2*z2^2*x1*beta^2
    FP_SUB(t1, P->X, t1);           // t1 = z2^2*x1*beta^2-Xfinal
    FP_MUL(t2, t1, t4);             // t4 = alpha.(z2^2*x1*beta^2-Xfinal)
    FP_MUL(P->Y, t3, t2);           // t2 = y1*beta^3
    FP_SUB(t4, t2, P->Y);           // Yfinal = alpha.(z2^2*x1*beta^2-Xfinal) - y1*beta^3
    
cleanup:
    FP_ZERO(t1);
    FP_ZERO(t2);
    FP_ZERO(t3);
    FP_ZERO(t4);

    return;
}


BOOL ECC_DBLMUL_W(POINT_WAFF *P_table, dig *k, POINT_WAFF Q, dig *l, POINT_WAFF R, MemType memory_use, PCurveStruct JacCurve)
{ // Wrapper for double-base scalar multiplication R = k.P + l.Q, where P = P_table
  // P is a fixed-base and Q is a variable-base
  // Weierstrass a=-3 curve
    unsigned int w_P;
    
    if (memory_use < 0 || memory_use >= MemTypeSize) {
        return FALSE;
    }

    if (memory_use == MEM_LARGE) {
        w_P = W_P_MEM_LARGE;
    } else if (memory_use == MEM_COMPACT) {
        w_P = W_P_MEM_COMPACT;
    } else {
        w_P = W_P_MEM_LARGE;
    }

    return ECC_DBLMUL_INTERNAL_W(P_table, k, Q, l, R, w_P, JacCurve);
}


BOOL ECC_DBLMUL_INTERNAL_W(POINT_WAFF *P_table, dig *k, POINT_WAFF Q, dig *l, POINT_WAFF R, unsigned int w_P, PCurveStruct JacCurve)
{ // Double-base scalar multiplication R = k.P + l.Q, where P = P_table[0], using wNAF with Interleaving
  // P is a fixed-base and Q is a variable-base
  // Weierstrass a=-3 curve
    unsigned int npoints, position;
    int i, digits_P[sizeof(FELM)*8 + 1]={0}, digits_Q[sizeof(FELM)*8 + 1]={0};         
    POINT_WJAC T; 
    POINT_PRECOMP_WCHU table[1 << (W_VARBASE-2)], S;
    POINT_WAFF SS;
    FELM t1, t2, t3;
    BOOL OK = FALSE;
#ifdef ML_COUNT
    nmul++; nsqr+=2; nadd+=5; 
#endif 
            
    // SECURITY NOTE: this function is intended for an operation not requiring constant-time execution such as signature verification. 
    
    /***** Input validation: *****/   
    if (P_table == NULL) {    // Full point validation for P is done during offline precomputation
        return FALSE;
    }  
    // Check if Q is the point at infinity (0,0)
    if (ECC_IS_INFINITY_WAFF(Q, JacCurve) == TRUE) {    
        return FALSE;   
    } 
    // Are scalars k, l in [1,r-1]?                
    if ((FP_ISZERO(k) == TRUE) || (MOD_EVAL(k, JacCurve->prime) == FALSE) || (FP_ISZERO(l) == TRUE) || (MOD_EVAL(l, JacCurve->prime) == FALSE)) {
        return FALSE;
    }
    // Are Q: (x,y) in [0,p-1]?
    if (MOD_EVAL(Q->x, JacCurve->prime) == FALSE || MOD_EVAL(Q->y, JacCurve->prime) == FALSE) {  
        return FALSE;
    }
    // Does Q lie on the curve?
    FP_SQR(Q->y, t1);           // y^2
    FP_SQR(Q->x, t2);
    FP_MUL(Q->x, t2, t3); 
    FP_ADD(t3, JacCurve->parameter2, t2);
    FP_ADD(Q->x, Q->x, t3); 
    FP_ADD(Q->x, t3, t3); 
    FP_SUB(t2, t3, t2);         // x^3 - 3x + b
    FP_SUB(t1, t2, t1); 
    if (FP_ISZERO(t1)==FALSE) {
        return FALSE;
    }
    /****************************/

    npoints = 1 << (W_VARBASE-2); 

    ECC_PRECOMP_WJAC(Q, table, npoints, JacCurve);               // Precomputation of points table[0],...,table[npoints-1]
    wNAF_recode(k, JacCurve->rbits, w_P, digits_P);              // Recode k and l to the wNAF representation
    wNAF_recode(l, JacCurve->rbits, W_VARBASE, digits_Q);
    FP_ZERO(T->X); FP_ZERO(T->Y); T->Y[0] = 1; FP_ZERO(T->Z);    // Initialize T as the point at infinity (0:1:0)  

    for (i = JacCurve->rbits; i >= 0; i--)
    {
        if (digits_Q[i] == 0) {
            ECCDOUBLE_WJAC(T, JacCurve);                         // Double (T_X:T_Y:T_Z) = 2(T_X:T_Y:T_Z)
        } else if (digits_Q[i] < 0) {
            position = (-digits_Q[i])/2;                         // Load S = (X_S:Y_S:Z_S:Z_S^2:Z_S^3) with a point from the precomputed table
            FP_COPY(table[position]->X, S->X);
            FP_COPY(table[position]->Y, S->Y);
            FP_COPY(table[position]->Z, S->Z);
            FP_COPY(table[position]->Z2, S->Z2);
            FP_COPY(table[position]->Z3, S->Z3);                
            FP_NEG(JacCurve->prime, S->Y);                      // Negate S
            ECCDOUBLEADD_CONDITIONALS_WJAC(S, T, JacCurve);     // Double-add (X_T:Y_T:Z_T) = 2(X_T:Y_T:Z_T) + (X_S:Y_S:Z_S:Z_S^2:Z_S^3)
#ifdef ML_COUNT
            nadd++; 
#endif
        } else if (digits_Q[i] > 0) {            
            position = (digits_Q[i])/2;                         // Load S = (X_S:Y_S:Z_S:Z_S^2:Z_S^3) with a point from the precomputed table
            FP_COPY(table[position]->X, S->X);
            FP_COPY(table[position]->Y, S->Y);
            FP_COPY(table[position]->Z, S->Z);
            FP_COPY(table[position]->Z2, S->Z2);
            FP_COPY(table[position]->Z3, S->Z3);
            ECCDOUBLEADD_CONDITIONALS_WJAC(S, T, JacCurve);     // Double-add (X_T:Y_T:Z_T) = 2(X_T:Y_T:Z_T) + (X_S:Y_S:Z_S:Z_S^2:Z_S^3)
        }

        if (digits_P[i] < 0) {                           
            position = (-digits_P[i])/2;                        // Load SS = (x_SS:y_SS) with a point from the precomputed table
            ECCCOPY_W(P_table[position], SS);
            FP_NEG(JacCurve->prime, SS->y);                     // Negate SS
            ECCMADD_CONDITIONALS_WJAC(SS, T, JacCurve);         // Mixed double-add (X_T:Y_T:Z_T) = 2(X_T:Y_T:Z_T) + (x_SS:y_SS)
#ifdef ML_COUNT
            nadd++; 
#endif
        } else if (digits_P[i] > 0) { 
            position = (digits_P[i])/2;                         // Load SS = (x_SS:y_SS) with a point from the precomputed table
            ECCCOPY_W(P_table[position], SS);
            ECCMADD_CONDITIONALS_WJAC(SS, T, JacCurve);         // Mixed double-add (X_T:Y_T:Z_T) = 2(X_T:Y_T:Z_T) + (x_SS:y_SS)
        }
    }
    OK = ECCNORM_W(T, R, JacCurve);                             // Output R = (x,y)
    
// cleanup
    for (i = 0; i < (sizeof(FELM)*8 + 1); i++) {
        digits_P[i] = 0;
        digits_Q[i] = 0;
    }
    ECCZERO_WJAC(T);
    for (i = 0; i < (1 << (W_VARBASE-2)); i++) {
        ECCZERO_WCHU(table[i]);
    }
    ECCZERO_WCHU(S);
    ECCZERO_WAFF(SS);
    FP_ZERO(t1);
    FP_ZERO(t2);
    FP_ZERO(t3);
    
    return OK;
}


POINT_WAFF* ECC_PRECOMP_DBLMUL_W(POINT_WAFF P, MemType memory_use, PCurveStruct JacCurve)
{ // Wrapper for precomputation scheme using affine coordinates for the fixed-base in double-scalar multiplication. Function returns NULL on error.
  // Weierstrass a=-3 curve
    unsigned int w_P;
    
    if (memory_use < 0 || memory_use >= MemTypeSize) {
        return FALSE;
    }

    if (memory_use == MEM_LARGE) { 
        w_P = W_P_MEM_LARGE;
    } else if (memory_use == MEM_COMPACT) {
        w_P = W_P_MEM_COMPACT;
    } else {
        w_P = W_P_MEM_LARGE;
    }                          

    return ECC_PRECOMP_DBLMUL_INTERNAL_W(P, w_P, JacCurve);
}



POINT_WAFF* ECC_PRECOMP_DBLMUL_INTERNAL_W(POINT_WAFF P, unsigned int w_P, PCurveStruct JacCurve)
{ // Precomputation scheme using affine coordinates for the fixed base in double-scalar multiplication. Function returns NULL on error.
  // Weierstrass a=-3 curve
    POINT_WAFF *T = NULL;
    POINT_PRECOMP_WCHU *T_CHU = NULL;
    unsigned int i, npoints; 
    FELM t1, t2, t3;  
                
    // SECURITY NOTE: precomputation for double-scalar multiplication uses public inputs. 

    /***** Input validation: *****/
    // Check if P is the point at infinity (0,0)
    if (ECC_IS_INFINITY_WAFF(P, JacCurve) == TRUE) { 
        return NULL;   
    } 
    // Are (x,y) in [0,p-1]?
    if (MOD_EVAL(P->x, JacCurve->prime) == FALSE || MOD_EVAL(P->y, JacCurve->prime) == FALSE) {  
        return NULL;
    }
    // Does P lie on the curve?
    FP_SQR(P->y, t1);           // y^2
    FP_SQR(P->x, t2);
    FP_MUL(P->x, t2, t3); 
    FP_ADD(t3, JacCurve->parameter2, t2);
    FP_ADD(P->x, P->x, t3); 
    FP_ADD(P->x, t3, t3); 
    FP_SUB(t2, t3, t2);         // x^3 - 3x + b
    FP_SUB(t1, t2, t1); 
    if (FP_ISZERO(t1)==FALSE) {
        T = NULL;
        goto cleanup;
    }
    /*****************************/

    npoints = (1 << (w_P-2));
    T_CHU = (POINT_PRECOMP_WCHU*)calloc(npoints, sizeof(POINT_PRECOMP_WCHU));    // Allocating space for temporary table
    if (T_CHU == NULL) {
        T = NULL;
        goto cleanup;
    }
    T = (POINT_WAFF*)calloc(npoints, sizeof(POINT_WAFF));                        // Allocating space for precomputed table
    if (T == NULL) {
        goto cleanup;
    }

    ECC_PRECOMP_WJAC(P, T_CHU, npoints, JacCurve);
    for (i=0; i<npoints; i++)
    {            
        FP_COPY(T_CHU[i]->Z, t1);  
        FP_INV(t1);                          // t1 = Z^-1
        FP_SQR(t1, t2);                      // t2 = Z^-2
        FP_MUL(T_CHU[i]->X, t2, T[i]->x);    // x = X/Z^2
        FP_MUL(t1, t2, t3);                  // t3 = Z^-3
        FP_MUL(T_CHU[i]->Y, t3, T[i]->y);    // y = Y/Z^3 
    }
    
cleanup:
    if (T_CHU != NULL) {
        free(T_CHU);
    }
    FP_ZERO(t1);
    FP_ZERO(t2);
    FP_ZERO(t3);

    return T;
}


BOOL ECC_DESTROY_PRECOMP_W(POINT_WAFF* T_fixed)
{ // Frees memory occupied by a precomputation table "T_fixed" used during fixed-base or double-scalar multiplications.
  // This function must be called once done using a table generated by ecc_precomp_fixed_Jacxxx or ecc_precomp_dblmul_Jacxxx. 
  // Weierstrass a=-3 curve
        
    if (T_fixed != NULL) {
        free(T_fixed); 
    }

    return TRUE;
}



/*******************************************************************************************/
/********** POINT/SCALAR FUNCTIONS FOR TWISTED EDWARDS a=-1 CURVE ("TED" CURVES) ***********/

void ECCSET_TE(POINT_TE P, PCurveStruct TedCurve)
{ // Set generator (x,y) on twisted Edwards a=-1 curve
    dig i, num_words = (TedCurve->pbits+ML_WORD-1)/ML_WORD;

    for (i = 0; i < num_words; i++) 
    {
        P->x[i] = TedCurve->generator_x[i]; 
        P->y[i] = TedCurve->generator_y[i]; 
    }
    return;
}


BOOL ECC_IS_NEUTRAL_AFF_TE(POINT_TE P, PCurveStruct TedCurve)
{ // Check if affine point P is the neutral point (0,1) 
    dig i, c, num_words = (TedCurve->pbits+ML_WORD-1)/ML_WORD;
    BOOL answer = TRUE;

    c = P->x[0] | (P->y[0] ^ 1);
    for (i = 1; i < num_words; i++)
    {
        c = c | P->x[i] | P->y[i]; 
    }
    answer = (BOOL)((~(0 - c) >> (ML_WORD-1)) & (~c >> (ML_WORD-1)));

    return answer;
}


BOOL ECC_IS_NEUTRAL_EXT_TE(POINT_EXT_TE P, PCurveStruct TedCurve)
{ // Check if projective point P is the neutral point (0:1:1) 
    dig i, c, num_words = (TedCurve->pbits+ML_WORD-1)/ML_WORD;
    BOOL answer = TRUE;
    
    c = P->X[0] | ((P->Y[0] | P->Z[0]) ^ 1);
    for (i = 1; i < num_words; i++)
    {
        c = c | P->X[i] | P->Y[i] | P->Z[i];
    }
    answer = (BOOL)((~(0 - c) >> (ML_WORD-1)) & (~c >> (ML_WORD-1)));

    return answer;
}


void ECCNORM_TE(POINT_EXT_TE Q, POINT_TE P, PCurveStruct TedCurve)
{ // Normalize a projective twisted Edwards point Q = (X:Y:Z) -> P = (x,y)
    FELM t1, t2;   
#ifdef ML_COUNT
    ninv++; nmul+=2; 
#endif
    UNREFERENCED_PARAMETER(TedCurve);
    
    FP_COPY(Q->Z, t1);  
    FP_INV(t1);                      // t1 = Z^-1
    FP_MUL(Q->X, t1, t2);            // t2 = X/Z
    FP_COPY(t2, P->x);               // x = X/Z
    FP_MUL(Q->Y, t1, t2);            // t2 = Y/Z
    FP_COPY(t2, P->y);               // y = Y/Z
    
// cleanup
    FP_ZERO(t1);
    FP_ZERO(t2);         

    return;
}


void ECCDOUBLE_EXT_TE(POINT_EXT_TE P, PCurveStruct TedCurve)
{ // Point doubling 2P
  // Twisted Edwards a=-1 curve
  // Input: P = (X1,Y1,Z1) in twisted Edwards coordinates
  // Output: 2P = (X2,Y2,Z2,Ta,Tb) in extended twisted Edwards coordinates, where T2 = Ta*Tb
    FELM t1, t2;    
#ifdef ML_COUNT
    nmul+=4; nsqr+=3; nadd+=5;
#endif
    UNREFERENCED_PARAMETER(TedCurve);
          
    // SECURITY NOTE: this function does not produce exceptions. For more information, refer to "Selecting elliptic curves for cryptography: 
    //                an efficiency and security analysis", http://eprint.iacr.org/2014/130
        
    FP_SQR(P->X, t1);               // t1 = x1^2
    FP_SQR(P->Y, t2);               // t2 = y1^2
    FP_ADD(P->Y, P->Y, P->Y);       // y = 2y1
    FP_ADD(t1, t2, P->Tb);          // Tb = x1^2+y1^2
    FP_SUB(t2, t1, t1);             // t1 = y1^2-x1^2
    FP_SQR(P->Z, t2);               // t2 = z1^2              
    FP_MUL(P->X, P->Y, P->Ta);      // Ta = 2x1y1
    FP_ADD(t2, t2, t2);             // t2 = 2z1^2
    FP_SUB(t2, t1, t2);             // t2 = 2z1^2-(y1^2-x1^2)
    FP_MUL(t1, P->Tb, P->Y);        // Yfinal = (x1^2+y1^2)(y1^2-x1^2)
    FP_MUL(t1, t2, P->Z);           // Zfinal = (y1^2-x1^2)[2z1^2-(y1^2-x1^2)]    
    FP_MUL(t2, P->Ta, P->X);        // Xfinal = 2x1y1[2z1^2-(y1^2-x1^2)]
    
// cleanup
    FP_ZERO(t1);
    FP_ZERO(t2);

    return;
}


void ECCUADD_EXT_TE(POINT_EXT_TE Q, POINT_EXT_TE P, PCurveStruct TedCurve)      
{ // Unified point addition P = P+Q or P = P+P
  // Twisted Edwards a=-1 curve
  // Inputs: P = (X1,Y1,Z1,T1a,T1b) in extended twisted Edwards coordinates, where T1 = T1a*T1b
  //         Q = (X2,Y2,Z2,T2a,T2b) in extended twisted Edwards coordinates, where T2 = T2a*T2b    
  // Output: P = (X1,Y1,Z1,T1a,T1b) in extended twisted Edwards coordinates, where T1 = T1a*T1b
    POINT_PRECOMP_EXT_TE QQ;
    FELM t1;    
#ifdef ML_COUNT
    nmul+=2; nadd+=4;
#endif
    UNREFERENCED_PARAMETER(TedCurve);

    // SECURITY NOTE: this addition function is complete (i.e., it works for any possible inputs, including the cases P!=Q, P=Q, P=-Q and P=neutral) on large prime-order twisted Edwards 
    //                subgroups. For more information, refer to "Selecting elliptic curves for cryptography: an efficiency and security analysis", http://eprint.iacr.org/2014/130
            
    FP_MUL(Q->Ta, TedCurve->parameter2, t1);    // t1 = d*Ta
    FP_ADD(Q->X, Q->Y, QQ->XY);                 // X_QQ = X2+Y2
    FP_SUB(Q->Y, Q->X, QQ->YX);                 // Y_QQ = Y2-X2
    FP_ADD(t1, t1, t1);                         // t1 = 2*d*Ta
    FP_ADD(Q->Z, Q->Z, QQ->Z2);                 // Z_QQ = 2*Z2
    FP_MUL(Q->Tb, t1, QQ->T2);                  // T_QQ = 2*d*T2

    ECCUADD_EXT_INTERNAL_TE(QQ, P, TedCurve);
    
// cleanup
    FP_ZERO(t1);
    ECCZERO_PRECOMP_EXT_TE(QQ);

    return;
}


void ECCUADD_EXT_INTERNAL_TE(POINT_PRECOMP_EXT_TE Q, POINT_EXT_TE P, PCurveStruct TedCurve)      
{ // Unified point addition P = P+Q or P = P+P
  // Twisted Edwards a=-1 curve
  // Inputs: P = (X1,Y1,Z1,Ta,Tb) in extended twisted Edwards coordinates, where T1 = Ta*Tb
  //         Q = (X2+Y2,Y2-X2,2Z2,2dT2) in extended twisted Edwards coordinates    
  // Output: P = (X1,Y1,Z1,Ta,Tb) in extended twisted Edwards coordinates, where T1 = Ta*Tb
    FELM t1, t2, t3;    
#ifdef ML_COUNT
    nmul+=8; nadd+=6;
#endif
    UNREFERENCED_PARAMETER(TedCurve);

    // SECURITY NOTE: this addition function is complete (i.e., it works for any possible inputs, including the cases P!=Q, P=Q, P=-Q and P=neutral) on large prime-order twisted Edwards 
    //                subgroups. For more information, refer to "Selecting elliptic curves for cryptography: an efficiency and security analysis", http://eprint.iacr.org/2014/130
        
    FP_MUL(P->Ta, P->Tb, t1);           // t1 = T1
    FP_COPY(t1, P->Ta);                 // Ta = T1
    FP_MUL(Q->Z2, P->Z, t1);            // t1 = 2z1.z2        
    FP_MUL(P->Ta, Q->T2, t2);           // t2 = 2dT1.T2 
    FP_COPY(t2, P->Ta);
    FP_ADD(P->X, P->Y, t2);             // t2 = (x1+y1)    
    FP_SUB(t1, P->Ta, t3);              // t3 = theta
    FP_ADD(t1, P->Ta, t1);              // t1 = alpha
    FP_SUB(P->Y, P->X, P->Z);           // z = (y1-x1)
    FP_MUL(Q->XY, t2, P->Ta);           // Ta = (x1+y1)(x2+y2)
    FP_MUL(Q->YX, P->Z, P->X);          // x1 = (y1-x1)(y2-x2)
    FP_SUB(P->Ta, P->X, P->Tb);         // Tbfinal = beta
    FP_ADD(P->Ta, P->X, P->Ta);         // Tafinal = omega
    FP_MUL(P->Tb, t3, P->X);            // Xfinal = beta.theta
    FP_MUL(t1, t3, P->Z);               // Zfinal = theta. alpha
    FP_MUL(P->Ta, t1, P->Y);            // Yfinal = alpha.omega
    
// cleanup
    FP_ZERO(t1);
    FP_ZERO(t2);
    FP_ZERO(t3);

    return;
}


// Functions for the precomputation (twisted Edwards a=-1 curves)

static void ECCADD_PRECOMP_TE(POINT_PRECOMP_EXT_TE P, POINT_PRECOMP_EXT_TE Q, POINT_PRECOMP_EXT_TE R, PCurveStruct TedCurve)      
{ // Special point addition R = P+Q for the precomputation
  // Twisted Edwards a=-1 curve
  // Inputs:  P = (X1+Y1,Y1-X1,Z2,T2) in modified extended coordinates 
  //          Q = (X2+Y2,Y2-X2,2Z2,2dT2) in extended twisted Edwards coordinates
  // Outputs: R = (X3+Y3,Y3-X3,2Z3,2dT3) in extended twisted Edwards coordinates
    FELM t1, t2, t3, t4;    
#ifdef ML_COUNT
    nmul+=9; nadd+=8;
#endif
     
    // SECURITY NOTE: this function does not produce exceptions in the context of variable-base precomputation. For more information,
    //                refer to "Selecting elliptic curves for cryptography: an efficiency and security analysis", http://eprint.iacr.org/2014/130
    
    FP_MUL(Q->Z2, P->Z2, t1);               // t1 = 2z1.z2        
    FP_MUL(P->T2, Q->T2, t2);               // t2 = 2dT1.T2   
    FP_SUB(t1, t2, t3);                     // t3 = theta
    FP_ADD(t1, t2, t1);                     // t1 = alpha
    FP_MUL(Q->XY, P->XY, t2);               // t2 = (x1+y1)(x2+y2)
    FP_MUL(Q->YX, P->YX, R->Z2);            // RZ2 = (y1-x1)(y2-x2)
    FP_SUB(t2, R->Z2, R->T2);               // RT2 = beta
    FP_ADD(t2, R->Z2, t2);                  // t2 = omega
    FP_MUL(R->T2, t3, R->XY);               // RX2 = beta.theta
    FP_MUL(t1, t3, R->Z2);                  // RZ2 = theta. alpha
    FP_MUL(R->T2, t2, t3);                  // t3 = beta.omega
    FP_MUL(t1, t2, t4);                     // t4 = alpha.omega
    FP_SUB(t4, R->XY, R->YX);               // RYXfinal = Y3-X3        
    FP_ADD(R->XY, t4, R->XY);               // RXYfinal = X3+Y3
    FP_ADD(R->Z2, R->Z2, R->Z2);            // RZ2final = 2Z3
    FP_ADD(t3, t3, t3);                     // RT2 = 2T3
    FP_MUL(TedCurve->parameter2, t3, R->T2);  // RT2final = 2dT3  
    
// cleanup
    FP_ZERO(t1);
    FP_ZERO(t2);
    FP_ZERO(t3);
    FP_ZERO(t4);
    
    return;
}


void ECC_PRECOMP_EXT_TE(POINT_EXT_TE P, POINT_PRECOMP_EXT_TE *T, unsigned int npoints, PCurveStruct TedCurve)
{ // Precomputation scheme using extended coordinates (X+Y,Y-X,2Z,2dT)
  // Twisted Edwards a=-1 curve
    POINT_PRECOMP_EXT_TE P2;
    POINT_EXT_TE Q;
    FELM t1; 
    unsigned int i;  
#ifdef ML_COUNT
    nmul+=3; nadd+=6;
#endif
    
    // SECURITY NOTE: this function does not produce exceptions in the context of variable-base scalar multiplication and double-scalar multiplication.
    //                For more information, refer to "Selecting elliptic curves for cryptography: an efficiency and security analysis", http://eprint.iacr.org/2014/130

    // Generating P2 = 2(X1,Y1,Z1,T1a,T1b) = (XP2,YP2,ZP2,TP2) and T[0] = P = (X1,Y1,Z1,T1) 
    ECCCOPY_EXT_TE(P, Q);
    FP_COPY(P->X, T[0]->XY);
    FP_COPY(P->Y, T[0]->YX);
    FP_COPY(P->Z, T[0]->Z2);
    FP_MUL(P->Ta, P->Tb, T[0]->T2);
    ECCDOUBLE_EXT_TE(Q, TedCurve);

    // Converting P2 to a modified extended coordinate system (XP2+YP2,Y2P-X2P,ZP2,TP2)
    FP_ADD(Q->X, Q->Y, P2->XY);                 // P2.xy = XP2+YP2
    FP_SUB(Q->Y, Q->X, P2->YX);                 // P2.yx = Y2P-X2P
    FP_COPY(Q->Z, P2->Z2);                      // P2.z2 = Z2P
    FP_MUL(Q->Ta, Q->Tb, P2->T2);               // P2.t2 = T2P
    
    // Converting P[0] to extended twisted Edwards representation for precomputed points (X1+Y1,Y1-X1,2*Z1,2*d*T1)
    FP_ADD(T[0]->XY, T[0]->YX, t1);             // t1 = X1+Y1
    FP_SUB(T[0]->YX, T[0]->XY, T[0]->YX);       // T[0].yx = Y1-X1
    FP_COPY(t1, T[0]->XY);                      // T[0].xy = X1+Y1
    FP_ADD(T[0]->Z2, T[0]->Z2, T[0]->Z2);       // T[0].z2 = 2Z1
    FP_ADD(T[0]->T2, T[0]->T2, t1);             // T[0].t2 = 2T1
    FP_MUL(TedCurve->parameter2, t1, T[0]->T2); // T[0].t2 = 2dT1  
    
    for (i=1; i<npoints; i++) {
        // T[i] = 2P+T[i-1] = (2*i+1)P = (XP2+YP2,Y2P-X2P,ZP2,TP2) + (X_(2*i-1)+Y_(2*i-1), Y_(2*i-1)-X_(2*i-1), 2Z_(2*i-1), 2dT_(2*i-1)) = (X_(2*i+1)+Y_(2*i+1), Y_(2*i+1)-X_(2*i+1), 2Z_(2*i+1), 2dT_(2*i+1))
        ECCADD_PRECOMP_TE(P2, T[i-1], T[i], TedCurve);
    }
    
// cleanup
    ECCZERO_PRECOMP_EXT_TE(P2);
    ECCZERO_EXT_TE(Q);
    FP_ZERO(t1);

    return;
}


BOOL ECC_MUL_TE(POINT_TE P, dig *k, POINT_TE Q, PCurveStruct TedCurve)
{ // Variable-base scalar multiplication Q = k.P using fixed-window method 
  // Twisted Edwards a=-1 curve
    unsigned int t = (TedCurve->rbits+(W_VARBASE-2))/(W_VARBASE-1);        // Fixed length of the fixed window representation 
    unsigned int npoints = 1 << (W_VARBASE-2);
    unsigned int num_words = (TedCurve->nbits+ML_WORD-1)/(ML_WORD);        // Number of words to represent field elements or elements in Z_r
    dig j, k_temp[((sizeof(FELM)*8)+ML_WORD-1)/ML_WORD] = {0};
    int digits[((sizeof(FELM)*8)+W_VARBASE-2)/(W_VARBASE-1) + 1]={0};
    sdig i, odd = 0;
    POINT_EXT_TE T; 
    POINT_PRECOMP_EXT_TE table[1 << (W_VARBASE-2)], R;
    FELM t1, t2, t3, t4;  
#ifdef ML_COUNT
    nmul+=2; nsqr+=2; nadd+=9; nlut+=(t+1);
#endif 
        
    // SECURITY NOTE: the crypto sensitive part of this function is protected against timing attacks and runs in constant-time on large prime-order twisted Edwards subgroups. 
    //                Conditional if-statements evaluate public data only and the number of iterations for all loops is public. Refer to "Selecting elliptic curves 
    //                for cryptography: an efficiency and security analysis", http://eprint.iacr.org/2014/130, for the full proof demonstrating exception-less, 
    //                constant-time execution of scalar multiplication.
    // DISCLAIMER:    the protocol designer is responsible for guarantying that early termination produced after detecting errors during input validation
    //                (of scalar k or base point P) does not leak any secret information. 

    
    /***** Input validation and elimination of small torsion: *****/
    // Check if P is the neutral point (0,1) 
    if (ECC_IS_NEUTRAL_AFF_TE(P, TedCurve) == TRUE) {                                     
        return FALSE;   
    } 
    // Is scalar k in [1,r-1]?                
    if ((FP_ISZERO(k) == TRUE) || (MOD_EVAL(k, TedCurve->prime) == FALSE)) {
        return FALSE;
    }
    // Are (x,y) in [0,p-1]?
    if (MOD_EVAL(P->x, TedCurve->prime) == FALSE || MOD_EVAL(P->y, TedCurve->prime) == FALSE) {  
        return FALSE;
    }
    // Does P lie on the curve?
    FP_SQR(P->y, t3);           
    FP_SQR(P->x, t2);
    FP_SUB(t3, t2, t1);             // -x^2 + y^2 
    FP_MUL(t2, t3, t4);
    FP_MUL(TedCurve->parameter2, t4, t3);
    FP_ZERO(t4);  t4[0] = 1;        // t4 = 1
    FP_ADD(t3, t4, t2);             // 1 + dx^2y^2
    FP_SUB(t1, t2, t1); 
    if (FP_ISZERO(t1) == FALSE) {     
        return FALSE;
    }
    // Eliminate small torsion (assuming co-factor 4)
    ECCCONVERT_AFF_TO_EXTPROJ_TE(P, T);
    ECCDOUBLE_EXT_TE(T, TedCurve);
    ECCDOUBLE_EXT_TE(T, TedCurve);                        // 4*P = T = (X_T:Y_T:Z_T:Ta_T:Tb_T)
    // Is the new point the neutral point (0:1:1)?
    if (ECC_IS_NEUTRAL_EXT_TE(T, TedCurve) == TRUE) {     // Assuming that input point is public so this does not leak info
        return FALSE;
    }
    /****************************/

    ECC_PRECOMP_EXT_TE(T, table, npoints, TedCurve);   // Precomputation of points T[0],...,T[npoints-1]
     
    odd = -((sdig)k[0] & 1);     
    FP_SUB(TedCurve->order, k, k_temp);                // Converting scalar to odd (r-k if even)
    for (j = 0; j < num_words; j++)                    // If (odd) then k = k_temp else k = k 
    {
        k_temp[j] = (odd & (k[j] ^ k_temp[j])) ^ k_temp[j];
    }

    fixed_window_recode(k_temp, TedCurve->rbits, W_VARBASE, digits); 

    LUT_EXT_TE(table, R, digits[t], npoints, TedCurve); 
    FP_SUB(R->XY, R->YX, T->X);                        // Initialize T = (X_T:Y_T:Z_T:Ta_T:Tb_T) with a point from the precomputed table
    FP_DIV2(T->X, T->X);
    FP_ADD(R->YX, T->X, T->Y); 
    FP_DIV2(R->Z2, T->Z); 
    FP_COPY(T->X, T->Ta); 
    FP_COPY(T->Y, T->Tb); 

    for (i = (t-1); i >= 0; i--)
    {
        for (j = 0; j < (W_VARBASE-1); j++)
        {
            ECCDOUBLE_EXT_TE(T, TedCurve);                      // Double (X_T:Y_T:Z_T:Ta_T:Tb_T) = 2(X_T:Y_T:Z_T:Ta_T:Tb_T)
        }
        LUT_EXT_TE(table, R, digits[i], npoints, TedCurve);     // Load R = (XY_R:YX_R:Z2_R:T2_R) with a point from the precomputed table
        ECCUADD_EXT_INTERNAL_TE(R, T, TedCurve);                // Complete addition (X_T:Y_T:Z_T:Ta_T:Tb_T) = (X_T:Y_T:Z_T:Ta_T:Tb_T) + (XY_R:YR_R:Z2_R:T2_R)   
    }

    FP_COPY(T->X, k_temp);
    FP_NEG(TedCurve->prime, k_temp);                            // Correcting scalar (-Tx if even)
    for (j = 0; j < num_words; j++)                             // If (even) then Tx = -Tx 
    {
        T->X[j] = (odd & (T->X[j] ^ k_temp[j])) ^ k_temp[j];
    }
    ECCNORM_TE(T, Q, TedCurve);                                 // Output Q = (x,y)
       
// cleanup
    for (j = 0; j < ((sizeof(FELM)*8)+ML_WORD-1)/ML_WORD; j++) {
        k_temp[j] = 0;
    }
    for (j = 0; j < ((sizeof(FELM)*8)+W_VARBASE-2)/(W_VARBASE-1) + 1; j++) {
        digits[j] = 0;
    }
    ECCZERO_EXT_TE(T);
    ECCZERO_PRECOMP_EXT_TE(R);
    for (j = 0; j < (1 << (W_VARBASE-2)); j++) {
        ECCZERO_PRECOMP_EXT_TE(table[j]);
    }
    FP_ZERO(t1);
    FP_ZERO(t2);
    FP_ZERO(t3);
    FP_ZERO(t4);

    return TRUE;
}


void ECCUMADD_EXT_TE(POINT_PRECOMP_EXTAFF_TE Q, POINT_EXT_TE P, PCurveStruct TedCurve)      
{ // Unified mixed point addition P = P+Q or P = P+P
  // Twisted Edwards a=-1 curve
  // Inputs: P = (X1,Y1,Z1,Ta,Tb) in extended twisted Edwards coordinates, where T1 = Ta*Tb
  //         Q = (x2+y2,y2-x2,2dt2) in extended affine coordinates    
  // Output: P = (X1,Y1,Z1,Ta,Tb) in extended twisted Edwards coordinates, where T1 = Ta*Tb
    FELM t1, t2, t3;   
#ifdef ML_COUNT
    nmul+=7; nadd+=7; 
#endif
    UNREFERENCED_PARAMETER(TedCurve);

    // SECURITY NOTE: this addition function is complete (i.e., it works for any possible inputs, including the cases P!=Q, P=Q, P=-Q and P=neutral) on large prime-order twisted Edwards
    //                subgroups. For more information, refer to "Selecting elliptic curves for cryptography: an efficiency and security analysis", http://eprint.iacr.org/2014/130
                
    FP_MUL(P->Ta, P->Tb, t1);           // t1 = T1
    FP_COPY(t1, P->Ta);                 // Ta = T1
    FP_ADD(P->Z, P->Z, t1);             // t1 = 2z1        
    FP_MUL(P->Ta, Q->t2, t2);           // t2 = 2dT1.T2 
    FP_COPY(t2, P->Ta);   
    FP_SUB(t1, P->Ta, t3);              // t3 = theta
    FP_ADD(t1, P->Ta, t1);              // t1 = alpha
    FP_ADD(P->X, P->Y, t2);             // t2 = (x1+y1) 
    FP_MUL(Q->xy, t2, P->Ta);           // Ta = (x1+y1)(x2+y2)
    FP_SUB(P->Y, P->X, t2);             // t2 = (y1-x1)
    FP_MUL(Q->yx, t2, P->X);            // x1 = (y1-x1)(y2-x2)
    FP_SUB(P->Ta, P->X, P->Tb);         // Tbfinal = beta
    FP_ADD(P->Ta, P->X, P->Ta);         // Tafinal = omega
    FP_MUL(P->Tb, t3, P->X);            // Xfinal = beta.theta
    FP_MUL(t1, t3, P->Z);               // Zfinal = theta. alpha
    FP_MUL(P->Ta, t1, P->Y);            // Yfinal = alpha.omega
    
// cleanup
    FP_ZERO(t1);
    FP_ZERO(t2);
    FP_ZERO(t3);

    return;


}BOOL ECC_MUL_FIXED_TE(POINT_PRECOMP_EXTAFF_TE *P_table, dig *k, POINT_TE Q, MemType memory_use, PCurveStruct TedCurve)
{ // Wrapper for fixed-base scalar multiplication Q = k.P, where P = P_table 
  // Twisted Edwards a=-1 curve
    unsigned int w, v;
    
    if (memory_use < 0 || memory_use >= MemTypeSize) {
        return FALSE;
    }

    if (memory_use == MEM_LARGE) {
        w = W_MEM_LARGE;
        v = V_MEM_LARGE;
    } else if (memory_use == MEM_COMPACT) {
        w = W_MEM_COMPACT_TE;
        v = V_MEM_COMPACT_TE;
    } else {
        w = W_MEM_LARGE;
        v = V_MEM_LARGE;
    }

    return ECC_MUL_FIXED_INTERNAL_TE(P_table, k, Q, w, v, TedCurve);
}



BOOL ECC_MUL_FIXED_INTERNAL_TE(POINT_PRECOMP_EXTAFF_TE *P_table, dig *k, POINT_TE Q, unsigned int w, unsigned int v, PCurveStruct TedCurve)
{ // Fixed-base scalar multiplication Q = k.P, where P = P_table, using the Modified LSB-set Comb method 
  // Twisted Edwards a=-1 curve
    unsigned int j, npoints, d, e, l, num_words = (TedCurve->nbits+ML_WORD-1)/ML_WORD;
    int i, ii, digits[2*MAXBITLENGTH] = {0};                         /////
    dig k_temp[((sizeof(FELM)*8)+ML_WORD-1)/ML_WORD] = {0};
    POINT_EXT_TE T; 
    POINT_PRECOMP_EXTAFF_TE R;
    sdig odd = 0;
    signed long digit = 0;
#ifdef ML_COUNT
    nadd+=4; nlut++;
#endif
        
    // SECURITY NOTE: the crypto sensitive part of this function is protected against timing attacks and runs in constant-time on large prime-order twisted Edwards subgroups. 
    //                Conditional if-statements evaluate public data only and the number of iterations for all loops is public. Refer to "Selecting elliptic curves 
    //                for cryptography: an efficiency and security analysis", http://eprint.iacr.org/2014/130, for the full proof demonstrating exception-less, 
    //                constant-time execution of scalar multiplication.
    // DISCLAIMER:    the protocol designer is responsible for guarantying that early termination produced after detecting errors during input validation
    //                (of scalar k) does not leak any secret information. 

    
    /***** Input validation: *****/                     
    // Is scalar k in [1,r-1]?                
    if ((FP_ISZERO(k) == TRUE) || (MOD_EVAL(k, TedCurve->prime) == FALSE)) {
        return FALSE;
    }
    if (P_table == NULL) {    // Full point validation is done during offline precomputation
        return FALSE;
    }
    /*****************************/

    e = (TedCurve->rbits+w*v-1)/(w*v);  
    if (TedCurve->rbits-e*w*v == 0) {
        return FALSE;                                  // This parameter selection is not allowed   
    }
    d = e*v;     
    l = d*w;                                           // Fixed length of the mLSB-set representation
    npoints = v*(1 << (w-1));

    odd = -((sdig)k[0] & 1);     
    FP_SUB(TedCurve->order, k, k_temp);                // Converting scalar to odd (r-k if even)
    for (j=0; j<num_words; j++)                        // If (odd) then k = k_temp else k = k 
    {
        k_temp[j] = (odd & (k[j] ^ k_temp[j])) ^ k_temp[j];
    }
    mLSB_set_recode(k_temp, TedCurve->rbits, l, d, digits); 

    // Extracting initial digit 
    digit = digits[w*d-1];
    for (i=(int)((w-1)*d-1); i>=(int)(2*d-1); i=i-d)           
    {
        digit = 2*digit + digits[i];
    }
    // Initialize T = (x+y,y-x,2t) with a point from the table
    LUT_EXTAFF_TE(P_table+(v-1)*(1 << (w-1)), R, digit, digits[d-1], 1 << (w-1), TedCurve);

    FP_SUB(R->xy, R->yx, T->X);                     
    FP_DIV2(T->X, T->X);
    FP_ADD(R->yx, T->X, T->Y); 
    FP_ZERO(T->Z); T->Z[0] = 1; 
    FP_COPY(T->X, T->Ta); 
    FP_COPY(T->Y, T->Tb);                            // Initial point T = (X:Y:1:Ta:Tb)

    for (j = 0; j < (v-1); j++)
    {
        digit = digits[w*d-(j+1)*e-1];
        for (i=(int)((w-1)*d-(j+1)*e-1); i>=(int)(2*d-(j+1)*e-1); i=i-d)           
        {
            digit = 2*digit + digits[i];
        }
        LUT_EXTAFF_TE(P_table+(v-j-2)*(1 << (w-1)), R, digit, digits[d-(j+1)*e-1], 1 << (w-1), TedCurve);  // Load R = (xy_R:yx_R:t2_R) with a point from the precomputed table
        ECCUMADD_EXT_TE(R, T, TedCurve);             // Complete mixed addition (X_T:Y_T:Z_T:Ta_T:Tb_T) = (X_T:Y_T:Z_T:Ta_T:Tb_T) + (xy_R:yx_R:t2_R)
#ifdef ML_COUNT
        nlut++;
#endif                      
    }

    for (ii = (e-2); ii >= 0; ii--)
    {
        ECCDOUBLE_EXT_TE(T, TedCurve); 
        for (j = 0; j <= (v-1); j++)
        {
            digit = digits[w*d-j*e+ii-e];

            for (i=(int)((w-1)*d-j*e+ii-e); i>=(int)(2*d-j*e+ii-e); i=i-d)           
            {
                digit = 2*digit + digits[i];
            }
            LUT_EXTAFF_TE(P_table+(v-j-1)*(1 << (w-1)), R, digit, digits[d-j*e+ii-e], 1 << (w-1), TedCurve);  // Load R = (xy_R:yx_R:t2_R) with a point from the precomputed table
            ECCUMADD_EXT_TE(R, T, TedCurve);         // Complete mixed addition (X_T:Y_T:Z_T:Ta_T:Tb_T) = (X_T:Y_T:Z_T:Ta_T:Tb_T) + (xy_R:yx_R:t2_R) 
#ifdef ML_COUNT
            nlut++;
#endif                                       
        }        
    } 

    FP_COPY(T->X, k_temp);
    FP_NEG(TedCurve->prime, k_temp);                             // Correcting scalar (-Tx if even)
    for (j=0; j<num_words; j++)                                  // If (even) then Tx = -Tx 
    {
        T->X[j] = (odd & (T->X[j] ^ k_temp[j])) ^ k_temp[j];
    }
    ECCNORM_TE(T, Q, TedCurve);                                  // Output Q = (x,y)
    
// cleanup
    for (j = 0; j < ((sizeof(FELM)*8)+ML_WORD-1)/ML_WORD; j++) {
        k_temp[j] = 0;
    }
    for (j = 0; j < (2*MAXBITLENGTH); j++) {
        digits[j] = 0;
    }
    ECCZERO_EXT_TE(T);
    ECCZERO_PRECOMP_EXTAFF_TE(R);    
    
    return TRUE;
}


POINT_PRECOMP_EXTAFF_TE* ECC_PRECOMP_FIXED_TE(POINT_TE P, MemType memory_use, PCurveStruct TedCurve)
{ // Wrapper for precomputation scheme using extended affine coordinates (x+y,y-x,2dt) for fixed-base scalar multiplication
  // Twisted Edwards a=-1 curve
    unsigned int w, v, d, e;
    
    if (memory_use < 0 || memory_use >= MemTypeSize) {
        return FALSE;
    }

    if (memory_use == MEM_LARGE) { 
        w = W_MEM_LARGE;
        v = V_MEM_LARGE;
    } else if (memory_use == MEM_COMPACT) {
        w = W_MEM_COMPACT_TE;
        v = V_MEM_COMPACT_TE;
    } else {
        w = W_MEM_LARGE;
        v = V_MEM_LARGE;
    }
    e = (TedCurve->rbits+w*v-1)/(w*v);    
    if (TedCurve->rbits-e*w*v == 0)     // This parameter selection is not allowed
        return FALSE;       
    d = e*v;                             

    return ECC_PRECOMP_FIXED_INTERNAL_TE(P, w, v, d, e, TedCurve);
}


POINT_PRECOMP_EXTAFF_TE* ECC_PRECOMP_FIXED_INTERNAL_TE(POINT_TE P, unsigned int w, unsigned int v, unsigned int d, unsigned int e, PCurveStruct TedCurve)
{ // Precomputation scheme using extended affine coordinates (x+y,y-x,2dt) for fixed-base scalar multiplication. Function returns NULL on error.
  // Twisted Edwards a=-1 curve
    POINT_TE A;       
	POINT_EXT_TE R, B, base[WMAX]; 
    POINT_PRECOMP_EXT_TE baseb[WMAX], RR;
    POINT_PRECOMP_EXTAFF_TE *T;
    unsigned int i, j, k, npoints, index;
    unsigned long index_group;
    FELM t1, t2, t3, t4;
                
    // SECURITY NOTE: precomputation for fixed-base scalar multiplication uses public inputs. 
    // DISCLAIMER:    the protocol designer is responsible for guarantying that early termination produced after detecting errors during input validation
    //                (of base point P) does not leak any secret information. 

    /***** Input validation and elimination of small torsion: *****/
    // Check if P is the neutral point (0,1)
    if (ECC_IS_NEUTRAL_AFF_TE(P, TedCurve) == TRUE) {                             
        return NULL;   
    } 
    // Are (x,y) in [0,p-1]?
    if (MOD_EVAL(P->x, TedCurve->prime) == FALSE || MOD_EVAL(P->y, TedCurve->prime) == FALSE) {  
        return NULL;
    }
    // Does P lie on the curve?
    FP_SQR(P->y, t3);           
    FP_SQR(P->x, t2);
    FP_SUB(t3, t2, t1);             // -x^2 + y^2 
    FP_MUL(t2, t3, t4);
    FP_MUL(TedCurve->parameter2, t4, t3);
    FP_ZERO(t4); t4[0] = 1;         // t4 = 1 
    FP_ADD(t3, t4, t2);             // 1 + dx^2y^2
    FP_SUB(t1, t2, t1); 
    if (FP_ISZERO(t1)==FALSE) {     
        T = NULL;
        goto cleanup;
    }
    // Eliminate small torsion (assuming co-factor 4)
    ECCCONVERT_AFF_TO_EXTPROJ_TE(P, B);
    ECCDOUBLE_EXT_TE(B, TedCurve);
    ECCDOUBLE_EXT_TE(B, TedCurve);                      // 4*P = T = (X_T:Y_T:Z_T:Ta_T:Tb_T)
    // Is the new point the neutral point (0:1:1)?
    if (ECC_IS_NEUTRAL_EXT_TE(B, TedCurve) == TRUE) {   // Assuming that input point is public so this does not leak info{     
        T = NULL;
        goto cleanup;
    }
    /*************************/
        
    if (TedCurve->rbits-e*w*v == 0) {
        return FALSE;                                   // This parameter selection is not allowed
    }

    npoints = v*(1 << (w-1));
    T = (POINT_PRECOMP_EXTAFF_TE*)calloc(npoints, sizeof(POINT_PRECOMP_EXTAFF_TE));   // Allocating space for table
    if (T == NULL) {
        goto cleanup;
    }
    
    ECCNORM_TE(B, A, TedCurve);
    ECCCONVERT_AFF_TO_EXTPROJ_TE(A, base[0]);           // base[0] = A = (X:Y:1:Ta:Tb)
    
    FP_COPY(TedCurve->parameter2, t1); 

    // Compute base point for each w (or row)
	for (i = 0; i < (w-1); i++) {
		ECCCOPY_EXT_TE(base[i], base[i+1]);
        FP_ADD(base[i]->X, base[i]->Y, baseb[i]->XY);
        FP_SUB(base[i]->Y, base[i]->X, baseb[i]->YX);
        FP_ADD(base[i]->Z, base[i]->Z, baseb[i]->Z2);
        FP_MUL(base[i]->Ta, t1, t2);
        FP_MUL(base[i]->Tb, t2, baseb[i]->T2);     
        FP_ADD(baseb[i]->T2, baseb[i]->T2, baseb[i]->T2);                  // baseb in coordinates (x+y,y-x,2z,2dt)
		for (j = 0; j < d; j++) ECCDOUBLE_EXT_TE(base[i+1], TedCurve);     // base[i+1] = 2^d base[i]
	}

    FP_ADD(base[w-1]->X, base[w-1]->Y, baseb[w-1]->XY);
    FP_SUB(base[w-1]->Y, base[w-1]->X, baseb[w-1]->YX);
    FP_ADD(base[w-1]->Z, base[w-1]->Z, baseb[w-1]->Z2);
    FP_MUL(base[w-1]->Ta, t1, t2);
    FP_MUL(base[w-1]->Tb, t2, baseb[w-1]->T2);
    FP_ADD(baseb[w-1]->T2, baseb[w-1]->T2, baseb[w-1]->T2);  // baseb in (x+y,y-x,2z,2dt)
        
    FP_COPY(A->x, T[0]->xy);                                 // T[0] = A in (x,y)
    FP_COPY(A->y, T[0]->yx); 
    
    // Compute precomputed points for the first table
    index = 0;
    index_group = 1;
    for (i = 0; i < (w-1); i++)                              // T[index] = (1 + u_0.2^d + ... + u_{w-2}.2^((w-1)d)) B
    {
        for (j = 0; j < index_group; j++)
        {
            FP_ADD(T[j]->xy, T[j]->yx, RR->XY); 
            FP_SUB(T[j]->yx, T[j]->xy, RR->YX); 
            FP_ZERO(RR->Z2); RR->Z2[0] = 1;                           
            FP_MUL(T[j]->xy, T[j]->yx, RR->T2); 
            ECCADD_PRECOMP_TE(RR, baseb[i+1], RR, TedCurve);
            index++;
            FP_SUB(RR->XY, RR->YX, R->X);   
            FP_DIV2(R->X, R->X);
            FP_ADD(RR->YX, R->X, R->Y); 
            FP_DIV2(RR->Z2, R->Z);                           // R in (x,y,z)  
            ECCNORM_TE(R, A, TedCurve);   
            ECCCONVERT_AFF_TO_EXTPROJ_TE(A, R);             
            FP_COPY(R->X, T[index]->xy);                             
            FP_COPY(R->Y, T[index]->yx);                     // T[] in (x,y)
        }
        index_group = 2*index_group;
    }
                
    // Compute precomputed points for the remaining tables
    index++;
    for (i = 0; i < (v-1); i++)                              // T[index] = 2^(ev) (1 + u_0.2^d + ... + u_{w-2}.2^((w-1)d)) B
    {
        for (j = 0; j < index; j++)
        {
            FP_COPY(T[i*index + j]->xy, R->X);                             
            FP_COPY(T[i*index + j]->yx, R->Y); 
            FP_ZERO(R->Z); R->Z[0] = 1; 
		    for (k = 0; k < e; k++) ECCDOUBLE_EXT_TE(R, TedCurve);     // 2^(ev) * X * B
            ECCNORM_TE(R, A, TedCurve);   
            ECCCONVERT_AFF_TO_EXTPROJ_TE(A, R);          
            FP_COPY(R->X, T[(i+1)*index + j]->xy);                             
            FP_COPY(R->Y, T[(i+1)*index + j]->yx);                     // Precomputed points T[] in (x,y)
        }
    }

    for (i = 0; i < npoints; i++)
    {
        FP_MUL(T[i]->xy, T[i]->yx, R->X);
        FP_MUL(TedCurve->parameter2, R->X, T[i]->t2);
        FP_ADD(T[i]->t2, T[i]->t2, T[i]->t2);
        FP_ADD(T[i]->xy, T[i]->yx, R->X);
        FP_SUB(T[i]->yx, T[i]->xy, T[i]->yx);
        FP_COPY(R->X, T[i]->xy);                             // Precomputed points T[] in coordinates (x+y,y-x,2dt)
    }
    
cleanup:
    ECCZERO_TE(A);                
    ECCZERO_EXT_TE(R);
    ECCZERO_EXT_TE(B);
    ECCZERO_PRECOMP_EXT_TE(RR);
    for (j = 0; j < WMAX; j++) {
        ECCZERO_EXT_TE(base[j]);
        ECCZERO_PRECOMP_EXT_TE(baseb[j]);
    }
    FP_ZERO(t1);
    FP_ZERO(t2);
    FP_ZERO(t3);   
    FP_ZERO(t4);    

    return T;
}


BOOL ECC_DBLMUL_TE(POINT_PRECOMP_EXTAFF_TE *P_table, dig *k, POINT_TE Q, dig *l, POINT_TE R, MemType memory_use, PCurveStruct TedCurve)
{ // Wrapper for double-base scalar multiplication R = k.P + l.Q, where P = P_table
  // P is a fixed-base and Q is a variable-base
  // Twisted Edwards a=-1 curve
    unsigned int w_P;
    
    if (memory_use < 0 || memory_use >= MemTypeSize) {
        return FALSE;
    }

    if (memory_use == MEM_LARGE) {
        w_P = W_P_MEM_LARGE;
    } else if (memory_use == MEM_COMPACT) {
        w_P = W_P_MEM_COMPACT;
    } else {
        w_P = W_P_MEM_LARGE;
    }

    return ECC_DBLMUL_INTERNAL_TE(P_table, k, Q, l, R, w_P, TedCurve);
}


BOOL ECC_DBLMUL_INTERNAL_TE(POINT_PRECOMP_EXTAFF_TE *P_table, dig *k, POINT_TE Q, dig *l, POINT_TE R, unsigned int w_P, PCurveStruct TedCurve)
{ // Double-base scalar multiplication R = k.P + l.Q using wNAF with Interleaving
  // P is a fixed-base and Q is a variable-base
  // Twisted Edwards a=-1 curve
    unsigned int npoints, position;
    int i, digits_P[sizeof(FELM)*8 + 1]={0}, digits_Q[sizeof(FELM)*8 + 1]={0};
    POINT_EXT_TE T; 
    POINT_PRECOMP_EXT_TE table[1 << (W_VARBASE-2)], S;
    POINT_PRECOMP_EXTAFF_TE SS;
    FELM t1, t2, t3, t4;    
#ifdef ML_COUNT
    nmul+=2; nsqr+=2; nadd+=5; 
#endif 
            
    // SECURITY NOTE: this function is intended for a non-constant-time operation such as signature verification. 
    
    /***** Input validation and small torsion elimination: *****/ 
    if (P_table == NULL) {    // Full point validation for P is done during offline precomputation
        return FALSE;
    } 
    // Check if Q is the neutral point (0,1)
    if (ECC_IS_NEUTRAL_AFF_TE(Q, TedCurve) == TRUE) {               
        return FALSE;   
    } 
    // Are scalars k, l in [1,r-1]?                
    if ((FP_ISZERO(k) == TRUE) || (MOD_EVAL(k, TedCurve->prime) == FALSE) || (FP_ISZERO(l) == TRUE) || (MOD_EVAL(l, TedCurve->prime) == FALSE)) {
        return FALSE;
    }
    // Are Q: (x,y) in [0,p-1]?
    if (MOD_EVAL(Q->x, TedCurve->prime) == FALSE || MOD_EVAL(Q->y, TedCurve->prime) == FALSE) {  
        return FALSE;
    }  
    // Does Q lie on the curve?
    FP_SQR(Q->y, t3);           
    FP_SQR(Q->x, t2);
    FP_SUB(t3, t2, t1);              // -x^2 + y^2 
    FP_MUL(t2, t3, t4);
    FP_MUL(TedCurve->parameter2, t4, t3);
    FP_ZERO(t4); t4[0] = 1;          // t4 = 1
    FP_ADD(t3, t4, t2);              // 1 + dx^2y^2
    FP_SUB(t1, t2, t1); 
    if (FP_ISZERO(t1)==FALSE) {     
        return FALSE;
    }
    // Eliminate small torsion (assuming co-factor 4)
    ECCCONVERT_AFF_TO_EXTPROJ_TE(Q, T);
    ECCDOUBLE_EXT_TE(T, TedCurve);
    ECCDOUBLE_EXT_TE(T, TedCurve);                      // 4*Q = T = (X_T:Y_T:Z_T:Ta_T:Tb_T)
    // Is the new point the neutral point (0:1:1)?
    if (ECC_IS_NEUTRAL_EXT_TE(T, TedCurve) == TRUE) {   // Assuming that input point is public so this does not leak info
        return FALSE;
    }
    /*****************************/

    npoints = 1 << (W_VARBASE-2); 

    ECC_PRECOMP_EXT_TE(T, table, npoints, TedCurve);    // Precomputation of points table[0],...,table[npoints-1]
    wNAF_recode(k, TedCurve->rbits, w_P, digits_P);     // Recode k and l to the wNAF representation
    wNAF_recode(l, TedCurve->rbits, W_VARBASE, digits_Q);
    FP_ZERO(T->X); FP_ZERO(T->Y); T->Y[0] = 1; FP_ZERO(T->Z); T->Z[0] = 1;  // Initialize T as the neutral point (0:1:1)   

    for (i = TedCurve->rbits; i >= 0; i--)
    {
        ECCDOUBLE_EXT_TE(T, TedCurve);                 // Double (X_T:Y_T:Z_T:Ta_T:Tb_T) = 2(X_T:Y_T:Z_T:Ta_T:Tb_T)
        if (digits_Q[i] < 0) {
            position = (-digits_Q[i])/2;               // Load S = (XY_S:YX_S:Z2_S:T2_S) = (X+Y,Y-X,2Z,2dT) from a point in the precomputed table
            FP_COPY(table[position]->XY, S->YX);
            FP_COPY(table[position]->YX, S->XY);
            FP_COPY(table[position]->Z2, S->Z2);
            FP_COPY(table[position]->T2, S->T2);
            FP_NEG(TedCurve->prime, S->T2);            // Negate S 
            ECCUADD_EXT_INTERNAL_TE(S, T, TedCurve);   // Complete addition (X_T:Y_T:Z_T:Ta_T:Tb_T) = (X_T:Y_T:Z_T:Ta_T:Tb_T) + (XY_S:YX_S:Z2_S:T2_S) 
#ifdef ML_COUNT
            nadd++; 
#endif                        
        } else if (digits_Q[i] > 0) {            
            position = (digits_Q[i])/2;                // Load S = (XY_S:YX_S:Z2_S:T2_S) = (X+Y,Y-X,2Z,2dT) from a point in the precomputed table
            FP_COPY(table[position]->XY, S->XY);
            FP_COPY(table[position]->YX, S->YX);
            FP_COPY(table[position]->Z2, S->Z2);
            FP_COPY(table[position]->T2, S->T2);
            ECCUADD_EXT_INTERNAL_TE(S, T, TedCurve);   // Complete addition (X_T:Y_T:Z_T:Ta_T:Tb_T) = (X_T:Y_T:Z_T:Ta_T:Tb_T) + (XY_S:YX_S:Z2_S:T2_S) 
        }

        if (digits_P[i] < 0) {                           
            position = (-digits_P[i])/2;               // Load SS = (xy_SS:yx_SS:t2_SS) = (x+y,y-x,2dt) from a point in the precomputed table
            FP_COPY(P_table[position]->xy, SS->yx);
            FP_COPY(P_table[position]->yx, SS->xy);
            FP_COPY(P_table[position]->t2, SS->t2);
            FP_NEG(TedCurve->prime, SS->t2);           // Negate SS
            ECCUMADD_EXT_TE(SS, T, TedCurve);          // Complete mixed addition (X_T:Y_T:Z_T:Ta_T:Tb_T) = (X_T:Y_T:Z_T:Ta_T:Tb_T) + (xy_SS:yx_SS:t2_SS)
#ifdef ML_COUNT
            nadd++; 
#endif                        
        } else if (digits_P[i] > 0) { 
            position = (digits_P[i])/2;                // Load SS = (xy_SS:yx_SS:t2_SS) = (x+y,y-x,2dt) from a point in the precomputed table
            FP_COPY(P_table[position]->xy, SS->xy);
            FP_COPY(P_table[position]->yx, SS->yx);
            FP_COPY(P_table[position]->t2, SS->t2);
            ECCUMADD_EXT_TE(SS, T, TedCurve);          // Complete mixed addition (X_T:Y_T:Z_T:Ta_T:Tb_T) = (X_T:Y_T:Z_T:Ta_T:Tb_T) + (xy_SS:yx_SS:t2_SS)
        }
    }
    ECCNORM_TE(T, R, TedCurve);                        // Output R = (x,y)
    
// cleanup
    for (i = 0; i < (sizeof(FELM)*8 + 1); i++) {
        digits_P[i] = 0;
        digits_Q[i] = 0;
    }
    ECCZERO_EXT_TE(T);
    for (i = 0; i < (1 << (W_VARBASE-2)); i++) {
        ECCZERO_PRECOMP_EXT_TE(table[i]);
    }
    ECCZERO_PRECOMP_EXT_TE(S);
    ECCZERO_PRECOMP_EXTAFF_TE(SS);
    FP_ZERO(t1);
    FP_ZERO(t2);
    FP_ZERO(t3);
    FP_ZERO(t4);
    
    return TRUE;
}


POINT_PRECOMP_EXTAFF_TE* ECC_PRECOMP_DBLMUL_TE(POINT_TE P, MemType memory_use, PCurveStruct TedCurve)
{ // Wrapper for precomputation scheme using extended affine coordinates (x+y,y-x,2dt) for the fixed-base in double-scalar multiplication
  // Twisted Edwards a=-1 curve
    unsigned int w_P;
    
    if (memory_use < 0 || memory_use >= MemTypeSize) {
        return FALSE;
    }

    if (memory_use == MEM_LARGE) { 
        w_P = W_P_MEM_LARGE;
    } else if (memory_use == MEM_COMPACT) {
        w_P = W_P_MEM_COMPACT;
    } else {
        w_P = W_P_MEM_LARGE;
    }                          

    return ECC_PRECOMP_DBLMUL_INTERNAL_TE(P, w_P, TedCurve);
}



POINT_PRECOMP_EXTAFF_TE* ECC_PRECOMP_DBLMUL_INTERNAL_TE(POINT_TE P, unsigned int w_P, PCurveStruct TedCurve)
{ // Precomputation scheme using extended affine coordinates (x+y,y-x,2dt) for the fixed base in double-scalar multiplication. Function returns NULL on error.
  // Twisted Edwards a=-1 curve
    POINT_PRECOMP_EXTAFF_TE *T = NULL;
    POINT_PRECOMP_EXT_TE *T_EXT = NULL;
    POINT_EXT_TE B;
    unsigned int i, npoints; 
    FELM t1, t2, t3, t4; 
                
    // SECURITY NOTE: precomputation for double-scalar multiplication uses public inputs. 

    /***** Input validation and elimination of small torsion: *****/
    // Check if P is the neutral point (0,1)
    if (ECC_IS_NEUTRAL_AFF_TE(P, TedCurve) == TRUE) {           
        return NULL;   
    } 
    // Are (x,y) in [0,p-1]?
    if (MOD_EVAL(P->x, TedCurve->prime) == FALSE || MOD_EVAL(P->y, TedCurve->prime) == FALSE) {  
        return NULL;
    }
    // Does P lie on the curve?
    FP_SQR(P->y, t3);           
    FP_SQR(P->x, t2);
    FP_SUB(t3, t2, t1);             // -x^2 + y^2 
    FP_MUL(t2, t3, t4);
    FP_MUL(TedCurve->parameter2, t4, t3);
    FP_ZERO(t4); t4[0] = 1;         // t4 = 1 
    FP_ADD(t3, t4, t2);             // 1 + dx^2y^2
    FP_SUB(t1, t2, t1); 
    if (FP_ISZERO(t1)==FALSE) {
        T = NULL;
        goto cleanup;
    }
    // Eliminate small torsion (assuming co-factor 4)
    ECCCONVERT_AFF_TO_EXTPROJ_TE(P, B);
    ECCDOUBLE_EXT_TE(B, TedCurve);
    ECCDOUBLE_EXT_TE(B, TedCurve);                      // 4*P = T = (X_T:Y_T:Z_T:Ta_T:Tb_T)
    // Is the new point the neutral point (0:1:1)?
    if (ECC_IS_NEUTRAL_EXT_TE(B, TedCurve) == TRUE) {   // Assuming that input point is public so this does not leak info
        T = NULL;
        goto cleanup;
    }
    /*************************/

    npoints = (1 << (w_P-2));
    T_EXT = (POINT_PRECOMP_EXT_TE*)calloc(npoints, sizeof(POINT_PRECOMP_EXT_TE));    // Allocating space for temporary table
    if (T_EXT == NULL) {
        T = NULL;
        goto cleanup;
    }
    T = (POINT_PRECOMP_EXTAFF_TE*)calloc(npoints, sizeof(POINT_PRECOMP_EXTAFF_TE));  // Allocating space for precomputed table
    if (T == NULL) {
        goto cleanup;
    }

    ECC_PRECOMP_EXT_TE(B, T_EXT, npoints, TedCurve);         
    for (i=0; i<npoints; i++)
    {
        FP_DIV2(T_EXT[i]->Z2, t1);
        FP_INV(t1);
        FP_MUL(T_EXT[i]->XY, t1, T[i]->xy);
        FP_MUL(T_EXT[i]->YX, t1, T[i]->yx);
        FP_MUL(T_EXT[i]->T2, t1, T[i]->t2);
    }
    
cleanup:
    if (T_EXT != NULL) {
        free(T_EXT);
    }
    ECCZERO_EXT_TE(B);
    FP_ZERO(t1);
    FP_ZERO(t2);
    FP_ZERO(t3);
    FP_ZERO(t4);

    return T;
}


BOOL ECC_DESTROY_PRECOMP_TE(POINT_PRECOMP_EXTAFF_TE* T_fixed)
{ // Frees memory occupied by a precomputation table "T_fixed" used during fixed-base or double-scalar multiplication.
  // This function must be called once done using a table generated by ecc_precomp_fixed_Tedxxx or ecc_precomp_dblmul_Tedxxx. 
  // Twisted Edwards a=-1 curve
        
    if (T_fixed != NULL) {
        free(T_fixed); 
    }

    return TRUE;
}