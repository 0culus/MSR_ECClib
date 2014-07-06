/**************************************************************************
* Suite for benchmarking/testing curve operations for MSR ECClib
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
* Abstract: benchmarking/testing curve operations
*
* This software is based on the article by Joppe Bos, Craig Costello, 
* Patrick Longa and Michael Naehrig, "Selecting elliptic curves for
* cryptography: an efficiency and security analysis", preprint available
* at http://eprint.iacr.org/2014/130.
***************************************************************************/  

#include <stdio.h>
#include <windows.h>
#include "msr_ecclib.h"
#include "tests.h"


#ifdef ECCURVES_256

BOOL ecc_test256_w(PCurveStruct JacCurve)
{ // Tests for curve Jac256
    BOOL OK = TRUE;
    dig n;
    BOOL passed;
    point_jac_Jac256 P, R;
    dig k[ML_WORDS256], kk[ML_WORDS256];
    point_Jac256 A, AA, B, BB, *T_fixed = NULL;
    
    printf("\n\nTESTING \n"); 
    printf("--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Curve arithmetic: Weierstrass a=-3 over GF(2^256-189) \n\n"); 

    // Point doubling (Weierstrass a=-3)
    passed = TRUE;
    eccset_Jac256(A, JacCurve); 
    eccconvert_aff_to_jac_Jac256(A, P);

    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        eccdouble_jac_Jac256(P, JacCurve);          // 2*P
        eccdouble_waff_256(A, JacCurve);            // 2*A
    }
    eccnorm_Jac256(P, AA, JacCurve);

    if (fpcompare256(A->x,AA->x)!=0 || fpcompare256(A->y,AA->y)!=0) passed=FALSE;
    if (passed==TRUE) printf("  Point doubling tests .................................................................... PASSED");
    else { printf("  Point doubling tests ... FAILED"); printf("\n"); return FALSE; }
    printf("\n");
    
    // "Complete" point addition (Weierstrass a=-3)
    passed = TRUE;
    eccset_Jac256(A, JacCurve);
    eccconvert_aff_to_jac_Jac256(A, P);
    ecccopy_jac_Jac256(P, R); 
    eccdouble_jac_Jac256(P, JacCurve);                  // P = 2P 
    eccset_Jac256(AA, JacCurve);
    eccdouble_waff_256(AA, JacCurve);                   // AA = 2A

    for (n=0; n<1; n++)
    {
        eccadd_jac_Jac256(R, P, JacCurve);              // P = P+Q
        eccadd_waff_256(A, AA, JacCurve);               // AA = AA+A
    }    
    eccnorm_Jac256(P, A, JacCurve);
    
    if (fpcompare256(A->x,AA->x)!=0 || fpcompare256(A->y,AA->y)!=0) passed=FALSE;
    if (passed==TRUE) printf("  (Complete) point addition tests ......................................................... PASSED");
    else { printf("  (Complete) point addition tests ... FAILED"); printf("\n"); return FALSE; }
    printf("\n");
    
    // Variable-base scalar multiplication (Weierstrass a=-3)
    eccset_Jac256(A, JacCurve); 
    
    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        random256(k, JacCurve->order); 

        ecc_mul_waff_256(A, k, AA, JacCurve);    
        ecc_scalar_mul_Jac256(A, k, B, JacCurve);
        
        if (fpcompare256(AA->x,B->x)!=0 || fpcompare256(AA->y,B->y)!=0) { passed=FALSE; break; }
    }    

    if (passed==TRUE) printf("  Variable-base scalar multiplication tests ............................................... PASSED");
    else { printf("  Variable-base scalar multiplication tests ... FAILED"); printf("\n"); return FALSE; }
    printf("\n");
        
    // Fixed-base scalar multiplication (Weierstrass a=-3)
    eccset_Jac256(A, JacCurve); 
    
    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        random256(k, JacCurve->order); 
        ecc_mul_waff_256(A, k, AA, JacCurve);  

        T_fixed = ecc_precomp_fixed_Jac256(A, MEM_LARGE, JacCurve);
        ecc_scalar_mul_fixed_Jac256(T_fixed, k, B, MEM_LARGE, JacCurve);
        ecc_destroy_precomp_Jac256(T_fixed); 
        
        if (fpcompare256(AA->x,B->x)!=0 || fpcompare256(AA->y,B->y)!=0) { passed=FALSE; break; }
    }    

    if (passed==TRUE) printf("  Fixed-base scalar multiplication tests .................................................. PASSED");
    else { printf("  Fixed-base scalar multiplication tests ... FAILED"); printf("\n"); return FALSE; }
    printf("\n");
    
    // Double-scalar multiplication (Weierstrass a=-3)
    passed = TRUE;
    eccset_Jac256(A, JacCurve); 
    random256(k, JacCurve->order); 
    ecc_mul_waff_256(A, k, B, JacCurve);    // Base points are A and B
    
    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        random256(k, JacCurve->order); 
        random256(kk, JacCurve->order); 

        ecc_mul_waff_256(A, k, AA, JacCurve);    
        ecc_mul_waff_256(B, kk, BB, JacCurve);
        eccadd_waff_256(BB, AA, JacCurve);

        T_fixed = ecc_precomp_dblmul_Jac256(A, MEM_LARGE, JacCurve);
        ecc_double_scalar_mul_Jac256(T_fixed, k, B, kk, BB, MEM_LARGE, JacCurve); 
        ecc_destroy_precomp_Jac256(T_fixed); 
        
        if (fpcompare256(AA->x,BB->x)!=0 || fpcompare256(AA->y,BB->y)!=0) { passed=FALSE; break; }
    }      

    if (passed==TRUE) printf("  Double-scalar multiplication tests ...................................................... PASSED");
    else { printf("  Double-scalar multiplication tests ... FAILED"); printf("\n"); return FALSE; }
    printf("\n");
    
    return OK;
}


BOOL ecc_test256_te(PCurveStruct TedCurve)
{ // Tests for curve Ted256
    BOOL OK = TRUE;
    dig n;
    BOOL passed;
    point_Ted256 A, AA, U, V;
    point_extproj_Ted256 P, R;
    point_Jac256 PP, QQ, RR, UU, VV;
    CurveStruct WeierstrassCurve = {0};
    dig256 d, t1, t2, t3;
    dig k[ML_WORDS256], kk[ML_WORDS256];
    point_extaff_precomp_Ted256 *T_fixed= NULL;
    
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Curve arithmetic: twisted Edwards a=-1 over GF(2^256-189) \n\n"); 

    // Set Weierstrass curve isomorphic to TedCurve
    fpcopy256(TedCurve->parameter2, d);
    fpzero256(t3); t3[0] = 14;
    fpmul256(t3, d, t1);                               // t1 = 14d
    fpsqr256(d, t2);                                   // t2 = d^2
    fpsub256(t1, t2, t2);                              // t2 = 14d-d^2
    fpzero256(t1); t1[0] = 1;    
    fpsub256(t2, t1, t2);                              // t2 = 14d-d^2-1
    fpzero256(t1); t1[0] = 48;
    fpinv256(t1);                                      // t1 = 1/48
    fpmul256(t1, t2, WeierstrassCurve.parameter1);     // a = (14d-d^2-1)/48
    WeierstrassCurve.nbits = TedCurve->nbits;
    WeierstrassCurve.rbits = TedCurve->rbits;
    WeierstrassCurve.pbits = TedCurve->pbits;

    // Point doubling (twisted Edwards a=-1)
    passed = TRUE;
    eccset_Ted256(A, TedCurve); 
    eccconvert_aff_to_extproj_Ted256(A, P);
    ecc_Ted256_to_weierstrass(A, PP, TedCurve);

    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        eccdouble_extproj_Ted256(P, TedCurve);        // 2*P
        eccdouble_waff_256(PP, &WeierstrassCurve);    // 2*PP
    }
    eccnorm_Ted256(P, A, TedCurve);
    ecc_weierstrass_to_Ted256(PP, AA, TedCurve);

    if (fpcompare256(A->x,AA->x)!=0 || fpcompare256(A->y,AA->y)!=0) passed=FALSE;
    if (passed==TRUE) printf("  Point doubling tests .................................................................... PASSED");
    else { printf("  Point doubling tests ... FAILED"); printf("\n"); return FALSE; }
    printf("\n");
    
    // Point addition (twisted Edwards a=-1)
    passed = TRUE;
    eccset_Ted256(A, TedCurve);
    eccconvert_aff_to_extproj_Ted256(A, P);
    ecccopy_extproj_Ted256(P, R);
    eccdouble_extproj_Ted256(P, TedCurve);             // P = 2P    
        
    eccset_Ted256(A, TedCurve); 
    ecc_Ted256_to_weierstrass(A, QQ, TedCurve);
    ecc_Ted256_to_weierstrass(A, PP, TedCurve);
    eccdouble_waff_256(PP, &WeierstrassCurve);         // PP = 2P

    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        eccadd_extproj_Ted256(R, P, TedCurve);         // P = P+Q
        eccadd_waff_256(QQ, PP, &WeierstrassCurve);    // PP = PP+QQ
    }    
    eccnorm_Ted256(P, A, TedCurve);
    ecc_weierstrass_to_Ted256(PP, AA, TedCurve);
    
    if (fpcompare256(A->x,AA->x)!=0 || fpcompare256(A->y,AA->y)!=0) passed=FALSE;
    if (passed==TRUE) printf("  Point addition tests .................................................................... PASSED");
    else { printf("  Point addition tests ... FAILED"); printf("\n"); return FALSE; }
    printf("\n");
    
    // Variable-base scalar multiplication (twisted Edwards a=-1)
    eccset_Ted256(A, TedCurve); 
    ecc_Ted256_to_weierstrass(A, PP, TedCurve);
    eccdouble_waff_256(PP, &WeierstrassCurve);
    eccdouble_waff_256(PP, &WeierstrassCurve);
        
    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        random256(k, TedCurve->order); 
        
        ecc_mul_waff_256(PP, k, RR, &WeierstrassCurve);    
        ecc_scalar_mul_Ted256(A, k, AA, TedCurve);
        ecc_Ted256_to_weierstrass(AA, QQ, TedCurve);
        
        if (fpcompare256(QQ->x,RR->x)!=0 || fpcompare256(QQ->y,RR->y)!=0) { passed=FALSE; break; }
    }    

    if (passed==TRUE) printf("  Variable-base scalar multiplication tests ............................................... PASSED");
    else { printf("  Variable-base scalar multiplication tests ... FAILED"); printf("\n"); return FALSE; }
    printf("\n");

    // Fixed-base scalar multiplication (twisted Edwards a=-1)
    eccset_Ted256(A, TedCurve); 
    ecc_Ted256_to_weierstrass(A, PP, TedCurve);
    eccdouble_waff_256(PP, &WeierstrassCurve);
    eccdouble_waff_256(PP, &WeierstrassCurve);
    
    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        random256(k, TedCurve->order); 
        ecc_mul_waff_256(PP, k, RR, &WeierstrassCurve);   

        T_fixed = ecc_precomp_fixed_Ted256(A, MEM_LARGE, TedCurve);
        ecc_scalar_mul_fixed_Ted256(T_fixed, k, AA, MEM_LARGE, TedCurve);
        ecc_Ted256_to_weierstrass(AA, QQ, TedCurve);
        ecc_destroy_precomp_Ted256(T_fixed); 
        
        if (fpcompare256(QQ->x,RR->x)!=0 || fpcompare256(QQ->y,RR->y)!=0) { passed=FALSE; break; }
    }    

    if (passed==TRUE) printf("  Fixed-base scalar multiplication tests .................................................. PASSED");
    else { printf("  Fixed-base scalar multiplication tests ... FAILED"); printf("\n"); return FALSE; }
    printf("\n");
    
    // Double-scalar multiplication (twisted Edwards a=-1)
    passed = TRUE;
    eccset_Ted256(A, TedCurve); 
    ecc_Ted256_to_weierstrass(A, PP, TedCurve);
    eccdouble_waff_256(PP, &WeierstrassCurve);
    eccdouble_waff_256(PP, &WeierstrassCurve);
    random256(k, TedCurve->order); 
    ecc_mul_waff_256(PP, k, RR, &WeierstrassCurve);    
    ecc_weierstrass_to_Ted256(RR, AA, TedCurve);    // Base points are (A, AA) in twisted Edwards
    eccdouble_waff_256(RR, &WeierstrassCurve);
    eccdouble_waff_256(RR, &WeierstrassCurve);      // Base points are 4*(PP, RR) in Weierstrass
    
    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        random256(k, TedCurve->order); 
        random256(kk, TedCurve->order); 

        ecc_mul_waff_256(PP, k, UU, &WeierstrassCurve);    
        ecc_mul_waff_256(RR, kk, VV, &WeierstrassCurve);
        eccadd_waff_256(VV, UU, &WeierstrassCurve);  
        ecc_weierstrass_to_Ted256(UU, U, TedCurve); 

        T_fixed = ecc_precomp_dblmul_Ted256(A, MEM_LARGE, TedCurve);
        ecc_double_scalar_mul_Ted256(T_fixed, k, AA, kk, V, MEM_LARGE, TedCurve); 
        ecc_destroy_precomp_Ted256(T_fixed); 
        
        if (fpcompare256(U->x,V->x)!=0 || fpcompare256(U->y,V->y)!=0) { passed=FALSE; break; }
    }      

    if (passed==TRUE) printf("  Double-scalar multiplication tests ...................................................... PASSED");
    else { printf("  Double-scalar multiplication tests ... FAILED"); printf("\n"); return FALSE; }
    printf("\n");

    return OK;
}

#endif


#ifdef ECCURVES_384

BOOL ecc_test384_w(PCurveStruct JacCurve)
{ // Tests for curve Jac384
    BOOL OK = TRUE;
    dig n;
    BOOL passed;
    point_jac_Jac384 P, R;
    dig384 k = {0}, kk = {0};
    point_Jac384 A, AA, B, BB, *T_fixed = NULL;
    
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Curve arithmetic: Weierstrass a=-3 over GF(2^384-317) \n\n"); 

    // Point doubling (Weierstrass a=-3)
    passed = TRUE;
    eccset_Jac384(A, JacCurve); 
    eccconvert_aff_to_jac_Jac384(A, P);

    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        eccdouble_jac_Jac384(P, JacCurve);          // 2*P
        eccdouble_waff_384(A, JacCurve);            // 2*A
    }
    eccnorm_Jac384(P, AA, JacCurve);

    if (fpcompare384(A->x,AA->x)!=0 || fpcompare384(A->y,AA->y)!=0) passed=FALSE;
    if (passed==TRUE) printf("  Point doubling tests .................................................................... PASSED");
    else { printf("  Point doubling tests ... FAILED"); printf("\n"); return FALSE; }
    printf("\n");
    
    // "Complete" point addition (Weierstrass a=-3)
    passed = TRUE;
    eccset_Jac384(A, JacCurve);
    eccconvert_aff_to_jac_Jac384(A, P);
    ecccopy_jac_Jac384(P, R); 
    eccdouble_jac_Jac384(P, JacCurve);                  // P = 2P 
    eccset_Jac384(AA, JacCurve);
    eccdouble_waff_384(AA, JacCurve);                   // AA = 2A

    for (n=0; n<1; n++)
    {
        eccadd_jac_Jac384(R, P, JacCurve);              // P = P+Q
        eccadd_waff_384(A, AA, JacCurve);               // AA = AA+A
    }    
    eccnorm_Jac384(P, A, JacCurve);
    
    if (fpcompare384(A->x,AA->x)!=0 || fpcompare384(A->y,AA->y)!=0) passed=FALSE;
    if (passed==TRUE) printf("  (Complete) point addition tests ......................................................... PASSED");
    else { printf("  (Complete) point addition tests ... FAILED"); printf("\n"); return FALSE; }
    printf("\n");
    
    // Variable-base scalar multiplication (Weierstrass a=-3)
    eccset_Jac384(A, JacCurve); 
    
    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        random384(k, JacCurve->order); 

        ecc_mul_waff_384(A, k, AA, JacCurve);    
        ecc_scalar_mul_Jac384(A, k, B, JacCurve);
        
        if (fpcompare384(AA->x,B->x)!=0 || fpcompare384(AA->y,B->y)!=0) { passed=FALSE; break; }
    }    

    if (passed==TRUE) printf("  Variable-base scalar multiplication tests ............................................... PASSED");
    else { printf("  Variable-base scalar multiplication tests ... FAILED"); printf("\n"); return FALSE; }
    printf("\n");
        
    // Fixed-base scalar multiplication (Weierstrass a=-3)
    eccset_Jac384(A, JacCurve); 
    
    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        random384(k, JacCurve->order); 
        ecc_mul_waff_384(A, k, AA, JacCurve);  

        T_fixed = ecc_precomp_fixed_Jac384(A, MEM_LARGE, JacCurve);
        ecc_scalar_mul_fixed_Jac384(T_fixed, k, B, MEM_LARGE, JacCurve);
        ecc_destroy_precomp_Jac384(T_fixed); 
        
        if (fpcompare384(AA->x,B->x)!=0 || fpcompare384(AA->y,B->y)!=0) { passed=FALSE; break; }
    }    

    if (passed==TRUE) printf("  Fixed-base scalar multiplication tests .................................................. PASSED");
    else { printf("  Fixed-base scalar multiplication tests ... FAILED"); printf("\n"); return FALSE; }
    printf("\n");
    
    // Double-scalar multiplication (Weierstrass a=-3)
    passed = TRUE;
    eccset_Jac384(A, JacCurve); 
    random384(k, JacCurve->order); 
    ecc_mul_waff_384(A, k, B, JacCurve);    // Base points are A and B
    
    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        random384(k, JacCurve->order); 
        random384(kk, JacCurve->order); 

        ecc_mul_waff_384(A, k, AA, JacCurve);    
        ecc_mul_waff_384(B, kk, BB, JacCurve);
        eccadd_waff_384(BB, AA, JacCurve);

        T_fixed = ecc_precomp_dblmul_Jac384(A, MEM_LARGE, JacCurve);
        ecc_double_scalar_mul_Jac384(T_fixed, k, B, kk, BB, MEM_LARGE, JacCurve); 
        ecc_destroy_precomp_Jac384(T_fixed); 
        
        if (fpcompare384(AA->x,BB->x)!=0 || fpcompare384(AA->y,BB->y)!=0) { passed=FALSE; break; }
    }      

    if (passed==TRUE) printf("  Double-scalar multiplication tests ...................................................... PASSED");
    else { printf("  Double-scalar multiplication tests ... FAILED"); printf("\n"); return FALSE; }
    printf("\n");
    
    return OK;
}


BOOL ecc_test384_te(PCurveStruct TedCurve)
{ // Tests for curve Ted384
    BOOL OK = TRUE;
    dig n;
    BOOL passed;
    point_Ted384 A, AA, U, V;
    point_extproj_Ted384 P, R;
    point_Jac384 PP, QQ, RR, UU, VV;
    CurveStruct WeierstrassCurve = {0};
    dig384 d, t1, t2, t3;
    dig k[ML_WORDS384], kk[ML_WORDS384];
    point_extaff_precomp_Ted384 *T_fixed= NULL;
    
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Curve arithmetic: twisted Edwards a=-1 over GF(2^384-317) \n\n"); 

    // Set Weierstrass curve isomorphic to TedCurve
    fpcopy384(TedCurve->parameter2, d);
    fpzero384(t3); t3[0] = 14;
    fpmul384(t3, d, t1);                               // t1 = 14d
    fpsqr384(d, t2);                                   // t2 = d^2
    fpsub384(t1, t2, t2);                              // t2 = 14d-d^2
    fpzero384(t1); t1[0] = 1;    
    fpsub384(t2, t1, t2);                              // t2 = 14d-d^2-1
    fpzero384(t1); t1[0] = 48;
    fpinv384(t1);                                      // t1 = 1/48
    fpmul384(t1, t2, WeierstrassCurve.parameter1);     // a = (14d-d^2-1)/48
    WeierstrassCurve.nbits = TedCurve->nbits;
    WeierstrassCurve.rbits = TedCurve->rbits;
    WeierstrassCurve.pbits = TedCurve->pbits;

    // Point doubling (twisted Edwards a=-1)
    passed = TRUE;
    eccset_Ted384(A, TedCurve); 
    eccconvert_aff_to_extproj_Ted384(A, P);
    ecc_Ted384_to_weierstrass(A, PP, TedCurve);

    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        eccdouble_extproj_Ted384(P, TedCurve);        // 2*P
        eccdouble_waff_384(PP, &WeierstrassCurve);    // 2*PP
    }
    eccnorm_Ted384(P, A, TedCurve);
    ecc_weierstrass_to_Ted384(PP, AA, TedCurve);

    if (fpcompare384(A->x,AA->x)!=0 || fpcompare384(A->y,AA->y)!=0) passed=FALSE;
    if (passed==TRUE) printf("  Point doubling tests .................................................................... PASSED");
    else { printf("  Point doubling tests ... FAILED"); printf("\n"); return FALSE; }
    printf("\n");
    
    // Point addition (twisted Edwards a=-1)
    passed = TRUE;
    eccset_Ted384(A, TedCurve);
    eccconvert_aff_to_extproj_Ted384(A, P);
    ecccopy_extproj_Ted384(P, R);
    eccdouble_extproj_Ted384(P, TedCurve);             // P = 2P    
        
    eccset_Ted384(A, TedCurve); 
    ecc_Ted384_to_weierstrass(A, QQ, TedCurve);
    ecc_Ted384_to_weierstrass(A, PP, TedCurve);
    eccdouble_waff_384(PP, &WeierstrassCurve);         // PP = 2P

    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        eccadd_extproj_Ted384(R, P, TedCurve);         // P = P+Q
        eccadd_waff_384(QQ, PP, &WeierstrassCurve);    // PP = PP+QQ
    }    
    eccnorm_Ted384(P, A, TedCurve);
    ecc_weierstrass_to_Ted384(PP, AA, TedCurve);
    
    if (fpcompare384(A->x,AA->x)!=0 || fpcompare384(A->y,AA->y)!=0) passed=FALSE;
    if (passed==TRUE) printf("  Point addition tests .................................................................... PASSED");
    else { printf("  Point addition tests ... FAILED"); printf("\n"); return FALSE; }
    printf("\n");
    
    // Variable-base scalar multiplication (twisted Edwards a=-1)
    eccset_Ted384(A, TedCurve); 
    ecc_Ted384_to_weierstrass(A, PP, TedCurve);
    eccdouble_waff_384(PP, &WeierstrassCurve);
    eccdouble_waff_384(PP, &WeierstrassCurve);
        
    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        random384(k, TedCurve->order); 
        
        ecc_mul_waff_384(PP, k, RR, &WeierstrassCurve);    
        ecc_scalar_mul_Ted384(A, k, AA, TedCurve);
        ecc_Ted384_to_weierstrass(AA, QQ, TedCurve);
        
        if (fpcompare384(QQ->x,RR->x)!=0 || fpcompare384(QQ->y,RR->y)!=0) { passed=FALSE; break; }
    }    

    if (passed==TRUE) printf("  Variable-base scalar multiplication tests ............................................... PASSED");
    else { printf("  Variable-base scalar multiplication tests ... FAILED"); printf("\n"); return FALSE; }
    printf("\n");

    // Fixed-base scalar multiplication (twisted Edwards a=-1)
    eccset_Ted384(A, TedCurve); 
    ecc_Ted384_to_weierstrass(A, PP, TedCurve);
    eccdouble_waff_384(PP, &WeierstrassCurve);
    eccdouble_waff_384(PP, &WeierstrassCurve);
    
    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        random384(k, TedCurve->order); 
        ecc_mul_waff_384(PP, k, RR, &WeierstrassCurve);   

        T_fixed = ecc_precomp_fixed_Ted384(A, MEM_LARGE, TedCurve);
        ecc_scalar_mul_fixed_Ted384(T_fixed, k, AA, MEM_LARGE, TedCurve);
        ecc_Ted384_to_weierstrass(AA, QQ, TedCurve);
        ecc_destroy_precomp_Ted384(T_fixed); 
        
        if (fpcompare384(QQ->x,RR->x)!=0 || fpcompare384(QQ->y,RR->y)!=0) { passed=FALSE; break; }
    }    

    if (passed==TRUE) printf("  Fixed-base scalar multiplication tests .................................................. PASSED");
    else { printf("  Fixed-base scalar multiplication tests ... FAILED"); printf("\n"); return FALSE; }
    printf("\n");
    
    // Double-scalar multiplication (twisted Edwards a=-1)
    passed = TRUE;
    eccset_Ted384(A, TedCurve); 
    ecc_Ted384_to_weierstrass(A, PP, TedCurve);
    eccdouble_waff_384(PP, &WeierstrassCurve);
    eccdouble_waff_384(PP, &WeierstrassCurve);
    random384(k, TedCurve->order); 
    ecc_mul_waff_384(PP, k, RR, &WeierstrassCurve);    
    ecc_weierstrass_to_Ted384(RR, AA, TedCurve);    // Base points are (A, AA) in twisted Edwards
    eccdouble_waff_384(RR, &WeierstrassCurve);
    eccdouble_waff_384(RR, &WeierstrassCurve);      // Base points are 4*(PP, RR) in Weierstrass
    
    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        random384(k, TedCurve->order); 
        random384(kk, TedCurve->order); 

        ecc_mul_waff_384(PP, k, UU, &WeierstrassCurve);    
        ecc_mul_waff_384(RR, kk, VV, &WeierstrassCurve);
        eccadd_waff_384(VV, UU, &WeierstrassCurve);  
        ecc_weierstrass_to_Ted384(UU, U, TedCurve); 

        T_fixed = ecc_precomp_dblmul_Ted384(A, MEM_LARGE, TedCurve);
        ecc_double_scalar_mul_Ted384(T_fixed, k, AA, kk, V, MEM_LARGE, TedCurve); 
        ecc_destroy_precomp_Ted384(T_fixed); 
        
        if (fpcompare384(U->x,V->x)!=0 || fpcompare384(U->y,V->y)!=0) { passed=FALSE; break; }
    }      

    if (passed==TRUE) printf("  Double-scalar multiplication tests ...................................................... PASSED");
    else { printf("  Double-scalar multiplication tests ... FAILED"); printf("\n"); return FALSE; }
    printf("\n");

    return OK;
}

#endif


#ifdef ECCURVES_512

BOOL ecc_test512_w(PCurveStruct JacCurve)
{ // Tests for curve Jac512
    BOOL OK = TRUE;
    dig n;
    BOOL passed;
    point_jac_Jac512 P, R;
    dig512 k = {0}, kk = {0};
    point_Jac512 A, AA, B, BB, *T_fixed = NULL;
    
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Curve arithmetic: Weierstrass a=-3 over GF(2^512-569) \n\n"); 

    // Point doubling (Weierstrass a=-3)
    passed = TRUE;
    eccset_Jac512(A, JacCurve); 
    eccconvert_aff_to_jac_Jac512(A, P);

    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        eccdouble_jac_Jac512(P, JacCurve);          // 2*P
        eccdouble_waff_512(A, JacCurve);            // 2*A
    }
    eccnorm_Jac512(P, AA, JacCurve);

    if (fpcompare512(A->x,AA->x)!=0 || fpcompare512(A->y,AA->y)!=0) passed=FALSE;
    if (passed==TRUE) printf("  Point doubling tests .................................................................... PASSED");
    else { printf("  Point doubling tests ... FAILED"); printf("\n"); return FALSE; }
    printf("\n");
        
    // "Complete" point addition (Weierstrass a=-3)
    passed = TRUE;
    eccset_Jac512(A, JacCurve);
    eccconvert_aff_to_jac_Jac512(A, P);
    ecccopy_jac_Jac512(P, R); 
    eccdouble_jac_Jac512(P, JacCurve);                  // P = 2P 
    eccset_Jac512(AA, JacCurve);
    eccdouble_waff_512(AA, JacCurve);                   // AA = 2A

    for (n=0; n<1; n++)
    {
        eccadd_jac_Jac512(R, P, JacCurve);              // P = P+Q
        eccadd_waff_512(A, AA, JacCurve);               // AA = AA+A
    }    
    eccnorm_Jac512(P, A, JacCurve);
    
    if (fpcompare512(A->x,AA->x)!=0 || fpcompare512(A->y,AA->y)!=0) passed=FALSE;
    if (passed==TRUE) printf("  (Complete) point addition tests ......................................................... PASSED");
    else { printf("  (Complete) point addition tests ... FAILED"); printf("\n"); return FALSE; }
    printf("\n");
    
    // Variable-base scalar multiplication (Weierstrass a=-3)
    eccset_Jac512(A, JacCurve); 
    
    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        random512(k, JacCurve->order); 

        ecc_mul_waff_512(A, k, AA, JacCurve);    
        ecc_scalar_mul_Jac512(A, k, B, JacCurve);
        
        if (fpcompare512(AA->x,B->x)!=0 || fpcompare512(AA->y,B->y)!=0) { passed=FALSE; break; }
    }    

    if (passed==TRUE) printf("  Variable-base scalar multiplication tests ............................................... PASSED");
    else { printf("  Variable-base scalar multiplication tests ... FAILED"); printf("\n"); return FALSE; }
    printf("\n");
        
    // Fixed-base scalar multiplication (Weierstrass a=-3)
    eccset_Jac512(A, JacCurve); 
    
    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        random512(k, JacCurve->order); 
        ecc_mul_waff_512(A, k, AA, JacCurve);  

        T_fixed = ecc_precomp_fixed_Jac512(A, MEM_LARGE, JacCurve);
        ecc_scalar_mul_fixed_Jac512(T_fixed, k, B, MEM_LARGE, JacCurve);
        ecc_destroy_precomp_Jac512(T_fixed); 
        
        if (fpcompare512(AA->x,B->x)!=0 || fpcompare512(AA->y,B->y)!=0) { passed=FALSE; break; }
    }    

    if (passed==TRUE) printf("  Fixed-base scalar multiplication tests .................................................. PASSED");
    else { printf("  Fixed-base scalar multiplication tests ... FAILED"); printf("\n"); return FALSE; }
    printf("\n");
    
    // Double-scalar multiplication (Weierstrass a=-3)
    passed = TRUE;
    eccset_Jac512(A, JacCurve); 
    random512(k, JacCurve->order); 
    ecc_mul_waff_512(A, k, B, JacCurve);    // Base points are A and B
    
    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        random512(k, JacCurve->order); 
        random512(kk, JacCurve->order); 

        ecc_mul_waff_512(A, k, AA, JacCurve);    
        ecc_mul_waff_512(B, kk, BB, JacCurve);
        eccadd_waff_512(BB, AA, JacCurve);

        T_fixed = ecc_precomp_dblmul_Jac512(A, MEM_LARGE, JacCurve);
        ecc_double_scalar_mul_Jac512(T_fixed, k, B, kk, BB, MEM_LARGE, JacCurve); 
        ecc_destroy_precomp_Jac512(T_fixed); 
        
        if (fpcompare512(AA->x,BB->x)!=0 || fpcompare512(AA->y,BB->y)!=0) { passed=FALSE; break; }
    }      

    if (passed==TRUE) printf("  Double-scalar multiplication tests ...................................................... PASSED");
    else { printf("  Double-scalar multiplication tests ... FAILED"); printf("\n"); return FALSE; }
    printf("\n");
    
    return OK;
}


BOOL ecc_test512_te(PCurveStruct TedCurve)
{ // Tests for curve Ted512
    BOOL OK = TRUE;
    dig n;
    BOOL passed;
    point_Ted512 A, AA, U, V;
    point_extproj_Ted512 P, R;
    point_Jac512 PP, QQ, RR, UU, VV;
    CurveStruct WeierstrassCurve = {0};
    dig512 d, t1, t2, t3;
    dig k[ML_WORDS512], kk[ML_WORDS512];
    point_extaff_precomp_Ted512 *T_fixed= NULL;
    
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Curve arithmetic: twisted Edwards a=-1 over GF(2^512-569) \n\n"); 

    // Set Weierstrass curve isomorphic to TedCurve
    fpcopy512(TedCurve->parameter2, d);
    fpzero512(t3); t3[0] = 14;
    fpmul512(t3, d, t1);                               // t1 = 14d
    fpsqr512(d, t2);                                   // t2 = d^2
    fpsub512(t1, t2, t2);                              // t2 = 14d-d^2
    fpzero512(t1); t1[0] = 1;    
    fpsub512(t2, t1, t2);                              // t2 = 14d-d^2-1
    fpzero512(t1); t1[0] = 48;
    fpinv512(t1);                                      // t1 = 1/48
    fpmul512(t1, t2, WeierstrassCurve.parameter1);     // a = (14d-d^2-1)/48
    WeierstrassCurve.nbits = TedCurve->nbits;
    WeierstrassCurve.rbits = TedCurve->rbits;
    WeierstrassCurve.pbits = TedCurve->pbits;

    // Point doubling (twisted Edwards a=-1)
    passed = TRUE;
    eccset_Ted512(A, TedCurve); 
    eccconvert_aff_to_extproj_Ted512(A, P);
    ecc_Ted512_to_weierstrass(A, PP, TedCurve);

    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        eccdouble_extproj_Ted512(P, TedCurve);        // 2*P
        eccdouble_waff_512(PP, &WeierstrassCurve);    // 2*PP
    }
    eccnorm_Ted512(P, A, TedCurve);
    ecc_weierstrass_to_Ted512(PP, AA, TedCurve);

    if (fpcompare512(A->x,AA->x)!=0 || fpcompare512(A->y,AA->y)!=0) passed=FALSE;
    if (passed==TRUE) printf("  Point doubling tests .................................................................... PASSED");
    else { printf("  Point doubling tests ... FAILED"); printf("\n"); return FALSE; }
    printf("\n");
    
    // Point addition (twisted Edwards a=-1)
    passed = TRUE;
    eccset_Ted512(A, TedCurve);
    eccconvert_aff_to_extproj_Ted512(A, P);
    ecccopy_extproj_Ted512(P, R);
    eccdouble_extproj_Ted512(P, TedCurve);             // P = 2P    
        
    eccset_Ted512(A, TedCurve); 
    ecc_Ted512_to_weierstrass(A, QQ, TedCurve);
    ecc_Ted512_to_weierstrass(A, PP, TedCurve);
    eccdouble_waff_512(PP, &WeierstrassCurve);         // PP = 2P

    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        eccadd_extproj_Ted512(R, P, TedCurve);         // P = P+Q
        eccadd_waff_512(QQ, PP, &WeierstrassCurve);    // PP = PP+QQ
    }    
    eccnorm_Ted512(P, A, TedCurve);
    ecc_weierstrass_to_Ted512(PP, AA, TedCurve);
    
    if (fpcompare512(A->x,AA->x)!=0 || fpcompare512(A->y,AA->y)!=0) passed=FALSE;
    if (passed==TRUE) printf("  Point addition tests .................................................................... PASSED");
    else { printf("  Point addition tests ... FAILED"); printf("\n"); return FALSE; }
    printf("\n");
    
    // Variable-base scalar multiplication (twisted Edwards a=-1)
    eccset_Ted512(A, TedCurve); 
    ecc_Ted512_to_weierstrass(A, PP, TedCurve);
    eccdouble_waff_512(PP, &WeierstrassCurve);
    eccdouble_waff_512(PP, &WeierstrassCurve);
        
    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        random512(k, TedCurve->order); 
        
        ecc_mul_waff_512(PP, k, RR, &WeierstrassCurve);    
        ecc_scalar_mul_Ted512(A, k, AA, TedCurve);
        ecc_Ted512_to_weierstrass(AA, QQ, TedCurve);
        
        if (fpcompare512(QQ->x,RR->x)!=0 || fpcompare512(QQ->y,RR->y)!=0) { passed=FALSE; break; }
    }    

    if (passed==TRUE) printf("  Variable-base scalar multiplication tests ............................................... PASSED");
    else { printf("  Variable-base scalar multiplication tests ... FAILED"); printf("\n"); return FALSE; }
    printf("\n");

    // Fixed-base scalar multiplication (twisted Edwards a=-1)
    eccset_Ted512(A, TedCurve); 
    ecc_Ted512_to_weierstrass(A, PP, TedCurve);
    eccdouble_waff_512(PP, &WeierstrassCurve);
    eccdouble_waff_512(PP, &WeierstrassCurve);
    
    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        random512(k, TedCurve->order); 
        ecc_mul_waff_512(PP, k, RR, &WeierstrassCurve);   

        T_fixed = ecc_precomp_fixed_Ted512(A, MEM_LARGE, TedCurve);
        ecc_scalar_mul_fixed_Ted512(T_fixed, k, AA, MEM_LARGE, TedCurve);
        ecc_Ted512_to_weierstrass(AA, QQ, TedCurve);
        ecc_destroy_precomp_Ted512(T_fixed); 
        
        if (fpcompare512(QQ->x,RR->x)!=0 || fpcompare512(QQ->y,RR->y)!=0) { passed=FALSE; break; }
    }    

    if (passed==TRUE) printf("  Fixed-base scalar multiplication tests .................................................. PASSED");
    else { printf("  Fixed-base scalar multiplication tests ... FAILED"); printf("\n"); return FALSE; }
    printf("\n");
    
    // Double-scalar multiplication (twisted Edwards a=-1)
    passed = TRUE;
    eccset_Ted512(A, TedCurve); 
    ecc_Ted512_to_weierstrass(A, PP, TedCurve);
    eccdouble_waff_512(PP, &WeierstrassCurve);
    eccdouble_waff_512(PP, &WeierstrassCurve);
    random512(k, TedCurve->order); 
    ecc_mul_waff_512(PP, k, RR, &WeierstrassCurve);    
    ecc_weierstrass_to_Ted512(RR, AA, TedCurve);    // Base points are (A, AA) in twisted Edwards
    eccdouble_waff_512(RR, &WeierstrassCurve);
    eccdouble_waff_512(RR, &WeierstrassCurve);      // Base points are 4*(PP, RR) in Weierstrass
    
    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        random512(k, TedCurve->order); 
        random512(kk, TedCurve->order); 

        ecc_mul_waff_512(PP, k, UU, &WeierstrassCurve);    
        ecc_mul_waff_512(RR, kk, VV, &WeierstrassCurve);
        eccadd_waff_512(VV, UU, &WeierstrassCurve);  
        ecc_weierstrass_to_Ted512(UU, U, TedCurve); 

        T_fixed = ecc_precomp_dblmul_Ted512(A, MEM_LARGE, TedCurve);
        ecc_double_scalar_mul_Ted512(T_fixed, k, AA, kk, V, MEM_LARGE, TedCurve); 
        ecc_destroy_precomp_Ted512(T_fixed); 
        
        if (fpcompare512(U->x,V->x)!=0 || fpcompare512(U->y,V->y)!=0) { passed=FALSE; break; }
    }      

    if (passed==TRUE) printf("  Double-scalar multiplication tests ...................................................... PASSED");
    else { printf("  Double-scalar multiplication tests ... FAILED"); printf("\n"); return FALSE; }
    printf("\n");

    return OK;
}

#endif



/****************** BENCHMARK TESTS *******************/
/******************************************************/

#ifdef ECCURVES_256

BOOL ecc_run256_w(PCurveStruct JacCurve)
{ // Benchmarking for curve Jac256
    BOOL OK = TRUE;
    unsigned int n;
    unsigned long long cycles, cycles1, cycles2;
    point_jac_Jac256 P, R;
    point_Jac256 A, AA, B;
    dig k[ML_WORDS256], kk[ML_WORDS256];
    point_Jac256 *T_fixed = NULL;

    SetThreadAffinityMask(GetCurrentThread(), 1);       // All threads are set to run in the same node
    SetThreadPriority(GetCurrentThread(), 2);           // Set to highest priority
        
    printf("\n\nBENCHMARKING \n"); 
    printf("--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Curve arithmetic: Weierstrass a=-3 over GF(2^256-189) \n\n"); 
    
    // Point doubling (Weierstrass a=-3)
    eccset_Jac256(A, JacCurve);
    eccconvert_aff_to_jac_Jac256(A, P);

    cycles = 0;
    for (n=0; n<ML_BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        eccdouble_jac_Jac256(P, JacCurve);
        eccdouble_jac_Jac256(P, JacCurve);
        eccdouble_jac_Jac256(P, JacCurve);
        eccdouble_jac_Jac256(P, JacCurve);
        eccdouble_jac_Jac256(P, JacCurve);
        eccdouble_jac_Jac256(P, JacCurve);
        eccdouble_jac_Jac256(P, JacCurve);
        eccdouble_jac_Jac256(P, JacCurve);
        eccdouble_jac_Jac256(P, JacCurve);
        eccdouble_jac_Jac256(P, JacCurve);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Point doubling runs in .......................................... %8lld cycles", cycles/(ML_BENCH_LOOPS*10));
    printf("\n"); 
    
    // "Complete" point addition (Weierstrass a=-3)    
    eccset_Jac256(A, JacCurve);
    eccconvert_aff_to_jac_Jac256(A, P);
    ecccopy_jac_Jac256(P, R); 
    eccdouble_jac_Jac256(P, JacCurve);

    cycles = 0;
    for (n=0; n<ML_BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        eccadd_jac_Jac256(R, P, JacCurve);
        eccadd_jac_Jac256(R, P, JacCurve);
        eccadd_jac_Jac256(R, P, JacCurve);
        eccadd_jac_Jac256(R, P, JacCurve);
        eccadd_jac_Jac256(R, P, JacCurve);
        eccadd_jac_Jac256(R, P, JacCurve);
        eccadd_jac_Jac256(R, P, JacCurve);
        eccadd_jac_Jac256(R, P, JacCurve);
        eccadd_jac_Jac256(R, P, JacCurve);
        eccadd_jac_Jac256(R, P, JacCurve);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  (Complete) point addition runs in ............................... %8lld cycles", cycles/(ML_BENCH_LOOPS*10));
    printf("\n"); 
    
    // Variable-base scalar multiplication (Weierstrass a=-3)
    eccset_Jac256(A, JacCurve);
    random256(k, JacCurve->order); 
    
    cycles = 0;
    for (n=0; n<ML_SHORT_BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        ecc_scalar_mul_Jac256(A, k, AA, JacCurve);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    } 
    printf("  Variable-base scalar mul runs in ................................ %8lld cycles", cycles/(ML_SHORT_BENCH_LOOPS));
#ifdef ML_COUNT
    ninv=0; nmul=0; nsqr=0; nadd=0; nlut=0; nclut=0; ncsel=0;
    ecc_scalar_mul_Jac256(A, k, AA, JacCurve);
    printf(" using %d inversions, %d muls, %d sqrs, %d adds, %d luts, %d 5-point luts and %d field selects", ninv, nmul, nsqr, nadd, nlut, nclut, ncsel);
#endif
    printf("\n"); 

    // Fixed-base scalar multiplication (Weierstrass a=-3)
    eccset_Jac256(A, JacCurve); 
    T_fixed = ecc_precomp_fixed_Jac256(A, MEM_LARGE, JacCurve);

    cycles = 0;
    for (n=0; n<ML_SHORT_BENCH_LOOPS; n++)
    {
        random256(k, JacCurve->order);  
        cycles1 = cpucycles();
        ecc_scalar_mul_fixed_Jac256(T_fixed, k, AA, MEM_LARGE, JacCurve);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    ecc_destroy_precomp_Jac256(T_fixed); 

    printf("  Fixed-base scalar mul (memory model=MEM_LARGE) runs in .......... %8lld cycles", cycles/(ML_SHORT_BENCH_LOOPS));
    printf("\n"); 
    
    // Double-scalar multiplication (Weierstrass a=-3)
    eccset_Jac256(A, JacCurve); 
    random256(k, JacCurve->order); 
    ecc_mul_waff_256(A, k, AA, JacCurve);    // Base points are A and AA

    T_fixed = ecc_precomp_dblmul_Jac256(AA, MEM_LARGE, JacCurve);
    cycles = 0;
    for (n=0; n<ML_SHORT_BENCH_LOOPS; n++)
    {
        random256(k, JacCurve->order); 
        random256(kk, JacCurve->order);
        cycles1 = cpucycles();
        ecc_double_scalar_mul_Jac256(T_fixed, k, A, kk, B, MEM_LARGE, JacCurve);  
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    ecc_destroy_precomp_Jac256(T_fixed); 

    printf("  Double-base scalar mul (memory model=MEM_LARGE) runs in ......... %8lld cycles", cycles/(ML_SHORT_BENCH_LOOPS));
    printf("\n"); 
    
    return OK;
}


BOOL ecc_run256_te(PCurveStruct TedCurve)
{ // Benchmarking for curve Ted256
    BOOL OK = TRUE;
    unsigned int n;
    unsigned long long cycles, cycles1, cycles2;
    point_Ted256 A, AA;
    point_extproj_Ted256 P, R;
    dig k[ML_WORDS256], kk[ML_WORDS256];
    point_extaff_precomp_Ted256 *T_fixed = NULL;

    SetThreadAffinityMask(GetCurrentThread(), 1);       // All threads are set to run in the same node
    SetThreadPriority(GetCurrentThread(), 2);           // Set to highest priority
        
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Curve arithmetic: twisted Edwards a=-1 over GF(2^256-189) \n\n"); 

    // Point doubling (twisted Edwards a=-1)
    eccset_Ted256(A, TedCurve);
    eccconvert_aff_to_extproj_Ted256(A, P);

    cycles = 0;
    for (n=0; n<ML_BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        eccdouble_extproj_Ted256(P, TedCurve);
        eccdouble_extproj_Ted256(P, TedCurve);
        eccdouble_extproj_Ted256(P, TedCurve);
        eccdouble_extproj_Ted256(P, TedCurve);
        eccdouble_extproj_Ted256(P, TedCurve);
        eccdouble_extproj_Ted256(P, TedCurve);
        eccdouble_extproj_Ted256(P, TedCurve);
        eccdouble_extproj_Ted256(P, TedCurve);
        eccdouble_extproj_Ted256(P, TedCurve);
        eccdouble_extproj_Ted256(P, TedCurve);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Point doubling runs in .......................................... %8lld cycles", cycles/(ML_BENCH_LOOPS*10));
    printf("\n"); 

    // Point addition (twisted Edwards a=-1)    
    eccset_Ted256(A, TedCurve);
    eccconvert_aff_to_extproj_Ted256(A, P);
    ecccopy_extproj_Ted256(P, R);
    eccdouble_extproj_Ted256(P, TedCurve);

    cycles = 0;
    for (n=0; n<ML_BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        eccadd_extproj_Ted256(R, P, TedCurve);
        eccadd_extproj_Ted256(R, P, TedCurve);
        eccadd_extproj_Ted256(R, P, TedCurve);
        eccadd_extproj_Ted256(R, P, TedCurve);
        eccadd_extproj_Ted256(R, P, TedCurve);
        eccadd_extproj_Ted256(R, P, TedCurve);
        eccadd_extproj_Ted256(R, P, TedCurve);
        eccadd_extproj_Ted256(R, P, TedCurve);
        eccadd_extproj_Ted256(R, P, TedCurve);
        eccadd_extproj_Ted256(R, P, TedCurve);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  (Complete) point addition runs in ............................... %8lld cycles", cycles/(ML_BENCH_LOOPS*10));
    printf("\n"); 

    // Variable-base scalar multiplication (twisted Edwards a=-1)
    eccset_Ted256(A, TedCurve);
    random256(k, TedCurve->order); 
    
    cycles = 0;
    for (n=0; n<ML_SHORT_BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        ecc_scalar_mul_Ted256(A, k, AA, TedCurve);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Variable-base scalar mul runs in ................................ %8lld cycles", cycles/(ML_SHORT_BENCH_LOOPS));
#ifdef ML_COUNT
    ninv=0; nmul=0; nsqr=0; nadd=0; nlut=0; 
    ecc_scalar_mul_Ted256(A, k, AA, TedCurve);
    printf(" using %d inversions, %d muls, %d sqrs, %d adds and %d luts", ninv, nmul, nsqr, nadd, nlut);
#endif
    printf("\n"); 

    // Fixed-base scalar multiplication (twisted Edwards a=-1)
    eccset_Ted256(A, TedCurve); 
    T_fixed = ecc_precomp_fixed_Ted256(A, MEM_LARGE, TedCurve);

    cycles = 0;
    for (n=0; n<ML_SHORT_BENCH_LOOPS; n++)
    {
        random256(k, TedCurve->order);  
        cycles1 = cpucycles();
        ecc_scalar_mul_fixed_Ted256(T_fixed, k, AA, MEM_LARGE, TedCurve);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    ecc_destroy_precomp_Ted256(T_fixed); 

    printf("  Fixed-base scalar mul (memory model=MEM_LARGE) runs in .......... %8lld cycles", cycles/(ML_SHORT_BENCH_LOOPS));
    printf("\n"); 
    
    // Double-scalar multiplication (twisted Edwards a=-1)
    eccset_Ted256(A, TedCurve); 
    random256(k, TedCurve->order); 
    ecc_scalar_mul_Ted256(A, k, AA, TedCurve);    // Base points are A and AA
            
    T_fixed = ecc_precomp_dblmul_Ted256(AA, MEM_LARGE, TedCurve);
    cycles = 0;
    for (n=0; n<ML_SHORT_BENCH_LOOPS; n++)
    {
        random256(k, TedCurve->order); 
        random256(kk, TedCurve->order);
        cycles1 = cpucycles();
        ecc_double_scalar_mul_Ted256(T_fixed, k, A, kk, AA, MEM_LARGE, TedCurve); 
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    ecc_destroy_precomp_Ted256(T_fixed); 

    printf("  Double-base scalar mul (memory model=MEM_LARGE) runs in ......... %8lld cycles", cycles/(ML_SHORT_BENCH_LOOPS));
    printf("\n"); 

    return OK;
}

#endif


#ifdef ECCURVES_384

BOOL ecc_run384_w(PCurveStruct JacCurve)
{ // Benchmarking for curve Jac384
    BOOL OK = TRUE;
    unsigned int n;
    unsigned long long cycles, cycles1, cycles2;
    point_Jac384 A, AA, B;
    point_jac_Jac384 P, R;
    dig384 k, kk;
    point_Jac384 *T_fixed = NULL;

    SetThreadAffinityMask(GetCurrentThread(), 1);       // All threads are set to run in the same node
    SetThreadPriority(GetCurrentThread(), 2);           // Set to highest priority
        
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Curve arithmetic: Weierstrass a=-3 over GF(2^384-317) \n\n"); 

    // Point doubling (Weierstrass a=-3)   
    eccset_Jac384(A, JacCurve);
    eccconvert_aff_to_jac_Jac384(A, P);

    cycles = 0;
    for (n=0; n<ML_BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        eccdouble_jac_Jac384(P, JacCurve);
        eccdouble_jac_Jac384(P, JacCurve);
        eccdouble_jac_Jac384(P, JacCurve);
        eccdouble_jac_Jac384(P, JacCurve);
        eccdouble_jac_Jac384(P, JacCurve);
        eccdouble_jac_Jac384(P, JacCurve);
        eccdouble_jac_Jac384(P, JacCurve);
        eccdouble_jac_Jac384(P, JacCurve);
        eccdouble_jac_Jac384(P, JacCurve);
        eccdouble_jac_Jac384(P, JacCurve);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Point doubling runs in .......................................... %8lld cycles", cycles/(ML_BENCH_LOOPS*10));
    printf("\n"); 
        
    // "Complete" point addition (Weierstrass a=-3)    
    eccset_Jac384(A, JacCurve);
    eccconvert_aff_to_jac_Jac384(A, P);
    ecccopy_jac_Jac384(P, R); 
    eccdouble_jac_Jac384(P, JacCurve);

    cycles = 0;
    for (n=0; n<ML_BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        eccadd_jac_Jac384(R, P, JacCurve);
        eccadd_jac_Jac384(R, P, JacCurve);
        eccadd_jac_Jac384(R, P, JacCurve);
        eccadd_jac_Jac384(R, P, JacCurve);
        eccadd_jac_Jac384(R, P, JacCurve);
        eccadd_jac_Jac384(R, P, JacCurve);
        eccadd_jac_Jac384(R, P, JacCurve);
        eccadd_jac_Jac384(R, P, JacCurve);
        eccadd_jac_Jac384(R, P, JacCurve);
        eccadd_jac_Jac384(R, P, JacCurve);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  (Complete) point addition runs in ............................... %8lld cycles", cycles/(ML_BENCH_LOOPS*10));
    printf("\n"); 

    // Variable-base scalar multiplication (Weierstrass a=-3)
    eccset_Jac384(A, JacCurve);
    random384(k, JacCurve->order);        
    
    cycles = 0;
    for (n=0; n<ML_SHORT_BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        ecc_scalar_mul_Jac384(A, k, AA, JacCurve);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Variable-base scalar mul runs in ................................ %8lld cycles", cycles/(ML_SHORT_BENCH_LOOPS));
#ifdef ML_COUNT
    ninv=0; nmul=0; nsqr=0; nadd=0; nlut=0; nclut=0; ncsel=0;
    ecc_scalar_mul_Jac384(A, k, AA, JacCurve);
    printf(" using %d inversions, %d muls, %d sqrs, %d adds, %d luts, %d 5-point luts and %d field selects", ninv, nmul, nsqr, nadd, nlut, nclut, ncsel);
#endif    
    printf("\n"); 
    
    // Fixed-base scalar multiplication (Weierstrass a=-3)
    eccset_Jac384(A, JacCurve); 
    T_fixed = ecc_precomp_fixed_Jac384(A, MEM_LARGE, JacCurve);

    cycles = 0;
    for (n=0; n<ML_SHORT_BENCH_LOOPS; n++)
    {
        random384(k, JacCurve->order);  
        cycles1 = cpucycles();
        ecc_scalar_mul_fixed_Jac384(T_fixed, k, AA, MEM_LARGE, JacCurve);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    ecc_destroy_precomp_Jac384(T_fixed); 

    printf("  Fixed-base scalar mul (memory model=MEM_LARGE) runs in .......... %8lld cycles", cycles/(ML_SHORT_BENCH_LOOPS));
    printf("\n"); 
    
    // Double-scalar multiplication (Weierstrass a=-3)
    eccset_Jac384(A, JacCurve); 
    random384(k, JacCurve->order); 
    ecc_mul_waff_384(A, k, AA, JacCurve);    // Base points are A and AA

    T_fixed = ecc_precomp_dblmul_Jac384(AA, MEM_LARGE, JacCurve);
    cycles = 0;
    for (n=0; n<ML_SHORT_BENCH_LOOPS; n++)
    {
        random384(k, JacCurve->order); 
        random384(kk, JacCurve->order);
        cycles1 = cpucycles();
        ecc_double_scalar_mul_Jac384(T_fixed, k, A, kk, B, MEM_LARGE, JacCurve);  
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    ecc_destroy_precomp_Jac384(T_fixed); 

    printf("  Double-base scalar mul (memory model=MEM_LARGE) runs in ......... %8lld cycles", cycles/(ML_SHORT_BENCH_LOOPS));
    printf("\n"); 
    
    return OK;
}


BOOL ecc_run384_te(PCurveStruct TedCurve)
{ // Benchmarking for curve Ted384
    BOOL OK = TRUE;
    unsigned int n;
    unsigned long long cycles, cycles1, cycles2;
    point_Ted384 A, AA;
    point_extproj_Ted384 P, R;
    dig k[ML_WORDS384], kk[ML_WORDS384];
    point_extaff_precomp_Ted384 *T_fixed = NULL;

    SetThreadAffinityMask(GetCurrentThread(), 1);       // All threads are set to run in the same node
    SetThreadPriority(GetCurrentThread(), 2);           // Set to highest priority
        
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Curve arithmetic: twisted Edwards a=-1 over GF(2^384-317) \n\n"); 

    // Point doubling (twisted Edwards a=-1)
    eccset_Ted384(A, TedCurve);
    eccconvert_aff_to_extproj_Ted384(A, P);

    cycles = 0;
    for (n=0; n<ML_BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        eccdouble_extproj_Ted384(P, TedCurve);
        eccdouble_extproj_Ted384(P, TedCurve);
        eccdouble_extproj_Ted384(P, TedCurve);
        eccdouble_extproj_Ted384(P, TedCurve);
        eccdouble_extproj_Ted384(P, TedCurve);
        eccdouble_extproj_Ted384(P, TedCurve);
        eccdouble_extproj_Ted384(P, TedCurve);
        eccdouble_extproj_Ted384(P, TedCurve);
        eccdouble_extproj_Ted384(P, TedCurve);
        eccdouble_extproj_Ted384(P, TedCurve);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Point doubling runs in .......................................... %8lld cycles", cycles/(ML_BENCH_LOOPS*10));
    printf("\n"); 

    // Point addition (twisted Edwards a=-1)    
    eccset_Ted384(A, TedCurve);
    eccconvert_aff_to_extproj_Ted384(A, P);
    ecccopy_extproj_Ted384(P, R);
    eccdouble_extproj_Ted384(P, TedCurve);

    cycles = 0;
    for (n=0; n<ML_BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        eccadd_extproj_Ted384(R, P, TedCurve);
        eccadd_extproj_Ted384(R, P, TedCurve);
        eccadd_extproj_Ted384(R, P, TedCurve);
        eccadd_extproj_Ted384(R, P, TedCurve);
        eccadd_extproj_Ted384(R, P, TedCurve);
        eccadd_extproj_Ted384(R, P, TedCurve);
        eccadd_extproj_Ted384(R, P, TedCurve);
        eccadd_extproj_Ted384(R, P, TedCurve);
        eccadd_extproj_Ted384(R, P, TedCurve);
        eccadd_extproj_Ted384(R, P, TedCurve);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  (Complete) point addition runs in ............................... %8lld cycles", cycles/(ML_BENCH_LOOPS*10));
    printf("\n"); 

    // Variable-base scalar multiplication (twisted Edwards a=-1)
    eccset_Ted384(A, TedCurve);
    random384(k, TedCurve->order); 
    
    cycles = 0;
    for (n=0; n<ML_SHORT_BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        ecc_scalar_mul_Ted384(A, k, AA, TedCurve);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Variable-base scalar mul runs in ................................ %8lld cycles", cycles/(ML_SHORT_BENCH_LOOPS));
#ifdef ML_COUNT
    ninv=0; nmul=0; nsqr=0; nadd=0; nlut=0; 
    ecc_scalar_mul_Ted384(A, k, AA, TedCurve);
    printf(" using %d inversions, %d muls, %d sqrs, %d adds and %d luts", ninv, nmul, nsqr, nadd, nlut);
#endif   
    printf("\n"); 

    // Fixed-base scalar multiplication (twisted Edwards a=-1)
    eccset_Ted384(A, TedCurve); 
    T_fixed = ecc_precomp_fixed_Ted384(A, MEM_LARGE, TedCurve);

    cycles = 0;
    for (n=0; n<ML_SHORT_BENCH_LOOPS; n++)
    {
        random384(k, TedCurve->order);  
        cycles1 = cpucycles();
        ecc_scalar_mul_fixed_Ted384(T_fixed, k, AA, MEM_LARGE, TedCurve);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    ecc_destroy_precomp_Ted384(T_fixed); 

    printf("  Fixed-base scalar mul (memory model=MEM_LARGE) runs in .......... %8lld cycles", cycles/(ML_SHORT_BENCH_LOOPS));
    printf("\n"); 
    
    // Double-scalar multiplication (twisted Edwards a=-1)
    eccset_Ted384(A, TedCurve); 
    random384(k, TedCurve->order); 
    ecc_scalar_mul_Ted384(A, k, AA, TedCurve);    // Base points are A and AA
            
    T_fixed = ecc_precomp_dblmul_Ted384(AA, MEM_LARGE, TedCurve);
    cycles = 0;
    for (n=0; n<ML_SHORT_BENCH_LOOPS; n++)
    {
        random384(k, TedCurve->order); 
        random384(kk, TedCurve->order);
        cycles1 = cpucycles();
        ecc_double_scalar_mul_Ted384(T_fixed, k, A, kk, AA, MEM_LARGE, TedCurve); 
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    ecc_destroy_precomp_Ted384(T_fixed); 

    printf("  Double-base scalar mul (memory model=MEM_LARGE) runs in ......... %8lld cycles", cycles/(ML_SHORT_BENCH_LOOPS));
    printf("\n"); 

    return OK;
}

#endif


#ifdef ECCURVES_512

BOOL ecc_run512_w(PCurveStruct JacCurve)
{ // Benchmarking for curve Jac512
    BOOL OK = TRUE;
    unsigned int n;
    unsigned long long cycles, cycles1, cycles2;
    point_Jac512 A, AA, B;
    point_jac_Jac512 P, R;
    dig512 k, kk;
    point_Jac512 *T_fixed;

    SetThreadAffinityMask(GetCurrentThread(), 1);       // All threads are set to run in the same node
    SetThreadPriority(GetCurrentThread(), 2);           // Set to highest priority
        
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Curve arithmetic: Weierstrass a=-3 over GF(2^512-569) \n\n"); 

    // Point doubling (Weierstrass a=-3) 
    eccset_Jac512(A, JacCurve);
    eccconvert_aff_to_jac_Jac512(A, P);

    cycles = 0;
    for (n=0; n<ML_BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        eccdouble_jac_Jac512(P, JacCurve);
        eccdouble_jac_Jac512(P, JacCurve);
        eccdouble_jac_Jac512(P, JacCurve);
        eccdouble_jac_Jac512(P, JacCurve);
        eccdouble_jac_Jac512(P, JacCurve);
        eccdouble_jac_Jac512(P, JacCurve);
        eccdouble_jac_Jac512(P, JacCurve);
        eccdouble_jac_Jac512(P, JacCurve);
        eccdouble_jac_Jac512(P, JacCurve);
        eccdouble_jac_Jac512(P, JacCurve);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Point doubling runs in .......................................... %8lld cycles", cycles/(ML_BENCH_LOOPS*10));
    printf("\n"); 
    
    // "Complete" point addition (Weierstrass a=-3)    
    eccset_Jac512(A, JacCurve);
    eccconvert_aff_to_jac_Jac512(A, P);
    ecccopy_jac_Jac512(P, R); 
    eccdouble_jac_Jac512(P, JacCurve);

    cycles = 0;
    for (n=0; n<ML_BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        eccadd_jac_Jac512(R, P, JacCurve);
        eccadd_jac_Jac512(R, P, JacCurve);
        eccadd_jac_Jac512(R, P, JacCurve);
        eccadd_jac_Jac512(R, P, JacCurve);
        eccadd_jac_Jac512(R, P, JacCurve);
        eccadd_jac_Jac512(R, P, JacCurve);
        eccadd_jac_Jac512(R, P, JacCurve);
        eccadd_jac_Jac512(R, P, JacCurve);
        eccadd_jac_Jac512(R, P, JacCurve);
        eccadd_jac_Jac512(R, P, JacCurve);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  (Complete) point addition runs in ............................... %8lld cycles", cycles/(ML_BENCH_LOOPS*10));
    printf("\n"); 
    
    // Variable-base scalar multiplication (Weierstrass a=-3)
    eccset_Jac512(A, JacCurve);
    random512(k, JacCurve->order);        
    
    cycles = 0;
    for (n=0; n<ML_SHORT_BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        ecc_scalar_mul_Jac512(A, k, AA, JacCurve);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Variable-base scalar mul runs in ................................ %8lld cycles", cycles/(ML_SHORT_BENCH_LOOPS));
#ifdef ML_COUNT
    ninv=0; nmul=0; nsqr=0; nadd=0; nlut=0; nclut=0; ncsel=0;
    ecc_scalar_mul_Jac512(A, k, AA, JacCurve);
    printf(" using %d inversions, %d muls, %d sqrs, %d adds, %d luts, %d 5-point luts and %d field selects", ninv, nmul, nsqr, nadd, nlut, nclut, ncsel);
#endif  
    printf("\n"); 
    
    // Fixed-base scalar multiplication (Weierstrass a=-3)
    eccset_Jac512(A, JacCurve); 
    T_fixed = ecc_precomp_fixed_Jac512(A, MEM_LARGE, JacCurve);

    cycles = 0;
    for (n=0; n<ML_SHORT_BENCH_LOOPS; n++)
    {
        random512(k, JacCurve->order);  
        cycles1 = cpucycles();
        ecc_scalar_mul_fixed_Jac512(T_fixed, k, AA, MEM_LARGE, JacCurve);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    ecc_destroy_precomp_Jac512(T_fixed); 

    printf("  Fixed-base scalar mul (memory model=MEM_LARGE) runs in .......... %8lld cycles", cycles/(ML_SHORT_BENCH_LOOPS));
    printf("\n"); 
    
    // Double-scalar multiplication (Weierstrass a=-3)
    eccset_Jac512(A, JacCurve); 
    random512(k, JacCurve->order); 
    ecc_mul_waff_512(A, k, AA, JacCurve);    // Base points are A and AA

    T_fixed = ecc_precomp_dblmul_Jac512(AA, MEM_LARGE, JacCurve);
    cycles = 0;
    for (n=0; n<ML_SHORT_BENCH_LOOPS; n++)
    {
        random512(k, JacCurve->order); 
        random512(kk, JacCurve->order);
        cycles1 = cpucycles();
        ecc_double_scalar_mul_Jac512(T_fixed, k, A, kk, B, MEM_LARGE, JacCurve);  
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    ecc_destroy_precomp_Jac512(T_fixed); 

    printf("  Double-base scalar mul (memory model=MEM_LARGE) runs in ......... %8lld cycles", cycles/(ML_SHORT_BENCH_LOOPS));
    printf("\n"); 
    
    return OK;
}


BOOL ecc_run512_te(PCurveStruct TedCurve)
{ // Benchmarking for curve Ted512
    BOOL OK = TRUE;
    unsigned int n;
    unsigned long long cycles, cycles1, cycles2;
    point_Ted512 A, AA;
    point_extproj_Ted512 P, R;
    dig k[ML_WORDS512], kk[ML_WORDS512];
    point_extaff_precomp_Ted512 *T_fixed = NULL;

    SetThreadAffinityMask(GetCurrentThread(), 1);       // All threads are set to run in the same node
    SetThreadPriority(GetCurrentThread(), 2);           // Set to highest priority
        
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Curve arithmetic: twisted Edwards a=-1 over GF(2^512-569) \n\n"); 

    // Point doubling (twisted Edwards a=-1)
    eccset_Ted512(A, TedCurve);
    eccconvert_aff_to_extproj_Ted512(A, P);

    cycles = 0;
    for (n=0; n<ML_BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        eccdouble_extproj_Ted512(P, TedCurve);
        eccdouble_extproj_Ted512(P, TedCurve);
        eccdouble_extproj_Ted512(P, TedCurve);
        eccdouble_extproj_Ted512(P, TedCurve);
        eccdouble_extproj_Ted512(P, TedCurve);
        eccdouble_extproj_Ted512(P, TedCurve);
        eccdouble_extproj_Ted512(P, TedCurve);
        eccdouble_extproj_Ted512(P, TedCurve);
        eccdouble_extproj_Ted512(P, TedCurve);
        eccdouble_extproj_Ted512(P, TedCurve);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Point doubling runs in .......................................... %8lld cycles", cycles/(ML_BENCH_LOOPS*10));
    printf("\n"); 

    // Point addition (twisted Edwards a=-1)    
    eccset_Ted512(A, TedCurve);
    eccconvert_aff_to_extproj_Ted512(A, P);
    ecccopy_extproj_Ted512(P, R);
    eccdouble_extproj_Ted512(P, TedCurve);

    cycles = 0;
    for (n=0; n<ML_BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        eccadd_extproj_Ted512(R, P, TedCurve);
        eccadd_extproj_Ted512(R, P, TedCurve);
        eccadd_extproj_Ted512(R, P, TedCurve);
        eccadd_extproj_Ted512(R, P, TedCurve);
        eccadd_extproj_Ted512(R, P, TedCurve);
        eccadd_extproj_Ted512(R, P, TedCurve);
        eccadd_extproj_Ted512(R, P, TedCurve);
        eccadd_extproj_Ted512(R, P, TedCurve);
        eccadd_extproj_Ted512(R, P, TedCurve);
        eccadd_extproj_Ted512(R, P, TedCurve);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  (Complete) point addition runs in ............................... %8lld cycles", cycles/(ML_BENCH_LOOPS*10));
    printf("\n"); 

    // Variable-base scalar multiplication (twisted Edwards a=-1)
    eccset_Ted512(A, TedCurve);
    random512(k, TedCurve->order); 
    
    cycles = 0;
    for (n=0; n<ML_SHORT_BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        ecc_scalar_mul_Ted512(A, k, AA, TedCurve);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Variable-base scalar mul runs in ................................ %8lld cycles", cycles/(ML_SHORT_BENCH_LOOPS));
#ifdef ML_COUNT
    ninv=0; nmul=0; nsqr=0; nadd=0; nlut=0;
    ecc_scalar_mul_Ted512(A, k, AA, TedCurve);
    printf(" using %d inversions, %d muls, %d sqrs, %d adds and %d luts", ninv, nmul, nsqr, nadd, nlut);
#endif  
    printf("\n"); 

    // Fixed-base scalar multiplication (twisted Edwards a=-1)    
    eccset_Ted512(A, TedCurve); 
    T_fixed = ecc_precomp_fixed_Ted512(A, MEM_LARGE, TedCurve);

    cycles = 0;
    for (n=0; n<ML_SHORT_BENCH_LOOPS; n++)
    {
        random512(k, TedCurve->order);  
        cycles1 = cpucycles();
        ecc_scalar_mul_fixed_Ted512(T_fixed, k, AA, MEM_LARGE, TedCurve);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    ecc_destroy_precomp_Ted512(T_fixed); 

    printf("  Fixed-base scalar mul (memory model=MEM_LARGE) runs in .......... %8lld cycles", cycles/(ML_SHORT_BENCH_LOOPS));
    printf("\n"); 
    
    // Double-scalar multiplication (twisted Edwards a=-1)
    eccset_Ted512(A, TedCurve); 
    random512(k, TedCurve->order); 
    ecc_scalar_mul_Ted512(A, k, AA, TedCurve);    // Base points are A and AA
            
    T_fixed = ecc_precomp_dblmul_Ted512(AA, MEM_LARGE, TedCurve);
    cycles = 0;
    for (n=0; n<ML_SHORT_BENCH_LOOPS; n++)
    {
        random512(k, TedCurve->order); 
        random512(kk, TedCurve->order);
        cycles1 = cpucycles();
        ecc_double_scalar_mul_Ted512(T_fixed, k, A, kk, AA, MEM_LARGE, TedCurve); 
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    ecc_destroy_precomp_Ted512(T_fixed); 

    printf("  Double-base scalar mul (memory model=MEM_LARGE) runs in ......... %8lld cycles", cycles/(ML_SHORT_BENCH_LOOPS));
    printf("\n"); 

    return OK;
}

#endif



int main()
{
    BOOL OK = TRUE;

    if (is_avx_supported() == FALSE) {
        printf("\n  UNSUPPORTED PLATFORM - support for AVX instructions cannot be detected \n");
        return FALSE;
    }
	
#ifdef ECCURVES_256
    OK = OK && ecc_test256_w(&curve_Jac256);       // Test "Jac256", Weierstrass a=-3 curve with p = 2^256-189
    OK = OK && ecc_test256_te(&curve_Ted256);      // Test "Ted256", twisted Edwards a=-1 curve with p = 2^256-189
#endif	
#ifdef ECCURVES_384
	OK = OK && ecc_test384_w(&curve_Jac384);       // Test "Jac384", Weierstrass a=-3 curve with p = 2^384-317
    OK = OK && ecc_test384_te(&curve_Ted384);      // Test "Ted384", twisted Edwards a=-1 curve with p = 2^384-317
#endif	
#ifdef ECCURVES_512
    OK = OK && ecc_test512_w(&curve_Jac512);       // Test "Jac512", Weierstrass a=-3 curve with p = 2^512-569
    OK = OK && ecc_test512_te(&curve_Ted512);      // Test "Ted512", twisted Edwards a=-1 curve with p = 2^512-569
#endif	
    
#ifdef ECCURVES_256
    OK = OK && ecc_run256_w(&curve_Jac256);        // Benchmark "Jac256", Weierstrass a=-3 curve with p = 2^256-189
    OK = OK && ecc_run256_te(&curve_Ted256);       // Benchmark "Ted256", twisted Edwards a=-1 curve with p = 2^256-189
#endif	
#ifdef ECCURVES_384
    OK = OK && ecc_run384_w(&curve_Jac384);        // Benchmark "Jac384", Weierstrass a=-3 curve with p = 2^384-317
    OK = OK && ecc_run384_te(&curve_Ted384);       // Benchmark "Ted384", twisted Edwards a=-1 curve with p = 2^384-317
#endif	
#ifdef ECCURVES_512
    OK = OK && ecc_run512_w(&curve_Jac512);        // Benchmark "Jac512", Weierstrass a=-3 curve with p = 2^512-569
    OK = OK && ecc_run512_te(&curve_Ted512);       // Benchmark "Ted512", twisted Edwards a=-1 curve with p = 2^512-569
#endif	
    
    return OK;
}