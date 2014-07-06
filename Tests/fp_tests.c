/**************************************************************************
* Suite for benchmarking/testing field operations for MSR ECClib
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
* Abstract: benchmarking/testing field operations
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

BOOL fp_test256(PCurveStruct PCurve)
{ // Tests of field arithmetic over GF(2^256-189)
    BOOL passed, OK = TRUE;
    dig n;
    dig256 a, b, c, d, e, f, p;

    fp_prime256(p, PCurve);
    printf("\n\nTESTING \n"); 
    printf("--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Field arithmetic over GF(2^256-189): \n\n"); 

    // Fp addition with p = 2^256-189
    passed = TRUE;
    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        random256(a, PCurve->prime); random256(b, PCurve->prime); random256(c, PCurve->prime); random256(d, PCurve->prime); random256(e, PCurve->prime); random256(f, PCurve->prime); 

        fpadd256(a, b, d); fpadd256(d, c, e);                 // e = (a+b)+c
        fpadd256(b, c, d); fpadd256(d, a, f);                 // f = a+(b+c)
        if (fpcompare256(e,f)!=0) { passed=FALSE; break; }

        fpadd256(a, b, d);                                     // d = a+b 
        fpadd256(b, a, e);                                     // e = b+a
        if (fpcompare256(d,e)!=0) { passed=FALSE; break; }

        fpzero256(b);
        fpadd256(a, b, d);                                     // d = a+0 
        if (fpcompare256(a,d)!=0) { passed=FALSE; break; }
        
        fpzero256(b);
        fpcopy256(a, d);     
        fpneg256(PCurve->prime, d);                      
        fpadd256(a, d, e);                                     // e = a+(-a)
        if (fpcompare256(e,b)!=0) { passed=FALSE; break; }
    }
    if (passed==TRUE) printf("  Addition tests (associativity, commutativity, identity, inverse)......................... PASSED");
    else { printf("  Addition tests... FAILED"); printf("\n"); return FALSE; }
    printf("\n");

    // Fp subtraction with p = 2^256-189
    passed = TRUE;
    for (n=0; n<ML_TEST_LOOPS; n++)
    {        
        random256(a, PCurve->prime); random256(b, PCurve->prime); random256(c, PCurve->prime); random256(d, PCurve->prime); random256(e, PCurve->prime); 

        fpsub256(a, b, d);                                    // d = a-b 
        fpsub256(b, a, e);                                    // e = b-a     
        fpneg256(PCurve->prime, e);     
        if (fpcompare256(d,e)!=0) { passed=FALSE; break; }

        fpzero256(b);
        fpsub256(a, b, d);                                    // d = a-0 
        if (fpcompare256(a,d)!=0) { passed=FALSE; break; }
                 
        fpsub256(a, a, d);                                    // e = a-(a)
        if (fp_iszero256(d) == FALSE) { passed=FALSE; break; }
    }
    if (passed==TRUE) printf("  Subtraction tests (anti-commutativity, identity, inverse)................................ PASSED");
    else { printf("  Subtraction tests... FAILED"); printf("\n"); return FALSE; }
    printf("\n");

    // Fp division by 2 with p = 2^256-189
    passed = TRUE;
    for (n=0; n<ML_TEST_LOOPS; n++)
    {        
        random256(a, PCurve->prime); random256(b, PCurve->prime); random256(c, PCurve->prime); random256(d, PCurve->prime); random256(e, PCurve->prime); 

        fpdiv2_256(a, c);                                  // c = a/2
        fpadd256(c, c, b);                                 // b = a 
        if (fpcompare256(a,b)!=0) { passed=FALSE; break; }

        fpdiv2_256(a, c);                                  // c = a/2
        fpzero256(b); b[0] = 2;
        fpmul256(c, b, d);                                 // d = a 
        if (fpcompare256(a,d)!=0) { passed=FALSE; break; }

        fpzero256(b);
        fpdiv2_256(b, c);                                  // 0 
        if (fpcompare256(c,b)!=0) { passed=FALSE; break; }
    }
    if (passed==TRUE) printf("  Division by 2 tests ..................................................................... PASSED");
    else { printf("  Division by 2 tests... FAILED"); printf("\n"); return FALSE; }
    printf("\n");

    // Fp negation with p = 2^256-189
    passed = TRUE;
    for (n=0; n<ML_TEST_LOOPS; n++)
    {        
        random256(a, PCurve->prime); random256(b, PCurve->prime); random256(c, PCurve->prime); random256(d, PCurve->prime); random256(e, PCurve->prime); 

        fpcopy256(a, c);
        fpneg256(PCurve->prime, a);                        // -a 
        fpneg256(PCurve->prime, a);                        // -(-a) 
        if (fpcompare256(a,c)!=0) { passed=FALSE; break; }

        fpsub256(a, b, c);                                 // c = a-b 
        fpneg256(PCurve->prime, b);                        // -b 
        fpadd256(a, b, d);                                 // d = a+(-b) 
        if (fpcompare256(c,d)!=0) { passed=FALSE; break; }

        fpzero256(b);
        fpneg256(PCurve->prime, b);                        // -0 
        if (fpcompare256(p,b)!=0) { passed=FALSE; break; }
    }
    if (passed==TRUE) printf("  Negation tests .......................................................................... PASSED");
    else { printf("  Negation tests... FAILED"); printf("\n"); return FALSE; }
    printf("\n");

    // Fp multiplication with p = 2^256-189
    passed = TRUE;
    for (n=0; n<ML_TEST_LOOPS; n++)
    {    
        random256(a, PCurve->prime); random256(b, PCurve->prime); random256(c, PCurve->prime); random256(d, PCurve->prime); random256(e, PCurve->prime); random256(f, PCurve->prime);

        fpmul256(a, b, d); fpmul256(d, c, e);                        // e = (a*b)*c
        fpmul256(b, c, d); fpmul256(d, a, f);                        // f = a*(b*c)
        if (fpcompare256(e,f)!=0) { passed=FALSE; break; }

        fpadd256(b, c, d); fpmul256(a, d, e);                        // e = a*(b+c)
        fpmul256(a, b, d); fpmul256(a, c, f); fpadd256(d, f, f);     // f = a*b+a*c
        if (fpcompare256(e,f)!=0) { passed=FALSE; break; }

        fpmul256(a, b, d);                                           // d = a*b 
        fpmul256(b, a, e);                                           // e = b*a 
        if (fpcompare256(d,e)!=0) { passed=FALSE; break; }

        fpzero256(b); b[0]=1; 
        fpmul256(a, b, d);                                           // d = a*1 
        if (fpcompare256(a,d)!=0) { passed=FALSE; break; }

        fpzero256(b);
        fpmul256(a, b, d);                                           // d = a*0 
        if (fpcompare256(b,d)!=0) { passed=FALSE; break; }
    }
    if (passed==TRUE) printf("  Multiplication tests (associativity, distributive, commutativity, identity, null)........ PASSED");
    else { printf("  Multiplication tests... FAILED"); printf("\n"); return FALSE; }
    printf("\n");

    // Fp squaring with p = 2^256-189
    passed = TRUE;
    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        random256(a, PCurve->prime); random256(b, PCurve->prime); random256(c, PCurve->prime); 

        fpsqr256(a, b);                                            // b = a^2
        fpmul256(a, a, c);                                         // c = a*a 
        if (fpcompare256(b,c)!=0) { passed=FALSE; break; }

        fpzero256(a);
        fpsqr256(a, d);                                            // d = 0^2 
        if (fpcompare256(a,d)!=0) { passed=FALSE; break; }
    }
    if (passed==TRUE) printf("  Squaring tests........................................................................... PASSED");
    else { printf("  Squaring tests... FAILED"); printf("\n"); return FALSE; }
    printf("\n");

    // Fp inversion with p = 2^256-189
    passed = TRUE;
    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        random256(a, PCurve->prime); random256(b, PCurve->prime); random256(c, PCurve->prime); random256(d, PCurve->prime);   

        fpzero256(d); d[0]=1; 
        fpcopy256(a, b);                            
        fpinv256(a);                                
        fpmul256(a, b, c);                                      // c = a*a^-1 
        if (fpcompare256(c,d)!=0) { passed=FALSE; break; }
    }
    if (passed==TRUE) printf("  Inversion tests.......................................................................... PASSED");
    else { printf("  Inversion tests... FAILED"); printf("\n"); return FALSE; }
    printf("\n");

    return OK;
}

#endif


#ifdef ECCURVES_384

BOOL fp_test384(PCurveStruct PCurve)
{ // Tests of field arithmetic over GF(2^384-317)
    BOOL passed, OK = TRUE;
    dig n;
    dig384 a, b, c, d, e, f, p;

    fp_prime384(p, PCurve);
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Field arithmetic over GF(2^384-317): \n\n"); 

    // Fp addition with p = 2^384-317
    passed = TRUE;
    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        random384(a, PCurve->prime); random384(b, PCurve->prime); random384(c, PCurve->prime); random384(d, PCurve->prime); random384(e, PCurve->prime); random384(f, PCurve->prime); 

        fpadd384(a, b, d); fpadd384(d, c, e);                  // e = (a+b)+c
        fpadd384(b, c, d); fpadd384(d, a, f);                  // f = a+(b+c)
        if (fpcompare384(e,f)!=0) { passed=FALSE; break; }

        fpadd384(a, b, d);                                     // d = a+b 
        fpadd384(b, a, e);                                     // e = b+a
        if (fpcompare384(d,e)!=0) { passed=FALSE; break; }

        fpzero384(b);
        fpadd384(a, b, d);                                     // d = a+0 
        if (fpcompare384(a,d)!=0) { passed=FALSE; break; }
        
        fpzero384(b);
        fpcopy384(a, d);     
        fpneg384(PCurve->prime, d);                      
        fpadd384(a, d, e);                                     // e = a+(-a)
        if (fpcompare384(e,b)!=0) { passed=FALSE; break; }
    }
    if (passed==TRUE) printf("  Addition tests (associativity, commutativity, identity, inverse)......................... PASSED");
    else { printf("Addition tests... FAILED"); printf("\n"); return FALSE; }
    printf("\n");

    // Fp subtraction with p = 2^384-317
    passed = TRUE;
    for (n=0; n<ML_TEST_LOOPS; n++)
    {        
        random384(a, PCurve->prime); random384(b, PCurve->prime); random384(c, PCurve->prime); random384(d, PCurve->prime); random384(e, PCurve->prime); 

        fpsub384(a, b, d);                                    // d = a-b 
        fpsub384(b, a, e);                                    // e = b-a     
        fpneg384(PCurve->prime, e);     
        if (fpcompare384(d,e)!=0) { passed=FALSE; break; }

        fpzero384(b);
        fpsub384(a, b, d);                                    // d = a-0 
        if (fpcompare384(a,d)!=0) { passed=FALSE; break; }
                 
        fpsub384(a, a, d);                                    // e = a-(a)
        if (fp_iszero384(d) == FALSE) { passed=FALSE; break; }
    }
    if (passed==TRUE) printf("  Subtraction tests (anti-commutativity, identity, inverse)................................ PASSED");
    else { printf("  Subtraction tests... FAILED"); printf("\n"); return FALSE; }
    printf("\n");

    // Fp division by 2 with p = 2^384-317
    passed = TRUE;
    for (n=0; n<ML_TEST_LOOPS; n++)
    {        
        random384(a, PCurve->prime); random384(b, PCurve->prime); random384(c, PCurve->prime); random384(d, PCurve->prime); random384(e, PCurve->prime); 

        fpdiv2_384(a, c);                                  // c = a/2
        fpadd384(c, c, b);                                 // b = a 
        if (fpcompare384(a,b)!=0) { passed=FALSE; break; }
        
        fpdiv2_384(a, c);                                  // c = a/2
        fpzero384(b); b[0] = 2;
        fpmul384(c, b, d);                                 // d = a 
        if (fpcompare384(a,d)!=0) { passed=FALSE; break; }
        
        fpzero384(b);
        fpdiv2_384(b, c);                                  // 0 
        if (fpcompare384(c,b)!=0) { passed=FALSE; break; }
    }
    if (passed==TRUE) printf("  Division by 2 tests ..................................................................... PASSED");
    else { printf("  Division by 2 tests... FAILED"); printf("\n"); return FALSE; }
    printf("\n");

    // Fp negation with p = 2^384-317
    passed = TRUE;
    for (n=0; n<ML_TEST_LOOPS; n++)
    {        
        random384(a, PCurve->prime); random384(b, PCurve->prime); random384(c, PCurve->prime); random384(d, PCurve->prime); random384(e, PCurve->prime); 

        fpcopy384(a, c);
        fpneg384(PCurve->prime, a);                        // -a 
        fpneg384(PCurve->prime, a);                        // -(-a) 
        if (fpcompare384(a,c)!=0) { passed=FALSE; break; }

        fpsub384(a, b, c);                                 // c = a-b 
        fpneg384(PCurve->prime, b);                        // -b 
        fpadd384(a, b, d);                                 // d = a+(-b) 
        if (fpcompare384(c,d)!=0) { passed=FALSE; break; }

        fpzero384(b);
        fpneg384(PCurve->prime, b);                        // -0 
        if (fpcompare384(p,b)!=0) { passed=FALSE; break; }
    }
    if (passed==TRUE) printf("  Negation tests .......................................................................... PASSED");
    else { printf("  Negation tests... FAILED"); printf("\n"); return FALSE; }
    printf("\n");
    
    // Fp multiplication with p = 2^384-317
    passed = TRUE;
    for (n=0; n<ML_TEST_LOOPS; n++)
    {    
        random384(a, PCurve->prime); random384(b, PCurve->prime); random384(c, PCurve->prime); random384(d, PCurve->prime); random384(e, PCurve->prime); random384(f, PCurve->prime);

        fpmul384(a, b, d); fpmul384(d, c, e);                     // e = (a*b)*c
        fpmul384(b, c, d); fpmul384(d, a, f);                     // f = a*(b*c)
        if (fpcompare384(e,f)!=0) { passed=FALSE; break; }

        fpadd384(b, c, d); fpmul384(a, d, e);                     // e = a*(b+c)
        fpmul384(a, b, d); fpmul384(a, c, f); fpadd384(d, f, f);  // f = a*b+a*c
        if (fpcompare384(e,f)!=0) { passed=FALSE; break; }

        fpmul384(a, b, d);                                        // d = a*b 
        fpmul384(b, a, e);                                        // e = b*a 
        if (fpcompare384(d,e)!=0) { passed=FALSE; break; }
        
        fpzero384(b); b[0] = 1;
        fpmul384(a, b, d);                                        // d = a*1 
        if (fpcompare384(a,d)!=0) { passed=FALSE; break; }

        fpzero384(b);
        fpmul384(a, b, d);                                        // d = a*0 
        if (fpcompare384(b,d)!=0) { passed=FALSE; break; }
    }
    if (passed==TRUE) printf("  Multiplication tests (associativity, distributive, commutativity, identity, null)........ PASSED");
    else { printf("  Multiplication tests... FAILED"); printf("\n"); return FALSE; }
    printf("\n");

    // Fp squaring with p = 2^384-317
    passed = TRUE;
    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        random384(a, PCurve->prime); random384(b, PCurve->prime); random384(c, PCurve->prime); 

        fpsqr384(a, b);                                            // b = a^2
        fpmul384(a, a, c);                                         // c = a*a 
        if (fpcompare384(b,c)!=0) { passed=FALSE; break; }
        
        fpzero384(a);
        fpsqr384(a, d);                                            // d = 0^2 
        if (fpcompare384(a,d)!=0) { passed=FALSE; break; }
    }
    if (passed==TRUE) printf("  Squaring tests........................................................................... PASSED");
    else { printf("  Squaring tests... FAILED"); printf("\n"); return FALSE; }
    printf("\n");
    
    // Fp inversion with p = 2^384-317
    passed = TRUE;
    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        random384(a, PCurve->prime); random384(b, PCurve->prime); random384(c, PCurve->prime); random384(d, PCurve->prime);   
                
        fpzero384(d); d[0]=1;  
        fpcopy384(a, b);                            
        fpinv384(a);                                
        fpmul384(a, b, c);                                        // c = a*a^-1 
        if (fpcompare384(c,d)!=0) { passed=FALSE; break; }
    }
    if (passed==TRUE) printf("  Inversion tests.......................................................................... PASSED");
    else { printf("  Inversion tests... FAILED"); printf("\n"); return FALSE; }
    printf("\n");
    
    return OK;
}

#endif


#ifdef ECCURVES_512

BOOL fp_test512(PCurveStruct PCurve)
{ // Tests of field arithmetic over GF(2^512-569)
    BOOL passed, OK = TRUE;
    dig n;
    dig512 a, b, c, d, e, f, p;

    fp_prime512(p, PCurve);
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Field arithmetic over GF(2^512-569): \n\n"); 

    // Fp addition with p = 2^512-569
    passed = TRUE;
    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        random512(a, PCurve->prime); random512(b, PCurve->prime); random512(c, PCurve->prime); random512(d, PCurve->prime); random512(e, PCurve->prime); random512(f, PCurve->prime); 

        fpadd512(a, b, d); fpadd512(d, c, e);                  // e = (a+b)+c
        fpadd512(b, c, d); fpadd512(d, a, f);                  // f = a+(b+c)
        if (fpcompare512(e,f)!=0) { passed=FALSE; break; }

        fpadd512(a, b, d);                                     // d = a+b 
        fpadd512(b, a, e);                                     // e = b+a
        if (fpcompare512(d,e)!=0) { passed=FALSE; break; }

        fpzero512(b);
        fpadd512(a, b, d);                                     // d = a+0 
        if (fpcompare512(a,d)!=0) { passed=FALSE; break; }
        
        fpzero512(b);
        fpcopy512(a, d);     
        fpneg512(PCurve->prime, d);                      
        fpadd512(a, d, e);                                     // e = a+(-a)
        if (fpcompare512(e,b)!=0) { passed=FALSE; break; }
    }
    if (passed==TRUE) printf("  Addition tests (associativity, commutativity, identity, inverse)......................... PASSED");
    else { printf("Addition tests... FAILED"); printf("\n"); return FALSE; }
    printf("\n");

    // Fp subtraction with p = 2^512-569
    passed = TRUE;
    for (n=0; n<ML_TEST_LOOPS; n++)
    {        
        random512(a, PCurve->prime); random512(b, PCurve->prime); random512(c, PCurve->prime); random512(d, PCurve->prime); random512(e, PCurve->prime); 

        fpsub512(a, b, d);                                    // d = a-b 
        fpsub512(b, a, e);                                    // e = b-a     
        fpneg512(PCurve->prime, e);     
        if (fpcompare512(d,e)!=0) { passed=FALSE; break; }

        fpzero512(b);
        fpsub512(a, b, d);                                    // d = a-0 
        if (fpcompare512(a,d)!=0) { passed=FALSE; break; }
                 
        fpsub512(a, a, d);                                    // e = a-(a)
        if (fp_iszero512(d) == FALSE) { passed=FALSE; break; }
    }
    if (passed==TRUE) printf("  Subtraction tests (anti-commutativity, identity, inverse)................................ PASSED");
    else { printf("  Subtraction tests... FAILED"); printf("\n"); return FALSE; }
    printf("\n");

    // Fp division by 2 with p = 2^512-569
    passed = TRUE;
    for (n=0; n<ML_TEST_LOOPS; n++)
    {        
        random512(a, PCurve->prime); random512(b, PCurve->prime); random512(c, PCurve->prime); random512(d, PCurve->prime); random512(e, PCurve->prime); 

        fpdiv2_512(a, c);                                  // c = a/2
        fpadd512(c, c, b);                                 // b = a 
        if (fpcompare512(a,b)!=0) { passed=FALSE; break; }
        
        fpdiv2_512(a, c);                                  // c = a/2
        fpzero512(b); b[0] = 2;
        fpmul512(c, b, d);                                 // d = a 
        if (fpcompare512(a,d)!=0) { passed=FALSE; break; }
        
        fpzero512(b);
        fpdiv2_512(b, c);                                  // 0 
        if (fpcompare512(c,b)!=0) { passed=FALSE; break; }
    }
    if (passed==TRUE) printf("  Division by 2 tests ..................................................................... PASSED");
    else { printf("  Division by 2 tests... FAILED"); printf("\n"); return FALSE; }
    printf("\n");

    // Fp negation with p = 2^512-569
    passed = TRUE;
    for (n=0; n<ML_TEST_LOOPS; n++)
    {        
        random512(a, PCurve->prime); random512(b, PCurve->prime); random512(c, PCurve->prime); random512(d, PCurve->prime); random512(e, PCurve->prime); 

        fpcopy512(a, c);
        fpneg512(PCurve->prime, a);                        // -a 
        fpneg512(PCurve->prime, a);                        // -(-a) 
        if (fpcompare512(a,c)!=0) { passed=FALSE; break; }

        fpsub512(a, b, c);                                 // c = a-b 
        fpneg512(PCurve->prime, b);                        // -b 
        fpadd512(a, b, d);                                 // d = a+(-b) 
        if (fpcompare512(c,d)!=0) { passed=FALSE; break; }

        fpzero512(b);
        fpneg512(PCurve->prime, b);                        // -0 
        if (fpcompare512(p,b)!=0) { passed=FALSE; break; }
    }
    if (passed==TRUE) printf("  Negation tests .......................................................................... PASSED");
    else { printf("  Negation tests... FAILED"); printf("\n"); return FALSE; }
    printf("\n");
    
    // Fp multiplication with p = 2^512-569
    passed = TRUE;
    for (n=0; n<ML_TEST_LOOPS; n++)
    {    
        random512(a, PCurve->prime); random512(b, PCurve->prime); random512(c, PCurve->prime); random512(d, PCurve->prime); random512(e, PCurve->prime); random512(f, PCurve->prime);

        fpmul512(a, b, d); fpmul512(d, c, e);                     // e = (a*b)*c
        fpmul512(b, c, d); fpmul512(d, a, f);                     // f = a*(b*c)
        if (fpcompare512(e,f)!=0) { passed=FALSE; break; }

        fpadd512(b, c, d); fpmul512(a, d, e);                     // e = a*(b+c)
        fpmul512(a, b, d); fpmul512(a, c, f); fpadd512(d, f, f);  // f = a*b+a*c
        if (fpcompare512(e,f)!=0) { passed=FALSE; break; }

        fpmul512(a, b, d);                                        // d = a*b 
        fpmul512(b, a, e);                                        // e = b*a 
        if (fpcompare512(d,e)!=0) { passed=FALSE; break; }
        
        fpzero512(b); b[0] = 1;
        fpmul512(a, b, d);                                        // d = a*1 
        if (fpcompare512(a,d)!=0) { passed=FALSE; break; }

        fpzero512(b);
        fpmul512(a, b, d);                                        // d = a*0 
        if (fpcompare512(b,d)!=0) { passed=FALSE; break; }
    }
    if (passed==TRUE) printf("  Multiplication tests (associativity, distributive, commutativity, identity, null)........ PASSED");
    else { printf("  Multiplication tests... FAILED"); printf("\n"); return FALSE; }
    printf("\n");
    
    // Fp squaring with p = 2^512-569
    passed = TRUE;
    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        random512(a, PCurve->prime); random512(b, PCurve->prime); random512(c, PCurve->prime); 

        fpsqr512(a, b);                                            // b = a^2
        fpmul512(a, a, c);                                         // c = a*a 
        if (fpcompare512(b,c)!=0) { passed=FALSE; break; }
        
        fpzero512(a);
        fpsqr512(a, d);                                            // d = 0^2 
        if (fpcompare512(a,d)!=0) { passed=FALSE; break; }
    }
    if (passed==TRUE) printf("  Squaring tests........................................................................... PASSED");
    else { printf("  Squaring tests... FAILED"); printf("\n"); return FALSE; }
    printf("\n");
        
    // Fp inversion with p = 2^512-569
    passed = TRUE;
    for (n=0; n<ML_TEST_LOOPS; n++)
    {
        random512(a, PCurve->prime); random512(b, PCurve->prime); random512(c, PCurve->prime); random512(d, PCurve->prime);   
                
        fpzero512(d); d[0]=1;  
        fpcopy512(a, b);                            
        fpinv512(a);                                
        fpmul512(a, b, c);                                        // c = a*a^-1 
        if (fpcompare512(c,d)!=0) { passed=FALSE; break; }
    }
    if (passed==TRUE) printf("  Inversion tests.......................................................................... PASSED");
    else { printf("  Inversion tests... FAILED"); printf("\n"); return FALSE; }
    printf("\n");
    
    return OK;
}

#endif


#ifdef ECCURVES_256

BOOL fp_run256(PCurveStruct PCurve)
{ // Benchmarking of field arithmetic over GF(2^256-189)
    BOOL OK = TRUE;
    dig n;
    unsigned long long cycles, cycles1, cycles2;
    dig256 a, b, c, d, e, f, p;

    SetThreadAffinityMask(GetCurrentThread(), 1);       // All threads are set to run in the same node
    SetThreadPriority(GetCurrentThread(), 2);           // Set to highest priority
        
    fp_prime256(p, PCurve);
    printf("\n\nBENCHMARKING \n"); 
    printf("--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Field arithmetic over GF(2^256-189): \n\n"); 

    // Fp addition with p = 2^256-189
    cycles = 0;
    for (n=0; n<ML_BENCH_LOOPS; n++)
    {
        random256(a, PCurve->prime); random256(b, PCurve->prime); random256(c, PCurve->prime); random256(d, PCurve->prime); random256(e, PCurve->prime); random256(f, PCurve->prime); 

        cycles1 = cpucycles();
        fpadd256(a, b, c);
        fpadd256(d, e, f);
        fpadd256(b, c, d);
        fpadd256(e, f, a);
        fpadd256(c, d, e);
        fpadd256(f, a, b);
        fpadd256(d, e, f);
        fpadd256(a, b, c);
        fpadd256(e, f, a);
        fpadd256(b, c, d);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Addition runs in ................... %lld cycles", cycles/(ML_BENCH_LOOPS*10));
    printf("\n"); 

    // Fp subtraction with p = 2^256-189
    cycles = 0;
    for (n=0; n<ML_BENCH_LOOPS; n++)
    {
        random256(a, PCurve->prime); random256(b, PCurve->prime); random256(c, PCurve->prime); random256(d, PCurve->prime); random256(e, PCurve->prime); random256(f, PCurve->prime); 

        cycles1 = cpucycles();
        fpsub256(a, b, c);
        fpsub256(d, e, f);
        fpsub256(b, c, d);
        fpsub256(e, f, a);
        fpsub256(c, d, e);
        fpsub256(f, a, b);
        fpsub256(d, e, f);
        fpsub256(a, b, c);
        fpsub256(e, f, a);
        fpsub256(b, c, d);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Subtraction runs in ................ %lld cycles", cycles/(ML_BENCH_LOOPS*10));
    printf("\n");

    // Fp division by 2 with p = 2^256-189
    cycles = 0;
    for (n=0; n<ML_BENCH_LOOPS; n++)
    {
        random256(a, PCurve->prime); random256(b, PCurve->prime); random256(c, PCurve->prime); random256(d, PCurve->prime); random256(e, PCurve->prime); random256(f, PCurve->prime); 

        cycles1 = cpucycles();
        fpdiv2_256(a, b);
        fpdiv2_256(c, d);
        fpdiv2_256(e, f);
        fpdiv2_256(b, c);
        fpdiv2_256(d, e);
        fpdiv2_256(f, a);
        fpdiv2_256(c, d);
        fpdiv2_256(e, f);
        fpdiv2_256(a, b);
        fpdiv2_256(c, d);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Division by 2 runs in .............. %lld cycles", cycles/(ML_BENCH_LOOPS*10));
    printf("\n");

    // Fp negation with p = 2^256-189
    cycles = 0;
    for (n=0; n<ML_BENCH_LOOPS; n++)
    {
        random256(a, PCurve->prime); random256(b, PCurve->prime); random256(c, PCurve->prime); random256(d, PCurve->prime); random256(e, PCurve->prime); random256(f, PCurve->prime); 

        cycles1 = cpucycles();
        fpneg256(PCurve->prime, a);
        fpneg256(PCurve->prime, b);
        fpneg256(PCurve->prime, c);
        fpneg256(PCurve->prime, d);
        fpneg256(PCurve->prime, e);
        fpneg256(PCurve->prime, f);
        fpneg256(PCurve->prime, a);
        fpneg256(PCurve->prime, b);
        fpneg256(PCurve->prime, c);
        fpneg256(PCurve->prime, d);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Negation runs in ................... %lld cycles", cycles/(ML_BENCH_LOOPS*10));
    printf("\n");

    // Fp squaring with p = 2^256-189
    cycles = 0;
    for (n=0; n<ML_BENCH_LOOPS; n++)
    {
        random256(a, PCurve->prime); random256(b, PCurve->prime); random256(c, PCurve->prime); random256(d, PCurve->prime); random256(e, PCurve->prime); random256(f, PCurve->prime); 

        cycles1 = cpucycles();
        fpsqr256(a, b);
        fpsqr256(c, d);
        fpsqr256(e, f);
        fpsqr256(b, c);
        fpsqr256(d, e);
        fpsqr256(f, a);
        fpsqr256(c, d);
        fpsqr256(e, f);
        fpsqr256(a, c);
        fpsqr256(d, e);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Squaring runs in ................... %lld cycles", cycles/(ML_BENCH_LOOPS*10));
    printf("\n");

    // Fp multiplication with p = 2^256-189
    cycles = 0;
    for (n=0; n<ML_BENCH_LOOPS; n++)
    {
        random256(a, PCurve->prime); random256(b, PCurve->prime); random256(c, PCurve->prime); random256(d, PCurve->prime); random256(e, PCurve->prime); random256(f, PCurve->prime); 

        cycles1 = cpucycles();
        fpmul256(a, b, c);
        fpmul256(d, e, f);
        fpmul256(b, c, d);
        fpmul256(e, f, a);
        fpmul256(c, d, e);
        fpmul256(f, a, b);
        fpmul256(d, e, f);
        fpmul256(a, b, c);
        fpmul256(e, f, a);
        fpmul256(b, c, d);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Multiplication runs in ............. %lld cycles", cycles/(ML_BENCH_LOOPS*10));
    printf("\n");

    // Fp inversion with p = 2^256-189
    cycles = 0;
    for (n=0; n<ML_SHORT_BENCH_LOOPS; n++)
    {
        random256(a, PCurve->prime); random256(b, PCurve->prime); random256(c, PCurve->prime); random256(d, PCurve->prime); random256(e, PCurve->prime); random256(f, PCurve->prime); 

        cycles1 = cpucycles();
        fpinv256(a);
        fpinv256(b);
        fpinv256(c);
        fpinv256(d);
        fpinv256(e);
        fpinv256(f);
        fpinv256(a);
        fpinv256(b);
        fpinv256(c);
        fpinv256(d);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Inversion runs in .................. %lld cycles", cycles/(ML_SHORT_BENCH_LOOPS*10));
    printf("\n");

    return OK;

}

#endif


#ifdef ECCURVES_384

BOOL fp_run384(PCurveStruct PCurve)
{ // Benchmarking of field arithmetic over GF(2^384-317)
    BOOL OK = TRUE;
    dig n;
    unsigned long long cycles, cycles1, cycles2;
    dig384 a, b, c, d, e, f, p;

    SetThreadAffinityMask(GetCurrentThread(), 1);       // All threads are set to run in the same node
    SetThreadPriority(GetCurrentThread(), 2);           // Set to highest priority
        
    fp_prime384(p, PCurve);
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Field arithmetic over GF(2^384-317): \n\n"); 

    // Fp addition with p = 2^384-317
    cycles = 0;
    for (n=0; n<ML_BENCH_LOOPS; n++)
    {
        random384(a, PCurve->prime); random384(b, PCurve->prime); random384(c, PCurve->prime); random384(d, PCurve->prime); random384(e, PCurve->prime); random384(f, PCurve->prime); 

        cycles1 = cpucycles();
        fpadd384(a, b, c);
        fpadd384(d, e, f);
        fpadd384(b, c, d);
        fpadd384(e, f, a);
        fpadd384(c, d, e);
        fpadd384(f, a, b);
        fpadd384(d, e, f);
        fpadd384(a, b, c);
        fpadd384(e, f, a);
        fpadd384(b, c, d);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Addition runs in ................... %lld cycles", cycles/(ML_BENCH_LOOPS*10));
    printf("\n"); 

    // Fp subtraction with p = 2^384-317
    cycles = 0;
    for (n=0; n<ML_BENCH_LOOPS; n++)
    {
        random384(a, PCurve->prime); random384(b, PCurve->prime); random384(c, PCurve->prime); random384(d, PCurve->prime); random384(e, PCurve->prime); random384(f, PCurve->prime); 

        cycles1 = cpucycles();
        fpsub384(a, b, c);
        fpsub384(d, e, f);
        fpsub384(b, c, d);
        fpsub384(e, f, a);
        fpsub384(c, d, e);
        fpsub384(f, a, b);
        fpsub384(d, e, f);
        fpsub384(a, b, c);
        fpsub384(e, f, a);
        fpsub384(b, c, d);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Subtraction runs in ................ %lld cycles", cycles/(ML_BENCH_LOOPS*10));
    printf("\n");

    // Fp division by 2 with p = 2^384-317
    cycles = 0;
    for (n=0; n<ML_BENCH_LOOPS; n++)
    {
        random384(a, PCurve->prime); random384(b, PCurve->prime); random384(c, PCurve->prime); random384(d, PCurve->prime); random384(e, PCurve->prime); random384(f, PCurve->prime); 

        cycles1 = cpucycles();
        fpdiv2_384(a, b);
        fpdiv2_384(c, d);
        fpdiv2_384(e, f);
        fpdiv2_384(b, c);
        fpdiv2_384(d, e);
        fpdiv2_384(f, a);
        fpdiv2_384(c, d);
        fpdiv2_384(e, f);
        fpdiv2_384(a, b);
        fpdiv2_384(c, d);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Division by 2 runs in .............. %lld cycles", cycles/(ML_BENCH_LOOPS*10));
    printf("\n");

    // Fp negation with p = 2^384-317
    cycles = 0;
    for (n=0; n<ML_BENCH_LOOPS; n++)
    {
        random384(a, PCurve->prime); random384(b, PCurve->prime); random384(c, PCurve->prime); random384(d, PCurve->prime); random384(e, PCurve->prime); random384(f, PCurve->prime); 

        cycles1 = cpucycles();
        fpneg384(PCurve->prime, a);
        fpneg384(PCurve->prime, b);
        fpneg384(PCurve->prime, c);
        fpneg384(PCurve->prime, d);
        fpneg384(PCurve->prime, e);
        fpneg384(PCurve->prime, f);
        fpneg384(PCurve->prime, a);
        fpneg384(PCurve->prime, b);
        fpneg384(PCurve->prime, c);
        fpneg384(PCurve->prime, d);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Negation runs in ................... %lld cycles", cycles/(ML_BENCH_LOOPS*10));
    printf("\n");
    
    // Fp squaring with p = 2^384-317
    cycles = 0;
    for (n=0; n<ML_BENCH_LOOPS; n++)
    {
        random384(a, PCurve->prime); random384(b, PCurve->prime); random384(c, PCurve->prime); random384(d, PCurve->prime); random384(e, PCurve->prime); random384(f, PCurve->prime); 

        cycles1 = cpucycles();
        fpsqr384(a, b);
        fpsqr384(c, d);
        fpsqr384(e, f);
        fpsqr384(b, c);
        fpsqr384(d, e);
        fpsqr384(f, a);
        fpsqr384(c, d);
        fpsqr384(e, f);
        fpsqr384(a, c);
        fpsqr384(d, e);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Squaring runs in ................... %lld cycles", cycles/(ML_BENCH_LOOPS*10));
    printf("\n");

    // Fp multiplication with p = 2^384-317
    cycles = 0;
    for (n=0; n<ML_BENCH_LOOPS; n++)
    {
        random384(a, PCurve->prime); random384(b, PCurve->prime); random384(c, PCurve->prime); random384(d, PCurve->prime); random384(e, PCurve->prime); random384(f, PCurve->prime); 

        cycles1 = cpucycles();
        fpmul384(a, b, c);
        fpmul384(d, e, f);
        fpmul384(b, c, d);
        fpmul384(e, f, a);
        fpmul384(c, d, e);
        fpmul384(f, a, b);
        fpmul384(d, e, f);
        fpmul384(a, b, c);
        fpmul384(e, f, a);
        fpmul384(b, c, d);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Multiplication runs in ............. %lld cycles", cycles/(ML_BENCH_LOOPS*10));
    printf("\n");

    // Fp inversion with p = 2^384-317
    cycles = 0;
    for (n=0; n<ML_SHORT_BENCH_LOOPS; n++)
    {
        random384(a, PCurve->prime); random384(b, PCurve->prime); random384(c, PCurve->prime); random384(d, PCurve->prime); random384(e, PCurve->prime); random384(f, PCurve->prime); 

        cycles1 = cpucycles();
        fpinv384(a);
        fpinv384(b);
        fpinv384(c);
        fpinv384(d);
        fpinv384(e);
        fpinv384(f);
        fpinv384(a);
        fpinv384(b);
        fpinv384(c);
        fpinv384(d);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Inversion runs in .................. %lld cycles", cycles/(ML_SHORT_BENCH_LOOPS*10));
    printf("\n");
    
    return OK;
}

#endif


#ifdef ECCURVES_512

BOOL fp_run512(PCurveStruct PCurve)
{ // Benchmarking of field arithmetic over GF(2^512-569)
    BOOL OK = TRUE;
    dig n;
    unsigned long long cycles, cycles1, cycles2;
    dig512 a, b, c, d, e, f, p;

    SetThreadAffinityMask(GetCurrentThread(), 1);       // All threads are set to run in the same node
    SetThreadPriority(GetCurrentThread(), 2);           // Set to highest priority
        
    fp_prime512(p, PCurve);
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Field arithmetic over GF(2^512-569): \n\n"); 

    // Fp addition with p = 2^512-569
    cycles = 0;
    for (n=0; n<ML_BENCH_LOOPS; n++)
    {
        random512(a, PCurve->prime); random512(b, PCurve->prime); random512(c, PCurve->prime); random512(d, PCurve->prime); random512(e, PCurve->prime); random512(f, PCurve->prime); 

        cycles1 = cpucycles();
        fpadd512(a, b, c);
        fpadd512(d, e, f);
        fpadd512(b, c, d);
        fpadd512(e, f, a);
        fpadd512(c, d, e);
        fpadd512(f, a, b);
        fpadd512(d, e, f);
        fpadd512(a, b, c);
        fpadd512(e, f, a);
        fpadd512(b, c, d);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Addition runs in ................... %lld cycles", cycles/(ML_BENCH_LOOPS*10));
    printf("\n"); 

    // Fp subtraction with p = 2^512-569
    cycles = 0;
    for (n=0; n<ML_BENCH_LOOPS; n++)
    {
        random512(a, PCurve->prime); random512(b, PCurve->prime); random512(c, PCurve->prime); random512(d, PCurve->prime); random512(e, PCurve->prime); random512(f, PCurve->prime); 

        cycles1 = cpucycles();
        fpsub512(a, b, c);
        fpsub512(d, e, f);
        fpsub512(b, c, d);
        fpsub512(e, f, a);
        fpsub512(c, d, e);
        fpsub512(f, a, b);
        fpsub512(d, e, f);
        fpsub512(a, b, c);
        fpsub512(e, f, a);
        fpsub512(b, c, d);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Subtraction runs in ................ %lld cycles", cycles/(ML_BENCH_LOOPS*10));
    printf("\n");

    // Fp division by 2 with p = 2^512-569
    cycles = 0;
    for (n=0; n<ML_BENCH_LOOPS; n++)
    {
        random512(a, PCurve->prime); random512(b, PCurve->prime); random512(c, PCurve->prime); random512(d, PCurve->prime); random512(e, PCurve->prime); random512(f, PCurve->prime); 

        cycles1 = cpucycles();
        fpdiv2_512(a, b);
        fpdiv2_512(c, d);
        fpdiv2_512(e, f);
        fpdiv2_512(b, c);
        fpdiv2_512(d, e);
        fpdiv2_512(f, a);
        fpdiv2_512(c, d);
        fpdiv2_512(e, f);
        fpdiv2_512(a, b);
        fpdiv2_512(c, d);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Division by 2 runs in .............. %lld cycles", cycles/(ML_BENCH_LOOPS*10));
    printf("\n");

    // Fp negation with p = 2^512-569
    cycles = 0;
    for (n=0; n<ML_BENCH_LOOPS; n++)
    {
        random512(a, PCurve->prime); random512(b, PCurve->prime); random512(c, PCurve->prime); random512(d, PCurve->prime); random512(e, PCurve->prime); random512(f, PCurve->prime); 

        cycles1 = cpucycles();
        fpneg512(PCurve->prime, a);
        fpneg512(PCurve->prime, b);
        fpneg512(PCurve->prime, c);
        fpneg512(PCurve->prime, d);
        fpneg512(PCurve->prime, e);
        fpneg512(PCurve->prime, f);
        fpneg512(PCurve->prime, a);
        fpneg512(PCurve->prime, b);
        fpneg512(PCurve->prime, c);
        fpneg512(PCurve->prime, d);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Negation runs in ................... %lld cycles", cycles/(ML_BENCH_LOOPS*10));
    printf("\n");
    
    // Fp squaring with p = 2^512-569
    cycles = 0;
    for (n=0; n<ML_BENCH_LOOPS; n++)
    {
        random512(a, PCurve->prime); random512(b, PCurve->prime); random512(c, PCurve->prime); random512(d, PCurve->prime); random512(e, PCurve->prime); random512(f, PCurve->prime); 

        cycles1 = cpucycles();
        fpsqr512(a, b);
        fpsqr512(c, d);
        fpsqr512(e, f);
        fpsqr512(b, c);
        fpsqr512(d, e);
        fpsqr512(f, a);
        fpsqr512(c, d);
        fpsqr512(e, f);
        fpsqr512(a, c);
        fpsqr512(d, e);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Squaring runs in ................... %lld cycles", cycles/(ML_BENCH_LOOPS*10));
    printf("\n");
    
    // Fp multiplication with p = 2^512-569
    cycles = 0;
    for (n=0; n<ML_BENCH_LOOPS; n++)
    {
        random512(a, PCurve->prime); random512(b, PCurve->prime); random512(c, PCurve->prime); random512(d, PCurve->prime); random512(e, PCurve->prime); random512(f, PCurve->prime); 

        cycles1 = cpucycles();
        fpmul512(a, b, c);
        fpmul512(d, e, f);
        fpmul512(b, c, d);
        fpmul512(e, f, a);
        fpmul512(c, d, e);
        fpmul512(f, a, b);
        fpmul512(d, e, f);
        fpmul512(a, b, c);
        fpmul512(e, f, a);
        fpmul512(b, c, d);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Multiplication runs in ............. %lld cycles", cycles/(ML_BENCH_LOOPS*10));
    printf("\n");

    // Fp inversion with p = 2^512-569
    cycles = 0;
    for (n=0; n<ML_SHORT_BENCH_LOOPS; n++)
    {
        random512(a, PCurve->prime); random512(b, PCurve->prime); random512(c, PCurve->prime); random512(d, PCurve->prime); random512(e, PCurve->prime); random512(f, PCurve->prime); 

        cycles1 = cpucycles();
        fpinv512(a);
        fpinv512(b);
        fpinv512(c);
        fpinv512(d);
        fpinv512(e);
        fpinv512(f);
        fpinv512(a);
        fpinv512(b);
        fpinv512(c);
        fpinv512(d);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Inversion runs in .................. %lld cycles", cycles/(ML_SHORT_BENCH_LOOPS*10));
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
    OK = OK && fp_test256(&curve_Jac256);      // Test field operations with p = 2^256-189
#endif
#ifdef ECCURVES_384
	OK = OK && fp_test384(&curve_Jac384);      // Test field operations with p = 2^384-317
#endif
#ifdef ECCURVES_512
    OK = OK && fp_test512(&curve_Jac512);      // Test field operations with p = 2^512-569
#endif

#ifdef ECCURVES_256
    OK = OK && fp_run256(&curve_Jac256);       // Benchmark field operations with p = 2^256-189
#endif
#ifdef ECCURVES_384
    OK = OK && fp_run384(&curve_Jac384);       // Benchmark field operations with p = 2^384-317
#endif
#ifdef ECCURVES_512
    OK = OK && fp_run512(&curve_Jac512);       // Benchmark field operations with p = 2^512-569
#endif
    
    return OK;
}