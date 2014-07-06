/**************************************************************************
* Sample code using MSR ECClib
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
* Abstract: sample code
*
* SECURITY NOTE: this sample code is for DEMONSTRATION PURPOSES ONLY.
*                DO NOT USE IN A PRODUCT AS IS
*
* This software is based on the article by Joppe Bos, Craig Costello, 
* Patrick Longa and Michael Naehrig, "Selecting elliptic curves for
* cryptography: an efficiency and security analysis", preprint available
* at http://eprint.iacr.org/2014/130.
***************************************************************************/  

#include <stdio.h>
#include <windows.h>
#include "msr_ecclib.h"
#include "../Tests/tests.h"


static void to_print(dig *value, int number)
{
    int i; 	  
    for (i=(number-1); i>=0; i--) {
        printf("%llx", *(value+i));
	}
    return;
}


#ifdef ECCURVES_256

BOOL ecc_dh_Jac256(PCurveStruct JacCurve)
{ // Runs ephemeral Diffie-Hellman key exchange using curve Jac256 (without hashing)
  // SECURITY NOTE: this function is for DEMONSTRATION PURPOSES ONLY. DO NOT USE AS IS
    dig a[ML_WORDS256], b[ML_WORDS256];
    point_Jac256 P, PK_A, PK_B, SH_A, SH_B, *P_fixed = NULL;
    
    printf("\n\nTESTING OF EPHEMERAL DIFFIE-HELLMAN KEY EXCHANGE (W/O HASHING), CURVE \"JAC256\" \n"); 
    printf("--------------------------------------------------------------------------------------------------------\n\n"); 

    // Offline calculation of the precomputed table to be used during key generation
    eccset_Jac256(P, JacCurve);                                  // Set generator
    P_fixed = ecc_precomp_fixed_Jac256(P, MEM_LARGE, JacCurve);  // Precompute table for generator; MEM_LARGE requires 3*2^5 = 96 points = 6KB of memory

    // Alice computes her public key
    random256(a, JacCurve->order);                               // Get some value as Alice's secret key (SECURITY NOTE: this function is for testing only)
    ecc_scalar_mul_fixed_Jac256(P_fixed, a, PK_A, MEM_LARGE, JacCurve);    // Compute Alice's public key PK_A = a*P
        
    printf(" Alice chooses a random integer \"a\" and computes her public key PK_A: \n");
    printf(" a = "); 
    to_print((dig *)a, 4);
    printf("\n");
    printf(" PK_A = a*P = (");
    to_print((dig *)PK_A->x, 4);
    printf(", ");
    to_print((dig *)PK_A->y, 4);
    printf(")\n\n");
    printf("    Alice sends PK_A to Bob. \n\n");

    // Bob computes his public key
    random256(b, JacCurve->order);                               // Get some value as Bob's secret key (SECURITY NOTE: this function is for testing only)
    ecc_scalar_mul_fixed_Jac256(P_fixed, b, PK_B, MEM_LARGE, JacCurve);    // Compute Bob's public key PK_B = b*P
        
    printf(" Bob chooses a random integer \"b\" and computes his public key PK_B: \n");
    printf(" b = "); 
    to_print((dig *)b, 4);
    printf("\n");
    printf(" PK_B = b*P = (");
    to_print((dig *)PK_B->x, 4);
    printf(", ");
    to_print((dig *)PK_B->y, 4);
    printf(")\n\n");
    printf("    Bob sends PK_B to Alice. \n\n");

    ecc_destroy_precomp_Jac256(P_fixed);                         // Destroy precomputed table once done using it

    // Alice computes her shared key using Bob's public key
    ecc_scalar_mul_Jac256(PK_B, a, SH_A, JacCurve);              // SH_A = a*PK_B = a*(b*P)

    printf(" Alice computes her shared key: \n");
    printf(" SHARED_A = a*PK_B = ("); 
    to_print((dig *)SH_A->x, 4);
    printf(", ");
    to_print((dig *)SH_A->y, 4);
    printf(")\n\n");

    // Bob computes his shared key using Alice's public key
    ecc_scalar_mul_Jac256(PK_A, b, SH_B, JacCurve);              // SH_B = b*PK_A = b*(a*P)

    printf(" Bob computes his shared key: \n");
    printf(" SHARED_B = b*PK_A = ("); 
    to_print((dig *)SH_B->x, 4);
    printf(", ");
    to_print((dig *)SH_B->y, 4);
    printf(")\n\n");
        
    if (fpcompare256(SH_A->x,SH_B->x)==0 && fpcompare256(SH_A->y,SH_B->y)==0) { 
        printf("    Shared keys matched (SUCCESS)\n\n");
    } else {
        printf("    Shared keys do not matched (FAILED)\n\n");
    }
  
    return TRUE;

 }


BOOL ecc_dh_Ted256(PCurveStruct TedCurve)
{ // Runs ephemeral Diffie-Hellman key exchange using curve Ted256 (without hashing)
  // SECURITY NOTE: this function is for DEMONSTRATION PURPOSES ONLY. DO NOT USE AS IS
    dig a[ML_WORDS256], b[ML_WORDS256];
    point_Ted256 P, PK_A, PK_B, SH_A, SH_B;
    point_extaff_precomp_Ted256 *P_fixed = NULL;
    
    printf("\n\nTESTING OF EPHEMERAL DIFFIE-HELLMAN KEY EXCHANGE (W/O HASHING), CURVE \"TED256\" \n"); 
    printf("--------------------------------------------------------------------------------------------------------\n\n"); 

    // Offline calculation of the precomputed table to be used during key generation
    eccset_Ted256(P, TedCurve);                                  // Set generator
    P_fixed = ecc_precomp_fixed_Ted256(P, MEM_LARGE, TedCurve);  // Precompute table for 4*generator; MEM_LARGE requires 3*2^5 = 96 points = 9KB of memory

    // Alice computes her public key
    random256(a, TedCurve->order);                               // Get some value as Alice's secret key (SECURITY NOTE: this function is for testing only)
    ecc_scalar_mul_fixed_Ted256(P_fixed, a, PK_A, MEM_LARGE, TedCurve);    // Compute Alice's public key PK_A = a*(4*P)
        
    printf(" Alice chooses a random integer \"a\" and computes her public key PK_A: \n");
    printf(" a = "); 
    to_print((dig *)a, 4);
    printf("\n");
    printf(" PK_A = a*(4*P) = (");
    to_print((dig *)PK_A->x, 4);
    printf(", ");
    to_print((dig *)PK_A->y, 4);
    printf(")\n\n");
    printf("    Alice sends PK_A to Bob. \n\n");

    // Bob computes his public key
    random256(b, TedCurve->order);                               // Get some value as Bob's secret key (SECURITY NOTE: this function is for testing only)
    ecc_scalar_mul_fixed_Ted256(P_fixed, b, PK_B, MEM_LARGE, TedCurve);    // Compute Bob's public key PK_B = b*(4*P)
        
    printf(" Bob chooses a random integer \"b\" and computes his public key PK_B: \n");
    printf(" b = "); 
    to_print((dig *)b, 4);
    printf("\n");
    printf(" PK_B = b*(4*P) = (");
    to_print((dig *)PK_B->x, 4);
    printf(", ");
    to_print((dig *)PK_B->y, 4);
    printf(")\n\n");
    printf("    Bob sends PK_B to Alice. \n\n");

    ecc_destroy_precomp_Ted256(P_fixed);                         // Destroy precomputed table once done using it

    // Alice computes her shared key using Bob's public key
    ecc_scalar_mul_Ted256(PK_B, a, SH_A, TedCurve);              // SH_A = a*PK_B = a*(b*4*P)

    printf(" Alice computes her shared key: \n");
    printf(" SHARED_A = a*PK_B = ("); 
    to_print((dig *)SH_A->x, 4);
    printf(", ");
    to_print((dig *)SH_A->y, 4);
    printf(")\n\n");

    // Bob computes his shared key using Alice's public key
    ecc_scalar_mul_Ted256(PK_A, b, SH_B, TedCurve);              // SH_B = b*PK_A = b*(a*4*P)

    printf(" Bob computes his shared key: \n");
    printf(" SHARED_B = b*PK_A = ("); 
    to_print((dig *)SH_B->x, 4);
    printf(", ");
    to_print((dig *)SH_B->y, 4);
    printf(")\n\n");
        
    if (fpcompare256(SH_A->x,SH_B->x)==0 && fpcompare256(SH_A->y,SH_B->y)==0) { 
        printf("    Shared keys matched (SUCCESS)\n\n");
    } else {
        printf("    Shared keys do not matched (FAILED)\n\n");
    }
  
    return TRUE;

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
    OK = OK && ecc_dh_Jac256(&curve_Jac256);       // Diffie-Hellman key exchange test using "Jac256"
    OK = OK && ecc_dh_Ted256(&curve_Ted256);       // Diffie-Hellman key exchange test using "Ted256"
#endif
    
    return OK;
}