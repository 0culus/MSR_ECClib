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
* Abstract: definitions of 512-bit functions for curves Jac512 and Ted512
*
* This software is based on the article by Joppe Bos, Craig Costello, 
* Patrick Longa and Michael Naehrig, "Selecting elliptic curves for
* cryptography: an efficiency and security analysis", preprint available
* at http://eprint.iacr.org/2014/130.
***************************************************************************/  

#include "msr_ecclib.h"


// Definition of 512-bit field elements and elements in Z_r
#define FELM     dig512



// Definition of point representation types, 512-bit     
#define POINT_WAFF              point_Jac512     
#define POINT_WJAC              point_jac_Jac512
#define POINT_PRECOMP_WCHU      point_chu_precomp_Jac512 
#define POINT_TE                point_Ted512
#define POINT_EXT_TE            point_extproj_Ted512   
#define POINT_PRECOMP_EXT_TE    point_extproj_precomp_Ted512
#define POINT_PRECOMP_EXTAFF_TE point_extaff_precomp_Ted512         

// Definition of field operations with p = 2^512-569
#define FP_COPY                 fpcopy512 
#define FP_ZERO                 fpzero512  
#define FP_ISZERO               fp_iszero512    
#define FP_ADD                  fpadd512   
#define FP_SUB                  fpsub512     
#define FP_NEG                  fpneg512         
#define MOD_EVAL                mod_eval512     
#define FP_DIV2                 fpdiv2_512       
#define FP_MULC                 fpmulc512     
#define FP_SQR                  fpsqr512      
#define FP_MUL                  fpmul512
#define FP_INV                  fpinv512                                       

// Selection of low-level field implementations with p = 2^512-569
#define FP_ZERO_LOW             fpzero512_a  
#define FP_ADD_LOW              fpadd512_a       
#define FP_SUB_LOW              fpsub512_a     
#define FP_NEG_LOW              fpneg512_a       
#define FP_DIV2_LOW             fpdiv2_512_a  
#define FP_MULC_LOW             fpmul512_a     
#define FP_SQR_LOW              fpsqr512_a        
#define FP_MUL_LOW              fpmul512_a  
#define FP_INV_LOW              fpinv512_fixedchain

// Definition of point operations with p = 2^512-569
#define ECCSET_W                        eccset_Jac512
#define ECCCOPY_W                       ecccopy_Jac512
#define ECCCOPY_WJAC                    ecccopy_jac_Jac512
#define ECCCONVERT_AFF_TO_JAC_W         eccconvert_aff_to_jac_Jac512
#define ECCZERO_WAFF                    ecczero_Jac512
#define ECCZERO_WJAC                    ecczero_jac_Jac512
#define ECCZERO_WCHU                    ecczero_chu_Jac512
#define ECC_IS_INFINITY_WAFF            ecc_is_infinity_Jac512
#define ECC_IS_INFINITY_WJAC            ecc_is_infinity_jac_Jac512
#define ECCNORM_W                       eccnorm_Jac512
#define ECCDOUBLE_WJAC                  eccdouble_jac_Jac512
#define ECCMADD_CONDITIONALS_WJAC       eccadd_mixed_jac_conditionals_Jac512
#define ECCUADD_WJAC                    eccadd_jac_Jac512
#define ECCUADD_NO_INIT_WJAC            eccadd_jac_no_init_Jac512
#define ECCUMADD_WJAC                   eccadd_mixed_jac_Jac512
#define ECCDOUBLEADD_WJAC               eccdoubleadd_jac_Jac512
#define ECCDOUBLEADD_CONDITIONALS_WJAC  eccdoubleadd_jac_conditionals_Jac512
#define ECCADD_PRECOMP_WJAC             eccadd_jac_precomp_Jac512
#define ECC_PRECOMP_WJAC                ecc_precomp_jac_Jac512
#define ECC_MUL_W                       ecc_scalar_mul_Jac512
#define LUT_WCHU                        lut_chu_Jac512
#define ECC_PRECOMP_FIXED_W             ecc_precomp_fixed_Jac512
#define ECC_PRECOMP_FIXED_INTERNAL_W    ecc_precomp_fixed_internal_Jac512
#define ECC_MUL_FIXED_W                 ecc_scalar_mul_fixed_Jac512
#define ECC_MUL_FIXED_INTERNAL_W        ecc_scalar_mul_fixed_internal_Jac512
#define LUT_WAFF                        lut_aff_Jac512
#define ECC_PRECOMP_DBLMUL_W            ecc_precomp_dblmul_Jac512
#define ECC_PRECOMP_DBLMUL_INTERNAL_W   ecc_precomp_dblmul_internal_Jac512
#define ECC_DESTROY_PRECOMP_W           ecc_destroy_precomp_Jac512
#define ECC_DBLMUL_W                    ecc_double_scalar_mul_Jac512
#define ECC_DBLMUL_INTERNAL_W           ecc_double_scalar_mul_internal_Jac512
#define COMPLETE_EVAL                   complete_eval_Jac512
#define COMPLETE_SELECT                 complete_select_Jac512
#define COMPLETE_LUT4                   complete_lut4_Jac512
#define COMPLETE_LUT5                   complete_lut5_Jac512

#define ECCSET_TE                       eccset_Ted512
#define ECCCOPY_EXT_TE                  ecccopy_extproj_Ted512
#define ECCCONVERT_AFF_TO_EXTPROJ_TE    eccconvert_aff_to_extproj_Ted512
#define ECCZERO_TE                      ecczero_Ted512
#define ECCZERO_EXT_TE                  ecczero_extproj_Ted512
#define ECCZERO_PRECOMP_EXT_TE          ecczero_extproj_precomp_Ted512
#define ECCZERO_PRECOMP_EXTAFF_TE       ecczero_extaff_precomp_Ted512
#define ECC_IS_NEUTRAL_AFF_TE           ecc_is_neutral_Ted512
#define ECC_IS_NEUTRAL_EXT_TE           ecc_is_neutral_extproj_Ted512
#define ECCNORM_TE                      eccnorm_Ted512
#define ECCDOUBLE_EXT_TE                eccdouble_extproj_Ted512
#define ECCUADD_EXT_TE                  eccadd_extproj_Ted512
#define ECCUADD_EXT_INTERNAL_TE         eccadd_extproj_internal_Ted512
#define ECCUMADD_EXT_TE                 eccadd_mixed_extproj_Ted512
#define ECCADD_TE_PRECOMP               eccadd_extproj_precomp_Ted512
#define ECC_PRECOMP_EXT_TE              ecc_precomp_extproj_Ted512
#define ECC_MUL_TE                      ecc_scalar_mul_Ted512
#define LUT_EXT_TE                      lut_extproj_Ted512
#define ECC_MUL_FIXED_TE                ecc_scalar_mul_fixed_Ted512
#define ECC_MUL_FIXED_INTERNAL_TE       ecc_scalar_mul_fixed_internal_Ted512
#define LUT_EXTAFF_TE                   lut_extaff_Ted512
#define ECC_PRECOMP_FIXED_TE            ecc_precomp_fixed_Ted512
#define ECC_PRECOMP_FIXED_INTERNAL_TE   ecc_precomp_fixed_internal_Ted512
#define ECC_PRECOMP_DBLMUL_TE           ecc_precomp_dblmul_Ted512
#define ECC_PRECOMP_DBLMUL_INTERNAL_TE  ecc_precomp_dblmul_internal_Ted512
#define ECC_DESTROY_PRECOMP_TE          ecc_destroy_precomp_Ted512
#define ECC_DBLMUL_TE                   ecc_double_scalar_mul_Ted512
#define ECC_DBLMUL_INTERNAL_TE          ecc_double_scalar_mul_internal_Ted512


#include "fp_template.c"
#include "ecc_template.c"


