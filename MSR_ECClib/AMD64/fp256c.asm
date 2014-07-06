;**************************************************************************
; MSR ECClib, an efficient and secure elliptic curve cryptographic library
;
;    Copyright (c) Microsoft Corporation. All rights reserved.
;    Licensed under the Apache License, Version 2.0 (the "License"); 
;    you may not use these files except in compliance with the License. 
;    You may obtain a copy of the License at
;                http://www.apache.org/licenses/LICENSE-2.0
;    Unless required by applicable law or agreed to in writing, software 
;    distributed under the License is distributed on an "AS IS" BASIS, 
;    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or 
;    implied. See the License for the specific language governing 
;    permissions and limitations under the License.
;
;
; Abstract: field operations over GF(2^256-189)
;
; This software is based on the article by Joppe Bos, Craig Costello, 
; Patrick Longa and Michael Naehrig, "Selecting elliptic curves for
; cryptography: an efficiency and security analysis", preprint available
; at http://eprint.iacr.org/2014/130.
;**************************************************************************

include macamd64.inc

P256_0 equ 18446744073709551427		  ; Prime p = 2^256-189
P256_1 equ 18446744073709551615
P256_c equ 189                        ; Value c in p = 2^256-c


.code
;****************************************************************************************
; (Constant-time) field multiplication using integer multiplication by product scanning
; Operation: c [r8] = a [rcx] * b [rdx] mod p, p = 2^256-189
; NOTE: input should have r8 != rcx and r8 != rdx 
; Inputs: a, b in [0, p-1]
; Output: c in [0, p-1] 
;****************************************************************************************
NESTED_ENTRY fpmul256_a, _TEXT00
  rex_push_reg rbx
  push_reg     rdi
  push_reg     r12
  push_reg     r13
  END_PROLOGUE

  mov  rbx, rdx  
  
  mov  rax, [rbx]
  mul qword ptr [rcx]    ; a0*b0
  mov  r13, rax          ; C0
  mov  r9, rdx           
  
  xor  r11, r11
  mov  rax, [rbx]
  mul qword ptr [rcx+8]  ; a1*b0
  add  r9, rax
  mov  r10, rdx
  adc  r10, 0

  mov  rax, [rbx+8]      ; a0*b1
  mul qword ptr [rcx] 
  add  r9, rax
  adc  r10, rdx
  ;adc  r11, 0
  setc r11b
  mov  [r8+8], r9        ; C1

  xor  r9, r9    
  mov  rax, [rbx]
  mul qword ptr [rcx+16] ; a2*b0
  add  r10, rax
  adc  r11, rdx
  ;adc  r9, 0

  mov  rax, [rbx+8]
  mul qword ptr [rcx+8]  ; a1*b1
  add  r10, rax
  adc  r11, rdx
  ;adc  r9, 0
  setc r9b

  mov  rax, [rbx+16]
  mul qword ptr [rcx]    ; a0*b2
  add  r10, rax
  mov  [r8+16], r10      ; C2
  adc  r11, rdx
  adc  r9, 0

  xor  r10, r10
  mov  rax, [rbx]
  mul qword ptr [rcx+24] ; a3*b0 
  add  r11, rax
  adc  r9, rdx
  ;adc  r10, 0
  setc r10b

  mov  rax,[rbx+8]
  mul qword ptr [rcx+16] ; a2*b1
  add  r11, rax
  adc  r9, rdx
  adc  r10, 0

  mov  rax,[rbx+16]
  mul qword ptr [rcx+8]  ; a1*b2
  add  r11, rax
  adc  r9, rdx
  adc  r10, 0

  mov  rax, [rbx+24]
  mul qword ptr [rcx]    ; a0*b3
  add  r11, rax
  mov  [r8+24], r11      ; C3
  adc  r9, rdx
  adc  r10, 0

  xor  r11, r11
  mov  rax, [rbx+8]
  mul qword ptr [rcx+24] ; a3*b1
  add  r9, rax
  adc  r10, rdx
  ;adc  r11, 0
  setc r11b

  mov  rax, [rbx+16]
  mul qword ptr [rcx+16] ; a2*b2
  add  r9, rax
  adc  r10, rdx
  adc  r11, 0

  mov  rax, [rbx+24]
  mul qword ptr [rcx+8]  ; a1*b3
  add  r9, rax
  mov  rdi, r9           ; rdi = C4
  adc  r10, rdx
  adc  r11, 0

  xor  r9, r9
  mov  rax, [rbx+16]
  mul qword ptr [rcx+24] ; a3*b2
  add  r10, rax
  adc  r11, rdx
  ;adc  r9, 0
  setc r9b

  mov  rax, [rbx+24]
  mul qword ptr [rcx+16] ; a2*b3
  add  r10, rax          ; r10 = C5
  adc  r11, rdx
  adc  r9, 0

  mov  rax, [rbx+24]
  mul qword ptr [rcx+24] ; a3*b3
  add  r11, rax          ; r11 = C6
  adc  r9, rdx           ; r9 = C7

; Reduction

  mov  rax, P256_c
  mul  rdi 
  add  r13, rax         ; r13 = partial0
  adc  rdx, 0    
  mov  rdi, rdx 

  xor  r12, r12
  mov  rax, P256_c
  mul  r10 
  add  rax, rdi 
  ;adc  r12, 0
  setc r12b 
  mov  r10, [r8+8]    
  add  r10, rax          ; r10 = partial1
  adc  r12, rdx 

  xor  rdi, rdi
  mov  rax, P256_c
  mul  r11 
  add  rax, r12  
  ;adc  rdi, 0 
  setb dil  
  mov  r11, [r8+16]    
  add  r11, rax          ; r11 = partial2
  adc  rdi, rdx  

  xor  r12, r12
  mov  rax, P256_c
  mul  r9 
  add  rax, rdi 
  adc  r12, 1       
  mov  r9, [r8+24]    
  add  r9, rax           ; r9 = partial3
  adc  rdx, r12	         ; rdx = partial4 + 1
  
  xor  r12, r12
  mov  rax, P256_c         
  mul  rdx   
  add  r13, rax			 ; r13 = partial0     
  adc  r10, 0			 ; r10 = partial1
  adc  r11, 0            ; r11 = partial2
  adc  r9, 0             ; r9 = partial3
  
  mov  rax, P256_c      ; final correction
  cmovc rax, r12
  sub  r13, rax
  mov  [r8], r13  
  sbb  r10, 0
  mov  [r8+8], r10 
  sbb  r11, 0
  mov  [r8+16], r11  
  sbb  r9, 0
  mov  [r8+24], r9   

ifdef BEGIN_EPILOGUE
    BEGIN_EPILOGUE
endif  
  pop  r13
  pop  r12
  pop  rdi
  pop  rbx
  ret
NESTED_END fpmul256_a, _TEXT00


;****************************************************************************************
; (Constant-time) field squaring using integer multiplication by product scanning
; Operation: c [rdx] = a [rcx]^2 mod p, p = 2^256-189
; NOTE: input should have rdx != rcx 
; Input:  a in [0, p-1]
; Output: c in [0, p-1] 
;****************************************************************************************
NESTED_ENTRY fpsqr256_a, _TEXT00
  rex_push_reg rbx
  push_reg     rdi
  push_reg     r12
  push_reg     r13
  END_PROLOGUE

  mov  rbx, rdx

  xor  r10, r10
  mov  rax, [rcx]
  mul qword ptr [rcx+8]   
  add  rax, rax
  mov  r8, rax
  adc  rdx, rdx
  mov  r9, rdx
  ;adc  r10, 0            ; 2*a0*a1
  setb r10b

  mov  rax, [rcx]
  mul qword ptr rax      ; a0^2
  mov  [rbx], rax        ; C0
  add  r8, rdx
  mov  [rbx+8], r8       ; C1
  adc  r9, 0

  xor  rdi, rdi
  mov  rax, [rcx]
  mul qword ptr [rcx+16] 
  add  rax, rax
  mov  r8, rax
  adc  rdx, rdx
  mov  r11, rdx
  ;adc  rdi, 0             ; 2*a0*a2
  setb dil

  mov  rax, [rcx+8]
  mul qword ptr [rcx+8]   ; a1^2
  add  r8, rax
  adc  r11, rdx
  adc  rdi, 0

  mov  rax, [rcx]
  mul qword ptr [rcx+24]  ; a0*a3
  add  r8, r9
  mov  [rbx+16], r8       ; C2
  adc  r11, r10
  adc  rdi, 0
  mov  r8, rax
  mov  r10, rdx
  
  xor  r9, r9 
  mov  rax, [rcx+8]
  mul qword ptr [rcx+16]   ; a1*a2
  add  r8, rax
  adc  r10, rdx
  ;adc  r9, 0
  setb r9b
  add  r8, r8
  adc  r10, r10
  adc  r9, r9

  mov  rax, [rcx+8]
  mul qword ptr [rcx+24]   ; a1*a3
  add  r8, r11
  mov  [rbx+24], r8        ; C3
  adc  r10, rdi
  adc  r9, 0
  xor  rdi, rdi
  add  rax, rax            ; 2*a1*a3
  mov  r8, rax
  adc  rdx, rdx
  mov  r11, rdx
  ;adc  rdi, 0
  setb dil
  
  mov  rax, [rcx+16]
  mul qword ptr [rcx+16]   ; a2^2
  add  r8, r10
  adc  r9, r11
  adc  rdi, 0
  add  r8, rax             ; r8 = C4
;  mov  [rbx+32], r8
  adc  r9, rdx
  adc  rdi, 0

  xor  r11, r11
  mov  rax, [rcx+16]
  mul qword ptr [rcx+24]
  add  rax, rax
  adc  rdx, rdx
  ;adc  r11, 0               ; 2*a2*a3
  setb r11b
  add  r9, rax              ; r9 = C5
;  mov  [rbx+40], r9
  adc  rdi, rdx
  adc  r11, 0

  mov  rax, [rcx+24]
  mul qword ptr rax         ; a3^2
  add  rdi, rax             ; rdi = C6 
;  mov  [rbx+48], rdi
  adc  r11, rdx             ; r11 = C7
;  mov  [rbx+56], r11

; Reduction

  mov  rax, P256_c
  mul qword ptr r8 
  mov  r8, [rbx]
  add  r8, rax              ; r8 = partial0
  adc  rdx, 0    
  mov  r10, rdx 

  xor  r12, r12
  mov  rax, P256_c
  mul qword ptr r9 
  add  rax, r10 
  ;adc  r12, 0
  setb r12b 
  mov  r9, [rbx+8]    
  add  r9, rax              ; r9 = partial1
  adc  r12, rdx 

  xor  rcx, rcx
  mov  rax, P256_c
  mul qword ptr rdi 
  add  rax, r12  
  ;adc  rcx, 0  
  setb cl 
  mov  rdi, [rbx+16]    
  add  rdi, rax             ; rdi = partial2
  adc  rcx, rdx  

  xor  r10, r10
  mov  rax, P256_c
  mul qword ptr r11 
  add  rax, rcx  
  adc  r10, 1       
  mov  r11, [rbx+24]    
  add  r11, rax             ; r11 = partial3
  adc  r10, rdx	            ; r10 = partial4 + 1
  
  xor  r12, r12
  mov  rax, P256_c         
  mul qword ptr r10 
  add  r8, rax			    ; r8 = partial0     
  adc  r9, 0			    ; r9 = partial1
  adc  rdi, 0               ; rdi = partial2
  adc  r11, 0               ; r11 = partial3
  
  mov  rax, P256_c         ; final correction
  cmovc rax, r12
  sub  r8, rax
  mov  [rbx], r8
  sbb  r9, 0
  mov  [rbx+8], r9 
  sbb  rdi, 0
  mov  [rbx+16], rdi  
  sbb  r11, 0
  mov  [rbx+24], r11  

ifdef BEGIN_EPILOGUE
    BEGIN_EPILOGUE
endif  
  pop  r13
  pop  r12
  pop  rdi
  pop  rbx
  ret
NESTED_END fpsqr256_a, _TEXT00


;****************************************************************************************
; (Constant-time) field addition 
; Operation: c [r8] = a [rcx] + b [rdx] mod p, p = 2^256-189
; Inputs:  a, b in [0, p-1]
; Output: c in [0, p-1]
;****************************************************************************************
LEAF_ENTRY fpadd256_a, _TEXT00  
  mov  r9, [rcx]      ; a + P256_c  
  add  r9, P256_c
  mov  r10, [rcx+8]
  adc  r10, 0
  mov  r11, [rcx+16]
  adc  r11, 0
  mov  rax, [rcx+24]
  adc  rax, 0
  
  mov  rcx, [rdx]     ; (a+P256_c) + b 
  add  r9, rcx
  mov  rcx, [rdx+8]
  adc  r10, rcx
  mov  rcx, [rdx+16]
  adc  r11, rcx
  mov  rcx, [rdx+24]
  adc  rax, rcx
  
  mov  rdx, 0          ; if (carry) then c = (a+P256_c) + b
  mov  rcx, P256_c     ; else c = (a+P256_c) + b - P256_c
  cmovc rcx, rdx
  sub  r9, rcx
  mov  [r8], r9
  sbb  r10, 0
  mov  [r8+8], r10
  sbb  r11, 0
  mov  [r8+16], r11
  sbb  rax, 0
  mov  [r8+24], rax
  ret
LEAF_END fpadd256_a, _TEXT00


;****************************************************************************************
; (Constant-time) field subtraction
; Operation: c [r8] = a [rcx] + b [rdx] mod p, p = 2^256-189
; Input:  a, b in [0, p-1]
; Output: c in [0, p-1]
;****************************************************************************************
NESTED_ENTRY fpsub256_a, _TEXT00
  rex_push_reg   r12
  END_PROLOGUE

  xor  rax, rax       ; a - b
  mov  r9, [rcx]
  sub  r9, [rdx]
  mov  r10, [rcx+8]
  sbb  r10, [rdx+8]
  mov  r11, [rcx+16]
  sbb  r11, [rdx+16]
  mov  r12, [rcx+24]
  sbb  r12, [rdx+24]

  mov  rcx, P256_c    ; if (carry) then c = (a-b) - P256_c  
  cmovnc rcx, rax     ; else c = a - b
  sub  r9, rcx
  mov  [r8], r9
  sbb  r10, 0
  mov  [r8+8], r10
  sbb  r11, 0
  mov  [r8+16], r11
  sbb  r12, 0
  mov  [r8+24], r12 

ifdef BEGIN_EPILOGUE
    BEGIN_EPILOGUE
endif  
  pop  r12
  ret
NESTED_END fpsub256_a, _TEXT00


;****************************************************************************************
; (Constant-time) field division by 2 
; Operation: a [rdx] = a [rcx]/2 mod p, p = 2^256-189
; Input:  a in [0, p-1]
; Output: c in [0, p-1]
;****************************************************************************************
LEAF_ENTRY fpdiv2_256_a, _TEXT00
  xor  r9, r9
  mov  r8, [rcx]
  bt   r8, 0
  mov  r11, P256_0
  cmovnc r11, r9
  mov  rax, P256_1
  cmovnc rax, r9
	
  add  r8, r11         ; if (a mod 2 = 1) then temp = a + p
  mov  r9, [rcx+8]     ; else temp = a + 0
  adc  r9, rax
  mov  r10, [rcx+16]
  adc  r10, rax
  mov  r11, [rcx+24]
  adc  r11, rax
  mov  rax, 0
  adc  rax, 0
  
  shrd r8, r9, 1       ; c = temp/2
  mov  [rdx], r8
  shrd r9, r10, 1
  mov  [rdx+8], r9
  shrd r10, r11, 1
  mov  [rdx+16], r10
  shrd r11, rax, 1
  mov  [rdx+24], r11
  ret
LEAF_END fpdiv2_256_a, _TEXT00


;****************************************************************************************
; (Constant-time) field negation and subtraction from a modulus
; Operation: a [rdx] = modulus [rcx] - a [rdx]
;            if modulus = p = 2^256-189, then this performs a field negation -a (mod p)
; Input:  a in [0, modulus-1]
; Output: a in [0, modulus-1], rax = 1 (TRUE) if a <= modulus
;****************************************************************************************
LEAF_ENTRY fpneg256_a, _TEXT00
  xor  rax, rax
  mov  r9, [rcx]      ; a = modulus - a
  sub  r9, [rdx]
  mov  [rdx], r9
  mov  r9, [rcx+8]      
  sbb  r9, [rdx+8]
  mov  [rdx+8], r9
  mov  r9, [rcx+16]      
  sbb  r9, [rdx+16]
  mov  [rdx+16], r9
  mov  r9, [rcx+24]      
  sbb  r9, [rdx+24]
  mov  [rdx+24], r9
  
  setnb  al
  ret
LEAF_END fpneg256_a, _TEXT00


;*********************************************************************************************************************
; (Constant-time) Evaluation for the complete addition
; Operation: if [rcx] = 0 (P=-Q) then index=0, if [rdx] = 0 (P infinity) then index=1, if [r8] = 0 (P=Q) then index=2, 
;            else index=3
; Output:    if index(rax)=3 then mask [r9] = 0xff...ff, else mask [r9] = 0  
;**********************************************************************************************************************
NESTED_ENTRY complete_eval_Jac256, _TEXT00
  rex_push_reg   r12
  END_PROLOGUE

  xor    rax, rax
  mov    r11, 3        ; index 3 (P+Q) 
  mov    r12, [rcx]
  mov    r10, [rcx+8]
  or     r12, r10
  mov    r10, [rcx+16]
  or     r12, r10
  mov    r10, [rcx+24]
  or     r12, r10
  cmovnz rax, r11      ; index 0 (P=-Q) if [rcx]=0
  
  mov    r11, 2         
  mov    r10, [r8]
  or     r12, r10
  mov    r10, [r8+8]
  or     r12, r10
  mov    r10, [r8+16]
  or     r12, r10
  mov    r10, [r8+24]
  or     r12, r10
  cmovz  rax, r11      ; index 2 (P=Q) if [rcx] & [r8]=0
  
  mov    r11, 1        
  mov    r12, [rdx]
  mov    r10, [rdx+8]
  or     r12, r10
  mov    r10, [rdx+16]
  or     r12, r10
  mov    r10, [rdx+24]
  or     r12, r10
  cmovz  rax, r11      ; index 1 (P infinity) if [rdx]=0

  xor    rcx, rcx
  mov    r10, 18446744073709551615
  mov    r11, rax
  sub    r11, 3
  cmovz  rcx, r10       ; mask = 0xff...f if index=3, else mask = 0
  mov    [r9], rcx

ifdef BEGIN_EPILOGUE
    BEGIN_EPILOGUE
endif  
  pop    r12
  ret
NESTED_END complete_eval_Jac256, _TEXT00


;*******************************************************************************************
; (Constant-time) Field element selection for the complete addition
; Operation: if (r9)=0 then c [r8] = a [rcx], else if (r9) = 0xff...ff then c [r8] = b [rdx]
;*******************************************************************************************
NESTED_ENTRY complete_select_Jac256, _TEXT00
  alloc_stack   8
  END_PROLOGUE

  xor           rax, rax
  ;sub           rsp, 8
  mov           [rsp], r9
  vbroadcastss  ymm0, DWORD PTR [rsp]
  vmovdqu       ymm1, YMMWORD PTR [rcx]     ; ymm1=a
  vmovdqu       ymm2, YMMWORD PTR [rdx]     ; ymm2=b
  vblendvpd     ymm3, ymm1, ymm2, ymm0      ; if ymm0=0 then ymm3=a else ymm3=b 
  vmovdqu       YMMWORD PTR [r8], ymm3
  mov           [rsp], rax                  ; Clean stack
  add           rsp, 8

ifdef BEGIN_EPILOGUE
    BEGIN_EPILOGUE
endif  
  ret
NESTED_END complete_select_Jac256, _TEXT00


;****************************************************************************************
; (Constant-time) Point extraction from 4-LUT for the complete mixed addition
; Operation: use index (rdx) to extract point from [rcx] and pass it to [r8]
;****************************************************************************************
NESTED_ENTRY complete_lut4_Jac256, _TEXT00
  alloc_stack   8
  END_PROLOGUE

  xor          rax, rax  
  mov          r11, 18446744073709551615  
  ;sub          rsp, 8
  
  vmovdqu      ymm0, YMMWORD PTR [rcx]        ; Load table[0]
  vmovdqu      ymm1, YMMWORD PTR [rcx+32]
  vmovdqu      ymm2, YMMWORD PTR [rcx+64]

  dec          rdx                            ; Pass over table[1]
  cmovnz       rax, r11
  mov          [rsp], rax
  vbroadcastss ymm4, DWORD PTR [rsp]
  vmovdqu      ymm5, YMMWORD PTR [rcx+96]
  vblendvpd    ymm0, ymm5, ymm0, ymm4
  vmovdqu      ymm5, YMMWORD PTR [rcx+128]
  vblendvpd    ymm1, ymm5, ymm1, ymm4
  vmovdqu      ymm5, YMMWORD PTR [rcx+160]
  vblendvpd    ymm2, ymm5, ymm2, ymm4
  
  xor          rax, rax                       ; Pass over table[2]
  dec          rdx
  cmovnz       rax, r11
  mov          [rsp], rax
  vbroadcastss ymm4, DWORD PTR [rsp]
  vmovdqu      ymm5, YMMWORD PTR [rcx+192]
  vblendvpd    ymm0, ymm5, ymm0, ymm4
  vmovdqu      ymm5, YMMWORD PTR [rcx+224]
  vblendvpd    ymm1, ymm5, ymm1, ymm4
  vmovdqu      ymm5, YMMWORD PTR [rcx+256]
  vblendvpd    ymm2, ymm5, ymm2, ymm4
  
  xor          rax, rax                       ; Pass over table[3]
  dec          rdx
  cmovz        r11, rax
  mov          [rsp], r11
  vbroadcastss ymm4, DWORD PTR [rsp]
  vmovdqu      ymm5, YMMWORD PTR [rcx+288]
  vblendvpd    ymm0, ymm5, ymm0, ymm4
  vmovdqu      ymm5, YMMWORD PTR [rcx+320]
  vblendvpd    ymm1, ymm5, ymm1, ymm4
  vmovdqu      ymm5, YMMWORD PTR [rcx+352]
  vblendvpd    ymm2, ymm5, ymm2, ymm4
    
  mov          [rsp], rax                     ; Clean stack
  vmovdqu      YMMWORD PTR [r8], ymm0         ; [r8] = table[index]
  vmovdqu      YMMWORD PTR [r8+32], ymm1
  vmovdqu      YMMWORD PTR [r8+64], ymm2
  add          rsp, 8

ifdef BEGIN_EPILOGUE
    BEGIN_EPILOGUE
endif  
  ret
NESTED_END complete_lut4_Jac256, _TEXT00


;****************************************************************************************
; (Constant-time) Point extraction from 5-LUT for the complete addition
; Operation: use index (rdx) to extract point from [rcx] and pass it to [r8]
;****************************************************************************************
NESTED_ENTRY complete_lut5_Jac256, _TEXT00
  alloc_stack   8
  END_PROLOGUE

  xor          rax, rax  
  mov          r11, 18446744073709551615  
  ;sub          rsp, 8
  
  vmovdqu      ymm0, YMMWORD PTR [rcx]        ; Load table[0]
  vmovdqu      ymm1, YMMWORD PTR [rcx+32]
  vmovdqu      ymm2, YMMWORD PTR [rcx+64]

  dec          rdx                            ; Pass over table[1]
  cmovnz       rax, r11
  mov          [rsp], rax
  vbroadcastss ymm4, DWORD PTR [rsp]
  vmovdqu      ymm5, YMMWORD PTR [rcx+96]
  vblendvpd    ymm0, ymm5, ymm0, ymm4
  vmovdqu      ymm5, YMMWORD PTR [rcx+128]
  vblendvpd    ymm1, ymm5, ymm1, ymm4
  vmovdqu      ymm5, YMMWORD PTR [rcx+160]
  vblendvpd    ymm2, ymm5, ymm2, ymm4
  
  xor          rax, rax                       ; Pass over table[2]
  dec          rdx
  cmovnz       rax, r11
  mov          [rsp], rax
  vbroadcastss ymm4, DWORD PTR [rsp]
  vmovdqu      ymm5, YMMWORD PTR [rcx+192]
  vblendvpd    ymm0, ymm5, ymm0, ymm4
  vmovdqu      ymm5, YMMWORD PTR [rcx+224]
  vblendvpd    ymm1, ymm5, ymm1, ymm4
  vmovdqu      ymm5, YMMWORD PTR [rcx+256]
  vblendvpd    ymm2, ymm5, ymm2, ymm4
  
  xor          rax, rax                       ; Pass over table[3]
  dec          rdx
  cmovnz       rax, r11
  mov          [rsp], rax
  vbroadcastss ymm4, DWORD PTR [rsp]
  vmovdqu      ymm5, YMMWORD PTR [rcx+288]
  vblendvpd    ymm0, ymm5, ymm0, ymm4
  vmovdqu      ymm5, YMMWORD PTR [rcx+320]
  vblendvpd    ymm1, ymm5, ymm1, ymm4
  vmovdqu      ymm5, YMMWORD PTR [rcx+352]
  vblendvpd    ymm2, ymm5, ymm2, ymm4
  
  xor          rax, rax                       ; Pass over table[4]
  dec          rdx
  cmovz        r11, rax
  mov          [rsp], r11
  vbroadcastss ymm4, DWORD PTR [rsp]
  vmovdqu      ymm5, YMMWORD PTR [rcx+384]
  vblendvpd    ymm0, ymm5, ymm0, ymm4
  vmovdqu      ymm5, YMMWORD PTR [rcx+416]
  vblendvpd    ymm1, ymm5, ymm1, ymm4
  vmovdqu      ymm5, YMMWORD PTR [rcx+448]
  vblendvpd    ymm2, ymm5, ymm2, ymm4
    
  mov          [rsp], rax                     ; Clean stack
  vmovdqu      YMMWORD PTR [r8], ymm0         ; [r8] = table[index]
  vmovdqu      YMMWORD PTR [r8+32], ymm1
  vmovdqu      YMMWORD PTR [r8+64], ymm2
  add          rsp, 8
  
ifdef BEGIN_EPILOGUE
    BEGIN_EPILOGUE
endif  
  ret
NESTED_END complete_lut5_Jac256, _TEXT00


;****************************************************************************************
; Zeroing field element
;****************************************************************************************
LEAF_ENTRY fpzero256_a, _TEXT00
  xor          rax, rax
  mov          [rcx], rax
  mov          [rcx+8], rax 
  mov          [rcx+16], rax
  mov          [rcx+24], rax 
  ret
LEAF_END fpzero256_a, _TEXT00


END