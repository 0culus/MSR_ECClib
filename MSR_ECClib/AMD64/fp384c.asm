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
; Abstract: field operations over GF(2^384-317)
;
; This software is based on the article by Joppe Bos, Craig Costello, 
; Patrick Longa and Michael Naehrig, "Selecting elliptic curves for
; cryptography: an efficiency and security analysis", preprint available
; at http://eprint.iacr.org/2014/130.
;**************************************************************************

include macamd64.inc

P384_0 equ 18446744073709551299		  ; Prime p = 2^384-317
P384_1 equ 18446744073709551615
P384_2 equ 18446744073709551615
P384_3 equ 18446744073709551615
P384_4 equ 18446744073709551615
P384_5 equ 18446744073709551615
P384_c equ 317                        ; Value c in p = 2^384-c


.code
;****************************************************************************************
; (Constant-time) field multiplication using integer multiplication by product scanning
; Operation: c [r8] = a [rcx] * b [rdx] mod p, p = 2^384-317
; NOTE: input should have r8 != rcx and r8 != rdx 
; Inputs: a, b in [0, p-1]
; Output: c in [0, p-1] 
;****************************************************************************************
NESTED_ENTRY fpmul384_a, _TEXT00
  rex_push_reg rbx
  push_reg     rdi
  push_reg     r12
  push_reg     r13
  push_reg     r14
  push_reg     r15
  END_PROLOGUE

  mov  rbx, rdx  
  
  mov  rax, [rbx]
  mul qword ptr [rcx]    ; a0*b0
  mov  [r8], rax         ; C0
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
  adc  r11, 0
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
  adc  r9, 0

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
  adc  r10, 0

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
  adc  r10,0

  xor  r11, r11
  mov  rax, [rbx+8]
  mul qword ptr [rcx+24] ; a3*b1
  add  r9, rax
  adc  r10, rdx
  adc  r11, 0

  mov  rax, [rbx+32]
  mul qword ptr [rcx]    ; a0*b4
  add  r9, rax
  adc  r10, rdx
  adc  r11, 0

  mov  rax, [rbx]
  mul qword ptr [rcx+32] ; a4*b0
  add  r9, rax
  adc  r10, rdx
  adc  r11, 0

  mov  rax, [rbx+16]
  mul qword ptr [rcx+16] ; a2*b2
  add  r9, rax
  adc  r10, rdx
  adc  r11, 0

  mov  rax, [rbx+24]
  mul qword ptr [rcx+8]  ; a1*b3
  add  r9, rax
  mov  [r8+32], r9       ; C4
  adc  r10, rdx
  adc  r11, 0

  xor  r9, r9
  mov  rax, [rbx+16]
  mul qword ptr [rcx+24] ; a3*b2
  add  r10, rax
  adc  r11, rdx
  adc  r9, 0

  mov  rax, [rbx+32]
  mul qword ptr [rcx+8]  ; a1*b4
  add  r10, rax
  adc  r11, rdx
  adc  r9, 0

  mov  rax, [rbx+8]
  mul qword ptr [rcx+32] ; a4*b1
  add  r10, rax
  adc  r11, rdx
  adc  r9, 0

  mov  rax, [rbx+40]
  mul qword ptr [rcx]    ; a0*b5
  add  r10, rax
  adc  r11, rdx
  adc  r9, 0

  mov  rax, [rbx]
  mul qword ptr [rcx+40] ; a5*b0
  add  r10, rax
  adc  r11, rdx
  adc  r9, 0

  mov  rax, [rbx+24]
  mul qword ptr [rcx+16] ; a2*b3
  add  r10, rax          
  mov  [r8+40], r10      ; C5
  adc  r11, rdx
  adc  r9, 0
  
  xor  r10, r10
  mov  rax, [rbx+24]
  mul qword ptr [rcx+24] ; a3*b3
  add  r11, rax          
  adc  r9, rdx        
  adc  r10, 0    
  
  mov  rax, [rbx+32]
  mul qword ptr [rcx+16] ; a2*b4
  add  r11, rax          
  adc  r9, rdx        
  adc  r10, 0     
  
  mov  rax, [rbx+16]
  mul qword ptr [rcx+32] ; a4*b2
  add  r11, rax          
  adc  r9, rdx        
  adc  r10, 0        
  
  mov  rax, [rbx+40]
  mul qword ptr [rcx+8] ; a1*b5
  add  r11, rax          
  adc  r9, rdx        
  adc  r10, 0           
  
  mov  rax, [rbx+8]
  mul qword ptr [rcx+40] ; a5*b1
  add  r11, rax          ; r11 = C6        
  adc  r9, rdx        
  adc  r10, 0   
  
  xor  r12, r12
  mov  rax, [rbx+40]
  mul qword ptr [rcx+16] ; a2*b5
  add  r9, rax          
  adc  r10, rdx        
  adc  r12, 0    

  mov  rax, [rbx+16]
  mul qword ptr [rcx+40] ; a5*b2
  add  r9, rax          
  adc  r10, rdx        
  adc  r12, 0       

  mov  rax, [rbx+32]
  mul qword ptr [rcx+24] ; a3*b4
  add  r9, rax          
  adc  r10, rdx        
  adc  r12, 0         

  mov  rax, [rbx+24]
  mul qword ptr [rcx+32] ; a4*b3
  add  r9, rax           ; r9 = C7         
  adc  r10, rdx        
  adc  r12, 0     
  
  xor  r13, r13
  mov  rax, [rbx+40]
  mul qword ptr [rcx+24] ; a3*b5
  add  r10, rax          
  adc  r12, rdx        
  adc  r13, 0     
  
  mov  rax, [rbx+24]
  mul qword ptr [rcx+40] ; a5*b3
  add  r10, rax          
  adc  r12, rdx        
  adc  r13, 0          

  mov  rax, [rbx+32]
  mul qword ptr [rcx+32] ; a4*b4
  add  r10, rax          ; r10 = C8         
  adc  r12, rdx        
  adc  r13, 0         
  
  xor  r14, r14
  mov  rax, [rbx+40]
  mul qword ptr [rcx+32] ; a4*b5
  add  r12, rax          
  adc  r13, rdx        
  adc  r14, 0            

  mov  rax, [rbx+32]
  mul qword ptr [rcx+40] ; a5*b4
  add  r12, rax          ; r12 = C9         
  adc  r13, rdx        
  adc  r14, 0               

  mov  rax, [rbx+40]
  mul qword ptr [rcx+40] ; a5*b5
  add  r13, rax          ; r13 = C10          
  adc  r14, rdx          ; r14 = C11                         

; Reduction

  mov  rax, P384_c
  mul qword ptr r11 
  mov  r11, [r8]
  add  r11, rax         ; r11 = partial0
  adc  rdx, 0    
  mov  rdi, rdx 

  xor  rbx, rbx
  mov  rax, P384_c
  mul qword ptr r9 
  add  rax, rdi 
  adc  rbx, 0 
  mov  r9, [r8+8]    
  add  r9, rax          ; r9 = partial1
  adc  rbx, rdx 

  xor  rdi, rdi
  mov  rax, P384_c
  mul qword ptr r10 
  add  rax, rbx  
  adc  rdi, 0   
  mov  r10, [r8+16]    
  add  r10, rax          ; r10 = partial2
  adc  rdi, rdx  

  xor  rbx, rbx
  mov  rax, P384_c
  mul qword ptr r12 
  add  rax, rdi 
  adc  rbx, 0 
  mov  r12, [r8+24]    
  add  r12, rax          ; r12 = partial3
  adc  rbx, rdx 

  xor  rdi, rdi
  mov  rax, P384_c
  mul qword ptr r13 
  add  rax, rbx  
  adc  rdi, 0   
  mov  r13, [r8+32]    
  add  r13, rax          ; r13 = partial4
  adc  rdi, rdx  

  xor  rbx, rbx
  mov  rax, P384_c
  mul qword ptr r14 
  add  rax, rdi 
  adc  rbx, 1 
  mov  r14, [r8+40]    
  add  r14, rax          ; r14 = partial5
  adc  rdx, rbx          ; rdx = partial4 + 1 
  
  xor  rbx, rbx
  mov  rax, P384_c         
  mul qword ptr rdx   
  add  r11, rax			 ; r11 = partial0     
  adc  r9, 0			 ; r9 = partial1
  adc  r10, 0            ; r10 = partial2
  adc  r12, 0            ; r12 = partial4
  adc  r13, 0            ; r13 = partial5
  adc  r14, 0            ; r14 = partial6
  
  mov  rax, P384_c      ; final correction
  cmovc rax, rbx
  sub  r11, rax
  mov  [r8], r11  
  sbb  r9, 0
  mov  [r8+8], r9 
  sbb  r10, 0
  mov  [r8+16], r10  
  sbb  r12, 0
  mov  [r8+24], r12 
  sbb  r13, 0
  mov  [r8+32], r13 
  sbb  r14, 0
  mov  [r8+40], r14 
  
ifdef BEGIN_EPILOGUE
    BEGIN_EPILOGUE
endif  
  pop  r15
  pop  r14
  pop  r13
  pop  r12
  pop  rdi
  pop  rbx
  ret
NESTED_END fpmul384_a, _TEXT00


;****************************************************************************************
; (Constant-time) field squaring using integer multiplication by product scanning
; Operation: c [rdx] = a [rcx]^2 mod p, p = 2^384-317
; NOTE: input should have rdx != rcx 
; Input:  a in [0, p-1]
; Output: c in [0, p-1] 
;****************************************************************************************
NESTED_ENTRY fpsqr384_a, _TEXT00
  rex_push_reg rbx
  push_reg     rdi
  push_reg     rsi
  push_reg     r12
  push_reg     r13
  push_reg     r14
  END_PROLOGUE

  mov  rbx, rdx

  xor  r10, r10
  mov  rax, [rcx]
  mul qword ptr [rcx+8]   
  add  rax, rax
  mov  r8, rax
  adc  rdx, rdx
  mov  r9, rdx
  adc  r10, 0            ; 2*a0*a1

  mov  rax, [rcx]
  mul qword ptr rax      ; a0^2
  mov  [rbx], rax        ; C0
  add  r8, rdx
  mov  [rbx+8], r8       ; C1
  adc  r9, 0

  xor  r11, r11
  mov  rax, [rcx]
  mul qword ptr [rcx+16] 
  add  rax, rax
  mov  r8, rax
  adc  rdx, rdx
  mov  r14, rdx
  adc  r11, 0             ; 2*a0*a2

  mov  rax, [rcx+8]
  mul qword ptr rax       ; a1^2
  add  r8, rax
  adc  r14, rdx
  adc  r11, 0

  mov  rax, [rcx]
  mul qword ptr [rcx+24]  ; a0*a3
  add  r8, r9
  mov  [rbx+16], r8       ; C2
  adc  r14, r10
  adc  r11, 0
  mov  r8, rax
  mov  r10, rdx
  
  xor  r9, r9 
  mov  rax, [rcx+8]
  mul qword ptr [rcx+16]   ; a1*a2
  add  r8, rax
  adc  r10, rdx
  adc  r9, 0
  add  r8, r8
  adc  r10, r10
  adc  r9, r9              ; 2(a0*a3 + a1*a2)

  mov  rax, [rcx+8]
  mul qword ptr [rcx+24]   ; a1*a3
  add  r8, r14
  mov  [rbx+24], r8        ; C3
  adc  r10, r11
  adc  r9, 0
  xor  r11, r11
  mov  r8, rax
  mov  r14, rdx

  mov  rax, [rcx+32]
  mul qword ptr [rcx]      ; a0*a4
  add  r8, rax
  adc  r14, rdx
  adc  r11, 0
  add  r8, r8              
  adc  r14, r14
  adc  r11, r11            ; 2(a1*a3 + a0*a4)

  mov  rax, [rcx+16]
  mul qword ptr [rcx+16]   ; a2^2
  add  r10, r8
  adc  r14, r9
  adc  r11, 0
  add  r10, rax             ; C4
  mov  [rbx+32], r10
  adc  r14, rdx
  adc  r11, 0

  mov  rax, [rcx+16]
  mul qword ptr [rcx+24]    ; a2*a3
  xor  r9, r9
  mov  r8, rax
  mov  r10, rdx

  mov  rax, [rcx+32]
  mul qword ptr [rcx+8]    ; a1*a4
  add  r8, rax
  adc  r10, rdx
  adc  r9, 0

  mov  rax, [rcx+40]
  mul qword ptr [rcx]      ; a0*a5
  add  r8, rax
  adc  r10, rdx
  adc  r9, 0
  add  r8, r8              
  adc  r10, r10
  adc  r9, r9              ; 2(a2*a3 + a1*a4 + a0*a5)

  mov  rax, [rcx+40]
  mul qword ptr [rcx+8]    ; a1*a5
  add  r14, r8             ; C5
  mov  [rbx+40], r14
  adc  r11, r10
  adc  r9, 0
  xor  r14, r14
  mov  r8, rax
  mov  r10, rdx

  mov  rax, [rcx+32]
  mul qword ptr [rcx+16]    ; a2*a4
  add  r8, rax
  adc  r10, rdx
  adc  r14, 0
  add  r8, r8              
  adc  r10, r10
  adc  r14, r14            ; 2(a2*a4 + a1*a5)

  mov  rax, [rcx+24]
  mul qword ptr rax        ; a3^2
  add  r11, r8             
  adc  r9, r10
  adc  r14, 0
  add  r11, rax            ; r11 = C6
  adc  r9, rdx
  adc  r14, 0

  mov  rax, [rcx+40]
  mul qword ptr [rcx+16]   ; a2*a5
  xor  r12, r12
  mov  r8, rax
  mov  r10, rdx

  mov  rax, [rcx+32]
  mul qword ptr [rcx+24]   ; a3*a4
  add  r8, rax             
  adc  r10, rdx
  adc  r12, 0
  add  r8, r8        
  adc  r10, r10
  adc  r12, r12            ; 2(a2*a5 + a3*a4)
    
  xor  r13, r13
  mov  rax, [rcx+40]
  mul qword ptr [rcx+24] 
  add  r9, r8              ; r9 = C7
  adc  r10, r14
  adc  r12, 0
  add  rax, rax  
  mov  r8, rax            
  adc  rdx, rdx
  mov  r14, rdx
  adc  r13, 0               ; 2a3*a5

  mov  rax, [rcx+32]
  mul qword ptr rax         ; a4^2
  add  r8, rax             
  adc  r14, rdx
  adc  r13, 0
  add  r10, r8              ; r10 = C8
  adc  r12, r14
  adc  r13, 0
   
  xor  r14, r14
  mov  rax, [rcx+40]
  mul qword ptr [rcx+32]    ; a4*a5
  add  rax, rax             
  adc  rdx, rdx
  adc  r14, 0
  add  r12, rax             ; r12 = C9
  adc  r13, rdx
  adc  r14, 0

  mov  rax, [rcx+40]
  mul qword ptr rax         ; a5^2
  add  r13, rax             ; r13 = C10 
  adc  r14, rdx             ; r14 = C11

; Reduction  

  mov  rax, P384_c
  mul qword ptr r11 
  mov  r11, [rbx]
  add  r11, rax         ; r11 = partial0
  adc  rdx, 0    
  mov  rdi, rdx 

  xor  rsi, rsi
  mov  rax, P384_c
  mul qword ptr r9 
  add  rax, rdi 
  adc  rsi, 0 
  mov  r9, [rbx+8]    
  add  r9, rax          ; r9 = partial1
  adc  rsi, rdx 

  xor  rdi, rdi
  mov  rax, P384_c
  mul qword ptr r10 
  add  rax, rsi  
  adc  rdi, 0   
  mov  r10, [rbx+16]    
  add  r10, rax          ; r10 = partial2
  adc  rdi, rdx  

  xor  rsi, rsi
  mov  rax, P384_c
  mul qword ptr r12 
  add  rax, rdi 
  adc  rsi, 0 
  mov  r12, [rbx+24]    
  add  r12, rax          ; r12 = partial3
  adc  rsi, rdx 

  xor  rdi, rdi
  mov  rax, P384_c
  mul qword ptr r13 
  add  rax, rsi  
  adc  rdi, 0   
  mov  r13, [rbx+32]    
  add  r13, rax          ; r13 = partial4
  adc  rdi, rdx  

  xor  rsi, rsi
  mov  rax, P384_c
  mul qword ptr r14 
  add  rax, rdi 
  adc  rsi, 1 
  mov  r14, [rbx+40]    
  add  r14, rax          ; r14 = partial5
  adc  rdx, rsi          ; rdx = partial6 + 1 
  
  xor  rsi, rsi
  mov  rax, P384_c         
  mul qword ptr rdx   
  add  r11, rax			 ; r11 = partial0     
  adc  r9, 0			 ; r9 = partial1
  adc  r10, 0            ; r10 = partial2
  adc  r12, 0            ; r12 = partial4
  adc  r13, 0            ; r13 = partial5
  adc  r14, 0            ; r14 = partial6
  
  mov  rax, P384_c      ; final correction
  cmovc rax, rsi
  sub  r11, rax
  mov  [rbx], r11  
  sbb  r9, 0
  mov  [rbx+8], r9 
  sbb  r10, 0
  mov  [rbx+16], r10  
  sbb  r12, 0
  mov  [rbx+24], r12 
  sbb  r13, 0
  mov  [rbx+32], r13 
  sbb  r14, 0
  mov  [rbx+40], r14 
  
ifdef BEGIN_EPILOGUE
    BEGIN_EPILOGUE
endif  
  pop  r14
  pop  r13
  pop  r12
  pop  rsi
  pop  rdi
  pop  rbx
  ret
NESTED_END fpsqr384_a, _TEXT00


;*********************************************************************
; (Constant-time) field addition 
; Operation: c [r8] = a [rcx] + b [rdx] mod p, p = 2^384-317
; Input:  a, b in [0, p-1]
; Output: c in [0, p-1]
;*********************************************************************
NESTED_ENTRY fpadd384_a, _TEXT00
  rex_push_reg r12
  push_reg     r13
  push_reg     r14
  END_PROLOGUE
  
  mov  r9, [rcx]         ; a + P384_c
  add  r9, P384_c
  mov  r10, [rcx+8]
  adc  r10, 0
  mov  r11, [rcx+16]
  adc  r11, 0
  mov  r12, [rcx+24]
  adc  r12, 0
  mov  r13, [rcx+32]
  adc  r13, 0
  mov  r14, [rcx+40]
  adc  r14, 0
  
  mov  rcx, [rdx]       ; (a+P384_c) + b 
  add  r9, rcx
  mov  rcx, [rdx+8]
  adc  r10, rcx
  mov  rcx, [rdx+16]
  adc  r11, rcx
  mov  rcx, [rdx+24]
  adc  r12, rcx
  mov  rcx, [rdx+32]
  adc  r13, rcx
  mov  rcx, [rdx+40]
  adc  r14, rcx
  
  mov  rdx, 0           ; if (carry) then c = (a+P384_c) + b
  mov  rcx, P384_c     ; else c = (a+P384_c) + b - P384_c
  cmovc rcx, rdx
  sub  r9, rcx
  mov  [r8], r9
  sbb  r10, 0
  mov  [r8+8], r10
  sbb  r11, 0
  mov  [r8+16], r11
  sbb  r12, 0
  mov  [r8+24], r12
  sbb  r13, 0
  mov  [r8+32], r13
  sbb  r14, 0
  mov  [r8+40], r14
  
ifdef BEGIN_EPILOGUE
    BEGIN_EPILOGUE
endif  
  pop  r14
  pop  r13
  pop  r12
  ret
NESTED_END fpadd384_a, _TEXT00


;*********************************************************************
; (Constant-time) field subtraction
; Operation: c [r8] = a [rcx] + b [rdx] mod p, p = 2^384-317
; Input:  a, b in [0, p-1]
; Output: c in [0, p-1]
;*********************************************************************
NESTED_ENTRY fpsub384_a, _TEXT00
  rex_push_reg r12
  push_reg     r13
  push_reg     r14
  END_PROLOGUE

  xor  rax, rax        ; a - b
  mov  r9, [rcx]
  sub  r9, [rdx]
  mov  r10, [rcx+8]
  sbb  r10, [rdx+8]
  mov  r11, [rcx+16]
  sbb  r11, [rdx+16]
  mov  r12, [rcx+24]
  sbb  r12, [rdx+24]
  mov  r13, [rcx+32]
  sbb  r13, [rdx+32]
  mov  r14, [rcx+40]
  sbb  r14, [rdx+40]

  mov  rcx, P384_0    ; if (carry) then c = (a-b) - P384_c 
  cmovnc rcx, rax      ; else c = a - b
  mov  rdx, P384_1 
  cmovnc rdx, rax
  add  r9, rcx
  mov  [r8], r9
  adc  r10, rdx
  mov  [r8+8], r10
  adc  r11, rdx
  mov  [r8+16], r11
  adc  r12, rdx
  mov  [r8+24], r12
  adc  r13, rdx
  mov  [r8+32], r13
  adc  r14, rdx
  mov  [r8+40], r14
  
ifdef BEGIN_EPILOGUE
    BEGIN_EPILOGUE
endif  
  pop  r14
  pop  r13
  pop  r12
  ret
NESTED_END fpsub384_a, _TEXT00


;*********************************************************************
; (Constant-time) field division by 2 
; Operation: a [rdx] = a [rcx]/2 mod p, p = 2^384-317
; Input:  a in [0, p-1]
; Output: c in [0, p-1]
;*********************************************************************
NESTED_ENTRY fpdiv2_384_a, _TEXT00 
    rex_push_reg    r12
    push_reg        r13
    push_reg        rbx
  END_PROLOGUE

    xor     rbx, rbx
    mov     rax, [rcx]
    bt      rax, 0
    mov     r8, P384_0 
    cmovnc  r8, rbx
    mov     r9, P384_1     ; P384_1 = ... = P384_5
    cmovnc  r9, rbx
    mov     r10, r9  
    mov     r11, r9 
    mov     r12, r9 
    mov     r13, r9 
    
    add     r8, rax         ; if (a mod 2 = 1) then temp = a + p
    mov     rax, [rcx+8]    ; else temp = a + 0
    adc     r9, rax
    mov     rax, [rcx+16]
    adc     r10, rax 
    mov     rax, [rcx+24]
    adc     r11, rax
    mov     rax, [rcx+32]
    adc     r12, rax
    mov     rax, [rcx+40]
    adc     r13, rax
	adc		rbx, 0
  
	shrd	r8, r9, 1       ; c = temp/2
	mov		[rdx], r8
	shrd	r9, r10, 1
	mov		[rdx+8], r9
	shrd	r10, r11, 1
	mov		[rdx+16], r10
	shrd	r11, r12, 1
	mov		[rdx+24], r11	
	shrd	r12, r13, 1
	mov		[rdx+32], r12	
	shrd	r13, rbx, 1
	mov		[rdx+40], r13	
	
ifdef BEGIN_EPILOGUE
    BEGIN_EPILOGUE
endif  
	pop     rbx
	pop     r13
	pop     r12
	ret
NESTED_END fpdiv2_384_a, _TEXT00


;****************************************************************************************
; (Constant-time) field negation and subtraction from a modulus
; Operation: a [rdx] = modulus [rcx] - a [rdx]
;            if modulus = p = 2^384-317, then this performs a field negation -a (mod p)
; Input:  a in [0, modulus-1]
; Output: a in [0, modulus-1], rax = 1 (TRUE) if a <= modulus
;****************************************************************************************
LEAF_ENTRY fpneg384_a, _TEXT00
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
  mov  r9, [rcx+32]      
  sbb  r9, [rdx+32]
  mov  [rdx+32], r9
  mov  r9, [rcx+40]      
  sbb  r9, [rdx+40]
  mov  [rdx+40], r9
  
  setnb  al
  ret
LEAF_END fpneg384_a, _TEXT00


;*********************************************************************************************************************
; (Constant-time) Evaluation for the complete addition
; Operation: if [rcx] = 0 (P=-Q) then index=0, if [rdx] = 0 (P infinity) then index=1, if [r8] = 0 (P=Q) then index=2, 
;            else index=3
; Output:    if index(rax)=3 then mask [r9] = 0xff...ff, else mask [r9] = 0  
;*********************************************************************************************************************
NESTED_ENTRY complete_eval_Jac384, _TEXT00
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
  mov    r10, [rcx+32]
  or     r12, r10
  mov    r10, [rcx+40]
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
  mov    r10, [r8+32]
  or     r12, r10
  mov    r10, [r8+40]
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
  mov    r10, [rdx+32]
  or     r12, r10
  mov    r10, [rdx+40]
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
NESTED_END complete_eval_Jac384, _TEXT00


;*******************************************************************************************
; (Constant-time) Field element selection for the complete addition
; Operation: if (r9)=0 then c [r8] = a [rcx], else if (r9) = 0xff...ff then c [r8] = b [rdx]
;*******************************************************************************************
NESTED_ENTRY complete_select_Jac384, _TEXT00
  alloc_stack   8
  END_PROLOGUE

  xor           rax, rax
  ;sub           rsp, 8
  mov           [rsp], r9
  vbroadcastss  ymm0, DWORD PTR [rsp]
  vmovdqu       ymm1, YMMWORD PTR [rcx]       ; ymm1=a
  vmovdqu       ymm2, YMMWORD PTR [rdx]       ; ymm2=b
  vblendvpd     ymm3, ymm1, ymm2, ymm0        ; if ymm0=0 then ymm3=a else ymm3=b
  vmovdqu       YMMWORD PTR [r8], ymm3

  vmovdqu       xmm1, XMMWORD PTR [rcx+32]    ; xmm1=a
  vmovdqu       xmm2, XMMWORD PTR [rdx+32]    ; xmm2=b
  vblendvpd     xmm3, xmm1, xmm2, xmm0        ; if xmm0=0 then xmm3=a else xmm3=b
  vmovdqu       XMMWORD PTR [r8+32], xmm3
  mov           [rsp], rax                    ; Clean stack
  add           rsp, 8

ifdef BEGIN_EPILOGUE
    BEGIN_EPILOGUE
endif  
  ret
NESTED_END complete_select_Jac384, _TEXT00


;****************************************************************************************
; (Constant-time) Point extraction from 4-LUT for the complete mixed addition
; Operation: use index (rdx) to extract point from [rcx] and pass it to [r8]
;****************************************************************************************
NESTED_ENTRY complete_lut4_Jac384, _TEXT00
  alloc_stack  40
  END_PROLOGUE
  
  vmovdqu      YMMWORD PTR [rsp+8], ymm6      ; Save ymm6 in stack 

  xor          rax, rax  
  mov          r11, 18446744073709551615  
  ;sub          rsp, 8
  
  vmovdqu      ymm0, YMMWORD PTR [rcx]        ; Load table[0]
  vmovdqu      ymm1, YMMWORD PTR [rcx+32]
  vmovdqu      ymm2, YMMWORD PTR [rcx+64]
  vmovdqu      ymm3, YMMWORD PTR [rcx+96]
  vmovdqu      xmm4, XMMWORD PTR [rcx+128]

  dec          rdx                            ; Pass over table[1]
  cmovnz       rax, r11
  mov          [rsp], rax
  vbroadcastss ymm5, DWORD PTR [rsp]
  vmovdqu      ymm6, YMMWORD PTR [rcx+144]
  vblendvpd    ymm0, ymm6, ymm0, ymm5
  vmovdqu      ymm6, YMMWORD PTR [rcx+176]
  vblendvpd    ymm1, ymm6, ymm1, ymm5
  vmovdqu      ymm6, YMMWORD PTR [rcx+208]
  vblendvpd    ymm2, ymm6, ymm2, ymm5
  vmovdqu      ymm6, YMMWORD PTR [rcx+240]
  vblendvpd    ymm3, ymm6, ymm3, ymm5
  vmovdqu      xmm6, XMMWORD PTR [rcx+272]
  vblendvpd    xmm4, xmm6, xmm4, xmm5
  
  xor          rax, rax                        ; Pass over table[2]
  dec          rdx
  cmovnz       rax, r11
  mov          [rsp], rax
  vbroadcastss ymm5, DWORD PTR [rsp]
  vmovdqu      ymm6, YMMWORD PTR [rcx+288]
  vblendvpd    ymm0, ymm6, ymm0, ymm5
  vmovdqu      ymm6, YMMWORD PTR [rcx+320]
  vblendvpd    ymm1, ymm6, ymm1, ymm5
  vmovdqu      ymm6, YMMWORD PTR [rcx+352]
  vblendvpd    ymm2, ymm6, ymm2, ymm5
  vmovdqu      ymm6, YMMWORD PTR [rcx+384]
  vblendvpd    ymm3, ymm6, ymm3, ymm5
  vmovdqu      xmm6, XMMWORD PTR [rcx+416]
  vblendvpd    xmm4, xmm6, xmm4, xmm5
  
  xor          rax, rax                         ; Pass over table[3]
  dec          rdx
  cmovz        r11, rax
  mov          [rsp], r11
  vbroadcastss ymm5, DWORD PTR [rsp]
  vmovdqu      ymm6, YMMWORD PTR [rcx+432]
  vblendvpd    ymm0, ymm6, ymm0, ymm5
  vmovdqu      ymm6, YMMWORD PTR [rcx+464]
  vblendvpd    ymm1, ymm6, ymm1, ymm5
  vmovdqu      ymm6, YMMWORD PTR [rcx+496]
  vblendvpd    ymm2, ymm6, ymm2, ymm5
  vmovdqu      ymm6, YMMWORD PTR [rcx+528]
  vblendvpd    ymm3, ymm6, ymm3, ymm5
  vmovdqu      xmm6, XMMWORD PTR [rcx+560]
  vblendvpd    xmm4, xmm6, xmm4, xmm5
    
  mov          [rsp], rax                        ; Clean stack
  vmovdqu      YMMWORD PTR [r8], ymm0            ; [r8] = table[index]
  vmovdqu      YMMWORD PTR [r8+32], ymm1
  vmovdqu      YMMWORD PTR [r8+64], ymm2
  vmovdqu      YMMWORD PTR [r8+96], ymm3
  vmovdqu      XMMWORD PTR [r8+128], xmm4
    
  vmovdqu      ymm6, YMMWORD PTR [rsp+8]         ; Restore ymm6
  add          rsp, 40

ifdef BEGIN_EPILOGUE
    BEGIN_EPILOGUE
endif  
  ret
NESTED_END complete_lut4_Jac384, _TEXT00


;****************************************************************************************
; (Constant-time) Point extraction from 5-LUT for the complete addition
; Operation: use index (rdx) to extract point from [rcx] and pass it to [r8]
;****************************************************************************************
NESTED_ENTRY complete_lut5_Jac384, _TEXT00
  alloc_stack   40
  END_PROLOGUE
  
  vmovdqu      YMMWORD PTR [rsp+8], ymm6      ; Save ymm6 in stack

  xor          rax, rax  
  mov          r11, 18446744073709551615  
  ;sub          rsp, 8
  
  vmovdqu      ymm0, YMMWORD PTR [rcx]        ; Load table[0]
  vmovdqu      ymm1, YMMWORD PTR [rcx+32]
  vmovdqu      ymm2, YMMWORD PTR [rcx+64]
  vmovdqu      ymm3, YMMWORD PTR [rcx+96]
  vmovdqu      xmm4, XMMWORD PTR [rcx+128]

  dec          rdx                            ; Pass over table[1]
  cmovnz       rax, r11
  mov          [rsp], rax
  vbroadcastss ymm5, DWORD PTR [rsp]
  vmovdqu      ymm6, YMMWORD PTR [rcx+144]
  vblendvpd    ymm0, ymm6, ymm0, ymm5
  vmovdqu      ymm6, YMMWORD PTR [rcx+176]
  vblendvpd    ymm1, ymm6, ymm1, ymm5
  vmovdqu      ymm6, YMMWORD PTR [rcx+208]
  vblendvpd    ymm2, ymm6, ymm2, ymm5
  vmovdqu      ymm6, YMMWORD PTR [rcx+240]
  vblendvpd    ymm3, ymm6, ymm3, ymm5
  vmovdqu      xmm6, XMMWORD PTR [rcx+272]
  vblendvpd    xmm4, xmm6, xmm4, xmm5
  
  xor          rax, rax                        ; Pass over table[2]
  dec          rdx
  cmovnz       rax, r11
  mov          [rsp], rax
  vbroadcastss ymm5, DWORD PTR [rsp]
  vmovdqu      ymm6, YMMWORD PTR [rcx+288]
  vblendvpd    ymm0, ymm6, ymm0, ymm5
  vmovdqu      ymm6, YMMWORD PTR [rcx+320]
  vblendvpd    ymm1, ymm6, ymm1, ymm5
  vmovdqu      ymm6, YMMWORD PTR [rcx+352]
  vblendvpd    ymm2, ymm6, ymm2, ymm5
  vmovdqu      ymm6, YMMWORD PTR [rcx+384]
  vblendvpd    ymm3, ymm6, ymm3, ymm5
  vmovdqu      xmm6, XMMWORD PTR [rcx+416]
  vblendvpd    xmm4, xmm6, xmm4, xmm5
  
  xor          rax, rax                         ; Pass over table[3]
  dec          rdx
  cmovnz       rax, r11
  mov          [rsp], rax
  vbroadcastss ymm5, DWORD PTR [rsp]
  vmovdqu      ymm6, YMMWORD PTR [rcx+432]
  vblendvpd    ymm0, ymm6, ymm0, ymm5
  vmovdqu      ymm6, YMMWORD PTR [rcx+464]
  vblendvpd    ymm1, ymm6, ymm1, ymm5
  vmovdqu      ymm6, YMMWORD PTR [rcx+496]
  vblendvpd    ymm2, ymm6, ymm2, ymm5
  vmovdqu      ymm6, YMMWORD PTR [rcx+528]
  vblendvpd    ymm3, ymm6, ymm3, ymm5
  vmovdqu      xmm6, XMMWORD PTR [rcx+560]
  vblendvpd    xmm4, xmm6, xmm4, xmm5
  
  xor          rax, rax                         ; Pass over table[4]
  dec          rdx
  cmovz        r11, rax
  mov          [rsp], r11
  vbroadcastss ymm5, DWORD PTR [rsp]
  vmovdqu      ymm6, YMMWORD PTR [rcx+592]
  vblendvpd    ymm0, ymm6, ymm0, ymm5
  vmovdqu      ymm6, YMMWORD PTR [rcx+624]
  vblendvpd    ymm1, ymm6, ymm1, ymm5
  vmovdqu      ymm6, YMMWORD PTR [rcx+656]
  vblendvpd    ymm2, ymm6, ymm2, ymm5
  vmovdqu      ymm6, YMMWORD PTR [rcx+688]
  vblendvpd    ymm3, ymm6, ymm3, ymm5
  vmovdqu      xmm6, XMMWORD PTR [rcx+720]
  vblendvpd    xmm4, xmm6, xmm4, xmm5
    
  mov          [rsp], rax                        ; Clean stack
  vmovdqu      YMMWORD PTR [r8], ymm0            ; [r8] = table[index]
  vmovdqu      YMMWORD PTR [r8+32], ymm1
  vmovdqu      YMMWORD PTR [r8+64], ymm2
  vmovdqu      YMMWORD PTR [r8+96], ymm3
  vmovdqu      XMMWORD PTR [r8+128], xmm4
    
  vmovdqu      ymm6, YMMWORD PTR [rsp+8]         ; Restore ymm6
  add          rsp, 40

ifdef BEGIN_EPILOGUE
    BEGIN_EPILOGUE
endif  
  ret
NESTED_END complete_lut5_Jac384, _TEXT00


;****************************************************************************************
; Zeroing field element
;****************************************************************************************
LEAF_ENTRY fpzero384_a, _TEXT00
  xor          rax, rax
  mov          [rcx], rax
  mov          [rcx+8], rax 
  mov          [rcx+16], rax
  mov          [rcx+24], rax 
  mov          [rcx+32], rax 
  mov          [rcx+40], rax 
  ret
LEAF_END fpzero384_a, _TEXT00


END