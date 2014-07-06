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
; Abstract: field operations over GF(2^512-569)
;
; This software is based on the article by Joppe Bos, Craig Costello, 
; Patrick Longa and Michael Naehrig, "Selecting elliptic curves for
; cryptography: an efficiency and security analysis", preprint available
; at http://eprint.iacr.org/2014/130.
;**************************************************************************

include macamd64.inc

P512_0 equ 18446744073709551047		  ; Prime p = 2^512-569
P512_1 equ 18446744073709551615
P512_2 equ 18446744073709551615
P512_3 equ 18446744073709551615
P512_4 equ 18446744073709551615
P512_5 equ 18446744073709551615
P512_6 equ 18446744073709551615
P512_7 equ 18446744073709551615
P512_c equ 569                        ; Value c in p = 2^512-c


.code
;****************************************************************************************
; (Constant-time) field multiplication using integer multiplication by product scanning
; Operation: c [r8] = a [rcx] * b [rdx] mod p, p = 2^512-569
; NOTE: input should have r8 != rcx and r8 != rdx 
; Inputs: a, b in [0, p-1]
; Output: c in [0, p-1] 
;****************************************************************************************
NESTED_ENTRY fpmul512_a, _TEXT00
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
  
  xor  rdi, rdi
  mov  rax, [rbx+24]
  mul qword ptr [rcx+24] ; a3*b3
  add  r11, rax          
  adc  r9, rdx        
  adc  rdi, 0    
  
  mov  rax, [rbx+32]
  mul qword ptr [rcx+16] ; a2*b4
  add  r11, rax          
  adc  r9, rdx        
  adc  rdi, 0     
  
  mov  rax, [rbx+16]
  mul qword ptr [rcx+32] ; a4*b2
  add  r11, rax          
  adc  r9, rdx        
  adc  rdi, 0         
  
  mov  rax, [rbx+48]
  mul qword ptr [rcx] ; a0*b6
  add  r11, rax          
  adc  r9, rdx        
  adc  rdi, 0      
  
  mov  rax, [rbx]
  mul qword ptr [rcx+48] ; a6*b0
  add  r11, rax          
  adc  r9, rdx        
  adc  rdi, 0         
  
  mov  rax, [rbx+40]
  mul qword ptr [rcx+8] ; a1*b5
  add  r11, rax          
  adc  r9, rdx        
  adc  rdi, 0           
  
  mov  rax, [rbx+8]
  mul qword ptr [rcx+40] ; a5*b1
  add  r11, rax          
  mov  [r8+48], r11      ; C6       
  adc  r9, rdx        
  adc  rdi, 0   
  
  xor  r12, r12
  mov  rax, [rbx+40]
  mul qword ptr [rcx+16] ; a2*b5
  add  r9, rax          
  adc  rdi, rdx        
  adc  r12, 0    

  mov  rax, [rbx+16]
  mul qword ptr [rcx+40] ; a5*b2
  add  r9, rax          
  adc  rdi, rdx        
  adc  r12, 0        

  mov  rax, [rbx+8]
  mul qword ptr [rcx+48] ; a6*b1
  add  r9, rax          
  adc  rdi, rdx        
  adc  r12, 0           

  mov  rax, [rbx+48]
  mul qword ptr [rcx+8] ; a1*b6
  add  r9, rax          
  adc  rdi, rdx        
  adc  r12, 0          

  mov  rax, [rbx]
  mul qword ptr [rcx+56] ; a7*b0
  add  r9, rax          
  adc  rdi, rdx        
  adc  r12, 0          

  mov  rax, [rbx+56]
  mul qword ptr [rcx]    ; a0*b7
  add  r9, rax          
  adc  rdi, rdx        
  adc  r12, 0              

  mov  rax, [rbx+32]
  mul qword ptr [rcx+24] ; a3*b4
  add  r9, rax          
  adc  rdi, rdx        
  adc  r12, 0         

  mov  rax, [rbx+24]
  mul qword ptr [rcx+32] ; a4*b3
  add  r9, rax           
  mov  [r8+56], r9       ; C7          
  adc  rdi, rdx        
  adc  r12, 0     
  
  xor  r13, r13
  mov  rax, [rbx+40]
  mul qword ptr [rcx+24] ; a3*b5
  add  rdi, rax          
  adc  r12, rdx        
  adc  r13, 0     
  
  mov  rax, [rbx+24]
  mul qword ptr [rcx+40] ; a5*b3
  add  rdi, rax          
  adc  r12, rdx        
  adc  r13, 0       
  
  mov  rax, [rbx+56]
  mul qword ptr [rcx+8] ; a1*b7
  add  rdi, rax          
  adc  r12, rdx        
  adc  r13, 0       
  
  mov  rax, [rbx+8]
  mul qword ptr [rcx+56] ; a7*b1
  add  rdi, rax          
  adc  r12, rdx        
  adc  r13, 0      
  
  mov  rax, [rbx+16]
  mul qword ptr [rcx+48] ; a6*b2
  add  rdi, rax          
  adc  r12, rdx        
  adc  r13, 0      
  
  mov  rax, [rbx+48]
  mul qword ptr [rcx+16] ; a2*b6
  add  rdi, rax          
  adc  r12, rdx        
  adc  r13, 0            

  mov  rax, [rbx+32]
  mul qword ptr [rcx+32] ; a4*b4
  add  rdi, rax          ; C8 = rdi             
  adc  r12, rdx        
  adc  r13, 0         
  
  xor  r14, r14
  mov  rax, [rbx+40]
  mul qword ptr [rcx+32] ; a4*b5
  add  r12, rax          
  adc  r13, rdx        
  adc  r14, 0       
  
  mov  rax, [rbx+56]
  mul qword ptr [rcx+16] ; a2*b7
  add  r12, rax          
  adc  r13, rdx        
  adc  r14, 0       
  
  mov  rax, [rbx+16]
  mul qword ptr [rcx+56] ; a7*b2
  add  r12, rax          
  adc  r13, rdx        
  adc  r14, 0      
  
  mov  rax, [rbx+48]
  mul qword ptr [rcx+24] ; a3*b6
  add  r12, rax          
  adc  r13, rdx        
  adc  r14, 0     
  
  mov  rax, [rbx+24]
  mul qword ptr [rcx+48] ; a6*b3
  add  r12, rax          
  adc  r13, rdx        
  adc  r14, 0                           

  mov  rax, [rbx+32]
  mul qword ptr [rcx+40] ; a5*b4
  add  r12, rax          ; C9 = r12        
  adc  r13, rdx        
  adc  r14, 0          
  
  xor  r15, r15
  mov  rax, [rbx+56]
  mul qword ptr [rcx+24] ; a3*b7
  add  r13, rax          
  adc  r14, rdx        
  adc  r15, 0    
  
  mov  rax, [rbx+24]
  mul qword ptr [rcx+56] ; a7*b3
  add  r13, rax          
  adc  r14, rdx        
  adc  r15, 0    

  mov  rax, [rbx+48]
  mul qword ptr [rcx+32] ; a4*b6
  add  r13, rax          
  adc  r14, rdx        
  adc  r15, 0      

  mov  rax, [rbx+32]
  mul qword ptr [rcx+48] ; a6*b4
  add  r13, rax          
  adc  r14, rdx        
  adc  r15, 0     

  mov  rax, [rbx+40]
  mul qword ptr [rcx+40] ; a5*b5
  add  r13, rax          ; C10 = r13        
  adc  r14, rdx        
  adc  r15, 0             
  
  xor  r9, r9
  mov  rax, [rbx+56]
  mul qword ptr [rcx+32] ; a4*b7
  add  r14, rax          
  adc  r15, rdx        
  adc  r9, 0    
  
  mov  rax, [rbx+32]
  mul qword ptr [rcx+56] ; a7*b4
  add  r14, rax          
  adc  r15, rdx        
  adc  r9, 0      
  
  mov  rax, [rbx+48]
  mul qword ptr [rcx+40] ; a5*b6
  add  r14, rax          
  adc  r15, rdx        
  adc  r9, 0     
  
  mov  rax, [rbx+40]
  mul qword ptr [rcx+48] ; a6*b5
  add  r14, rax          ; C11 = r14
  adc  r15, rdx        
  adc  r9, 0               
  
  xor  r10, r10
  mov  rax, [rbx+56]
  mul qword ptr [rcx+40] ; a5*b7
  add  r15, rax          
  adc  r9, rdx        
  adc  r10, 0  
    
  mov  rax, [rbx+40]
  mul qword ptr [rcx+56] ; a7*b5
  add  r15, rax          
  adc  r9, rdx        
  adc  r10, 0    
    
  mov  rax, [rbx+48]
  mul qword ptr [rcx+48] ; a6*b6
  add  r15, rax          ; C12 = r15       
  adc  r9, rdx        
  adc  r10, 0                
  
  xor  r11, r11
  mov  rax, [rbx+56]
  mul qword ptr [rcx+48] ; a6*b7
  add  r9, rax          
  adc  r10, rdx        
  adc  r11, 0    
  
  mov  rax, [rbx+48]
  mul qword ptr [rcx+56] ; a7*b6
  add  r9, rax           ; C13 = r9       
  adc  r10, rdx        
  adc  r11, 0                                                 

  mov  rax, [rbx+56]
  mul qword ptr [rcx+56] ; a7*b7
  add  r10, rax          ; C14 = r10         
  adc  r11, rdx          ; C15 = r11                      

; Reduction

  mov  rax, P512_c
  mul qword ptr rdi 
  mov  rbx, [r8]
  add  rbx, rax         ; rbx = partial0
  adc  rdx, 0    
  mov  rcx, rdx 

  xor  rdi, rdi
  mov  rax, P512_c
  mul qword ptr r12
  add  rax, rcx  
  adc  rdi, 0   
  mov  r12, [r8+8]    
  add  r12, rax          ; r12 = partial1
  adc  rdi, rdx  

  xor  rcx, rcx
  mov  rax, P512_c
  mul qword ptr r13
  add  rax, rdi 
  adc  rcx, 0 
  mov  r13, [r8+16]    
  add  r13, rax          ; r13 = partial2
  adc  rcx, rdx 

  xor  rdi, rdi
  mov  rax, P512_c
  mul qword ptr r14
  add  rax, rcx  
  adc  rdi, 0   
  mov  r14, [r8+24]    
  add  r14, rax          ; r14 = partial3
  adc  rdi, rdx  

  xor  rcx, rcx
  mov  rax, P512_c
  mul qword ptr r15 
  add  rax, rdi 
  adc  rcx, 0 
  mov  r15, [r8+32]    
  add  r15, rax          ; r15 = partial4
  adc  rcx, rdx 

  xor  rdi, rdi
  mov  rax, P512_c
  mul qword ptr r9 
  add  rax, rcx  
  adc  rdi, 0   
  mov  r9, [r8+40]    
  add  r9, rax          ; r9 = partial5
  adc  rdi, rdx   

  xor  rcx, rcx
  mov  rax, P512_c
  mul qword ptr r10 
  add  rax, rdi 
  adc  rcx, 0 
  mov  r10, [r8+48]    
  add  r10, rax          ; r10 = partial6
  adc  rcx, rdx 

  xor  rdi, rdi
  mov  rax, P512_c
  mul qword ptr r11
  add  rax, rcx 
  adc  rdi, 1   
  mov  r11, [r8+56]    
  add  r11, rax          ; r11 = partial7
  adc  rdx, rdi          ; rdx = partial8 + 1
  
  xor  rdi, rdi
  mov  rax, P512_c         
  mul qword ptr rdx   
  add  rbx, rax			 ; rbx = partial0 
  adc  r12, 0			 ; r12 = partial1
  adc  r13, 0            ; r13 = partial2
  adc  r14, 0            ; r14 = partial3
  adc  r15, 0            ; r15 = partial4
  adc  r9, 0             ; r9 = partial5
  adc  r10, 0            ; r10 = partial6
  adc  r11, 0            ; r11 = partial7
    
  mov  rax, P512_c      ; final correction
  cmovc rax, rdi
  sub  rbx, rax
  mov  [r8], rbx  
  sbb  r12, 0
  mov  [r8+8], r12 
  sbb  r13, 0
  mov  [r8+16], r13  
  sbb  r14, 0
  mov  [r8+24], r14 
  sbb  r15, 0
  mov  [r8+32], r15 
  sbb  r9, 0
  mov  [r8+40], r9 
  sbb  r10, 0
  mov  [r8+48], r10 
  sbb  r11, 0
  mov  [r8+56], r11 
  
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
NESTED_END fpmul512_a, _TEXT00


;****************************************************************************************
; (Constant-time) field squaring using integer multiplication by product scanning
; Operation: c [rdx] = a [rcx]^2 mod p, p = 2^512-569
; NOTE: input should have rdx != rcx 
; Input:  a in [0, p-1]
; Output: c in [0, p-1] 
;****************************************************************************************
NESTED_ENTRY fpsqr512_a, _TEXT00
  rex_push_reg rbx
  push_reg     rsi
  push_reg     rdi
  push_reg     r12
  push_reg     r13
  push_reg     r14
  push_reg     r15
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

  mov  rax, [rcx+48]
  mul qword ptr [rcx]       ; a0*a6
  add  r8, rax
  adc  r10, rdx
  adc  r14, 0

  mov  rax, [rcx+32]
  mul qword ptr [rcx+16]    ; a2*a4
  add  r8, rax
  adc  r10, rdx
  adc  r14, 0
  add  r8, r8              
  adc  r10, r10
  adc  r14, r14            ; 2(a2*a4 + a1*a5 + a0*a6)

  mov  rax, [rcx+24]
  mul qword ptr rax        ; a3^2
  add  r11, r8             
  adc  r9, r10
  adc  r14, 0
  add  r11, rax            ; C6
  mov  [rbx+48], r11
  adc  r9, rdx
  adc  r14, 0

  mov  rax, [rcx+40]
  mul qword ptr [rcx+16]   ; a2*a5
  xor  r12, r12
  mov  r8, rax
  mov  r10, rdx

  mov  rax, [rcx+56]
  mul qword ptr [rcx]      ; a0*a7
  add  r8, rax
  adc  r10, rdx
  adc  r12, 0

  mov  rax, [rcx+48]
  mul qword ptr [rcx+8]    ; a1*a6
  add  r8, rax
  adc  r10, rdx
  adc  r12, 0

  mov  rax, [rcx+32]
  mul qword ptr [rcx+24]   ; a3*a4
  add  r8, rax             
  adc  r10, rdx
  adc  r12, 0
  add  r8, r8        
  adc  r10, r10
  adc  r12, r12            ; 2(a2*a5 + a3*a4 + a0*a7 + a1*a6)
    
  xor  r13, r13
  mov  rax, [rcx+40]
  mul qword ptr [rcx+24]   ; a3*a5
  add  r9, r8              ; C7
  mov  [rbx+56], r9
  adc  r10, r14
  adc  r12, 0
  xor  r13, r13
  mov  r8, rax
  mov  r14, rdx

  mov  rax, [rcx+56]
  mul qword ptr [rcx+8]     ; a1*a7
  add  r8, rax
  adc  r14, rdx
  adc  r13, 0

  mov  rax, [rcx+48]
  mul qword ptr [rcx+16]    ; a2*a6
  add  r8, rax
  adc  r14, rdx
  adc  r13, 0
  add  r8, r8              
  adc  r14, r14
  adc  r13, r13            ; 2(a2*a6 + a1*a7 + a3*a5)

  mov  rax, [rcx+32]
  mul qword ptr rax         ; a4^2
  add  r8, rax             
  adc  r14, rdx
  adc  r13, 0
  add  r10, r8              ; r10 = C8
  adc  r12, r14
  adc  r13, 0  

  mov  rax, [rcx+40]
  mul qword ptr [rcx+32]    ; a4*a5
  xor  r14, r14
  mov  r15, rax
  mov  r9, rdx

  mov  rax, [rcx+56]
  mul qword ptr [rcx+16]    ; a2*a7
  add  r15, rax
  adc  r9, rdx
  adc  r14, 0

  mov  rax, [rcx+48]
  mul qword ptr [rcx+24]   ; a3*a6
  add  r15, rax
  adc  r9, rdx
  adc  r14, 0
  add  r15, r15              
  adc  r9, r9
  adc  r14, r14            ; 2(a4*a5 + a2*a7 + a3*a6)

  mov  rax, [rcx+56]
  mul qword ptr [rcx+24]   ; a3*a7
  add  r12, r15            ; r12 = C9
  adc  r13, r9
  adc  r14, 0
  xor  r15, r15
  mov  r9, rax
  mov  rdi, rdx  

  mov  rax, [rcx+48]
  mul qword ptr [rcx+32]   ; a4*a6
  add  r9, rax
  adc  rdi, rdx
  adc  r15, 0
  add  r9, r9              
  adc  rdi, rdi
  adc  r15, r15            ; 2(a3*a7 + a4*a6)

  mov  rax, [rcx+40]
  mul qword ptr [rcx+40]   ; a5^2
  add  r13, r9
  adc  r14, rdi
  adc  r15, 0 
  mov  rdi, rax
  mov  r11, rdx
   
  mov  rax, [rcx+56]
  mul qword ptr [rcx+32]   ; a4*a7
  add  r13, rdi            ; r13 = C10
  adc  r14, r11
  adc  r15, 0 
  xor  r9, r9
  mov  rdi, rax
  mov  r11, rdx
   
  mov  rax, [rcx+48]
  mul qword ptr [rcx+40]   ; a5*a6
  add  rdi, rax
  adc  r11, rdx
  adc  r9, 0
  add  rdi, rdi
  adc  r11, r11
  adc  r9, r9              ; 2(a4*a7 + a5*a6)
    
  mov  rax, [rcx+56]
  mul qword ptr [rcx+40] 
  add  r14, rdi              ; r14 = C11
  adc  r15, r11
  adc  r9, 0
  xor  rdi, rdi
  add  rax, rax  
  mov  rsi, rax            
  adc  rdx, rdx
  mov  r11, rdx
  adc  rdi, 0               ; 2a5*a7

  mov  rax, [rcx+48]
  mul qword ptr rax         ; a6^2
  add  rsi, rax             
  adc  r11, rdx
  adc  rdi, 0
  add  r15, rsi              ; r15 = C12
  adc  r9, r11
  adc  rdi, 0
   
  xor  r11, r11
  mov  rax, [rcx+56]
  mul qword ptr [rcx+48]    ; a6*a7
  add  rax, rax             
  adc  rdx, rdx
  adc  r11, 0
  add  r9, rax             ; r9 = C13
  adc  rdi, rdx
  adc  r11, 0

  mov  rax, [rcx+56]
  mul qword ptr rax         ; a7^2
  add  rdi, rax             ; rdi = C14 
  adc  r11, rdx             ; r11 = C15                      

; Reduction

  mov  rax, P512_c
  mul qword ptr r10 
  mov  rsi, [rbx]
  add  rsi, rax         ; rsi = partial0
  adc  rdx, 0    
  mov  rcx, rdx 

  xor  r10, r10
  mov  rax, P512_c
  mul qword ptr r12
  add  rax, rcx  
  adc  r10, 0   
  mov  r12, [rbx+8]    
  add  r12, rax          ; r12 = partial1
  adc  r10, rdx  

  xor  rcx, rcx
  mov  rax, P512_c
  mul qword ptr r13
  add  rax, r10 
  adc  rcx, 0 
  mov  r13, [rbx+16]    
  add  r13, rax          ; r13 = partial2
  adc  rcx, rdx 

  xor  r10, r10
  mov  rax, P512_c
  mul qword ptr r14
  add  rax, rcx  
  adc  r10, 0   
  mov  r14, [rbx+24]    
  add  r14, rax          ; r14 = partial3
  adc  r10, rdx  

  xor  rcx, rcx
  mov  rax, P512_c
  mul qword ptr r15 
  add  rax, r10 
  adc  rcx, 0 
  mov  r15, [rbx+32]    
  add  r15, rax          ; r15 = partial4
  adc  rcx, rdx 

  xor  r10, r10
  mov  rax, P512_c
  mul qword ptr r9 
  add  rax, rcx  
  adc  r10, 0   
  mov  r9, [rbx+40]    
  add  r9, rax          ; r9 = partial5
  adc  r10, rdx   

  xor  rcx, rcx
  mov  rax, P512_c
  mul qword ptr rdi 
  add  rax, r10 
  adc  rcx, 0 
  mov  rdi, [rbx+48]    
  add  rdi, rax          ; rdi = partial6
  adc  rcx, rdx 

  xor  r10, r10
  mov  rax, P512_c
  mul qword ptr r11
  add  rax, rcx 
  adc  r10, 1   
  mov  r11, [rbx+56]    
  add  r11, rax          ; r11 = partial7
  adc  rdx, r10          ; rdx = partial8 + 1
  
  xor  r10, r10
  mov  rax, P512_c         
  mul qword ptr rdx   
  add  rsi, rax			 ; rsi = partial0 
  adc  r12, 0			 ; r12 = partial1
  adc  r13, 0            ; r13 = partial2
  adc  r14, 0            ; r14 = partial3
  adc  r15, 0            ; r15 = partial4
  adc  r9, 0             ; r9 = partial5
  adc  rdi, 0            ; rdi = partial6
  adc  r11, 0            ; r11 = partial7
    
  mov  rax, P512_c      ; final correction
  cmovc rax, r10
  sub  rsi, rax
  mov  [rbx], rsi  
  sbb  r12, 0
  mov  [rbx+8], r12 
  sbb  r13, 0
  mov  [rbx+16], r13  
  sbb  r14, 0
  mov  [rbx+24], r14 
  sbb  r15, 0
  mov  [rbx+32], r15 
  sbb  r9, 0
  mov  [rbx+40], r9 
  sbb  rdi, 0
  mov  [rbx+48], rdi 
  sbb  r11, 0
  mov  [rbx+56], r11 
  
ifdef BEGIN_EPILOGUE
    BEGIN_EPILOGUE
endif  
  pop  r15
  pop  r14
  pop  r13
  pop  r12
  pop  rdi
  pop  rsi
  pop  rbx
  ret
NESTED_END fpsqr512_a, _TEXT00


;*********************************************************************
; (Constant-time) field addition 
; Operation: c [r8] = a [rcx] + b [rdx] mod p, p = 2^512-569
; Input:  a, b in [0, p-1]
; Output: c in [0, p-1]
;*********************************************************************
NESTED_ENTRY fpadd512_a, _TEXT00
  rex_push_reg r12
  push_reg     r13
  push_reg     r14
  push_reg     r15
  push_reg     rbx
  END_PROLOGUE
  
  mov  r9, [rcx]        ; a + P512_c
  add  r9, P512_c
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
  mov  r15, [rcx+48]
  adc  r15, 0
  mov  rbx, [rcx+56]
  adc  rbx, 0
  
  mov  rcx, [rdx]       ; (a+P512_c) + b 
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
  mov  rcx, [rdx+48]
  adc  r15, rcx
  mov  rcx, [rdx+56]
  adc  rbx, rcx
  
  mov  rdx, 0           ; if (carry) then c = (a+P512_c) + b
  mov  rcx, P512_c      ; else c = (a+P512_c) + b - P512_c
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
  sbb  r15, 0
  mov  [r8+48], r15
  sbb  rbx, 0
  mov  [r8+56], rbx
  
ifdef BEGIN_EPILOGUE
    BEGIN_EPILOGUE
endif  
  pop  rbx
  pop  r15
  pop  r14
  pop  r13
  pop  r12
  ret
NESTED_END fpadd512_a, _TEXT00


;*********************************************************************
; (Constant-time) field subtraction
; Operation: c [r8] = a [rcx] + b [rdx] mod p, p = 2^512-569
; Input:  a, b in [0, p-1]
; Output: c in [0, p-1]
;*********************************************************************
NESTED_ENTRY fpsub512_a, _TEXT00
  rex_push_reg r12
  push_reg     r13
  push_reg     r14
  push_reg     r15
  push_reg     rbx
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
  mov  r15, [rcx+48]
  sbb  r15, [rdx+48]
  mov  rbx, [rcx+56]
  sbb  rbx, [rdx+56]

  mov  rcx, P512_0      ; if (carry) then c = (a-b) - P512_c 
  cmovnc rcx, rax       ; else c = a - b
  mov  rdx, P512_1 
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
  adc  r15, rdx
  mov  [r8+48], r15
  adc  rbx, rdx
  mov  [r8+56], rbx
  
ifdef BEGIN_EPILOGUE
    BEGIN_EPILOGUE
endif  
  pop  rbx
  pop  r15
  pop  r14
  pop  r13
  pop  r12
  ret
NESTED_END fpsub512_a, _TEXT00


;*********************************************************************
; (Constant-time) field division by 2 
; Operation: a [rdx] = a [rcx]/2 mod p, p = 2^512-569
; Input:  a in [0, p-1]
; Output: c in [0, p-1]
;*********************************************************************
NESTED_ENTRY fpdiv2_512_a, _TEXT00 
    rex_push_reg    r12
    push_reg        r13
    push_reg        r14
    push_reg        r15
    push_reg        rbx
  END_PROLOGUE

    xor     rbx, rbx
    mov     rax, [rcx]
    bt      rax, 0
    mov     r8, P512_0 
    cmovnc  r8, rbx
    mov     r9, P512_1     ; P512_1 = ... = P512_7
    cmovnc  r9, rbx
    mov     r10, r9  
    mov     r11, r9 
    mov     r12, r9 
    mov     r13, r9 
    mov     r14, r9 
    mov     r15, r9 
    
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
    mov     rax, [rcx+48]
    adc     r14, rax
    mov     rax, [rcx+56]
    adc     r15, rax
	adc		rbx, 0
  
	shrd	r8, r9, 1        ; c = temp/2
	mov		[rdx], r8
	shrd	r9, r10, 1
	mov		[rdx+8], r9
	shrd	r10, r11, 1
	mov		[rdx+16], r10
	shrd	r11, r12, 1
	mov		[rdx+24], r11	
	shrd	r12, r13, 1
	mov		[rdx+32], r12		
	shrd	r13, r14, 1
	mov		[rdx+40], r13	
	shrd	r14, r15, 1
	mov		[rdx+48], r14
	shrd	r15, rbx, 1
	mov		[rdx+56], r15	
	
ifdef BEGIN_EPILOGUE
    BEGIN_EPILOGUE
endif  
	pop     rbx
	pop     r15
	pop     r14
	pop     r13
	pop     r12
	ret
NESTED_END fpdiv2_512_a, _TEXT00


;****************************************************************************************
; (Constant-time) field negation and subtraction from a modulus
; Operation: a [rdx] = modulus [rcx] - a [rdx]
;            if modulus = p = 2^512-569, then this performs a field negation -a (mod p)
; Input:  a in [0, modulus-1]
; Output: a in [0, modulus-1], rax = 1 (TRUE) if a <= modulus
;****************************************************************************************
LEAF_ENTRY fpneg512_a, _TEXT00
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
  mov  r9, [rcx+48]      
  sbb  r9, [rdx+48]
  mov  [rdx+48], r9
  mov  r9, [rcx+56]      
  sbb  r9, [rdx+56]
  mov  [rdx+56], r9
  
  setnb  al
  ret
LEAF_END fpneg512_a, _TEXT00


;**********************************************************************************************************************
; (Constant-time) Evaluation for the complete addition
; Operation: if [rcx] = 0 (P=-Q) then index=0, if [rdx] = 0 (P infinity) then index=1, if [r8] = 0 (P=Q) then index=2, 
;            else index=3
; Output:    if index(rax)=3 then mask [r9] = 0xff...ff, else mask [r9] = 0  
;**********************************************************************************************************************
NESTED_ENTRY complete_eval_Jac512, _TEXT00
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
  mov    r10, [rcx+48]
  or     r12, r10
  mov    r10, [rcx+56]
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
  mov    r10, [r8+48]
  or     r12, r10
  mov    r10, [r8+56]
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
  mov    r10, [rdx+48]
  or     r12, r10
  mov    r10, [rdx+56]
  or     r12, r10
  cmovz  rax, r11      ; index 1 (P inf) if [rdx]=0

  xor    rcx, rcx
  mov    r10, 18446744073709551615
  mov    r11, rax
  sub    r11, 3
  cmovz  rcx, r10       ; mask = ff...f if index=3, else mask = 0
  mov    [r9], rcx
  
ifdef BEGIN_EPILOGUE
    BEGIN_EPILOGUE
endif  
  pop    r12
  ret
NESTED_END complete_eval_Jac512, _TEXT00


;*******************************************************************************************
; (Constant-time) Field element selection for the complete addition
; Operation: if (r9)=0 then c [r8] = a [rcx], else if (r9) = 0xff...ff then c [r8] = b [rdx]
;*******************************************************************************************
NESTED_ENTRY complete_select_Jac512, _TEXT00
  alloc_stack   8
  END_PROLOGUE

  xor           rax, rax
  ;sub           rsp, 8
  mov           [rsp], r9
  vbroadcastss  ymm0, DWORD PTR [rsp]
  vmovdqu       ymm1, YMMWORD PTR [rcx]         ; ymm1=a
  vmovdqu       ymm2, YMMWORD PTR [rdx]         ; ymm2=b
  vblendvpd     ymm3, ymm1, ymm2, ymm0          ; if ymm0=0 then ymm3=a else ymm3=b
  vmovdqu       YMMWORD PTR [r8], ymm3
  
  vmovdqu       ymm1, YMMWORD PTR [rcx+32]      ; ymm1=a
  vmovdqu       ymm2, YMMWORD PTR [rdx+32]      ; ymm2=b
  vblendvpd     ymm3, ymm1, ymm2, ymm0          ; if ymm0=0 then ymm3=a else ymm3=b
  vmovdqu       YMMWORD PTR [r8+32], ymm3  
  mov           [rsp], rax                      ; Clean stack     
  add           rsp, 8

ifdef BEGIN_EPILOGUE
    BEGIN_EPILOGUE
endif  
  ret
NESTED_END complete_select_Jac512, _TEXT00


;****************************************************************************************
; (Constant-time) Point extraction from 4-LUT for the complete mixed addition
; Operation: use index (rdx) to extract point from [rcx] and pass it to [r8]
;****************************************************************************************
NESTED_ENTRY complete_lut4_Jac512, _TEXT00
  alloc_stack  72
  END_PROLOGUE
  
  vmovdqu      YMMWORD PTR [rsp+8], ymm6      ; Save ymm6,7 in stack 
  vmovdqu      YMMWORD PTR [rsp+40], ymm7 

  xor          rax, rax  
  mov          r11, 18446744073709551615  
  ;sub          rsp, 8
  
  vmovdqu      ymm0, YMMWORD PTR [rcx]        ; Load table[0]
  vmovdqu      ymm1, YMMWORD PTR [rcx+32]
  vmovdqu      ymm2, YMMWORD PTR [rcx+64]
  vmovdqu      ymm3, YMMWORD PTR [rcx+96]        
  vmovdqu      ymm4, YMMWORD PTR [rcx+128]
  vmovdqu      ymm5, YMMWORD PTR [rcx+160]

  dec          rdx                            ; Pass over table[1]
  cmovnz       rax, r11
  mov          [rsp], rax
  vbroadcastss ymm6, DWORD PTR [rsp]
  vmovdqu      ymm7, YMMWORD PTR [rcx+192]
  vblendvpd    ymm0, ymm7, ymm0, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+224]
  vblendvpd    ymm1, ymm7, ymm1, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+256]
  vblendvpd    ymm2, ymm7, ymm2, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+288]
  vblendvpd    ymm3, ymm7, ymm3, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+320]
  vblendvpd    ymm4, ymm7, ymm4, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+352]
  vblendvpd    ymm5, ymm7, ymm5, ymm6
  
  xor          rax, rax                       ; Pass over table[2]
  dec          rdx
  cmovnz       rax, r11
  mov          [rsp], rax
  vbroadcastss ymm6, DWORD PTR [rsp]
  vmovdqu      ymm7, YMMWORD PTR [rcx+384]
  vblendvpd    ymm0, ymm7, ymm0, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+416]
  vblendvpd    ymm1, ymm7, ymm1, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+448]
  vblendvpd    ymm2, ymm7, ymm2, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+480]
  vblendvpd    ymm3, ymm7, ymm3, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+512]
  vblendvpd    ymm4, ymm7, ymm4, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+544]
  vblendvpd    ymm5, ymm7, ymm5, ymm6
  
  xor          rax, rax                       ; Pass over table[3]
  dec          rdx
  cmovz        r11, rax
  mov          [rsp], r11
  vbroadcastss ymm6, DWORD PTR [rsp]
  vmovdqu      ymm7, YMMWORD PTR [rcx+576]
  vblendvpd    ymm0, ymm7, ymm0, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+608]
  vblendvpd    ymm1, ymm7, ymm1, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+640]
  vblendvpd    ymm2, ymm7, ymm2, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+672]
  vblendvpd    ymm3, ymm7, ymm3, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+704]
  vblendvpd    ymm4, ymm7, ymm4, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+736]
  vblendvpd    ymm5, ymm7, ymm5, ymm6
    
  mov          [rsp], rax                      ; Clean stack
  vmovdqu      YMMWORD PTR [r8], ymm0          ; [r8] = table[index]
  vmovdqu      YMMWORD PTR [r8+32], ymm1
  vmovdqu      YMMWORD PTR [r8+64], ymm2
  vmovdqu      YMMWORD PTR [r8+96], ymm3
  vmovdqu      YMMWORD PTR [r8+128], ymm4
  vmovdqu      YMMWORD PTR [r8+160], ymm5
    
  vmovdqu      ymm6, YMMWORD PTR [rsp+8]       ; Restore ymm6,7
  vmovdqu      ymm7, YMMWORD PTR [rsp+40] 
  add          rsp, 72

ifdef BEGIN_EPILOGUE
    BEGIN_EPILOGUE
endif  
  ret
NESTED_END complete_lut4_Jac512, _TEXT00


;****************************************************************************************
; (Constant-time) Point extraction from 5-LUT for the complete addition
; Operation: use index (rdx) to extract point from [rcx] and pass it to [r8]
;****************************************************************************************
NESTED_ENTRY complete_lut5_Jac512, _TEXT00
  alloc_stack  72
  END_PROLOGUE
  
  vmovdqu      YMMWORD PTR [rsp+8], ymm6      ; Save ymm6,7 in stack 
  vmovdqu      YMMWORD PTR [rsp+40], ymm7 

  xor          rax, rax  
  mov          r11, 18446744073709551615  
  ;sub          rsp, 8
  
  vmovdqu      ymm0, YMMWORD PTR [rcx]        ; Load table[0]
  vmovdqu      ymm1, YMMWORD PTR [rcx+32]
  vmovdqu      ymm2, YMMWORD PTR [rcx+64]
  vmovdqu      ymm3, YMMWORD PTR [rcx+96]        
  vmovdqu      ymm4, YMMWORD PTR [rcx+128]
  vmovdqu      ymm5, YMMWORD PTR [rcx+160]

  dec          rdx                            ; Pass over table[1]
  cmovnz       rax, r11
  mov          [rsp], rax
  vbroadcastss ymm6, DWORD PTR [rsp]
  vmovdqu      ymm7, YMMWORD PTR [rcx+192]
  vblendvpd    ymm0, ymm7, ymm0, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+224]
  vblendvpd    ymm1, ymm7, ymm1, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+256]
  vblendvpd    ymm2, ymm7, ymm2, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+288]
  vblendvpd    ymm3, ymm7, ymm3, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+320]
  vblendvpd    ymm4, ymm7, ymm4, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+352]
  vblendvpd    ymm5, ymm7, ymm5, ymm6
  
  xor          rax, rax                       ; Pass over table[2]
  dec          rdx
  cmovnz       rax, r11
  mov          [rsp], rax
  vbroadcastss ymm6, DWORD PTR [rsp]
  vmovdqu      ymm7, YMMWORD PTR [rcx+384]
  vblendvpd    ymm0, ymm7, ymm0, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+416]
  vblendvpd    ymm1, ymm7, ymm1, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+448]
  vblendvpd    ymm2, ymm7, ymm2, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+480]
  vblendvpd    ymm3, ymm7, ymm3, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+512]
  vblendvpd    ymm4, ymm7, ymm4, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+544]
  vblendvpd    ymm5, ymm7, ymm5, ymm6
  
  xor          rax, rax                       ; Pass over table[3]
  dec          rdx
  cmovnz       rax, r11
  mov          [rsp], rax
  vbroadcastss ymm6, DWORD PTR [rsp]
  vmovdqu      ymm7, YMMWORD PTR [rcx+576]
  vblendvpd    ymm0, ymm7, ymm0, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+608]
  vblendvpd    ymm1, ymm7, ymm1, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+640]
  vblendvpd    ymm2, ymm7, ymm2, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+672]
  vblendvpd    ymm3, ymm7, ymm3, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+704]
  vblendvpd    ymm4, ymm7, ymm4, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+736]
  vblendvpd    ymm5, ymm7, ymm5, ymm6
  
  xor          rax, rax                       ; Pass over table[3]
  dec          rdx
  cmovz        r11, rax
  mov          [rsp], r11
  vbroadcastss ymm6, DWORD PTR [rsp]
  vmovdqu      ymm7, YMMWORD PTR [rcx+768]
  vblendvpd    ymm0, ymm7, ymm0, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+800]
  vblendvpd    ymm1, ymm7, ymm1, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+832]
  vblendvpd    ymm2, ymm7, ymm2, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+864]
  vblendvpd    ymm3, ymm7, ymm3, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+896]
  vblendvpd    ymm4, ymm7, ymm4, ymm6
  vmovdqu      ymm7, YMMWORD PTR [rcx+928]
  vblendvpd    ymm5, ymm7, ymm5, ymm6
    
  mov          [rsp], rax                      ; Clean stack
  vmovdqu      YMMWORD PTR [r8], ymm0          ; [r8] = table[index]
  vmovdqu      YMMWORD PTR [r8+32], ymm1
  vmovdqu      YMMWORD PTR [r8+64], ymm2
  vmovdqu      YMMWORD PTR [r8+96], ymm3
  vmovdqu      YMMWORD PTR [r8+128], ymm4
  vmovdqu      YMMWORD PTR [r8+160], ymm5
    
  vmovdqu      ymm6, YMMWORD PTR [rsp+8]       ; Restore ymm6,7
  vmovdqu      ymm7, YMMWORD PTR [rsp+40] 
  add          rsp, 72

ifdef BEGIN_EPILOGUE
    BEGIN_EPILOGUE
endif  
  ret
NESTED_END complete_lut5_Jac512, _TEXT00


;****************************************************************************************
; Zeroing field element
;****************************************************************************************
LEAF_ENTRY fpzero512_a, _TEXT00
  xor          rax, rax
  mov          [rcx], rax
  mov          [rcx+8], rax 
  mov          [rcx+16], rax
  mov          [rcx+24], rax 
  mov          [rcx+32], rax 
  mov          [rcx+40], rax 
  mov          [rcx+48], rax 
  mov          [rcx+56], rax 
  ret
LEAF_END fpzero512_a, _TEXT00


END