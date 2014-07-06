                           MSR ECCLib V1.0 (C Edition)
                           ===========================

The MSR ECCLib V1.0 library (C Edition) implements essential elliptic curve functions 
supporting a selection of the elliptic curves proposed in [ECC14]. This library was 
developed by Microsoft Research for experimentation purposes. All functions evaluating 
secret data have regular, constant-time execution, protecting against timing and cache attacks.

The library is made available under the the Apache License, Version 2.0 (the "License"); 
you may not use these files except in compliance with the License. You may obtain a copy of 
the License at http://www.apache.org/licenses/LICENSE-2.0. Unless required by applicable law 
or agreed to in writing, software distributed under the License is distributed on an "AS IS" 
BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the 
License for the specific language governing permissions and limitations under the License.


CONTENTS:
---------

MSR_ECClib_v1.0.sln    - Visual Studio 2012 solution file
MSR_ECClib/            - Library project directory
Sample/                - Sample project 
Tests/                 - Test projects


SUPPORTED PLATFORMS:
-------------------

MSR ECCLib V1.0 is currently supported in x64 platforms with AVX (Intel's Sandy Bridge and
more recent architectures).


BUILDING THE LIBRARY:
--------------------

Open the solution file (MSR_ECClib_v1.0.sln) in Visual Studio 2012, select "x64" as Platform, 
select "Release" as configuration option and select "Build Solution" from the "Build" menu.


RUNNING THE TESTS:
-----------------

After building the solution file, run ecc_tests.exe and fp_tests.exe generated at 
<LibraryPath>\x64\Release\ecc_tests\ and <LibraryPath>\x64\Release\fp_tests\, respectively, 
from the command prompt.


RUNNING THE SAMPLE CODE:
-----------------------

After building the solution file, run sample.exe generated at <LibraryPath>\x64\Release\sample\ 
from the command prompt. Follow the sample project and the code in <LibraryPath>\Sample\sample.c 
as an example of the use of the library.


USING THE LIBRARY:
-----------------

After building the solution file, add the MSR_ECClib_v1.0.lib file generated at 
<LibraryPath>\x64\Release\ to the set of References for a project, and add msr_ecclib.h 
located at <LibraryPath>\MSR_ECClib\ to the list of Header Files of a project.


REFERENCES:
-----------

[ECC14]   Joppe W. Bos and Craig Costello and Patrick Longa and Michael Naehrig. 
          Selecting Elliptic Curves for Cryptography: An Efficiency and Security Analysis.
          Cryptology ePrint Archive: Report 2014/130, February 2014. http://eprint.iacr.org/2014/130.