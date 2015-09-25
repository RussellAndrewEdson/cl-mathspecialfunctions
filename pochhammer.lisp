;;;; Copyright (c) 2015 Russell Andrew Edson
;;;; 
;;;; Permission is hereby granted, free of charge, to any person obtaining a
;;;; copy of this software and associated documentation files (the "Software"),
;;;; to deal in the Software without restriction, including without limitation
;;;; the rights to use, copy, modify, merge, publish, distribute, sublicense,
;;;; and/or sell copies of the Software, and to permit persons to whom the
;;;; Software is furnished to do so, subject to the following conditions:
;;;; 
;;;; The above copyright notice and this permission notice shall be included
;;;; in all copies or substantial portions of the Software.
;;;; 
;;;; THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
;;;; OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
;;;; FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
;;;; AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
;;;; LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
;;;; FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
;;;; DEALINGS IN THE SOFTWARE.


;;;; Code for approximations of the Pochhammer symbol (rising factorial).
;;;; Date: 25/09/2015

(in-package :cl-mathspecialfunctions)

;;; The Pochhammer symbol (rising factorial version) can be defined simply
;;; as a ratio of gamma functions:
(defun pochhammer-gamma-identity (a n)
  "Returns the Pochhammer symbol as a ratio of gamma functions."
  (/ (gamma (+ a n))
     (gamma a)))

;;; The Pochhammer symbol (rising factorial version). Here we use
;;; the gamma ratio definition for the computation.
;;;
;;; Code Usage Examples:
;;;
;;;   Computing pochhammer(3,4) = 360
;;;   > (pochhammer 3 4)
;;;   360.00092684735586d0
;;;
;;;   Computing pochhammer(17,17) = 4.150171972903141d23
;;;   > (pochhammer 17 17)
;;;   4.150178572068358d23
;;;
;;;   Computing pochhammer(-3,-4) = 1/840 = 0.0011904762
;;;   > (pochhammer -3 -4)
;;;   #C(0.001190473055379472d0 -0.0d0)
(defun pochhammer (a n)
  "Returns an approximation for the Pochhammer symbol (a)_n."
  (pochhammer-gamma-identity a n))
