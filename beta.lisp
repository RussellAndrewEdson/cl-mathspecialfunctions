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


;;;; Code for approximations of the beta function.
;;;; Date: 22/09/2015

(in-package :cl-mathspecialfunctions)

;;; The beta function is related closely to the gamma function,
;;; with the identity B(x,y) = G(x)G(y)/G(x+y).
(defun beta-gamma-identity (x y)
  "Returns the beta function approximation using the Gamma function."
  (/ (* (gamma x) (gamma y))
     (gamma (+ x y))))

;;; The beta function can be approximated using an asymptotic formula
;;; by James Stirling, particularly when x and y are large in magnitude.
(defun stirling-approximation (x y)
  "Approximates the beta function using Stirling's approximation."
  (* (sqrt (* 2 pi))
     (expt (/ x (+ x y)) x)
     (expt (/ y (+ x y)) y)
     (expt (/ (* x y) (+ x y)) -1/2)))

;;; The beta function. Here we use the gamma function identity for
;;; 'small' values of x and y, and use the Stirling asymptotic
;;; formula for large x, y.
;;;
;;; Code Usage Examples:
;;;
;;;   Computing beta(1,1/2) = 2
;;;   > (beta 1 1/2)
;;;   1.9888264803704336d0
;;;
;;;
;;;   Computing beta(80,65) = 2.0452135320392d-44
;;;   > (beta 80 65)
;;;   2.041639548717661d-44
;;;
;;;   Computing beta(200,201) = 4.856608623805591d-122
;;;   > (beta 200 201)
;;;   4.853581858408593d-122
(defun beta (x y)
  "Returns an approximation to beta(x,y), with the gamma/Stirling formulae."
  (let ((bignum 15))
    (if (and (< x bignum) (< y bignum))
        (beta-gamma-identity x y)
        (stirling-approximation x y))))
