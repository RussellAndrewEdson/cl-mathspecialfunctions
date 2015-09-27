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


;;;; Code for approximations of the Gauss hypergeometric function.
;;;; TODO
;;;; Date: 27/09/2015

(in-package :cl-mathspecialfunctions)

;;; Computes a partial summation of the hypergeometric series. The
;;; number of terms can be specified (default is 20 terms.)
;;;
;;; Note: This series only converges for |z| < 1.
(defun hypergeometric-series (a b c z &optional (num-terms 20))
  ;; Each term in the series is a ratio of factorials in a, b and c,
  ;; multiplied by a power of z and divided by a factorial. We can
  ;; keep track of termwise results and calculate the partial sum
  ;; imperatively to speed things up.
  (let ((partial-sum 1)
        (term 1))
    (loop for j from 1 to num-terms do
         (setf term (* term
                       (/ (* (+ a j -1) (+ b j -1)) (+ c j -1))
                       z
                       (/ 1 j)))
         (setf partial-sum (+ partial-sum term)))
    partial-sum))

;;; The Gauss Hypergeometric function, 2F1(a,b,c,z). Takes as an optional
;;; parameter a function that computes an approximation according to
;;; some scheme.
;;;
;;; Code Usage Examples:
;;;
;;;
(defun hypergeometric-gauss (a b c z
                             &optional
                               (approx-scheme #'hypergeometric-series))
  (funcall approx-scheme a b c z))
