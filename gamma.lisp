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


;;;; Code for approximations of the gamma function.
;;;; Date: 22/09/2015

(in-package :cl-mathspecialfunctions)

;;; An approximation to the gamma function proposed by Gergo Nemes in 2007,
;;; based on the Stirling approximation.
(defun nemes-approximation (z)
  "Approximates the gamma function using a scheme by Gergo Nemes."
  (* (sqrt (/ (* 2 pi) z))
     (expt (* (/ 1 (exp 1))
              (+ z (/ 1 (- (* 12 z) (/ 1 (* 10 z))))))
           z)))

;;; A gamma function approximation devised by John Spouge in 1994.
;;; It uses a parameter 'a' to truncate the coefficient series: the relative
;;; error is bounded by a^(-1/2)*(2*pi)^(-a-1/2).
;;;
;;; (Note: The coefficients become very large as a increases, so machine
;;; round-off error quickly becomes a factor for a larger than about 13
;;; for the usual double-precision floating-point arithmetic.)
(defun spouge-approximation (z &optional (a 7))
  "John Spouge's approximation for the gamma function, using a coefficients."
  (labels ((factorial (n)
             (if (<= n 1)
                 1
                 (reduce #'* (loop for i from 1 to n collect i)))))
    (let ((c0 (sqrt (* 2 pi)))
          (ck (mapcar
               (lambda (k) (* (/ (expt -1 (1- k)) (factorial (1- k)))
                              (expt (+ (- k) a) (- k 1/2))
                              (exp (+ (- k) a))))
               (loop for k from 1 to (1- a) collect k))))
      (* (expt (+ (1- z) a) (+ (1- z) 1/2))
         (exp (- (+ (1- z) a)))
         (+ c0
            (reduce #'+
                    (mapcar (lambda (k) (/ (elt ck (1- k)) (+ (1- z) k)))
                            (loop for k from 1 to (1- a) collect k))))))))

;;; A gamma function approximation by Robert Windschitl (2002), based on the
;;; Stirling approximation. 
(defun windschitl-approximation (z)
  "Computes the gamma function using Robert Windschitl's approximation."
  (* (sqrt (/ (* 2 pi) z))
     (expt (* (/ z (exp 1))
              (sqrt (+ (* z (sinh (/ 1 z))) (/ 1 (* 810 (expt z 6))))))
           z)))

;;; The gamma function -- takes as an optional parameter a function that
;;; computes an approximation according to some scheme.
;;;
;;; Code Usage Examples:
;;;
;;;   Computing gamma(5) = factorial(4) = 24
;;;   > (gamma 5)
;;;   23.999999197528098d0
;;;   > (gamma 5 #'windschitl-approximation)
;;;   24.000014164469363d0
;;;   > (gamma 5 #'spouge-approximation)
;;;   23.99999350818005d0
;;;
;;;   Computing gamma(1/2) = sqrt(pi) = 1.7724538509055159d0
;;;   > (gamma 1/2 #'nemes-approximation)
;;;   1.7630961081013576d0
;;;   > (gamma 1/2 #'windschitl-approximation)
;;;   1.7831932384034148d0
;;;   > (gamma 1/2 (lambda (z) (spouge-approximation z 13)))
;;;   1.7724538055449004d0
(defun gamma (z &optional (approx-scheme #'nemes-approximation))
  "Returns an approximation to gamma(z), using the Nemes' scheme by default."
  (funcall approx-scheme z))
