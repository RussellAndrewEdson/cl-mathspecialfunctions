;;; Copyright (c) 2015 Russell Andrew Edson
;;; 
;;; Permission is hereby granted, free of charge, to any person obtaining a
;;; copy of this software and associated documentation files (the "Software"),
;;; to deal in the Software without restriction, including without limitation
;;; the rights to use, copy, modify, merge, publish, distribute, sublicense,
;;; and/or sell copies of the Software, and to permit persons to whom the
;;; Software is furnished to do so, subject to the following conditions:
;;; 
;;; The above copyright notice and this permission notice shall be included
;;; in all copies or substantial portions of the Software.
;;; 
;;; THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
;;; OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
;;; FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
;;; AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
;;; LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
;;; FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
;;; DEALINGS IN THE SOFTWARE.


;;; Code for approximations of the gamma function.

(in-package :cl-mathspecialfunctions)

(defun nemes-approximation (z)
  "Approximates the gamma function using a scheme by Gergo Nemes."
  (* (sqrt (/ (* 2 pi) z))
     (expt (* (/ 1 (exp 1))
              (+ z (/ 1 (- (* 12 z) (/ 1 (* 10 z))))))
           z)))

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

(defun stirling-approximation (z)
  "Computes the gamma function using James Stirling's approximation."
  (* (sqrt (/ (* 2 pi) z))
     (expt (* (/ z (exp 1))
              (sqrt (+ (* z (sinh (/ 1 z))) (/ 1 (* 810 (expt z 6))))))
           z)))

(defun gamma (z &optional (approx-scheme #'nemes-approximation))
  "Returns an approximation to gamma(z), using the Nemes' scheme by default."
  (funcall approx-scheme z))
