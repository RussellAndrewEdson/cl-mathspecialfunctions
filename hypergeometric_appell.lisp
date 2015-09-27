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


;;;; Code for approximations of the Appell hypergeometric function.
;;;; TODO
;;;; Date: 27/09/2015

(in-package :cl-mathspecialfunctions)

;;; The Appell Hypergeometric function, F1(a,b,bb,c,x,y). We can compute
;;; this as a truncation of an infinite sum of Gauss hypergeometric
;;; functions -- the given parameter num-terms specifies how many terms
;;; of the summation to use (default is 20).
;;;
;;; Code Usage Examples:
;;;
;;;
(defun hypergeometric-appell (a b bb c x y &optional (num-terms 20))
  "Returns Appell's Hypergeometric function, F1(a,b,bb,c,x,y)."
  ;; As each term in the series contains a ratio of factorials in a, b
  ;; and c, we can keep track of the previous term each time and update
  ;; it, instead of computing factorials each time.
  (let ((partial-sum (hypergeometric-gauss a bb c y))
        (abcx-term 1))
    (loop for m from 1 to num-terms do
         (setf abcx-term (* abcx-term
                            (/ (* (+ a m -1) (+ b m -1)) (+ c m -1))
                            x
                            (/ 1 m)))
         (setf partial-sum
               (+ partial-sum
                  (* abcx-term
                     (hypergeometric-gauss (+ a m) bb (+ c m) y)))))
    partial-sum))
         
         
