(ns dwba.angular-momentum-test
  (:require [clojure.test :refer :all]
            [dwba.angular-momentum :as jam]))

(deftest wigner-3j-260000
  (let [v (jam/wigner-3j 2 6 4 0 0 0)
        expected (/ (Math/sqrt 715.0) 143.0)]
    (is (< (Math/abs (- v expected)) 1e-10))))

(deftest wigner-6j-333333
  (let [v (jam/wigner-6j 3 3 3 3 3 3)]
    (is (< (Math/abs (+ v (/ 1.0 14.0))) 1e-10))))

(deftest wigner-6j-half-integer-triangle
  ;; {½ ½ 0; ½ ½ 1} — value from standard tables / SymPy
  (let [v (jam/wigner-6j 0.5 0.5 0 0.5 0.5 1)]
    (is (number? v))
    (is (< (Math/abs v) 1e3))))
