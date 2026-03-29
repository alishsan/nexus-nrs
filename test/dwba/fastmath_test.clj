(ns dwba.fastmath-test
  (:require [clojure.test :refer :all]
            [fastmath.core :as fm]
            [fastmath.special :as special]
            [fastmath.complex :as cplx]))

(deftest tricomi-U-complex-test
  (testing "tricomis-U-complex U(a,b,z) with complex arguments"
    (let [z (cplx/complex 2 -1)  ; 2 - i
          a 4.0
          b cplx/I              ; imaginary unit
          result (special/tricomis-U-complex a b z)]
      (println "=== Tricomis U(a,b,z) Test ===")
      (println "Input:")
      (println "  z = 2 - i = " z)
      (println "  a = " a)
      (println "  b = cplx/I = " b)
      (println "Result:")
      (println "  special/tricomis-U-complex(4.0, cplx/I, 2-i) = " result)
      (println "  Real part = " (cplx/re result))
      (println "  Imaginary part = " (cplx/im result))
      (println "  Magnitude = " (cplx/abs result))
      (println "  Phase = " (cplx/arg result) " radians")
      (is (number? (cplx/re result)))
      (is (number? (cplx/im result)))
      (is (not (Double/isNaN (cplx/re result))))
      (is (not (Double/isNaN (cplx/im result)))))))

(deftest fastmath-basic-test
  (testing "Basic fastmath functions work"
    (is (= (fm/sin 0) 0.0))
    (is (= (fm/cos 0) 1.0))
    (is (= (fm/exp 0) 1.0))))

(deftest fastmath-complex-test
  (testing "Fastmath complex number operations"
    (let [z1 (cplx/complex 1 1)
          z2 (cplx/complex 2 2)
          sum (cplx/add z1 z2)]
      (is (= (cplx/re sum) 3.0))
      (is (= (cplx/im sum) 3.0)))))
