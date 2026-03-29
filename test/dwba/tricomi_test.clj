(ns dwba.tricomi-test
  (:require [clojure.test :refer :all]
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

(deftest tricomi-U-complex-additional-test
  (testing "tricomis-U-complex with different arguments"
    (let [z (cplx/complex 1 0)    ; 1 + 0i
          a 2.0
          b cplx/I               ; i
          result (special/tricomis-U-complex a b z)]
      (println "\n=== Additional Test ===")
      (println "special/tricomis-U-complex(2.0, cplx/I, 1) = " result)
      (println "  Real part = " (cplx/re result))
      (println "  Imaginary part = " (cplx/im result))
      (is (number? (cplx/re result)))
      (is (number? (cplx/im result))))))
