(ns dwba.tricomis-test
  (:require [clojure.test :refer :all]
            [fastmath.special :as special]
            [fastmath.complex :as cplx]))

(deftest tricomis-U-complex-test
  (testing "tricomis-U-complex function with complex arguments"
    (let [z (cplx/complex 2 -1)  ; 2 - i
          a 2.0
          b cplx/I               ; i
          result (special/tricomis-U-complex a b z)]
      (println "tricomis-U-complex(2.0, i, 2-i) = U(a,b,z) " result)
      (is (number? (cplx/re result)))
      (is (number? (cplx/im result)))
      (is (not (Double/isNaN (cplx/re result))))
      (is (not (Double/isNaN (cplx/im result)))))))

(deftest tricomis-U-complex-additional-test
  (testing "tricomis-U-complex with different arguments"
    (let [z (cplx/complex 1 0)    ; 1 + 0i
          a 1.0
          b cplx/I               ; i
          result (special/tricomis-U-complex a b z)]
      (println "tricomis-U-complex(1.0, i, 1) = U(a,b,z) " result)
      (is (number? (cplx/re result)))
      (is (number? (cplx/im result))))))
