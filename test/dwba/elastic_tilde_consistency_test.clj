(ns dwba.elastic-tilde-consistency-test
  "Verify **|f_C + f_N| = |f̃_C + e^{−2iσ_0} f_N|** (same dσ) — T&N / `doc/transfer_method.tex` §total elastic amplitude."
  (:require [clojure.test :refer [deftest is testing]]
            [functions :as f :refer [mass-factor-from-mu Z1Z2ee elastic-nuclear-amplitude-fn
                                     coulomb-scattering-amplitude-thompson-nunes-eq-3181
                                     coulomb-amplitude-tilde channel-sommerfeld-eta coulomb-sigma-L]]
            [complex :as c]))

(defn- mag-add [a b]
  (c/mag (c/add a b)))

(defn- tilde-vs-nontilde-mag-diff
  [E-CM V Lcut eta k e-2s0 angle-rad]
  (let [fc2 (coulomb-scattering-amplitude-thompson-nunes-eq-3181 angle-rad eta k)
        fct2 (coulomb-amplitude-tilde angle-rad eta k)
        fn2 (elastic-nuclear-amplitude-fn E-CM V angle-rad Lcut)
        fnt2 (c/mul e-2s0 fn2)]
    (Math/abs (- (mag-add fc2 fn2) (mag-add fct2 fnt2)))))

(deftest elastic-tilde-matches-nontilde-at-pi-sm148
  (let [m-alpha 3727.379
        A 148
        m-sm (* 931.5 A)
        Z-sm 62
        E-lab 50.0
        E-CM (* E-lab (/ m-sm (+ m-alpha m-sm)))
        mu (/ (* m-alpha m-sm) (+ m-alpha m-sm))
        z12 (* 2 Z-sm 1.44)
        V [65.0 7.5 0.67]
        imag [3.0 7.5 0.67]
        Lcut 45]
    (binding [f/mass-factor (mass-factor-from-mu mu)
              f/Z1Z2ee z12
              f/*elastic-imag-ws-params* imag]
      (let [k (Math/sqrt (* f/mass-factor E-CM))
            eta (channel-sommerfeld-eta E-CM)
            theta Math/PI
            fc (coulomb-scattering-amplitude-thompson-nunes-eq-3181 theta eta k)
            fct (coulomb-amplitude-tilde theta eta k)
            fn (elastic-nuclear-amplitude-fn E-CM V theta Lcut)
            s0 (coulomb-sigma-L 0 eta)
            e-2s0 (c/complex-polar (* -2.0 s0) 1.0)
            fnt (c/mul e-2s0 fn)
            m1 (mag-add fc fn)
            m2 (mag-add fct fnt)
            ;; also spot-check 90° and 30°
            ang90 (/ Math/PI 2.0)
            ang30 (/ Math/PI 6.0)
            chk (partial tilde-vs-nontilde-mag-diff E-CM V Lcut eta k e-2s0)]
        (is (< (Math/abs (- m1 m2)) 1e-9)
            (str "|f_C+f_N| vs |tilde| at pi: " m1 " vs " m2))
        (is (< (chk ang90) 1e-9) "tilde consistency at 90°")
        (is (< (chk ang30) 1e-9) "tilde consistency at 30°")))))
