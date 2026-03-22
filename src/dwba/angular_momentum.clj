(ns dwba.angular-momentum
  "Exact Wigner 3j / 6j symbols for integer and half-integer angular momenta.
  Algorithm ported from SymPy physics/wigner.py (Wei 1999; Rasch & Yu 2003) — same
  binomial alternating-sum structure, evaluated as double for use in DWBA factors.

  Reference: sympy.physics.wigner (BSD-licensed; Sage origin)."
  )

(defn- doubled-int
  "2×j as long; j must be integer or half-integer within tolerance."
  [x]
  (let [d (Math/round (* 2.0 (double x)))]
    (when (> (Math/abs (- (* 2.0 (double x)) (double d))) 1.0e-9)
      (throw (ex-info "expecting integer or half-integer angular momentum" {:value x})))
    (long d)))

(defn- dj-couple? [^long dj1 ^long dj2 ^long dj3]
  (and (zero? (mod (+ dj1 dj2 dj3) 2))
       (<= dj3 (+ dj1 dj2))
       (<= dj1 (+ dj2 dj3))
       (<= dj2 (+ dj1 dj3))))

(defn- comb ^clojure.lang.BigInt [^long n ^long k]
  (cond (< k 0) 0N
        (> k n) 0N
        (zero? k) 1N
        :else (let [k (long (min k (- n k)))]
                (reduce (fn [^clojure.lang.BigInt acc ^long i]
                          (/ (*' acc (-' n (-' k i))) (bigint i)))
                        1N (range 1 (inc k))))))

(defn wigner-3j
  "Wigner 3j symbol (j1 j2 j3; m1 m2 m3). Half-integers allowed. Returns double."
  [j1 j2 j3 m1 m2 m3]
  (let [dj1 (doubled-int j1) dj2 (doubled-int j2) dj3 (doubled-int j3)
        dm1 (doubled-int m1) dm2 (doubled-int m2) dm3 (doubled-int m3)]
    (cond (not= (+ dm1 dm2 dm3) 0) 0.0
          (not (dj-couple? dj1 dj2 dj3)) 0.0
          (> (Math/abs dm1) dj1) 0.0
          (> (Math/abs dm2) dj2) 0.0
          (> (Math/abs dm3) dj3) 0.0
          (not (and (zero? (mod (- dj1 dm1) 2))
                    (zero? (mod (- dj2 dm2) 2))
                    (zero? (mod (- dj3 dm3) 2)))) 0.0
          :else
          (let [sumj (quot (+ dj1 dj2 dj3) 2)
                jm1 (- sumj dj1)
                jm2 (- sumj dj2)
                jm3 (- sumj dj3)
                j1mm1 (quot (- dj1 dm1) 2)
                j2mm2 (quot (- dj2 dm2) 2)
                j3mm3 (quot (- dj3 dm3) 2)
                j1pm1 (quot (+ dj1 dm1) 2)
                imin (long (max 0 (- j1pm1 jm2) (- j2mm2 jm1)))
                imax (long (min jm3 j1pm1 j2mm2))
                sumres (loop [ii (long imin) acc 0N]
                         (if (> ii imax)
                           acc
                           (let [ti (*' (comb jm3 ii)
                                       (*' (comb jm2 (- j1pm1 ii))
                                           (comb jm1 (- j2mm2 ii))))]
                             (recur (inc ii) (-' ti acc)))))
                res-num (*' (comb dj1 jm2) (comb dj2 jm1))
                res-den (*' (comb sumj jm3) (comb dj1 j1mm1) (comb dj2 j2mm2) (comb dj3 j3mm3) (bigint (inc sumj)))]
            (if (zero? sumres)
              0.0
              (let [ressqrt (Math/sqrt (/ (double res-num) (double res-den)))
                    phase (if (even? (+ dj1 (quot (+ dj3 dm3) 2) imax)) 1.0 -1.0)]
                (* phase ressqrt (double sumres))))))))

(defn- tm->m [^long tm]
  (/ (double tm) 2.0))

(defn wigner-6j
  "Wigner 6j symbol {j1 j2 j3; j4 j5 j6}. Half-integers allowed. Returns double.
  Evaluated as Σ_m (-1)^{Σ(j-m)} × (four Wigner 3j) (Wikipedia definition).
  Matches SymPy e.g. wigner_6j(3,3,3,3,3,3) = -1/14; note {1,1,1;1,1,1} = 1/6."
  [j1 j2 j3 j4 j5 j6]
  (let [dj1 (doubled-int j1) dj2 (doubled-int j2) dj3 (doubled-int j3)
        dj4 (doubled-int j4) dj5 (doubled-int j5) dj6 (doubled-int j6)
        rng (fn [^long dj] (range (- dj) (+ dj 1) 2))]
    (cond (not (dj-couple? dj1 dj2 dj3)) 0.0
          (not (dj-couple? dj1 dj5 dj6)) 0.0
          (not (dj-couple? dj4 dj2 dj6)) 0.0
          (not (dj-couple? dj4 dj5 dj3)) 0.0
          :else
          (double
           (reduce
            +
            (for [tm1 (rng dj1) tm2 (rng dj2) tm3 (rng dj3)
                  tm4 (rng dj4) tm5 (rng dj5) tm6 (rng dj6)
                  :when (and (zero? (+ tm1 tm2 tm3))
                             (zero? (+ tm1 (- tm5) tm6))
                             (zero? (- (+ tm4 tm2) tm6))
                             (zero? (+ (- tm4) tm5 tm3)))]
              (let [m1 (tm->m tm1) m2 (tm->m tm2) m3 (tm->m tm3)
                    m4 (tm->m tm4) m5 (tm->m tm5) m6 (tm->m tm6)
                    ph (if (odd? (quot (- (+ dj1 dj2 dj3 dj4 dj5 dj6)
                                        (+ tm1 tm2 tm3 tm4 tm5 tm6))
                                     2))
                         -1.0
                         1.0)
                    t1 (wigner-3j j1 j2 j3 (- m1) (- m2) (- m3))
                    t2 (wigner-3j j1 j5 j6 m1 (- m5) m6)
                    t3 (wigner-3j j4 j2 j6 m4 m2 (- m6))
                    t4 (wigner-3j j4 j5 j3 (- m4) m5 m3)]
                (* ph t1 t2 t3 t4))))))))

(defn clebsch-gordan-exact
  "Clebsch–Gordan <j1 m1 j2 m2 | j3 m3> via 3j (Edmonds convention)."
  [j1 m1 j2 m2 j3 m3]
  (let [twice (long (+ (* 2.0 (double j1)) (* -2.0 (double j2)) (* 2.0 (double m3))))]
    (* (if (even? twice) 1.0 -1.0)
       (Math/sqrt (inc (* 2.0 (double j3))))
       (wigner-3j j1 j2 j3 m1 m2 (- m3)))))
