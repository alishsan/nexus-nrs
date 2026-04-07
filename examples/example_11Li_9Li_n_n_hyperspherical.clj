#!/usr/bin/env clojure
;; 11Li treated as 9Li + n + n in a three-body hyperspherical setup.
;;
;; What this script demonstrates:
;; 1) Mass-scaled Jacobi coordinates (x, y) for the 9Li+n+n partition.
;; 2) Hyperspherical variables (rho, alpha) and the 6D volume element Jacobian.
;; 3) Jacobian determinant from particle coordinates -> (Rcm, x, y).
;; 4) A simple hyperradial effective potential estimate for K=0.
;;
;; Run:
;;   lein run -m clojure.main examples/example_11Li_9Li_n_n_hyperspherical.clj

(ns examples.example-11Li-9Li-n-n-hyperspherical)

(def ^:private hbarc 197.3269804) ; MeV*fm

;; Masses in MeV/c^2 (close to common values; good for an educational example).
(def ^:private m-n 939.5654205)
(def ^:private m-9Li (+ (* 9.0 931.49410242) 24.954)) ; 9Li atomic-mass-based approximation

(def ^:private m1 m-9Li) ; core (9Li)
(def ^:private m2 m-n)   ; neutron 1
(def ^:private m3 m-n)   ; neutron 2
(def ^:private M (+ m1 m2 m3))

;; Jacobi reduced masses for the (23)-1 partition:
;; x ~ (r3-r2), y ~ (r1 - (m2 r2 + m3 r3)/(m2+m3))
(def ^:private mu-x (/ (* m2 m3) (+ m2 m3)))
(def ^:private mu-y (/ (* m1 (+ m2 m3)) M))

;; Mass scale m for dimensionless scaling factors; choose nucleon mass.
(def ^:private m-scale m-n)

(defn- v+ [a b] (mapv + a b))
(defn- v- [a b] (mapv - a b))
(defn- v* [s v] (mapv #(* s %) v))
(defn- dot [a b] (reduce + (map * a b)))
(defn- norm [v] (Math/sqrt (double (dot v v))))

(defn jacobi-mass-scaled
  "Mass-scaled Jacobi coordinates for 9Li+n+n partition:
   x = sqrt(mu_x/m) * (r3-r2)
   y = sqrt(mu_y/m) * (r1 - (m2 r2 + m3 r3)/(m2+m3))"
  [r1 r2 r3]
  (let [s-x (Math/sqrt (/ mu-x m-scale))
        s-y (Math/sqrt (/ mu-y m-scale))
        r23cm (v* (/ 1.0 (+ m2 m3)) (v+ (v* m2 r2) (v* m3 r3)))
        x (v* s-x (v- r3 r2))
        y (v* s-y (v- r1 r23cm))]
    {:x x :y y}))

(defn hyperspherical
  "Given mass-scaled Jacobi vectors x,y:
   rho = sqrt(|x|^2 + |y|^2), alpha = atan2(|x|, |y|)."
  [{:keys [x y]}]
  (let [xnorm (norm x)
        ynorm (norm y)
        rho (Math/sqrt (+ (* xnorm xnorm) (* ynorm ynorm)))
        alpha (Math/atan2 xnorm ynorm)]
    {:rho rho :alpha alpha :xnorm xnorm :ynorm ynorm}))

(defn hyperangular-jacobian
  "Jacobian factor in 6D hyperspherical decomposition:
   d^3x d^3y = rho^5 d rho * sin^2(alpha) cos^2(alpha) d alpha * dOmega_x * dOmega_y."
  [rho alpha]
  (* (Math/pow rho 5.0)
     (Math/pow (Math/sin alpha) 2.0)
     (Math/pow (Math/cos alpha) 2.0)))

(defn- det3
  "Determinant of 3x3 matrix represented as [[a b c] [d e f] [g h i]]."
  [[[a b c] [d e f] [g h i]]]
  (+ (* a (- (* e i) (* f h)))
     (* -1.0 b (- (* d i) (* f g)))
     (* c (- (* d h) (* e g)))))

(defn particle->Rxy-det-1d
  "1D Jacobian determinant for (r1,r2,r3)->(Rcm,x,y) mapping, one Cartesian component.
   Full 3D Jacobian is this determinant cubed (same linear map on x,y,z components)."
  []
  (let [sx (Math/sqrt (/ mu-x m-scale))
        sy (Math/sqrt (/ mu-y m-scale))
        denom23 (+ m2 m3)
        a11 (/ m1 M), a12 (/ m2 M), a13 (/ m3 M)
        ;; x = sx*(r3-r2)
        a21 0.0, a22 (- sx), a23 sx
        ;; y = sy*(r1 - (m2 r2 + m3 r3)/(m2+m3))
        a31 sy
        a32 (* -1.0 sy (/ m2 denom23))
        a33 (* -1.0 sy (/ m3 denom23))]
    (det3 [[a11 a12 a13]
           [a21 a22 a23]
           [a31 a32 a33]])))

(defn simpson
  "Simple Simpson integrator for even n."
  [f a b n]
  (let [n (if (even? n) n (inc n))
        h (/ (- b a) (double n))
        xs (map #(+ a (* h %)) (range (inc n)))
        ys (map f xs)
        s0 (+ (first ys) (last ys))
        s1 (reduce + (map #(nth ys %) (range 1 n 2)))
        s2 (reduce + (map #(nth ys %) (range 2 n 2)))]
    (* (/ h 3.0) (+ s0 (* 4.0 s1) (* 2.0 s2)))))

(defn verify-hyper-jacobian
  "Checks integral of exp(-rho^2) over R^6 using hyperspherical Jacobian.
   Exact value: pi^3."
  []
  (let [rho-int (simpson (fn [rho] (* (Math/pow rho 5.0) (Math/exp (- (* rho rho)))))
                         0.0 8.0 2000)
        alpha-int (simpson (fn [a] (* (Math/pow (Math/sin a) 2.0)
                                      (Math/pow (Math/cos a) 2.0)))
                           0.0 (/ Math/PI 2.0) 2000)
        omega (* 16.0 Math/PI Math/PI)
        val (* rho-int alpha-int omega)
        exact (Math/pow Math/PI 3.0)]
    {:numeric val
     :exact exact
     :rel-error (/ (Math/abs (- val exact)) exact)}))

(defn hyperradial-centrifugal
  "Effective hypercentrifugal term:
   V_hc(rho) = hbar^2/(2m) * (K+3/2)(K+5/2)/rho^2  [MeV], m = m-scale."
  [rho K]
  (let [pref (/ (* hbarc hbarc) (* 2.0 m-scale))
        kfac (* (+ K 1.5) (+ K 2.5))]
    (* pref (/ kfac (* rho rho)))))

(defn demo-configuration []
  ;; Example geometry in fm:
  ;; core at origin, neutrons opposite along z with slight asymmetry.
  (let [r1 [0.0 0.0 0.0]
        r2 [0.0 0.0 3.8]
        r3 [0.0 0.0 -4.2]
        jac (jacobi-mass-scaled r1 r2 r3)
        hs (hyperspherical jac)]
    {:r1 r1 :r2 r2 :r3 r3 :jac jac :hs hs}))

(defn -main [& _]
  (println "=== 11Li as 9Li + n + n (hyperspherical example) ===\n")

  (println "Mass setup (MeV/c^2):")
  (println (format "  m(9Li)=%.3f, m(n)=%.6f, total M=%.3f" m1 m-n M))
  (println (format "  mu_x (n-n)=%.6f, mu_y (core-(nn))=%.6f" mu-x mu-y))
  (println)

  (let [cfg (demo-configuration)
        {:keys [x y]} (:jac cfg)
        {:keys [rho alpha xnorm ynorm]} (:hs cfg)]
    (println "Sample geometry (fm):")
    (println (format "  r_core=%s" (pr-str (:r1 cfg))))
    (println (format "  r_n1  =%s" (pr-str (:r2 cfg))))
    (println (format "  r_n2  =%s" (pr-str (:r3 cfg))))
    (println)
    (println "Mass-scaled Jacobi vectors:")
    (println (format "  x (n2-n1)      = %s  |x|=%.5f" (pr-str x) xnorm))
    (println (format "  y (core-(nn)CM)= %s  |y|=%.5f" (pr-str y) ynorm))
    (println)
    (println "Hyperspherical variables:")
    (println (format "  rho   = %.6f fm" rho))
    (println (format "  alpha = %.6f rad (%.3f deg)" alpha (* 180.0 (/ alpha Math/PI))))
    (println (format "  J_hyper(rho,alpha)=rho^5 sin^2(a) cos^2(a)=%.6f"
                     (hyperangular-jacobian rho alpha)))
    (println))

  (let [j1 (particle->Rxy-det-1d)
        j3 (Math/pow j1 3.0)]
    (println "Jacobian of linear coordinate transform:")
    (println (format "  det_1D[(r1,r2,r3)->(Rcm,x,y)] = %.8f" j1))
    (println (format "  det_3D = det_1D^3 = %.8f" j3))
    (println "  (Constant as expected for a linear transformation.)")
    (println))

  (let [{:keys [numeric exact rel-error]} (verify-hyper-jacobian)]
    (println "Hyperspherical Jacobian validation with test integral:")
    (println "  I = ∫ d^3x d^3y exp(-rho^2) should equal pi^3")
    (println (format "  numeric=%.8f, exact=%.8f, rel.error=%.3e" numeric exact rel-error))
    (println))

  (println "Simple K=0 hyperradial centrifugal term estimates:")
  (doseq [rho [2.0 4.0 6.0 8.0]]
    (println (format "  rho=%.1f fm -> V_hc=%.4f MeV" rho (hyperradial-centrifugal rho 0))))

  (println "\nDone."))

(-main)
