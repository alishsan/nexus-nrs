(ns functions
(:require
[fastmath.core :as m]
[fastmath.polynomials :as poly]
 [fastmath.special.hypergeometric :as hg]
 [fastmath.special :as spec]
 [complex :as c]
   [fastmath.vector :as v]
   [dwba.finite-well :as fw :refer [j-l k-l j-l-deriv k-l-deriv j-ratio k-ratio
                                     finite-well-matching-error solver-step
                                     find-bound-state-finite-well find-bound-state-finite-well-nlz
                                     find-bound-state-finite-well-nlz-delta-so
                                     find-bound-state-finite-well-nlz-delta-so-over-v0
                                     find-discontinuities
                                     find-all-bound-states
                                     l-dot-s-spin-half
                                     delta-surface-thomas-so-shift-over-v0
                                     delta-surface-thomas-so-shift-MeV
                                     finite-well-delta-so-splitting-j-doublet-over-v0
                                     finite-well-delta-so-splitting-j-doublet-MeV
                                     lambda-pi-squared-fm2-default]]
))

(declare f-func f-func-deriv g-func g-func-deriv hankel0+ hankel0- phase-shift0 cm-to-lab-angle)

(def hbarc 197.7) ; MeV·fm (ℏc)

;(def mu 745) ;MeV/c^2; alpha+n
;(def mu 869.4) ; 14C+n
;(def mu 884.3); 16O + n
(def ^:dynamic *default-mu* 469.46) ; p+n reduced mass (MeV/c²)

;; Reaction-dependent: 2μ/(ℏc)²; k² = mass-factor * E. Override with (binding [mass-factor ...] ...) for other reactions.
(def ^:dynamic mass-factor (/ (* 2 *default-mu*) hbarc hbarc))

;; Reaction-dependent: Z₁Z₂e² in MeV·fm. Override with (binding [Z1Z2ee ...] ...) for other reactions.
(def ^:dynamic Z1Z2ee (* 2 1.44)) ; default: alpha + proton

(defn mass-factor-from-mu [mu-reduced]
  "Return mass-factor = 2μ/(ℏc)² for the given reduced mass μ (MeV/c²). Use with (binding [mass-factor (mass-factor-from-mu mu)] ...)."
  (/ (* 2 mu-reduced) hbarc hbarc))

;; Kinematic transformation functions
(defn lab-to-cm-energy [E-lab m1 m2]
  "Convert laboratory energy to center-of-mass energy"
  (* E-lab (/ m2 (+ m1 m2))))

(defn lab-to-cm-angle [theta-lab m1 m2]
  "Convert laboratory angle to center-of-mass angle - User's cteformula"
  (let [target (double theta-lab)
        f (fn [theta-cm] (- (cm-to-lab-angle theta-cm m1 m2) target))
        tol 1e-12
        max-iters 200
        samples 2000
        step (/ Math/PI samples)
        bisect-root
        (fn [a b]
          (loop [low a
                 high b
                 f-low (f a)
                 iter 0]
            (let [mid (* 0.5 (+ low high))
                  f-mid (f mid)]
              (if (or (>= iter max-iters)
                      (< (Math/abs f-mid) tol)
                      (< (Math/abs (- high low)) tol))
                mid
                (if (or (zero? f-mid)
                        (not= (m/signum f-low) (m/signum f-mid)))
                  (recur low mid f-low (inc iter))
                  (recur mid high f-mid (inc iter)))))))
        roots
        (reduce (fn [acc i]
                  (let [a (* i step)
                        b (* (inc i) step)
                        fa (f a)
                        fb (f b)]
                    (if (or (zero? fa)
                            (zero? fb)
                            (not= (m/signum fa) (m/signum fb)))
                      (conj acc (bisect-root a b))
                      acc)))
                []
                (range samples))
        fallback
        (reduce (fn [best i]
                  (let [x (* i step)
                        err (Math/abs (f x))]
                    (if (< err (:err best))
                      {:x x :err err}
                      best)))
                {:x 0.0 :err Double/POSITIVE_INFINITY}
                (range (inc samples)))]
    (if (seq roots)
      (let [best-root (apply min-key #(Math/abs (f %)) roots)]
        best-root)
      (:x fallback))))

(defn cm-to-lab-angle [theta-cm m1 m2]
  "Convert center-of-mass angle to laboratory angle - inverse of lab-to-cm-angle"
  (let [ratio (/ m1 m2)
        theta (Math/atan2 (Math/sin theta-cm)
                          (+ (Math/cos theta-cm) ratio))]
    (if (neg? theta) (+ theta Math/PI) theta)))

(defn jacobian-lab-to-cm [theta-lab m1 m2]
  "Calculate Jacobian for lab-to-CM transformation"
  (let [ratio (/ m1 (+ m1 m2))
        cos-theta-lab (Math/cos theta-lab)
        sin-theta-lab (Math/sin theta-lab)
        denominator (+ (* cos-theta-lab cos-theta-lab) 
                       (* ratio ratio) 
                       (- (* 2 ratio cos-theta-lab)))]
    (/ (* ratio ratio) denominator)))


(defn deriv
  ([fn1 x dx]
   (/ (- (fn1 (+ x dx)) (fn1 x)) dx))
  ([fn1 L x dx] ;for hankel functions with L dependence and complex numbers
   (c/div (c/subt2 (fn1 L (+ x dx)) (fn1 L x)) dx))

  ([fn1 L eta x dx] ;for hankel functions with L dependence and complex numbers
   (c/div (c/subt2 (fn1 L eta (+ x dx)) (fn1 L eta x)) dx))
)
(defn subtract-second [a b] [(first a) (- (second a) (second b))])

(defn to-vec2 [x] (v/vec2 (c/re x) (c/im x)))

(defn distorted-wave [k r]
  "Returns the distorted wave function at distance r for wavenumber k."
  (c/complex-polar (* k r) (* k r)))


(defn WS ;  Woods-Saxon potential V(R) = -V0/(1+exp ((r-R0)/a0))
  [r [V0 R0 a0]]
  ( / (* -1.0 V0) (+ 1.0 (Math/exp (/ (- r R0) a0))))
  )

(defn WS-complex ;  Woods-Saxon potential V(R) = -V0/(1+exp ((r-R0)/a0))
  [r [V0 R0 a0]]
  ( c/div (c/mul -1.0 V0) (+ 1.0 (Math/exp (/ (- r R0) a0))))
  )

(defn WS-complex-full [r [V0 R0 a0 W0]]
  "Complex Woods-Saxon potential: V(r) = -(V0 + iW0)/(1+exp((r-R0)/a0))"
  (let [real-part (/ (* -1.0 V0) (+ 1.0 (Math/exp (/ (- r R0) a0))))
        imag-part (/ (* -1.0 W0) (+ 1.0 (Math/exp (/ (- r R0) a0))))]
    (c/complex-cartesian real-part imag-part)))

;; --- Numerov Integration Functions ---

(defn woods-saxon-numerov [r v0 rad diff]
  "Woods-Saxon potential for Numerov method: V(r) = -V0/(1+exp((r-R0)/a0))"
  (/ (- v0) (+ 1.0 (Math/exp (/ (- r rad) diff)))))

(defn f-r-numerov ([r e l v0 rad diff mass-factor];mass-factor is different from the name-space mass factor
  "Effective potential function for Numerov integration"
  (if (zero? r)
    ;; At r=0, centrifugal term dominates: l(l+1)/r^2 -> infinity
    ;; But we never actually use r=0 in Numerov (starts at r=h)
    Double/POSITIVE_INFINITY
    (let [v-eff (+ (woods-saxon-numerov r v0 rad diff)
                   (/ (* l (inc l)) (* mass-factor r r)))]
      (* mass-factor (- v-eff e)))))

 ([r e l v0 rad diff] ; use the name-space mass-factor
  "Effective potential function for Numerov integration"
  (if (zero? r)
    ;; At r=0, centrifugal term dominates: l(l+1)/r^2 -> infinity
    ;; But we never actually use r=0 in Numerov (starts at r=h)
    Double/POSITIVE_INFINITY
    (let [v-eff (+ (woods-saxon-numerov r v0 rad diff)
                   (/ (* l (inc l)) (* mass-factor r r)))]
      (* mass-factor (- v-eff e)))))
  )

(defn thomas-spin-orbit-central-MeV
  "Thomas term matching **`dwba.transfer/optical-potential-woods-saxon`** (real part):
   **V_so · (l·s) · (1/r) · df_so/dr** with **f_so = (1 + exp((r−R_so)/a_so))^{-1}**.
   Returns MeV to add to the central Woods–Saxon **before** the reduced-mass centrifugal term."
  [r V-so R-so a-so l-dot-s]
  (if (or (< r 1e-14) (zero? l-dot-s))
    0.0
    (let [f-so (/ 1.0 (+ 1.0 (Math/exp (/ (- r R-so) a-so))))
          df-so-dr (/ (* f-so (- 1.0 f-so)) a-so)]
      (* V-so l-dot-s df-so-dr (/ 1.0 r)))))

(defn f-r-numerov-spin-orbit
  "Centrifugal + Woods–Saxon + Thomas spin-orbit; same Numerov **f** scaling as **`f-r-numerov`**.
   (No primitive type hints on arity >4 — Clojure compile limit.)"
  [r e l v0 rad diff mass-factor V-so R-so a-so l-dot-s]
  (if (zero? r)
    Double/POSITIVE_INFINITY
    (let [v-ws (woods-saxon-numerov r v0 rad diff)
          v-so (thomas-spin-orbit-central-MeV r V-so R-so a-so l-dot-s)
          v-cent (+ v-ws v-so)
          v-eff (+ v-cent (/ (* l (inc l)) (* mass-factor r r)))]
      (* mass-factor (- v-eff e)))))


(defn plot-function [f start end step & y];;"plots" function f vs. the first variable
(mapv (fn [x] [x (apply f x y)] ) (range start end step ))
  )

(defn bessel-start-l1 [r q]
  "Power series expansion for Riccati-Bessel function F1 near r=0 for l=1"
  ;; Use power series to avoid numerical underflow near origin
  ;; F1(qr) = sin(qr)/(qr) - cos(qr) ≈ (qr)²/3 - (qr)⁴/30 + ...
  (let [z (* q r)]
    ;; Power series: F1(z) ≈ z²/3 - z⁴/30 (accurate for small z, avoids underflow)
    (- (/ (* z z)3.0) 
       (/ (* z z z z) 30.0))))

(defn solve-numerov [e l v0 rad diff h r-max]
  "Solve the radial Schrödinger equation using the Numerov method"
  (let [steps (int (/ r-max h))
        q (Math/sqrt (* mass-factor (+ e v0)))
        ;; Initialize with Bessel Start
        u0 0.0
        u1 (bessel-start-l1 h q)
        
        ;; Pre-calculate f(r) values
        ;; fs[i] corresponds to f(i*h)
        ;; f(0) is infinite for l>0, but u(0)=0, so f(0)*u(0) = 0
        ;; We set f(0)=0 to avoid NaN from infinity*0
        fs (mapv (fn [r] 
                   (if (zero? r)
                     0.0  ; f(0) is infinite, but u(0)=0, so f(0)*u(0)=0 anyway
                     (f-r-numerov r e l v0 rad diff)))
                 (take (+ steps 2) (iterate #(+ % h) 0.0)))
        h2-12 (/ (* h h) 12.0)]
    
    (loop [n 1
           results [u0 u1]]
      (if (>= n (dec steps))
        results
        (let [un (get results n)        ; u at r = n*h
              un-1 (get results (dec n)) ; u at r = (n-1)*h = 0 when n=1
              ;; Numerov formula uses: f_{n-1}, f_n, f_{n+1}
              ;; where f_n = f(n*h)
              ;; When n=1: f[0]*u[0] = f(0)*0 = 0, so fs[0] value doesn't matter
              fn-1 (get fs (dec n))  ; f at r = (n-1)*h
              fn (get fs n)          ; f at r = n*h
              fn+1 (get fs (inc n))  ; f at r = (n+1)*h
              
              ;; Numerov Step:
              ;; un+1 (1 - h^2/12 fn+1) = 2un - un-1 + h^2/12 (10fn un + fn-1 un-1)
              numerator (+ (* 2.0 un) 
                           (- un-1) 
                           (* h2-12 (+ (* 10.0 fn un) (* fn-1 un-1))))
              denominator (- 1.0 (* h2-12 fn+1))
              un+1 (/ numerator denominator)]
          (recur (inc n) (conj results un+1)))))))

(defn check-wronskian [u e l v0 rad diff h]
  "Check Wronskian conservation for Numerov integration"
  ;; The paper's discrete Wronskian formula is: W_n = (1 - h²/12 f_{n+1}) u_{n+1} v_n - (1 - h²/12 f_n) u_n v_{n+1}
  ;; This is for two independent solutions u and v.
  ;;
  ;; For a single solution, we can't compute a true Wronskian (W(u,u) = 0).
  ;; However, the paper mentions checking ΔW = W_1 - W_0, which suggests they're checking
  ;; the conservation of a quantity related to the Numerov step structure.
  ;;
  ;; The Numerov algorithm preserves a symplectic structure. For a single solution,
  ;; we can check the conservation of a quantity that measures the "phase-space volume"
  ;; or the deviation from perfect conservation.
  ;;
  ;; A practical approach: Check if the quantity (1 - h²/12 f_n) u_n remains consistent,
  ;; or check the conservation of a related discrete invariant.
  ;;
  ;; Actually, let's check the paper's formula with v = u, but understand that this
  ;; measures something related to the step structure, not a true Wronskian:
  ;; W_n = (1 - h²/12 f_{n+1}) u_{n+1} u_n - (1 - h²/12 f_n) u_n u_{n+1}
  ;;     = (h²/12)(f_n - f_{n+1}) u_n u_{n+1}
  ;;
  ;; This quantity should be approximately constant (conserved) for the Numerov algorithm.
  ;; We check its drift from the first value.
  (let [h2-12 (/ (* h h) 12.0)
        ;; fs[i] corresponds to f(i*h)
        num-points (count u)
        fs (mapv #(f-r-numerov % e l v0 rad diff) 
                 (take num-points (iterate #(+ % h) 0.0)))]
    ;; Start from n=1 to avoid n=0 where u_0 = 0
    (for [n (range 1 (dec (count u)))]
      (let [un (get u n)
            un+1 (get u (inc n))
            fn (get fs n)
            fn+1 (get fs (inc n))
            ;; Discrete Wronskian-like quantity: (h²/12)(f_n - f_{n+1}) u_n u_{n+1}
            ;; This should be conserved (approximately constant)
            w-n (* h2-12 (- fn fn+1) un un+1)]
        w-n))))

(defn calculate-error [u-true u-test]
  "Calculate absolute errors between true and test wavefunction solutions"
  (mapv (fn [t tst] (Math/abs (- t tst))) u-true u-test))

(defn numerov-convergence-test [e l v0 rad diff h-fine h-test r-max]
  "Test Numerov convergence by comparing fine-grid solution with test solution"
  (let [u-true (solve-numerov e l v0 rad diff h-fine r-max)
        u-test (solve-numerov e l v0 rad diff h-test r-max)
        ;; Downsample fine solution to match test grid
        downsample-factor (int (/ h-test h-fine))
        u-true-downsampled (take-nth downsample-factor u-true)
        ;; Ensure same length (may differ by 1 due to rounding)
        min-len (min (count u-true-downsampled) (count u-test))
        errors (calculate-error (take min-len u-true-downsampled) 
                                (take min-len u-test))]
    {:max-error (apply max errors)
     :mean-error (/ (reduce + errors) (count errors))
     :errors errors
     :u-true-count (count u-true)
     :u-test-count (count u-test)
     :downsampled-count min-len}))

(defn naive-start-l1 [h]
  "Naive power series start: u(r) ≈ r^(l+1) for l=1"
  (Math/pow h 2))

(defn solve-numerov-naive [e l v0 rad diff h r-max]
  "Solve Numerov with naive r^(l+1) power series start"
  (let [steps (int (/ r-max h))
        ;; Initialize with naive start
        u0 0.0
        u1 (Math/pow h (inc l))
        
        ;; Pre-calculate f(r) values
        ;; fs[i] corresponds to f(i*h)
        ;; f(0) is infinite for l>0, but u(0)=0, so f(0)*u(0) = 0
        ;; We set f(0)=0 to avoid NaN from infinity*0
        fs (mapv (fn [r] 
                   (if (zero? r)
                     0.0  ; f(0) is infinite, but u(0)=0, so f(0)*u(0)=0 anyway
                     (f-r-numerov r e l v0 rad diff)))
                 (take (+ steps 2) (iterate #(+ % h) 0.0)))
        h2-12 (/ (* h h) 12.0)]
    
    (loop [n 1
           results [u0 u1]]
      (if (>= n (dec steps))
        results
        (let [un (get results n)        ; u at r = n*h
              un-1 (get results (dec n)) ; u at r = (n-1)*h
              ;; Numerov formula uses: f_{n-1}, f_n, f_{n+1}
              ;; where f_n = f(n*h)
              ;; When n=1: f[0]*u[0] = f(0)*0 = 0, so fs[0] value doesn't matter
              fn-1 (get fs (dec n))  ; f at r = (n-1)*h
              fn (get fs n)          ; f at r = n*h
              fn+1 (get fs (inc n))  ; f at r = (n+1)*h
              
              ;; Numerov Step:
              ;; un+1 (1 - h^2/12 fn+1) = 2un - un-1 + h^2/12 (10fn un + fn-1 un-1)
              numerator (+ (* 2.0 un) 
                           (- un-1) 
                           (* h2-12 (+ (* 10.0 fn un) (* fn-1 un-1))))
              denominator (- 1.0 (* h2-12 fn+1))
              un+1 (/ numerator denominator)]
          (recur (inc n) (conj results un+1)))))))

(defn calculate-stability-data [e l v0 rad diff h r-max]
  "Calculate Wronskian stability comparison between Bessel start and naive start"
  (let [;; Bessel Start Integration
        u-bessel (solve-numerov e l v0 rad diff h r-max)
        
        ;; Naive Power Start Integration (u1 = h^(l+1))
        u-naive (solve-numerov-naive e l v0 rad diff h r-max)
        
        ;; Discrete Wronskians
        w-bessel (check-wronskian u-bessel e l v0 rad diff h)
        w-naive (check-wronskian u-naive e l v0 rad diff h)]
    
    {:bessel-w-drift (if (seq w-bessel)
                       (apply max (map #(Math/abs (- % (first w-bessel))) w-bessel))
                       0.0)
     :naive-w-drift (if (seq w-naive)
                      (apply max (map #(Math/abs (- % (first w-naive))) w-naive))
                      0.0)
     :bessel-wronskian w-bessel
     :naive-wronskian w-naive
     :bessel-w-initial (first w-bessel)
     :naive-w-initial (first w-naive)}))

;; --- Phase Shift Extraction from Numerov Solution ---

(defn extract-wavefunction-at-boundary [u h r-boundary]
  "Extract wavefunction value and derivative at boundary from Numerov solution"
  (let [idx (int (/ r-boundary h))
        idx (min idx (- (count u) 2))  ; Ensure we can calculate derivative
        ;; Use the actual r value at this index, not r-boundary
        r-actual (* idx h)
        u-a (get u idx)
        ;; Central difference for derivative (more accurate)
        u-prime-a (if (and (> idx 0) (< idx (dec (count u))))
                    (/ (- (get u (inc idx)) (get u (dec idx))) (* 2 h))
                    ;; Forward/backward difference at boundaries
                    (if (zero? idx)
                      (/ (- (get u (inc idx)) u-a) h)
                      (/ (- u-a (get u (dec idx))) h)))]
    {:u u-a
     :u-prime u-prime-a
     :r r-actual  ; Use actual r value, not requested boundary
     :index idx}))

(defn r-matrix-from-numerov [u h r-boundary]
  "Calculate R-matrix from Numerov wavefunction solution"
  (let [{:keys [u u-prime r]} (extract-wavefunction-at-boundary u h r-boundary)]
    ;; R = u(a) / (a * u'(a))
    ;; Use actual r value from extraction, not requested boundary
    (if (or (zero? u-prime) (Double/isNaN u-prime) (Double/isInfinite u-prime))
      Double/NaN
      (/ u (* r u-prime)))))

(defn interpolate [wave-function r r-max]
  "Interpolate wave function value at radius r from discrete wave-function values.
   
   Parameters:
   - wave-function: Vector/list of wave-function values at discrete points
   - r: Radius at which to interpolate (must be between 0 and r-max)
   - r-max: Maximum radius corresponding to the last wave-function value
   
   Returns: Interpolated wave-function value at radius r
   
   Uses linear interpolation between the two nearest grid points.
   Assumes the first value corresponds to r=0 and the last to r=r-max."
  (let [n-points (count wave-function)]
    (cond
      ;; Edge case: empty or single point
      (<= n-points 1)
      (if (empty? wave-function)
        0.0
        (first wave-function))
      ;; If r is exactly at a grid point or beyond r-max
      (<= r 0.0)
      (first wave-function)
      (>= r r-max)
      (last wave-function)
      ;; Linear interpolation between two nearest points
      :else
      (let [;; Calculate step size: h = r-max / (n-points - 1)
            ;; since we have n-points values from r=0 to r=r-max
            h (/ r-max (dec n-points))
            ;; Find the index where r falls: idx = r / h
            idx (/ r h)
            idx-floor (int (Math/floor idx))
            idx-ceil (min (inc idx-floor) (dec n-points))
            ;; Handle edge cases
            idx-floor (max 0 (min idx-floor (dec n-points)))
            r-floor (* idx-floor h)
            r-ceil (* idx-ceil h)
            u-floor (get wave-function idx-floor)
            u-ceil (get wave-function idx-ceil)
            ;; Linear interpolation: u(r) = u_floor + (u_ceil - u_floor) * (r - r_floor) / (r_ceil - r_floor)
            weight (if (zero? (- r-ceil r-floor))
                     0.0
                     (/ (- r r-floor) (- r-ceil r-floor)))]
        (+ u-floor (* weight (- u-ceil u-floor)))))))

(defn phase-shift-from-numerov [u h r-boundary e l]
  "Extract phase shift from Numerov wavefunction solution using S-matrix method (same as phase-shift0)"
  (let [{:keys [u u-prime r]} (extract-wavefunction-at-boundary u h r-boundary)
        ;; Calculate Ra = u/u' (R-matrix times a, same as r-matrix-a returns)
        ;; r-matrix-a returns: Ra = u/dudr where dudr = u'
        Ra (/ u u-prime)
        k (m/sqrt (* mass-factor e))
        rho (* k r)
        ;; Use S-matrix method (same as s-matrix0 and phase-shift0)
        ;; S = (H- - Ra*k*H-') / (H+ - Ra*k*H+')
        ;; where H+ and H- are outgoing/incoming Hankel functions
        hankel-minus (hankel0- l rho)
        hankel-plus (hankel0+ l rho)
        hankel-minus-prime (deriv hankel0- l rho 0.000001)
        hankel-plus-prime (deriv hankel0+ l rho 0.000001)
        numerator (c/subt2 hankel-minus (c/mul Ra k hankel-minus-prime))
        denominator (c/subt2 hankel-plus (c/mul Ra k hankel-plus-prime))
        s-matrix (c/div numerator denominator)]
    (if (or (Double/isNaN Ra) (Double/isInfinite Ra)
            (Double/isNaN (c/arg s-matrix)))
      Double/NaN
      (/ (c/arg s-matrix) 2))))

(defn exact-phase-shift-numerov [e l v0 rad diff r-boundary]
  "Calculate 'exact' phase shift using very fine Numerov integration with Bessel start"
  ;; Use very fine grid with Bessel start as reference (matches paper methodology)
  ;; Paper uses h=0.001 or finer as "exact" reference for convergence comparison
  ;; This measures how quickly coarse-grid solutions converge to the fine-grid solution
  (let [h-fine 0.0001  ; Very fine grid for reference (finer than paper's 0.001 for better accuracy)
        u (solve-numerov e l v0 rad diff h-fine r-boundary)]
    (phase-shift-from-numerov u h-fine r-boundary e l)))

(defn phase-shift-convergence-table [e l v0 rad diff h-values r-boundary]
  "Generate phase shift convergence table matching paper format"
  ;; Use R-matrix method as exact reference (independent of Numerov initialization)
  (let [V [v0 rad diff]
        delta-exact (phase-shift0 e V r-boundary l)]
    (mapv (fn [h]
            (let [u-bessel (solve-numerov e l v0 rad diff h r-boundary)
                  u-naive (solve-numerov-naive e l v0 rad diff h r-boundary)
                  delta-bessel (phase-shift-from-numerov u-bessel h r-boundary e l)
                  delta-naive (phase-shift-from-numerov u-naive h r-boundary e l)]
              {:h h
               :naive-error (Math/abs (- delta-naive delta-exact))
               :bessel-error (Math/abs (- delta-bessel delta-exact))
               ;; Direct difference between the two Numerov schemes
               :scheme-diff (Math/abs (- delta-naive delta-bessel))
               ;; Ratio of naive to Bessel error (may be NaN/Inf if bessel-error ~ 0)
               :error-ratio (if (and (pos? (Math/abs (- delta-bessel delta-exact)))
                                     (Double/isFinite (Math/abs (- delta-naive delta-exact))))
                              (/ (Math/abs (- delta-naive delta-exact))
                                 (Math/abs (- delta-bessel delta-exact)))
                              Double/NaN)
               :naive-phase-shift delta-naive
               :bessel-phase-shift delta-bessel
               :exact-phase-shift delta-exact}))
          h-values)))

(defn print-convergence-table [table]
  "Print phase shift convergence table in paper format"
  (println "Phase Shift Error Convergence |δ_calc - δ_exact|")
  (println (apply str (repeat 60 "-")))
  (println (format "%-8s %-22s %-22s %-22s %-12s"
                   "h (fm)"
                   "Naive Start Error"
                   "Bessel-Start Error"
                   "|δ_naive - δ_bessel|"
                   "Naive/Bessel"))
  (println (apply str (repeat 110 "-")))
  (doseq [row table]
    (println
     (format "%-8.2f %-22.12e %-22.12e %-22.12e %-12.4f"
             (:h row)
             (:naive-error row)
             (:bessel-error row)
             (:scheme-diff row)
             (let [r (double (or (:error-ratio row) Double/NaN))]
               (if (Double/isNaN r) 0.0 r)))))
  (println ""))



;; Example of a single Numerov step with complex numbers
(defn numerov-step-complex [u-n u-prev f-prev f-n f-next h]
  (let [h2-12 (/ (* h h) 12.0)
        ;; Numerator: 2*u_n - u_{n-1} + (h^2/12)*(10*f_n*u_n + f_{n-1}*u_{n-1})
        term1 (c/subt2 (c/mul 2.0 u-n) u-prev)
        term2 (c/mul h2-12 (c/add (c/mul 10.0 (c/mul f-n u-n)) 
                                   (c/mul f-prev u-prev)))
        numerator (c/add term1 term2)
        ;; Denominator: 1 - (h^2/12)*f_{n+1}
        denominator (c/subt2 1.0 (c/mul h2-12 f-next))]
    (c/div numerator denominator)))

(defn- numerov-step-complex-guarded
  "Like **`numerov-step-complex`**, but nudges a vanishing denominator so high-**L** / large-|f| grids do not halt **`r-matrix-complex-imag-ws`** with **divide-by-zero**."
  [u-n u-prev f-prev f-n f-next h]
  (let [h2-12 (/ (* h h) 12.0)
        term1 (c/subt2 (c/mul 2.0 u-n) u-prev)
        term2 (c/mul h2-12 (c/add (c/mul 10.0 (c/mul f-n u-n))
                                  (c/mul f-prev u-prev)))
        numerator (c/add term1 term2)
        denominator (c/subt2 1.0 (c/mul h2-12 f-next))
        mag (double (c/mag denominator))
        denominator' (if (< mag 1e-14)
                       (c/add denominator (c/complex-cartesian 1e-14 0.0))
                       denominator)]
    (c/div numerator denominator')))

(def ^:private numerov-renorm-mag-lo
  "If **max |u|** on the integrated samples falls below this, multiply **every** **u** by the same factor (**R = u/(a u′)** unchanged)."

  1e-200)

(def ^:private numerov-renorm-mag-hi
  "If **max |u|** exceeds this, scale the whole solution vector down jointly."

  1e200)

(defn numerov-append-and-renormalize
  "Append **u-next** to the Numerov **u** vector and keep **max |u|** inside **[lo, hi]** by multiplying
  **all** samples by one real factor when needed. The reduced equation **u″ = f(r) u** is linear in **u**
  (Coulomb only in **V**), so **u/(a u′)** is unchanged by a joint rescale.

  Use only where **S** (or **R**) is the goal: repeated rescales accumulate a different factor per integration, which
  breaks meaningful **|**u(r)**|** / cross-**L** scaling for DWBA **`:raw`** waves. **`distorted-wave-optical`** omits this on purpose.

  **magmx** is the running **max |u|** before **u-next** is appended. Returns **`[u-vec′ magmx′]`** after any rescale."
  ([u-vec u-next magmx]
   (numerov-append-and-renormalize u-vec u-next numerov-renorm-mag-lo numerov-renorm-mag-hi magmx))
  ([u-vec u-next lo hi magmx]
   (let [lo (double lo)
         hi (double hi)
         magmx (double magmx)
         u-vec' (conj u-vec u-next)
         magmx' (Math/max magmx (double (c/mag u-next)))]
     (cond
       (or (Double/isNaN magmx') (Double/isInfinite magmx'))
       [u-vec' magmx']

       (zero? magmx')
       [u-vec' magmx']

       (> magmx' hi)
       (let [fac (/ hi magmx')]
         [(mapv #(c/mul (c/complex-cartesian fac 0.0) %) u-vec')
          hi])

       (< magmx' lo)
       (let [fac (/ lo magmx')]
         [(mapv #(c/mul (c/complex-cartesian fac 0.0) %) u-vec')
          lo])

       :else
       [u-vec' magmx']))))

(defn xi ;h-bar and speed of light c are set to 1
  [^double E V  ^double a ^long L]  ;no coulomb ;construct R-matrix * a depending on 1D Woods-Saxon potential V(R) = -V0/(1+exp ((r-R0)/a0)) V = [V0, R0, a0]
;choose u(0) = 0, u'(0) = 1 for initial conditions
(let [
dr 0.0000100]
  (loop [x dr pot 0 d2udr2 -0.1 
dudr 1
 ur (* dudr dr) accum []]
(if (> x a)
  (take-nth 1000 accum)
(recur  (+ x dr) (WS x V) (*  (+ (/ (* L (inc L)) (* x x)) (* mass-factor (-  pot E))) ur) (+ dudr (* d2udr2 dr))  (+ ur (* dudr dr)) (into accum [[x ur ]])) ) 
)))

(defn r-matrix-a 
  [^double E V  ^double a ^long L]  ;no coulomb ;construct R-matrix * a depending on 1D Woods-Saxon potential V(R) = -V0/(1+exp ((r-R0)/a0)) V = [V0, R0, a0]
;choose u(0) = 0, u'(0) = 1 for initial conditions
(let [
dr 0.001]
(loop [x dr pot 0 d2udr2 (/ 1. dr) dudr 1 ur dr]
(if (> x a)
(/ ur dudr) ;Ra = u/dudr  (dudr = d2udr2 * dr)
(recur  (+ x dr) (WS x V) (*  (+ (/ (* L (inc L)) (* x x)) (* mass-factor (-  pot E))) ur) (+ dudr (* d2udr2 dr))  (+ ur (* dudr dr))) ) 
)))

(defn r-matrix-step 
  [^double E ^double V ^double a ^long L]  ;no coulomb ;construct R-matrix 
  (let [  k (m/sqrt (*  mass-factor (- E V)))
rho (* k a) ]
    (/
     (f-func L rho)
     (f-func-deriv L rho)
     rho
)
)
)

(defn g-matrix-step ;ratio of G to G' , divided by a, same as G/rho G dot
  [^double E  ^double a ^long L]  ;no coulomb ;construct R-matrix 
  (let [  k (m/sqrt (*  mass-factor  E ))
rho (* k a) ]
    (/
     (g-func L rho)
     (g-func-deriv L rho)
     rho
)
)
)


(defn Coulomb-pot ([ r r0] ; potential of a uniformly charged sphere of charge e and radius r0
                   (if (> r r0) (/  Z1Z2ee r) (* r (/ Z1Z2ee r0 r0))))
   ([Ze r r0] ; potential of a uniformly charged sphere of charge e and radius r0
  (if (> r r0) (/  Ze r) (* r (/ Ze r0 r0))))
  ) 
  


(defn f-func [L rho] (* rho (spec/spherical-bessel-j L rho)))
(defn g-func [L rho] (* -1 rho (spec/spherical-bessel-y L rho)))

;; --- Finite Well Functions ---
;; These functions have been moved to dwba.finite-well namespace.
;; They are re-exported here for backward compatibility.

(defn bisection

  ([f [low high]]
(bisection f [low high] 1.e20 100)
   )
  ([f [low high] tolerance max-iters]
  "Generic bisection root-finding algorithm.
   
   Parameters:
   - f: Function to find root of (should change sign between low and high)
   - low: Lower bound of search interval
   - high: Upper bound of search interval
   - tolerance: Convergence tolerance (stop when |high - low| < tolerance)
   - max-iters: Maximum number of iterations
   
   Returns: {:root, :value, :iterations, :converged?}
   - root: The root found
   - value: f(root) - should be close to 0 (check |value| < tolerance for precision)
   - iterations: Number of iterations used
   - converged?: Whether convergence was achieved"
  (let [f-low (f low)
        f-high (f high)
        ;; Keep interval tolerance bounded so legacy callers that pass a huge tolerance
        ;; (to focus on function-value checks) still iterate meaningfully.
        interval-tol (min 1.0e-5 (Math/abs tolerance))]
    (if (= (m/signum f-low) (m/signum f-high))
      {:root low
       :value f-low
       :iterations 0
       :converged? false
       :error "Function has same sign at both endpoints"}
      (loop [low low
             high high
             iter 0]
        (if (or (>= iter max-iters)
                (< (Math/abs (- high low)) interval-tol))
          (let [mid (/ (+ low high) 2.0)
                f-mid (f mid)]
            {:root mid
             :value f-mid
             :iterations iter
             :converged? (or (< (Math/abs f-mid) tolerance)
                             (< (Math/abs (- high low)) interval-tol))})
          (let [mid (/ (+ low high) 2.0)
                f-mid (f mid)
                f-low (f low)]
            (if (= (m/signum f-low) (m/signum f-mid))
              (recur mid high (inc iter))
              (recur low mid (inc iter))))))))))

(defn secant [f x0 x1 tolerance max-iters]
  "Secant method for root finding.
   
   The secant method is faster than bisection for smooth functions, but doesn't
   guarantee convergence. It uses linear interpolation between two points to
   approximate the root.
   
   Parameters:
   - f: Function to find root of
   - x0: First initial guess
   - x1: Second initial guess (should be different from x0)
   - tolerance: Convergence tolerance (stop when |f(x)| < tolerance or |x_n - x_{n-1}| < tolerance)
   - max-iters: Maximum number of iterations
   
   Returns: {:root, :value, :iterations, :converged?}
   - root: The root found
   - value: f(root) - should be close to 0
   - iterations: Number of iterations used
   - converged?: Whether convergence was achieved
   
   Algorithm:
   x_{n+1} = x_n - f(x_n) * (x_n - x_{n-1}) / (f(x_n) - f(x_{n-1}))
   
   Example:
   (secant #(- (* % %) 4) 1.0 3.0 1e-10 100)
   => Finds root of x² - 4 = 0 (should be x = 2)"
  (let [f-x0 (f x0)
        f-x1 (f x1)]
    ;; Check for nil values
    (if (or (nil? f-x0) (nil? f-x1))
      {:root x0
       :value f-x0
       :iterations 0
       :converged? false
       :error "Function returned nil"}
      ;; Check if we already have a root
      (if (< (Math/abs f-x0) tolerance)
        {:root x0
         :value f-x0
         :iterations 0
         :converged? true}
        (if (< (Math/abs f-x1) tolerance)
          {:root x1
           :value f-x1
           :iterations 0
           :converged? true}
          ;; Check if denominator would be zero (function values are equal)
          (if (< (Math/abs (- f-x1 f-x0)) 1e-15)
            {:root x1
             :value f-x1
             :iterations 0
             :converged? false
             :error "Function values at initial guesses are too close"}
            ;; Iterate using secant method
            (loop [x-prev x0
                   x-curr x1
                   f-prev f-x0
                   f-curr f-x1
                   iter 0]
              (if (or (>= iter max-iters)
                      (< (Math/abs f-curr) tolerance)
                      (< (Math/abs (- x-curr x-prev)) tolerance))
                {:root x-curr
                 :value f-curr
                 :iterations iter
                 :converged? (< (Math/abs f-curr) tolerance)}
                (let [;; Secant method formula
                      denominator (- f-curr f-prev)
                      ;; Avoid division by zero
                      x-next (if (< (Math/abs denominator) 1e-15)
                               x-curr  ; If denominator too small, don't update
                               (- x-curr (* f-curr (/ (- x-curr x-prev) denominator))))
                      f-next (f x-next)]
                  (if (nil? f-next)
                    {:root x-curr
                     :value f-curr
                     :iterations iter
                     :converged? false
                     :error "Function returned nil during iteration"}
                    (recur x-curr x-next f-curr f-next (inc iter))))))))))))

;; Finite-well functions removed - now in dwba.finite-well namespace

(defn hankel0+ [L rho]
  (c/complex-cartesian (g-func L rho ) (f-func L rho))
  )

(defn hankel0- [L rho]
  (c/complex-cartesian (g-func L rho ) (* -1.0 (f-func L rho)))
  )

(defn f-func-deriv [^long L ^double rho] 
  (* 0.5 (+ (f-func (dec L) rho) 
(/ (f-func L rho) rho)
 (* -1 (f-func (inc L) rho))
)))

(defn g-func-deriv [^long L ^double rho] 
  (* 0.5 (+ (g-func (dec L) rho) 
(/ (g-func L rho) rho)
 (* -1 (g-func (inc L) rho))
)))


(defn hankel0+-deriv [L rho]
  (c/complex-cartesian (g-func-deriv L rho ) (f-func-deriv L rho))
  )

(defn hankel0--deriv [L rho]
  (c/complex-cartesian (g-func-deriv L rho ) (* -1.0 (f-func-deriv L rho)))
  )


(defn CL [L ^double eta] 
(->   (* Math/PI eta -0.5)
(Math/exp)
(* (Math/pow 2 L))
(c/mul (c/mag (c/gamma-complex   (c/complex-cartesian (inc L) eta))))
(c/div (m/factorial (inc L)))
)
)
;(defn hypergeometric-1F1)

(defn delta1 [k1 k2 R]; phase shift for L=1, square well
 (let [rho1 (* k1 R) rho2 (* k2 R)]   (- (m/acot (+ (* (/ rho2 rho1) (m/cot rho1 )) (/ 1. rho2) (/ (* -1. rho2) rho1 rho1))) rho2)))

(defn factorial [n]
(if (= n 1) 1 (* (factorial (- n 1)))
))

(defn s-matrix0 [^double E V  ^double a ^long L]
  (let [ra (r-matrix-a E V a L)
         k (m/sqrt (*  mass-factor E))
        rho (* k a)]
    (c/div (c/subt2 (hankel0- L rho) (c/mul ra k (deriv hankel0- L rho 0.000001)) )
         (c/subt2 (hankel0+ L rho) (c/mul ra k (deriv hankel0+ L rho 0.000001)) ))
 ))

(defn s-matrix-step [^double E ^double V  ^double a ^long L]
  (let [rm (r-matrix-step E V a L)
         k (m/sqrt (*  mass-factor E))
        rho (* k a)]
    (c/div (c/subt2 (hankel0- L rho) (c/mul rm rho (hankel0--deriv L rho )) )
         (c/subt2 (hankel0+ L rho) (c/mul rm rho (hankel0+-deriv L rho )) ))
 ))

(defn phase-shift-step  [^double E ^double V  ^double a ^long L ]
  (let [s (s-matrix-step E V a L)]
    (/ (c/arg s) 2)
    ))

(defn k-matrix0 [^double E V  ^double a ^long L]
  (let [ra (r-matrix-a E V a L)]
    (c/div (c/subt2  (c/mul ra (deriv f-func L a 0.0000001))  (f-func L a))
         (c/subt2 (g-func L a) (c/mul ra (deriv g-func L a 0.0000001)) ))
 ))

(defn phase-shift0  [^double E V  ^double a ^long L ]
  (let [s (s-matrix0 E V a L)]
    (/ (c/arg s) 2)
    ))

(defn sigma-L0 [E V a L]
(* (/ 2 E) Math/PI (+ (* 2 L) 1) (Math/pow (c/mag (c/subt2 1. (s-matrix0 E V a L)))  2) )
)


          

(defn rising-factorial
  "Rising (Pochhammer) factorial."
  [n x]
    (c/div (c/gamma-complex (c/add x n))
             (c/gamma-complex x)))


(defn pocn
  [ac b z n]
   (c/div (c/mul (c/npow z n) (rising-factorial n ac)) (c/mul (rising-factorial n b) (m/factorial n)))
)

(defn hypergeometric-complex-1F1
  "Kummer's (confluent hypergeometric, 1F1) function for compex arguments."
  [ac b z]
 (->> (range 20) (map #(pocn ac b z %)) (reduce c/add))   
)


(defn pocn2F0 ;used for hypergeometric-2F0
  [a1 a2 z n]
  (c/div  (c/mul  (c/npow z n)  (rising-factorial n a1) (rising-factorial n a2)) (m/factorial n))
)

(defn hypergeometric-complex-2F0
[a1 a2 z]
 ;(->> (range 20) (map #(pocn2F0 a1 a2 z %)) (reduce add))   
  (spec/hypergeometric-pFq-complex [(to-vec2 a1) (to-vec2 a2)] [] (to-vec2 z))
  )

(defn coulomb-F [L eta r]
(c/mul (CL L eta) (m/pow r (inc L)) (c/complex-polar (- r) 1)
     (hypergeometric-complex-1F1 (c/complex-cartesian (inc L) (- eta)) (* 2 (inc L)) (c/complex-cartesian 0 (* 2 r) )
                                 )
     ))

(defn hypergeometric-complex-U
  [a b z]
  (c/div 
   (apply c/complex-cartesian  (hypergeometric-complex-2F0 a ( c/subt  a b -1.)  (c/div -1. z) ))
    (c/cpowc  z a))    
  )

(defn Hankel+ [L, eta, rho]
  (let [sigmal (c/arg (apply c/complex-cartesian (spec/gamma-complex (v/vec2 (inc L) eta))))
        theta (+ rho (* L Math/PI -0.5) sigmal (* eta -1.0 (Math/log (* 2 rho))))
  a (c/complex-cartesian (inc L)  eta)
        ]
    (c/mul (c/complex-polar theta 1)
            (c/cpowc (c/complex-cartesian 0 (* -2 rho)) a) 
  ( hypergeometric-complex-U a (* 2 (inc L)) (c/complex-cartesian 0 (* -2.0 rho)))
)))

(defn Hankel- [L, eta, rho]
   (let [sigmal (c/arg (apply c/complex-cartesian (spec/gamma-complex (v/vec2 (inc L) eta))))
         theta (+ rho (* L Math/PI -0.5) sigmal (* eta -1.0 (Math/log (* 2 rho))))
         a (c/complex-cartesian (inc L) (* -1.0 eta))]
     (c/mul (c/complex-polar (* -1.0 theta) 1)
          (c/cpowc (c/complex-cartesian 0 (* 2 rho)) a) 
  ( hypergeometric-complex-U a (* 2 (inc L)) (c/complex-cartesian 0 (* 2.0 rho)))
)))

(defn coulomb-sigma-L
  "Coulomb partial-wave phase **σ_L(η) = arg Γ(L+1+iη)** (radians), same ingredient as **`Hankel+`** / **`Hankel-`**.
  Use for Austern **Eq. (5.6)** **e^{i(σ_{αL_α}+σ_{βL_β})}** when the **σ** are taken as **Coulomb** phases.

  **η** must match **`s-matrix`**: with **k = √(mass-factor·E)**, use **`channel-sommerfeld-eta`**. **L** ≥ 0."
  ^double [^long L ^double eta]
  (double (c/arg (apply c/complex-cartesian (spec/gamma-complex (v/vec2 (inc (double L)) eta))))))

(defn coulomb-phase-diff
  "**e^{2i(σ_L − σ_0)}** as a product of **L** simple complex ratios — **no Gamma evaluations**.

  **σ_L − σ_0 = Σ_{k=1}^{L} arg(k+iη)**, so **e^{2i(σ_L−σ_0)} = ∏_{k=1}^{L} (k+iη)/(k−iη)** (each factor is a pure phase of magnitude 1). Returns **1+0i** for **L = 0**.

  Note: the conjugate **∏(k−iη)/(k+iη)** = **e^{2i(σ_0−σ_L)}** arises in some references that define **σ** with the opposite sign or treat attractive Coulomb (**η < 0**) as the repulsive baseline."
  [^long L ^double eta]
  (loop [acc (c/complex-cartesian 1.0 0.0)
         k   (long 1)]
    (if (> k (long L))
      acc
      (let [norm-sq (+ (* (double k) (double k)) (* eta eta))]
        (recur (c/mul acc (c/complex-cartesian (/ (- (* (double k) (double k)) (* eta eta)) norm-sq)
                                               (/ (* 2.0 (double k) eta) norm-sq)))
               (inc k))))))

(defn coulomb-amplitude-tilde
  "**f̃_C(θ) = e^{−2iσ_0} f_C(θ)** (**fm**): Coulomb amplitude with the global **σ_0** phase removed.

  **f̃_C = −η/(2k sin²(θ/2)) · exp(−iη ln sin²(θ/2))** — no Gamma evaluations.

  Use together with **`coulomb-phase-diff`** to form **|f_C + f_N|² = |f̃_C + f̃_N|²**, where
  **f̃_N = Σ_{L=0}^{L_max} (−i/2k)(2L+1) P_L · e^{2i(σ_L−σ_0)} · (S^n_L − 1)**
  and **e^{2i(σ_L−σ_0)}** = **`coulomb-phase-diff`** (pure product, no Gamma)."
  [^double theta-rad ^double eta ^double k]
  (let [s2    (Math/sin (* 0.5 theta-rad))
        sinsq (max (* s2 s2) 1e-300)
        mag   (/ (Math/abs eta) (* 2.0 k sinsq))
        phase (* -1.0 eta (Math/log sinsq))]
    (c/complex-polar (+ phase Math/PI) mag)))

(defn channel-sommerfeld-eta
  "**η(E)** for current **`mass-factor`** (2μ/ℏ² in 1/(MeV·fm²)) and **`Z1Z2ee`** (MeV·fm), matching **`s-matrix`**:
  **k = √(mass-factor·E)**, **η = Z₁Z₂e² · mass-factor / (2k)**. **E** = CM energy (MeV).

  Bind **`mass-factor`** and **`Z1Z2ee`** separately for entrance (α) and exit (β) channels if they differ."
  ^double [^double E]
  (let [k (m/sqrt (* mass-factor E))]
    (* Z1Z2ee mass-factor (/ 1.0 (* 2.0 k)))))

(defn- partial-wave-exp2sigma-Sn-minus-one
  "Thompson & Nunes **(3.1.88)** nuclear factor **e^{2iσ_L}(S_L^n − 1)**.

  **`s-matrix`** is the same R-matrix / outgoing–ingoing Hankel ratio as **`s-matrix0`** (neutral), with **Hankel±** for **η ≠ 0** and **e^{2iσ_L}** used **only** here (not to define **`s-matrix`**).

  Pass **S_L^n = `s-matrix`** and **e^{2iσ_L}** = **`c/complex-polar (* 2 σ_L) 1`**, **σ_L = coulomb-sigma-L L η**."
  [S-L-n e2is]
  (c/mul e2is (c/subt S-L-n 1.0)))

(defn coulomb-scattering-amplitude-thompson-nunes-eq-3181
  "Point-charge Coulomb amplitude **f_C(θ)** (**fm**), Ian J. Thompson & Filomena M. Nunes,
  *Nuclear Reactions for Astrophysics* (Cambridge), **Eq. (3.181)** — non-relativistic CM, bare **η**.
  **η** = **`channel-sommerfeld-eta`** of **E_cm**; **k** = **`√(mass-factor·E_cm)`** (**fm⁻¹**);
  **σ_0** = **`coulomb-sigma-L`** **0** **η** = **arg Γ(1+iη)**.

  Implemented as **f_C = −(η/(2k sin²(θ/2))) exp(i(−η ln sin²(θ/2) + 2σ_0))** (equivalently **|η|/(2k sin²)** with the **π** phase fold for the leading minus)."
  [^double theta-rad eta-val ^double k]
  (let [eta (double eta-val)
        sigma-0 (coulomb-sigma-L 0 eta)
        s2 (Math/sin (* 0.5 theta-rad))
        sinsq (max (* s2 s2) 1e-300)
        mag (/ eta (* 2.0 k sinsq))
        phase (+ (* -1.0 eta (Math/log sinsq)) (* 2.0 sigma-0))]
    (c/complex-polar (+ phase Math/PI) mag)))

(def ^:dynamic *elastic-imag-ws-params*
  "When bound to [W0 R_W a_W] (W0 > 0), elastic s-matrix uses V(r)=V_C+V_real WS - i W0 f_W(r).
   Nil = real potential only (default)."
  nil)

(def ^:dynamic *elastic-match-radius-fm*
  "When set to positive **a** (fm), elastic **`s-matrix`** matches Coulomb **Hankel±** at **r = a**
   instead of **2(R_C + a₀)** from **V**. Default nil → **2(R_C + a₀)**."
  nil)

(def ^:dynamic *r-matrix-numerov-dr-fm*
  "Positive radial step **dr** (fm) for **`r-matrix-complex-imag-ws`** and Coulomb **`r-matrix`**
   in **`s-matrix`**. Default nil → **0.001** fm."
  nil)

(defn r-matrix-complex-imag-ws
  "Complex **R = u/(a u′)** at radius **a** with **V = V_C + V_WS,re − i W₀ f_W(r)**.

  Older explicit stepping on **u″ = f(r) u** spoiled **u** at large **L** for typical **dr** (**NaN `s-matrix`**
  beyond **L ≈ 45**), so elastic codes silently dropped high partial waves and **dσ/dΩ** plots looked **jagged**.
  This uses complex **Numerov** (**`numerov-step-complex`**) on a grid whose last **u** step lands at **a**,
  with **`numerov-append-and-renormalize`** so **|**u**|** does not under/overflow along **r** (**R** invariant).

  Optional **`dr-ref`** (fm): step size; default **0.001**. Use **0.1** with large **a** (e.g. **300** fm) for experiments."
  ([E V a L w-vec] (r-matrix-complex-imag-ws E V a L w-vec 0.001))
  ([E V a L w-vec dr-ref]
  (let [E (double E)
        a (double a)
        L (long L)
        [W0 RW aW] w-vec
        W0 (double W0)
        RW (double RW)
        aW (double aW)
        R0 (second V)
        dr-ref (double dr-ref)
        _ (when (not (pos? dr-ref))
            (throw (ex-info "r-matrix-complex-imag-ws: dr-ref must be positive" {:dr-ref dr-ref})))
        n-intervals (max 2 (long (Math/ceil (/ a dr-ref))))
        h (/ a (double n-intervals))
        imag-at (fn [^double r]
                  (* -1.0 W0 (/ 1.0 (+ 1.0 (Math/exp (/ (- r RW) aW))))))
        ;; f_0 … f_{N+1} with **N = n-intervals** so Numerov can reach u_N at r = a.
        rs (vec (take (+ n-intervals 2) (iterate #(+ ^double % h) 0.0)))
        fs (mapv (fn [^double r]
                   (if (< r 1e-14)
                     (c/complex-cartesian 0.0 0.0)
                     (let [v-re (+ (Coulomb-pot r R0) (WS r V))
                           v-im (imag-at r)
                           pot (c/complex-cartesian v-re v-im)
                           cent (/ (* L (inc L)) (* r r))
                           v-minus-e (c/subt pot (c/complex-cartesian E 0.0))]
                       (c/add (c/complex-cartesian cent 0.0) (c/mul mass-factor v-minus-e)))))
                 rs)
        u0 (c/complex-cartesian 0.0 0.0)
        u1 (c/complex-cartesian (Math/pow h (inc L)) 0.0)
        magmx0 (Math/max (double (c/mag u0)) (double (c/mag u1)))
        ;; Build **u_0 … u_N** on **r = 0, h, …, N h = a** (**N = n-intervals**).
        results (loop [res [u0 u1]
                       mag-mx magmx0]
                  (if (= (count res) (inc n-intervals))
                    res
                    (let [n (dec (count res))
                          un (peek res)
                          un-1 (get res (- n 1))
                          fn-1 (get fs (dec n))
                          fn (get fs n)
                          fn+1 (get fs (inc n))
                          unp1 (numerov-step-complex-guarded un un-1 fn-1 fn fn+1 h)
                          [res' mag-mx'] (numerov-append-and-renormalize res unp1 mag-mx)]
                      (recur res' mag-mx'))))
        u-at-a (peek results)
        u-prev (get results (- (count results) 2))
        u-2 (get results (- (count results) 3))
        ;; First-order forward difference; if **u′** is tiny (node at **a**), use 3-point backward stencil.
        dudr-fwd (c/div (c/subt2 u-at-a u-prev) (c/complex-cartesian h 0.0))
        dudr (if (< ^double (c/mag dudr-fwd) 1e-12)
               (c/div (c/add (c/subt2 (c/mul 3.0 u-at-a) (c/mul 4.0 u-prev)) u-2)
                      (c/complex-cartesian (* 2.0 h) 0.0))
               dudr-fwd)
        denom (c/mul (c/complex-cartesian a 0.0) dudr)
        denom' (if (< ^double (c/mag denom) 1e-22)
                   ;; Node of **u** near **a** → tiny **u′**; nudge denominator so **R** stays finite.
                   (c/mul (c/complex-cartesian a 0.0)
                          (c/add dudr (c/complex-cartesian 1e-22 0.0)))
                   denom)]
    (c/div u-at-a denom'))))

(defn r-matrix ([^double E V ^long L ]  ;with coulomb ;construct R-matrix * a depending on 1D Coulomb + Woods-Saxon potential V(R) = -V0/(1+exp ((r-R0)/a0)) V = [V0, R0, a0]
                                        ;choose u(0) = 0, u'(0) = 1 for initial conditions
                (let [
                      R0 (second V)
                      a0 (last V)
                      a (* 3 (+ R0 a0)) ; a outside of the nuclear range
                      dr 0.001]
                  (loop [x dr pot 0 d2udr2 (/ 1. dr) dudr 1 ur dr]
                    (if (> x a)
                      (/ ur dudr a); R = ur/ (a dudr) 

                      (recur  (+ x dr) (+ (Coulomb-pot x R0) (WS x V)) (*  (+ (/ (* L (inc L)) (* x x)) (* mass-factor (-  pot E))) ur) (+ dudr (* d2udr2 dr))  (+ ur (* dudr dr))) ) 


                    )))
([^double E V a ^long L ] (r-matrix E V a L 0.001))
([E V a L dr]  ;construct R-matrix * a depending on 1D Coulomb + Woods-Saxon potential V(R) = -V0/(1+exp ((r-R0)/a0)) V = [V0, R0, a0]
                                        ;choose u(0) = 0, u'(0) = 1 for initial conditions
                (let [
                      R0 (second V)
                      a0 (last V)
                      dr (double dr)]
                  (when (not (pos? dr))
                    (throw (ex-info "r-matrix [E V a L dr]: dr must be positive" {:dr dr})))
                  (loop [x dr pot 0 d2udr2 (/ 1. dr) dudr 1 ur dr]
                    (if (> x a)
                      (/ ur dudr a); R = ur/ (a dudr) 

                      (recur  (+ x dr) (+ (Coulomb-pot x R0) (WS x V)) (*  (+ (/ (* L (inc L)) (* x x)) (* mass-factor (-  pot E))) ur) (+ dudr (* d2udr2 dr))  (+ ur (* dudr dr))) ) 


                    )))
)

(defn- hankel-rho-deriv-step
  "Finite-difference step dρ for ∂H/∂ρ at matching radius ρ. A fixed 1e-7 fails for large L (Hankel
  varies too fast / loses precision). Step scales ~0.01–2% of ρ for heavy-ion matching."
  ^double [^double rho]
  (double (max 1.0e-8 (min 2.0e-2 (* rho 2.0e-4)))))

(defn- complex-finite?
  "True if z has finite real and imaginary parts (avoids NaN-poisoned partial-wave sums)."
  [z]
  (try
    (let [r (double (c/re z))
          im (double (c/im z))]
      (and (Double/isFinite r) (Double/isFinite im)))
    (catch Exception _ false)))

(defn- s-matrix-3-impl [E V L match-a radial-dr]
  (let [V-vec (vec V)
        V0 (double (first V-vec))
        Rc (double (second V-vec))
        k (m/sqrt (* mass-factor E))
        wv *elastic-imag-ws-params*
        imag-on? (and (sequential? wv) (pos? (double (first wv))))
        a (double match-a)
        radial-dr (double radial-dr)
        ;; Same matching ratio as **`s-matrix0`**: **(Ō⁻ − ρ R ∂Ō⁻) / (Ō⁺ − ρ R ∂Ō⁺)** with **Ō = Hankel±** and **η = Z₁Z₂e²·mass/(2k)**.
        ;; **No σ** in **`s-matrix`** — **σ** enters only **`partial-wave-exp2sigma-Sn-minus-one`** for **f_N**.
        ;;
        ;; **`pure-point-coulomb?`** is a **misleading** name: it does **not** mean “point Coulomb charge, but
        ;; nuclear physics can still be on elsewhere.” In this code the **entrance** optical channel is **only**
        ;; **`V = [V₀ R_C a]`** (real WS + **`Coulomb-pot r R_C`**) plus optional imag WS. So this flag really means:
        ;; **V₀ = 0** (no real Woods–Saxon well), **R_C = 0** (point **Z₁Z₂e²/r** in **`Coulomb-pot`**), and **no**
        ;; imaginary WS ⇒ there is **no** short-range optical potential to generate **S_L^n ≠ 1**; we short-circuit
        ;; to **S^n ≡ 1** (pure Coulomb pole, no nuclear/optical scattering in **V**).
        pure-point-coulomb? (and (zero? V0) (zero? Rc) (not imag-on?))]
    (if pure-point-coulomb?
      (c/complex-cartesian 1.0 0.0)
      (let [R (if imag-on?
                (r-matrix-complex-imag-ws E V-vec a L wv radial-dr)
                (c/complex-cartesian (r-matrix E V-vec a L radial-dr) 0.0))
            eta (* Z1Z2ee (/ mass-factor k 2))
            rho (* k a)
            rho-c (c/complex-cartesian rho 0.0)
            d-rho (hankel-rho-deriv-step rho)]
        (c/div (c/subt2 (Hankel- L eta rho) (c/mul rho-c R (deriv Hankel- L eta rho d-rho)))
               (c/subt2 (Hankel+ L eta rho) (c/mul rho-c R (deriv Hankel+ L eta rho d-rho))))))))

;; Cache must include mass-factor, Z1Z2ee, optional imaginary WS, match radius **a**, and **dr**: all affect **S**.
(def ^:private s-matrix-3-memo
  (memoize (fn [[E V L mf z12 wkey ma dr]]
             (binding [mass-factor (double mf)
                       Z1Z2ee (double z12)
                       *elastic-imag-ws-params* wkey]
               (s-matrix-3-impl E V L (double ma) (double dr))))))

(defn s-matrix
  "**S_L^n** — R-matrix match to Coulomb **Hankel±**, **same quotient as `s-matrix0`** with **hankel0± → Hankel±** and **η = Z₁Z₂e²·mass/(2k)**.

  **σ** does **not** appear in **`s-matrix`**. Coulomb **σ** is used only in **`partial-wave-exp2sigma-Sn-minus-one`** (**e^{2iσ}(S^n−1)**) and **f_C**.

  **Pure point Coulomb** (**V₀ = R_C = 0**, no imag WS): **1 + 0i**.

  Two arities: memoized **[E V L]** (optional **`*elastic-match-radius-fm*`**, **`*r-matrix-numerov-dr-fm*`**) vs explicit **a** in **[E V a L]**."
  ([E V L]
                (let [E (double E)
                      L (long L)
                      V-vec (vec V)
                      match-a (if *elastic-match-radius-fm*
                                (double *elastic-match-radius-fm*)
                                (* 2 (+ (double (second V-vec)) (double (last V-vec)))))
                      radial-dr (let [d *r-matrix-numerov-dr-fm*]
                                  (if (and d (pos? (double d)))
                                    (double d)
                                    0.001))]
                  (when (not (pos? match-a))
                    (throw (ex-info "s-matrix: *elastic-match-radius-fm* must be positive when set"
                                    {:match-a match-a})))
                  (s-matrix-3-memo [E V-vec L mass-factor Z1Z2ee *elastic-imag-ws-params* match-a radial-dr])))

  ([E V a L]
                (let [E (double E)
                      a (double a)
                      L (long L)
                      k (m/sqrt (*  mass-factor E))
                      R (/ (r-matrix-a E V a L) a)
                      eta (* Z1Z2ee mass-factor (/ 1. k 2))
                      rho (* k a)
                      d-rho (hankel-rho-deriv-step rho)]
                  (c/div (c/subt2 (Hankel- L eta rho) (c/mul rho R (deriv Hankel- L eta rho d-rho)))
                         (c/subt2 (Hankel+ L eta rho) (c/mul rho R (deriv Hankel+ L eta rho d-rho)))))))

;; When bound to **`(fn [E-cm V-params L] ...)`** returning complex **S_L^n**, **`differential-cross-section`** /
;; **`elastic-nuclear-amplitude-fn`** / **`total-cross-section`** / **`ftheta-L`** use it instead of **`s-matrix`**.
;; Nil = default **`s-matrix`**. (Plain **IFn** only — **`with-redefs`** on **`s-matrix`** fails: primitive **`IFn`**.)
(def ^:dynamic *partial-wave-s-matrix-fn* nil)

(defn- s-matrix-for-partial-wave-sum
  [E-cm V L]
  (if *partial-wave-s-matrix-fn*
    (*partial-wave-s-matrix-fn* E-cm V L)
    (s-matrix E-cm V L)))

(defn phase-shift
  "Same as **`phase-shift0`**: **½ arg(S_L^n)** with **S_L^n = `s-matrix`** (no added **σ**)."
  ([^double E V ^long L]
   (/ (double (c/arg (s-matrix E V L))) 2.0))
  ([^double E V a ^long L]
   (/ (double (c/arg (s-matrix E V a L))) 2.0)))


;R-matrix method
 (defn epsn0 [n V a] ;eigenenergies of trial R-matrix eigenfunctions with beta = 0, and L=1
   (+ (Math/pow (* Math/PI (/ (+ n 0.5) (Math/sqrt mass-factor) a) ) 2) V))

(defn rm-omega0 [n a r] (* (m/sqrt (/ 2 a)) (Math/sin (* Math/PI (+ n 0.5) (/ r a) ))))

(defn rm0-N [E V a N] ; R in the R-matrix method, for L= 0
  (reduce +  (map (fn [n] (/ (rm-omega0 n a a) (epsn0 n V a) mass-factor a)) (range 1 N)))
)
;end R-matrix method

(defn ftheta-L [^double E V ^long L theta]
  (let [k (m/sqrt (* mass-factor E))
        eta (channel-sommerfeld-eta E)
        sig (coulomb-sigma-L L eta)
        e2is (c/complex-polar (* 2.0 sig) 1.0)
        S-L-n (s-matrix-for-partial-wave-sum E V L)
        bracket (partial-wave-exp2sigma-Sn-minus-one S-L-n e2is)]
    (c/mul (c/div (c/complex-cartesian 0 -1) (* 2.0 k))
           (inc (* 2 L))
           (poly/eval-legendre-P L (m/cos theta))
           bracket)))

(defn hypergeometric-complex-U2
[a b z]
  (c/mul (c/div Math/PI (Math/sin (* Math/PI b)))  (c/subt (c/div (hypergeometric-complex-1F1 a b z) (c/mul
                                                                                  (c/gamma-complex (c/subt (c/add a 1) b)) (c/gamma-complex b))) 
                                         (c/mul (c/cpowc z (dec b)) (c/div (hypergeometric-complex-1F1  (c/subt (c/add a 1) b)  (- b 2) z) (c/mul (c/gamma-complex  a)  (c/gamma-complex (c/subt b 2))) ))))

;  (c/div   (hypergeometric-complex-1F1 a (- b 2) z)      (c/mul (gamma-complex (c/subt b 2)) (gamma-complex a))) 
  )

;; DWBA Differential Cross-Section Functions
;;
;; **Units:** Kinematics give scattering amplitude **f** in **fm** (k in fm⁻¹), so **|f|²** is **fm²/sr**.
;; The **×10** on **dσ** below is the usual nuclear cross-section convention **1 fm² = 10 mb**, i.e.
;; **dσ(mb/sr) = 10 × dσ(fm²/sr)** — not an extra physical factor.
(defn differential-cross-section
  "Differential cross-section |f|² (**mb/sr**) from partial-wave sum (|f|² in **fm²/sr**, then **×10**).
  Coulomb elastic: each wave uses **−i/(2k) · (2L+1) P_L · e^{2iσ_L}(S_L^n−1)** (T&N **(3.1.88)**), **S_L^n = `s-matrix`**.

  Fourth arg is either:
  - a long/int L-max: sum L = 0 … L-max (legacy);
  - a collection of non-negative longs: sum only those L (e.g. main-page angular momenta).
  Optional **`binding [*`partial-wave-s-matrix-fn*` (fn [E V L] ...)]`** supplies **S_L^n** (e.g. Numerov **R**
  → **S**) instead of **`s-matrix`**.

  Partial waves whose **`s-matrix`** / override or amplitude is non-finite (NaN/Inf from Coulomb U/Hankel at
  extreme L, η, ρ) are omitted so the sum stays finite; extend L only when numerics are stable."
  [E-cm ws-params theta-cm L-spec]
  (let [k (m/sqrt (* mass-factor E-cm))
        eta (channel-sommerfeld-eta E-cm)
        L-seq (if (number? L-spec)
                (range 0 (inc (long L-spec)))
                (sort (distinct (filter (fn [x] (>= (long x) 0)) (map long L-spec)))))
        amp-for-L
        (fn [^long L]
          (let [sig (coulomb-sigma-L L eta)
                e2is (c/complex-polar (* 2.0 sig) 1.0)
                S-L-n (s-matrix-for-partial-wave-sum E-cm ws-params L)
                bracket (partial-wave-exp2sigma-Sn-minus-one S-L-n e2is)
                f-L (c/mul (c/div (c/complex-cartesian 0 -1) (* 2.0 k))
                           (inc (* 2 L))
                           (poly/eval-legendre-P L (m/cos theta-cm))
                           bracket)]
            (when (and (complex-finite? S-L-n) (complex-finite? f-L))
              f-L)))
        amps (keep amp-for-L L-seq)
        total-amplitude
        (if (empty? amps)
          (c/complex-cartesian 0.0 0.0)
          (reduce (fn [acc a]
                    (if (complex-finite? a)
                      (c/add acc a)
                      acc))
                  (c/complex-cartesian 0.0 0.0)
                  amps))]
    (let [tr (double (c/re total-amplitude))
          ti (double (c/im total-amplitude))
          ;; |f|² is fm²/sr; ×10 ⇒ mb/sr (1 fm² = 10 mb).
          dsig (* 10.0 (+ (* tr tr) (* ti ti)))]
      (c/complex-cartesian (if (Double/isFinite dsig) dsig 0.0) 0.0))))

(defn elastic-nuclear-amplitude-fn
  "Elastic **f_N(θ)** (**fm**) — Thompson–Nunes **(3.1.88)** nuclear part, sum **L = 0 … L_cut**.
  Each term is **−i/(2k) · (2L+1) P_L(cos θ) · e^{2iσ_L}(S_L^n − 1)** with **S_L^n = `s-matrix`** (**(3.1.84)**).
  When **S_L^n = 1**, **f_N = 0**.

  **Numerics:** for charged projectiles on light targets, **S^n** from Coulomb–Hankel matching can grow the sum
  erratically past moderate **L_cut**; treat **L_cut** as a convergence knob (dashboard default 22) until high‑**L**
  matching is hardened.

  Analytic **f_C** is **not** included — add **`coulomb-scattering-amplitude-thompson-nunes-eq-3181`** for **f_C + f_N**."
  [E-cm ws-params theta-cm L-cut]
  (let [L-upper (long L-cut)
        k (m/sqrt (* mass-factor E-cm))
        eta (channel-sommerfeld-eta E-cm)
        L-seq (range 0 (inc L-upper))
        amp-for-L
        (fn [^long L]
          (let [sig (coulomb-sigma-L L eta)
                e2is (c/complex-polar (* 2.0 sig) 1.0)
                S-L-n (s-matrix-for-partial-wave-sum E-cm ws-params L)
                bracket (partial-wave-exp2sigma-Sn-minus-one S-L-n e2is)
                f-L (c/mul (c/div (c/complex-cartesian 0 -1) (* 2.0 k))
                           (inc (* 2 L))
                           (poly/eval-legendre-P L (m/cos theta-cm))
                           bracket)]
            (when (and (complex-finite? S-L-n) (complex-finite? f-L))
              f-L)))
        amps (keep amp-for-L L-seq)]
    (if (empty? amps)
      (c/complex-cartesian 0.0 0.0)
      (reduce (fn [acc a]
                (if (complex-finite? a)
                  (c/add acc a)
                  acc))
              (c/complex-cartesian 0.0 0.0)
              amps))))

(defn elastic-nuclear-amplitude-tilde-fn
  "Elastic **f̃_N(θ)** (**fm**): partial-wave sum **L = 0 … L_cut** with the **σ_0-stripped** Coulomb phase.

  Each term is **−i/(2k) · (2L+1) P_L(cos θ) · e^{2i(σ_L−σ_0)} · (S_L^n − 1)** where **e^{2i(σ_L−σ_0)}** =
  **`coulomb-phase-diff`** (product form, no Gamma) and **S_L^n** = **`s-matrix-for-partial-wave-sum`**.

  Same geometry as **`elastic-nuclear-amplitude-fn`**, but **(S^n−1)** is multiplied by **e^{2i(σ_L−σ_0)}** instead of
  using **e^{2iσ_L}(S^n−1)** in one bracket.  With **f̃_C** = **`coulomb-amplitude-tilde`**, **|f̃_C + f̃_N| = |f_C + f_N|**.

  Optional **`binding [*`partial-wave-s-matrix-fn*` (fn [E V L] …)]`** supplies **S_L^n** (e.g. Numerov **R** → **S**
  in **`example_16Odp`**); default is **`s-matrix`** on **ws-params**.

  **5-arg** **`[E-cm ws-params theta-cm L-cut η]`**: use this **η** (radians, dimensionless Sommerfeld) in **`coulomb-phase-diff`**
  instead of **`channel-sommerfeld-eta`** — for scripts that already fixed **η** from the entrance channel (must still
  **`binding`** **`mass-factor`** / **`Z1Z2ee`** so **k** matches that channel)."
  ([E-cm ws-params theta-cm L-cut]
   (elastic-nuclear-amplitude-tilde-fn E-cm ws-params theta-cm L-cut (channel-sommerfeld-eta E-cm)))
  ([E-cm ws-params theta-cm L-cut eta]
   (let [eta-d (double eta)
         L-upper (long L-cut)
         k (m/sqrt (* mass-factor E-cm))
         L-seq (range 0 (inc L-upper))
         amp-for-L
         (fn [^long L]
           (let [ph-prod (coulomb-phase-diff L eta-d)
                 S-L-n (s-matrix-for-partial-wave-sum E-cm ws-params L)
                 bracket (c/mul ph-prod (c/subt2 S-L-n 1.0))
                 f-L (c/mul (c/div (c/complex-cartesian 0 -1) (* 2.0 k))
                            (inc (* 2 L))
                            (poly/eval-legendre-P L (m/cos theta-cm))
                            bracket)]
             (when (and (complex-finite? S-L-n) (complex-finite? ph-prod) (complex-finite? f-L))
               f-L)))
         amps (keep amp-for-L L-seq)]
     (if (empty? amps)
       (c/complex-cartesian 0.0 0.0)
       (reduce (fn [acc a]
                 (if (complex-finite? a)
                   (c/add acc a)
                   acc))
               (c/complex-cartesian 0.0 0.0)
               amps)))))

(defn differential-cross-section-nuclear-cut
  "Elastic **dσ/dΩ** in **mb/sr**: **10 × |f_C(θ) + f_N(θ)|²** with **|·|²** in **fm²/sr** (**×10** = **1 fm² = 10 mb**).

  - **f_C** — **`coulomb-scattering-amplitude-thompson-nunes-eq-3181`** (Thompson & Nunes **Eq. (3.181)**).
  - **f_N** — **`elastic-nuclear-amplitude-fn`** with the same **L_cut** (**L = 0 … L_cut**).

  For **pure** point Coulomb (**`s-matrix` = 1**), **f_N = 0** and **|f_C|²** is Rutherford **dσ/dΩ**."
  [E-cm ws-params theta-cm L-cut]
  (let [k (m/sqrt (* mass-factor E-cm))
        eta (channel-sommerfeld-eta E-cm)
        f-c (coulomb-scattering-amplitude-thompson-nunes-eq-3181 theta-cm eta k)
        f-n (elastic-nuclear-amplitude-fn E-cm ws-params theta-cm L-cut)
        total-amplitude (c/add f-c f-n)]
    (let [tr (double (c/re total-amplitude))
          ti (double (c/im total-amplitude))
          ;; |f_C+f_N|² is fm²/sr; ×10 ⇒ mb/sr (1 fm² = 10 mb), same as **`differential-cross-section`**.
          dsig (* 10.0 (+ (* tr tr) (* ti ti)))]
      (c/complex-cartesian (if (Double/isFinite dsig) dsig 0.0) 0.0))))

(defn total-cross-section [E-cm ws-params L-max]
  "Total cross-section σ (integrated); returns **mb** (partial-wave sum × 10, 1 fm² = 10 mb).
  Each **L** uses **|1 − S_L^n|²** with **S_L^n = `s-matrix`**."
  (let [;; Sum over partial waves (kinematic area in fm² before ×10)
        total-sigma
        (reduce +
          (for [L (range 0 (inc L-max))]
            (let [S-L-n (s-matrix-for-partial-wave-sum E-cm ws-params L)
                  ;; Cross-section contribution for this L
                  sigma-L (* (/ 2 E-cm) Math/PI (inc (* 2 L))
                             (Math/pow (c/mag (c/subt 1.0 S-L-n)) 2))]
              sigma-L)))]
    (* 10.0 total-sigma)))
