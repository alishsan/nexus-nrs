(ns dwba.finite-well
  "Functions for solving bound states in a finite square well potential.
   
   This namespace provides functions for:
   - Spherical Bessel functions (j-l, k-l) and their derivatives
   - Matching condition calculations
   - Finding bound state energies using Newton-Raphson method
   - Finding all bound states for a given well depth
   - Selecting the **n**-th partial-wave state via **`find-bound-state-finite-well-nlz`**
   - First-order **Thomas spin–orbit** in the **sharp-surface** limit
     (**a_Woods–Saxon → 0** → **δ(r − R₀)**): **`delta-surface-thomas-so-shift-over-v0`** (**ΔE_so / V₀**,
     dimensionless) and **`delta-surface-thomas-so-shift-MeV`**, plus doublet splits and **nlz** helpers."
  (:require [fastmath.core :as m]))

(def ^:const lambda-pi-squared-fm2-default
  "Pion Compton **λ_π²** (**fm²**) for δ-surface SO estimate; ~(ℏ/(m_π c))² ≈ 2.0."
  2.0)

;; --- Spherical Bessel Function Approximations for low L (finite well) ---

(defn- j-l-double-fact-odd-denom
  "Denominator **(2l+1)!!** in **j_l(x) ∼ x^l / (2l+1)!!** as **x → 0**."
  [l]
  (reduce * 1 (range 1 (+ 2 (* 2 l)) 2)))

(defn j-l [l x]
  "Spherical Bessel function j_l(x) using recurrence relations.
   More efficient for low l values than general implementations."
  (let [x (double x)]
    (case l
      -1 (/ (m/cos x) x)
      (if (and (>= l 0) (< (m/abs x) 1.0e-14))
        (/ (m/pow x l) (double (j-l-double-fact-odd-denom l)))
        (case l
          0  (/ (m/sin x) x)
          1  (- (/ (m/sin x) (* x x)) (/ (m/cos x) x))
          ;; General recurrence: j_{l} = ((2l-1)/x)j_{l-1} - j_{l-2}
          (let [j-m2 (j-l (- l 2) x)
                j-m1 (j-l (- l 1) x)]
            (- (* (/ (dec (* 2 l)) x) j-m1) j-m2)))))))

(defn k-l [l x]
  "Modified Spherical Bessel function k_l(x) for bound states.
   Used for exponential decay outside the well."
  (case l
    -1 (/ (m/exp (- x)) x)  ;; Related to k_1 via recurrence
    0  (/ (m/exp (- x)) x)
    1  (* (/ (m/exp (- x)) x) (+ 1 (/ 1 x)))
    ;; General recurrence: k_{l+1} = k_{l-1} + (2l+1)/x * k_l
    (let [km2 (k-l (- l 2) x)
          km1 (k-l (- l 1) x)]
      (+ km2 (* (/ (dec (* 2 l)) x) km1)))))

;; --- Finite Well Bound State Functions ---

(defn j-l-deriv [l x]
  "Derivative of spherical Bessel function j_l(x).
   Uses the recurrence relation: j_l'(x) = (l/x) * j_l(x) - j_{l+1}(x)"
  (if (zero? x)
    (case l
      0 0.0
      1 (/ 1.0 3.0)
      Double/NaN)
    (- (* (/ l x) (j-l l x)) (j-l (inc l) x))))

(defn k-l-deriv [l x]
  "Derivative of modified spherical Bessel function k_l(x).
   Uses the recurrence relation: k_l'(x) = -(l/x) * k_l(x) - k_{l+1}(x)"
  (if (zero? x)
    Double/POSITIVE_INFINITY
    (- (* (/ l x) (k-l l x)) (k-l (inc l) x))))

(defn j-ratio [l x]
  "Ratio j_{l-1}(x) / j_l(x) for spherical Bessel functions."
  (if (zero? l)
    ;; For l=0, j_{-1} doesn't exist, use derivative form
    (/ (j-l-deriv 0 x) (j-l 0 x))
    (/ (j-l (dec l) x) (j-l l x))))

(defn k-ratio [l x]
  "Ratio k_{l-1}(x) / k_l(x) for modified spherical Bessel functions."
  (if (zero? l)
    ;; For l=0, k_{-1} doesn't exist, use derivative form
    (/ (k-l-deriv 0 x) (k-l 0 x))
    (/ (k-l (dec l) x) (k-l l x))))

(defn finite-well-matching-error [xi eta l]
  "Returns the mismatch in the log-derivative matching condition at the well boundary.
   
   Parameters:
   - xi: Dimensionless wave number inside well (xi = ka, where k = sqrt(2m(E+V0))/hbar)
   - eta: Dimensionless decay parameter outside well (eta = kappa*a, where kappa = sqrt(2m|E|)/hbar)
   - l: Orbital angular momentum quantum number
   
   Returns: Error in matching condition (should be 0 for bound state).
   
   Uses the recurrence relation form of the matching condition:
   xi * j_{l-1}(xi) / j_l(xi) + eta * k_{l-1}(eta) / k_l(eta) = 0
   
   Or equivalently:
   xi * j_{l-1}(xi) / j_l(xi) = -eta * k_{l-1}(eta) / k_l(eta)
   
   This is equivalent to the log-derivative form but more numerically stable.
   The sign convention: we want left - right = 0, so error = left + right (since right has negative sign)."
  (if (zero? l)
    ;; For l=0, j_{-1} doesn't exist, so use direct derivative form
    ;; At the same physical radius r=a, we need to match:
    ;; k * j_0'(xi)/j_0(xi) = kappa * k_0'(eta)/k_0(eta)
    ;; Since xi = k*a and eta = kappa*a, this becomes:
    ;; xi * j_0'(xi)/j_0(xi) = eta * k_0'(eta)/k_0(eta)
    (let [j0-xi (j-l 0 xi)
          j0-prime-xi (j-l-deriv 0 xi)
          k0-eta (k-l 0 eta)
          k0-prime-eta (k-l-deriv 0 eta)
          log-deriv-inside (* xi (/ j0-prime-xi j0-xi))  ; xi * j_0'(xi)/j_0(xi)
          log-deriv-outside (* eta (/ k0-prime-eta k0-eta))]  ; eta * k_0'(eta)/k_0(eta)
      (- log-deriv-inside log-deriv-outside))
    ;; For l > 0, use recurrence relation form
    (let [left (* xi (j-ratio l xi))
          right (* eta (k-ratio l eta))]
      ;; Matching condition: left + right = 0 (since right should be negative of left)
      (+ left right))))

(defn solver-step [e-ratio l z0]
  "Returns a map of {f(E), f'(E)} for the energy ratio e = |E|/V0.
   
   Parameters:
   - l: Orbital angular momentum quantum number
   - z0: Dimensionless well depth parameter
   - e-ratio: |E|/V0 (dimensionless energy, 0 < e-ratio < 1)
   
   Returns: {:f f-val, :f-prime f-prime}
   - f: Function value (matching error)
   - f-prime: Derivative of function w.r.t. e-ratio
   
   Uses Newton-Raphson method with analytical derivatives.
   
   Note: Returns NaN if e-ratio is outside valid range [0, 1]."
  ;; Check if e-ratio is in valid range (0, 1)
  ;; For e-ratio >= 1: xi = z0*sqrt(1-e) becomes imaginary (NaN)
  ;; For e-ratio <= 0: eta = z0*sqrt(e) becomes imaginary or zero (NaN)
  ;; Return NaN map if outside valid range (for plotting purposes)
  (if (or (<= e-ratio 0) (>= e-ratio 1.0))
    {:f Double/NaN :f-prime Double/NaN}
    (let [eta (* z0 (m/sqrt e-ratio))
          xi (* z0 (m/sqrt (- 1 e-ratio)))
          ;; Use finite-well-matching-error for consistency
          f-val (finite-well-matching-error xi eta l)
          ;; For derivative calculation, we need the intermediate values
          ;; For l=0: j-ratio = j_0'/j_0, for l>0: j-ratio = j_{l-1}/j_l
          left-val (* xi (j-ratio l xi))
          right-val (* eta (k-ratio l eta))
          ;; Raw derivatives w.r.t xi and eta
          ;; For l=0: j-ratio = j_0'/j_0, so left-val = xi * j_0'/j_0
          ;; For l>0: j-ratio = j_{l-1}/j_l, so left-val = xi * j_{l-1}/j_l
          [d-left-dxi d-right-deta]
          (if (zero? l)
            ;; For l=0: left-val = xi * j_0'/j_0, where j_0' = -j_1
            ;; d/dxi [xi * j_0'/j_0] = j_0'/j_0 + xi * d/dxi[j_0'/j_0]
            ;; Using j_0' = -j_1, we have j_0'/j_0 = -j_1/j_0
            ;; d/dxi[-j_1/j_0] = [-(j_1'*j_0 - j_1*j_0')] / j_0^2
            ;; Using j_1' = j_0/x - j_1 and j_0' = -j_1:
            ;; = [-(j_0/x - j_1)*j_0 - j_1*(-j_1)] / j_0^2
            ;; = [-(j_0^2/x - j_1*j_0) + j_1^2] / j_0^2
            ;; = -j_0/(x*j_0) + j_1/j_0 + j_1^2/j_0^2 = -1/x + j_1/j_0 + (j_1/j_0)^2
            (let [j0-xi (j-l 0 xi)
                  j0-prime-xi (j-l-deriv 0 xi)  ; = -j_1
                  j1-xi (j-l 1 xi)
                  j-ratio-val (/ j0-prime-xi j0-xi)  ; = -j_1/j_0
                  j1-over-j0 (/ j1-xi j0-xi)
                  d-j-ratio-dxi (+ (- (/ 1.0 xi)) j1-over-j0 (* j1-over-j0 j1-over-j0))
                  d-left-dxi-val (+ j-ratio-val (* xi d-j-ratio-dxi))
                  
                  ;; For k_0: k_0' = -k_1 (from k-l-deriv)
                  k0-eta (k-l 0 eta)
                  k0-prime-eta (k-l-deriv 0 eta)  ; = -k_1
                  k1-eta (k-l 1 eta)
                  k-ratio-val (/ k0-prime-eta k0-eta)  ; = -k_1/k_0
                  ;; d/deta[k_0'/k_0] = d/deta[-k_1/k_0]
                  ;; k_1' = -k_0/eta - k_1 (from k-l-deriv for l=1)
                  ;; d/deta[-k_1/k_0] = [-(k_1'*k_0 - k_1*k_0')] / k_0^2
                  ;; = [-(-k_0/eta - k_1)*k_0 - k_1*(-k_1)] / k_0^2
                  ;; = [(k_0^2/eta + k_1*k_0) + k_1^2] / k_0^2
                  ;; = k_0/(eta*k_0) + k_1/k_0 + (k_1/k_0)^2 = 1/eta + k_1/k_0 + (k_1/k_0)^2
                  k1-over-k0 (/ k1-eta k0-eta)
                  d-k-ratio-deta (+ (/ 1.0 eta) k1-over-k0 (* k1-over-k0 k1-over-k0))
                  ;; For f = xi * j_0'/j_0 - eta * k_0'/k_0
                  ;; d/deta[f] = -d/deta[eta * k_0'/k_0]
                  ;; d/deta[eta * k_0'/k_0] = k_0'/k_0 + eta * d/deta[k_0'/k_0]
                  ;; So d/deta[f] = -[k_0'/k_0 + eta * d/deta[k_0'/k_0]]
                  ;; For f = xi * j_0'/j_0 - eta * k_0'/k_0
                  ;; d/deta[f] = -d/deta[eta * k_0'/k_0]
                  ;; d/deta[eta * k_0'/k_0] = k_0'/k_0 + eta * d/deta[k_0'/k_0]
                  ;; Since f-prime > 0 everywhere, d/deta[f] > 0, so d/deta[eta*k_0'/k_0] < 0
                  ;; We store d-right-deta as d/deta[eta*k_0'/k_0] (which should be negative)
                  ;; Note: k-ratio-val = k_0'/k_0 is typically negative, making d-right-deta negative
                  d-right-deta-val (+ k-ratio-val (* eta d-k-ratio-deta))]
              [d-left-dxi-val d-right-deta-val])
            ;; For l>0: use recurrence relation
            (let [d-left-dxi-val (+ 1 (* (- (/ (dec (* 2 l)) xi)) left-val) (* left-val left-val))
                  d-right-deta-val (+ 1 (* (/ (dec (* 2 l)) eta) right-val) (* right-val right-val))]
              [d-left-dxi-val d-right-deta-val]))
          ;; Chain rule: d/de-ratio = (d/dxi * dxi/de) + (d/deta * deta/de)
          ;; dxi/de = -z0 / (2 * sqrt(1-e)) [negative: as e increases, xi decreases]
          ;; deta/de = z0 / (2 * sqrt(e)) [positive: as e increases, eta increases]
          dxi-de (/ (- z0) (* 2.0 (m/sqrt (- 1.0 e-ratio))))
          deta-de (/ z0 (* 2.0 (m/sqrt e-ratio)))
          ;; For l=0: f = xi*j_0'/j_0 - eta*k_0'/k_0
          ;; d/dxi[f] = d/dxi[xi*j_0'/j_0] = d-left-dxi
          ;; d/deta[f] = -d/deta[eta*k_0'/k_0] = d-right-deta (which is negative)
          ;; f-prime = d-left-dxi * dxi/de + d-right-deta * deta/de
          ;; Since dxi/de is negative and deta/de is positive:
          ;; f-prime = d-left-dxi * (negative) + (negative) * (positive)
          ;; If f is increasing, we need the second term to dominate positively
          ;; This suggests d-right-deta should actually be positive, not negative!
          ;; Let me check: if f = left - right, then d/deta[f] = -d/deta[right]
          ;; But d-right-deta is calculated as d/deta[right], so we need to negate it
          f-prime (if (zero? l)
                    (- (+ (* d-left-dxi dxi-de) (* (- d-right-deta) deta-de)))  ; Negate to fix sign - f is increasing so f-prime should be positive
                    (+ (* d-left-dxi dxi-de) (* d-right-deta deta-de)))]  ; d-right-deta is correct for l>0
      {:f f-val :f-prime f-prime})))

(defn- j-l-zero-in-xi-bracket
  "Bisection locate **ξ** in **(ξ-lo, ξ-hi)** with **j_l(ξ)=0** (**l ≥ 1**)."
  [l ^double xi-lo ^double xi-hi]
  (when (and (pos? l) (< xi-lo xi-hi))
    (let [fl (j-l l xi-lo)
          fh (j-l l xi-hi)]
      (when (and (not (Double/isNaN fl)) (not (Double/isNaN fh))
                 (< (* fl fh) 0.0))
        (loop [a xi-lo fa fl b xi-hi fb fh n 0]
          (if (>= n 72)
            (* 0.5 (+ a b))
            (let [m (* 0.5 (+ a b))
                  fm (j-l l m)]
              (if (< (m/abs fm) 1e-14)
                m
                (if (<= (* fa fm) 0.0)
                  (recur a fa m fm (inc n))
                  (recur m fm b fb (inc n)))))))))))

(defn- e-ratios-at-j-l-xi-poles
  "**e = |E|/V₀** where **ξ = z₀√(1−e)** hits a zero of **j_l** (denominator of **j_{l-1}/j_l**)."
  [l z0]
  (when (pos? l)
    (let [z (double z0)
          step 0.04]
      (loop [xa 0.02
             acc []]
        (if (>= xa z)
          (->> acc (filter some?) distinct sort vec)
          (let [xb (min z (+ xa step))]
            (recur xb
                   (if-let [xi0 (j-l-zero-in-xi-bracket l xa xb)]
                     (let [e (- 1.0 (/ (* xi0 xi0) (* z z)))]
                       (if (and (> e 0.01) (< e 0.99))
                         (conj acc e)
                         acc))
                     acc))))))))

(defn find-discontinuities [l z0]
  "Finds analytical discontinuities in the matching error function.
   
   For l=0: j_0'/j_0 = cot(xi) - 1/xi, which is discontinuous when sin(xi) = 0,
   i.e., when xi = n*π for n = 1, 2, 3, ...
   
   Since xi = z0 * sqrt(1 - e_ratio), we have:
   z0 * sqrt(1 - e_ratio) = n*π
   1 - e_ratio = (n*π/z0)^2
   e_ratio = 1 - (n*π/z0)^2
   
   For **l > 0**, includes **e** where **j_l(ξ)=0** (poles of **j_{l-1}/j_l**).
   
   Returns: Vector of e-ratio values where discontinuities occur, sorted."
  (if (zero? l)
    (let [max-n (int (/ z0 Math/PI))
          discontinuities (for [n (range 1 (inc max-n))
                                :let [n-pi-over-z0 (/ (* n Math/PI) z0)
                                      e-ratio (- 1.0 (* n-pi-over-z0 n-pi-over-z0))]]
                            (when (and (> e-ratio 0.01) (< e-ratio 0.99))
                              e-ratio))]
      (->> discontinuities
           (filter some?)
           sort
           vec))
    (e-ratios-at-j-l-xi-poles l z0)))

(defn scan-energy-range-dimensionless [l z0 step-size]
  "Scans the energy range and returns function values.
   
   Returns: Vector of {:e-ratio, :f} for each valid point."
  (let [test-points (range 0.01 0.99 step-size)]
    (->> test-points
         (map (fn [e-ratio]
                (let [{:keys [f]} (solver-step e-ratio l z0)]
                  (when (not (Double/isNaN f))
                    {:e-ratio e-ratio :f f}))))
         (filter some?)
         (vec))))

(defn find-numerical-discontinuities [results]
  "Finds numerical discontinuities (large jumps) in the function values.
   
   Returns: Vector of e-ratio values where discontinuities occur."
  (->> (partition 2 1 results)
       (keep (fn [[curr next]]
               (when (and curr next)
                 (let [jump-size (m/abs (- (:f next) (:f curr)))]
                   (when (> jump-size 1000.0)
                     (/ (+ (:e-ratio curr) (:e-ratio next)) 2.0))))))
       (vec)))

(defn find-sign-change-candidates [results all-disc-e-ratios]
  "Finds sign change candidates, excluding those near discontinuities.
   
   Returns: Vector of {:e1, :e2, :f1, :f2} for each sign change."
  (->> (partition 2 1 results)
       (keep (fn [[curr next]]
               (when (and curr next)
                 (let [f1 (:f curr)
                       f2 (:f next)
                       s1 (m/signum f1)
                       s2 (m/signum f2)
                       e1 (:e-ratio curr)
                       e2 (:e-ratio next)
                       e-mid (/ (+ e1 e2) 2.0)
                       near-disc? (some (fn [disc-e]
                                          (< (m/abs (- e-mid disc-e)) 0.02))
                                        all-disc-e-ratios)
                       is-sign-change? (and (not (zero? s1))
                                            (not (zero? s2))
                                            (not= s1 s2)
                                            (not near-disc?))]
                   (when is-sign-change?
                     {:e1 e1 :e2 e2 :f1 f1 :f2 f2})))))
       (vec)))

(defn find-local-minima-candidates [results]
  "Finds local minima in |f| as potential roots.
   
   Returns: Vector of {:e1, :e2, :f1, :f2} for each local minimum."
  (->> (partition 3 1 results)
       (keep (fn [[prev curr next]]
               (when (and prev curr next)
                 (let [f-prev (m/abs (:f prev))
                       f-curr (m/abs (:f curr))
                       f-next (m/abs (:f next))
                       is-minimum? (and (< f-curr f-prev)
                                       (< f-curr f-next)
                                       (< f-curr 10.0))]
                   (when is-minimum?
                     {:e1 (max (:e-ratio prev) 0.01)
                      :e2 (min (:e-ratio next) 0.99)
                      :f1 (:f curr)
                      :f2 (:f curr)})))))
       (vec)))

(defn- bisect-matching-error-root
  "**e = |E|/V₀** with **f(e)=0**, given **f(e-lo) f(e-hi) ≤ 0**."
  [e-lo e-hi l z0]
  (let [f-lo (:f (solver-step e-lo l z0))
        f-hi (:f (solver-step e-hi l z0))]
    (when (and (not (Double/isNaN f-lo)) (not (Double/isNaN f-hi))
               (<= (* (double f-lo) (double f-hi)) 0.0))
      (loop [a (double e-lo) fa f-lo b (double e-hi) fb f-hi n 0]
        (let [m (* 0.5 (+ a b))
              fm (:f (solver-step m l z0))
              tol 1e-14]
          (cond
            (>= n 96) [m (inc n)]
            (or (Double/isNaN fm) (Double/isInfinite fm))
            [(* 0.5 (+ a b)) (inc n)]
            (< (m/abs fm) tol) [m (inc n)]
            (<= (* fa fm) 0.0) (recur a fa m fm (inc n))
            :else (recur m fm b fb (inc n))))))))

(defn- newton-polish-matching-e
  "A few guarded Newton steps from **e0** when bisection did not apply."
  [e0 l z0 max-steps]
  (loop [e-ratio (double e0) k 0]
    (if (>= k (long max-steps))
      [e-ratio k]
      (let [{:keys [f f-prime]} (solver-step e-ratio l z0)
            ok? (and (not (Double/isNaN f)) (not (Double/isInfinite f))
                     (not (Double/isNaN f-prime)) (not (Double/isInfinite f-prime))
                     (> (m/abs f-prime) 1e-15))
            f-abs (m/abs f)]
        (if (and ok? (< f-abs 1e-12))
          [e-ratio (inc k)]
          (if-not ok?
            [e-ratio k]
            (let [step (min 0.2 (max -0.2 (/ f f-prime)))
                  e2 (max 0.01 (min 0.99 (- e-ratio step)))]
              (recur e2 (inc k)))))))))

(defn refine-root [candidate l z0]
  "Refines **f(e)=0** using **bisection** on the scan bracket when possible, else midpoint **+**
   a short Newton polish. Newton alone was fragile when **f′** and early-exit **stuck?** logic
   stopped before **|f|** was small."
  (let [{:keys [e1 e2 f1 f2]} candidate
        f1v (double (if (some? f1) f1 (:f (solver-step e1 l z0))))
        f2v (double (if (some? f2) f2 (:f (solver-step e2 l z0))))
        e-mid (* 0.5 (+ (double e1) (double e2)))]
    (if-let [[e-bis iters-b]
             (when (and (not (Double/isNaN f1v)) (not (Double/isNaN f2v))
                        (<= (* f1v f2v) 0.0))
               (bisect-matching-error-root e1 e2 l z0))]
      (let [e0 (double e-bis)
            [e-nw iters-n] (newton-polish-matching-e e0 l z0 12)
            xi (* z0 (m/sqrt (- 1.0 e-nw)))
            eta (* z0 (m/sqrt e-nw))
            ferr (finite-well-matching-error xi eta l)]
        {:e-ratio e-nw
         :xi xi
         :eta eta
         :energy e-nw
         :converged? (< (m/abs ferr) 1e-6)
         :matching-error ferr
         :iterations (+ iters-b iters-n)})
      (let [[e-nw iters-n] (newton-polish-matching-e e-mid l z0 40)
            xi (* z0 (m/sqrt (- 1.0 e-nw)))
            eta (* z0 (m/sqrt e-nw))
            ferr (finite-well-matching-error xi eta l)]
        {:e-ratio e-nw
         :xi xi
         :eta eta
         :energy e-nw
         :converged? (< (m/abs ferr) 1e-6)
         :matching-error ferr
         :iterations iters-n}))))

(defn deduplicate-states [states min-separation]
  "Removes duplicate states that are too close to each other.
   
   Parameters:
   - states: Vector of bound state maps, sorted by e-ratio
   - min-separation: Minimum e-ratio difference to consider states distinct
   
   Returns: Deduplicated vector of states."
  (loop [remaining states
         result []]
    (if (empty? remaining)
      result
      (let [current (first remaining)
            rest-states (rest remaining)
            too-close? (some (fn [prev-state]
                              (< (m/abs (- (:e-ratio current) (:e-ratio prev-state))) min-separation))
                            result)]
        (if too-close?
          (recur rest-states result)
          (recur rest-states (conj result current)))))))

(defn find-all-bound-states [l z0]
  "Finds all bound states for a given l and z0.
   
   Parameters:
   - l: Orbital angular momentum quantum number
   - z0: Dimensionless well depth parameter
   
   Returns: Vector of bound states, each with {:e-ratio, :xi, :eta, :energy, :converged?, :matching-error, :iterations}
   sorted by e-ratio (lowest energy first).
   
   Uses a fine grid scan to find sign changes, then refines brackets with **bisection** and a short Newton polish."
  (let [;; Find analytical discontinuities first
        analytical-discs (find-discontinuities l z0)
        ;; Scan energy range
        results (scan-energy-range-dimensionless l z0 0.001)
        ;; Find numerical discontinuities
        numerical-discs (find-numerical-discontinuities results)
        ;; Combine all discontinuities
        all-disc-e-ratios (sort (concat analytical-discs numerical-discs))
        ;; Find sign change candidates
        sign-changes (find-sign-change-candidates results all-disc-e-ratios)
        ;; Find local minima candidates
        local-minima (find-local-minima-candidates results)
        ;; Combine all candidates
        all-candidates (concat sign-changes local-minima)
        ;; Refine each candidate
        bound-states (mapv #(refine-root % l z0) all-candidates)
        ;; Sort by e-ratio
        sorted-states (sort-by :e-ratio bound-states)
        ;; Deduplicate
        deduplicated (deduplicate-states sorted-states 0.001)]
    deduplicated))

(defn- valid-finite-well-state?
  "Filter predicate: converged root with small matching error and physical e-ratio."
  [state tolerance]
  (and (:converged? state)
       (< (m/abs (:matching-error state)) tolerance)
       (>= (:e-ratio state) 0.01)
       (<= (:e-ratio state) 0.99)))

(defn find-bound-state-finite-well-nlz
  "Return the **n**-th bound state (1-based) for partial wave **l** and dimensionless depth **z0**.

   **n** is the spectroscopic index **for fixed l**: **n = 1** is the **most bound** state (largest
   **e-ratio** = |E|/V₀), **n = 2** is the first excited state in that partial wave, etc.

   Implementation: **`find-all-bound-states`**, then keep well-converged roots, sort by **e-ratio**
   **descending**, and take the **n**-th entry.

   Returns the same keys as **`find-bound-state-finite-well`** plus **:n**, **:l**, **:z0**.
   Returns **nil** if **n** is not a positive integer or exceeds the number of physical states.

   **z0** is the same dimensionless depth as elsewhere: **z0 = a √(2mV₀)/ℏ**."
  [n l z0]
  (when (pos-int? n)
    (let [all (find-all-bound-states l z0)
          valid (->> all
                     (filter #(valid-finite-well-state? % 1e-3))
                     (sort-by :e-ratio)
                     reverse
                     vec)]
      (when (<= n (count valid))
        (merge (nth valid (dec n))
               {:n n :l l :z0 z0})))))

(defn find-bound-state-finite-well [l z0]
  "Most bound state for partial wave **l** (**n = 1** in **`find-bound-state-finite-well-nlz`**).

   Returns **{:e-ratio, :xi, :eta, :energy, :converged?, :matching-error, :iterations}**;
   for multiple roots use **`find-all-bound-states`** or **`find-bound-state-finite-well-nlz`** with **n ≥ 2**."
  (when-let [st (find-bound-state-finite-well-nlz 1 l z0)]
    (dissoc st :n :l :z0)))

;; --- Sharp-surface (δ) Thomas spin–orbit first-order shift ---------------------------------

(defn l-dot-s-spin-half
  "Expectation **⟨l·s⟩** for spin **½**: **(j(j+1) − l(l+1) − 3/4)/2**.

   Same convention as **`dwba.transfer/l-dot-s-nucleon`** (this ns cannot require **transfer**)."
  [l j]
  (/ (- (* (double j) (+ (double j) 1.0))
        (* (double l) (+ (double l) 1.0))
        0.75)
     2.0))

(defn- square-well-u-matched-unnorm
  "Reduced radial **u(r)** (unnormalized) for spherical square well radius **a** **fm**:
   inside **r ≤ a**: **k_l(η) j_l(ξ r/a)**; outside: **j_l(ξ) k_l(η r/a)** — continuous at **r = a**."
  [l xi eta a r]
  (let [l (int l)
        xi (double xi)
        eta (double eta)
        a (double a)
        r (double r)]
    (if (<= r a)
      (* (k-l l eta) (j-l l (* (/ xi a) r)))
      (* (j-l l xi) (k-l l (* (/ eta a) r))))))

(defn- integrate-u-squared-trapezoid
  "∫₀^{r-max} u(r)² dr with uniform trapezoid rule (**steps** segments)."
  [u-fn r-max steps]
  (let [n (long steps)
        h (/ (double r-max) n)
        u2 (fn [i] (let [r (* h (long i))]
                     (m/pow (u-fn r) 2)))]
    (* h (+ (* 0.5 (+ (u2 0) (u2 n)))
            (reduce (fn [acc i] (+ acc (u2 i))) 0.0 (range 1 n))))))

(defn finite-well-u-squared-norm-integral-unnorm
  "∫₀^∞ u_unnorm(r)² dr for the matched square-well **u** (dimensionless **ξ, η**, radius **a** fm).
   Exterior is truncated at **r_max = a (1 + tail_factor/η)**."
  ([l xi eta a]
   (finite-well-u-squared-norm-integral-unnorm l xi eta a 40.0 512))
  ([l xi eta a tail-factor steps]
   (let [eta (max (double eta) 1.0e-6)
         a (double a)
         tf (double tail-factor)
         r-max (* a (+ 1.0 (/ tf eta)))
         ufn (fn [r] (square-well-u-matched-unnorm l xi eta a r))]
     (integrate-u-squared-trapezoid ufn r-max (long steps)))))

(defn finite-well-u-at-radius-normalized
  "Normalized reduced radial **u(a)** at well radius **a** (same **u** with ∫ u² dr = 1)."
  [l xi eta a]
  (let [I (finite-well-u-squared-norm-integral-unnorm l xi eta a)
        ua (square-well-u-matched-unnorm l xi eta a a)
        inv-norm (/ 1.0 (m/sqrt I))]
    (* ua inv-norm)))

(defn- thomas-so-delta-radial-prefactor-dimensionless
  "Pure number **(λ_π²/a²)·(a·u_norm(a)²) = λ_π² u_norm² / a** from matched square well (**ξ, η**, **a** **fm**)."
  [lambda-pi-sq-over-a-sq l xi eta a-fm]
  (let [un (finite-well-u-at-radius-normalized l xi eta a-fm)
        a (double a-fm)]
    (* (double lambda-pi-sq-over-a-sq) a (* un un))))

(defn delta-surface-thomas-so-shift-over-v0
  "First-order δ-surface Thomas SO shift **relative to well depth V₀**:

   **ΔE_so / V₀ = − (V_so⁽⁰⁾/V₀) · (λ_π²/a²) · (a·u_norm(a)²) · ⟨l·s⟩**,

   with **u** reduced radial (**∫ u² dr = 1**), **a = R₀** (**fm** only enters **u** normalization and
   **λ_π²/a²**). Inputs **V_so⁽⁰⁾/V₀** and **λ_π²/a²** are dimensionless; **`a-fm`** is needed to evaluate **u**.

   **ΔE_so (MeV) = V₀ · (ΔE_so/V₀)** when **V₀** is in **MeV**."
  [v-so-over-v0 lambda-pi-sq-over-a-sq a-fm l j xi eta]
  (* -1.0 (double v-so-over-v0)
     (thomas-so-delta-radial-prefactor-dimensionless lambda-pi-sq-over-a-sq l xi eta a-fm)
     (l-dot-s-spin-half l j)))

(defn delta-surface-thomas-so-shift-MeV
  "Same δ-surface formula as **`delta-surface-thomas-so-shift-over-v0`**, with **V_so⁽⁰⁾** and **λ_π²** in
   **MeV** / **fm²** — equals **V₀ × delta-surface-thomas-so-shift-over-v0** with **v = V_so⁽⁰⁾/V₀** and
   **λ_π²/a² = (λ_π² in fm²) / a²**.

   **ΔE_so = − V_so^(0) λ_π² (u_norm(R₀)² / R₀) ⟨l·s⟩**."
  [V-so-0 lambda-pi-sq-fm2 a-fm l j xi eta]
  (let [un (finite-well-u-at-radius-normalized l xi eta a-fm)
        a (double a-fm)]
    (* -1.0 (double V-so-0) (double lambda-pi-sq-fm2)
       (/ (* un un) a)
       (l-dot-s-spin-half l j))))

(defn finite-well-delta-so-splitting-j-doublet-over-v0
  "**(ΔE_so(j=l+½) − ΔE_so(j=l−½)) / V₀** from the same radial prefactor as **`delta-surface-thomas-so-shift-over-v0**.
   **nil** for **l = 0**."
  [v-so-over-v0 lambda-pi-sq-over-a-sq a-fm l xi eta]
  (when (pos? l)
    (let [pref (* -1.0 (double v-so-over-v0)
                  (thomas-so-delta-radial-prefactor-dimensionless lambda-pi-sq-over-a-sq l xi eta a-fm))
          j+ (+ (double l) 0.5)
          j- (- (double l) 0.5)]
      (* pref (- (l-dot-s-spin-half l j+) (l-dot-s-spin-half l j-))))))

(defn finite-well-delta-so-splitting-j-doublet-MeV
  "SO energy **split** **E(j=l+½) − E(j=l−½)** from the same δ-surface first-order formula
   (both built from normalized **u** at **a**). For **l = 0** returns **nil** (no doublet).

   Equals **V₀ × finite-well-delta-so-splitting-j-doublet-over-v0**."
  [V-so-0 lambda-pi-sq-fm2 a-fm l xi eta]
  (when (pos? l)
    (let [un (finite-well-u-at-radius-normalized l xi eta a-fm)
          pref (* -1.0 (double V-so-0) (double lambda-pi-sq-fm2)
                  (/ (* un un) (double a-fm)))
          j+ (+ (double l) 0.5)
          j- (- (double l) 0.5)]
      (* pref (- (l-dot-s-spin-half l j+) (l-dot-s-spin-half l j-))))))

(defn find-bound-state-finite-well-nlz-delta-so
  "Like **`find-bound-state-finite-well-nlz`** plus first-order δ-surface SO for given **j** (**l ± ½**).

   7-arg: **{:delta-e-so-MeV, :j}** merged into the state.

   8-arg: also **:e-central-MeV** = **−V₀·|E|/V₀**, **:e-total-MeV** = **:e-central-MeV + :delta-e-so-MeV**,
   and **:e-central-over-v0** / **:e-total-over-v0** (same convention as **`-delta-so-over-v0`**)."
  ([n l z0 a-fm V-so-0 lambda-pi-sq-fm2 j]
   (when-let [st (find-bound-state-finite-well-nlz n l z0)]
     (let [xi (:xi st)
           eta (:eta st)
           de (delta-surface-thomas-so-shift-MeV V-so-0 lambda-pi-sq-fm2 a-fm l j xi eta)]
       (merge st {:delta-e-so-MeV de :j j}))))
  ([n l z0 a-fm V-so-0 lambda-pi-sq-fm2 j V0-MeV]
   (when-let [m (find-bound-state-finite-well-nlz-delta-so n l z0 a-fm V-so-0 lambda-pi-sq-fm2 j)]
     (let [v0 (double V0-MeV)
           er (:e-ratio m)
           dec (:delta-e-so-MeV m)
           ecm (* (- v0) er)
           etm (+ ecm dec)
           deov (/ dec v0)]
       (merge m {:e-central-MeV ecm :e-total-MeV etm
                 :e-central-over-v0 (- er)
                 :e-total-over-v0 (+ (- er) deov)})))))

(defn find-bound-state-finite-well-nlz-delta-so-over-v0
  "Same idea as **`find-bound-state-finite-well-nlz-delta-so`**, with **:delta-e-so-over-v0**.

   **Energies (zero at threshold):** **:e-central-over-v0** = **−|E|/V₀**, **:e-total-over-v0** =
   **:e-central-over-v0 + :delta-e-so-over-v0**.  Sort by **:e-total-over-v0 ascending** = most bound first."
  [n l z0 a-fm v-so-over-v0 lambda-pi-sq-over-a-sq j]
  (when-let [st (find-bound-state-finite-well-nlz n l z0)]
    (let [xi (:xi st)
          eta (:eta st)
          er (:e-ratio st)
          de (delta-surface-thomas-so-shift-over-v0 v-so-over-v0 lambda-pi-sq-over-a-sq a-fm l j xi eta)
          ec (- er)
          et (+ ec de)]
      (merge st {:delta-e-so-over-v0 de :j j
                 :e-central-over-v0 ec
                 :e-total-over-v0 et}))))
