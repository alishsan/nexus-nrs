(ns functions.rmatrix-lagrange
  "Calculable R-matrix on a Lagrange mesh (Baye / Descouvemont–Baye, Rep. Prog. Phys. 2010).

  Single-channel scattering on **(0, a)** with a local central potential. The default
  **`:gauss`** path uses Gauss–Legendre mesh points **r_i = a x_i**, overlap **δ_ij**, local
  potential at mesh points, and analytical kinetic + Bloch **L(0)** (Appendix C, C10–C11).

  Phase shifts use the same neutral **Hankel±** matching as **`functions/phase-shift0`**."
  (:require [functions :as f]
            [complex :as c]))

;; --- Gauss–Legendre on (0, 1) ---

(defn- legendre-p
  "Legendre polynomial **P_n(x)** on **[-1, 1]**."
  [n x]
  (cond
    (zero? n) 1.0
    (== n 1) (double x)
    :else
    (loop [k 1 p0 1.0 p1 (double x)]
      (if (= k n)
        p1
        (let [pk (/ (- (* (+ (* 2 k) 1) x p1) (* k p0)) (inc k))]
          (recur (inc k) p1 pk))))))

(defn- legendre-dp
  "Derivative **P_n'(x)**."
  [n x]
  (if (zero? n)
    0.0
    (let [pn (legendre-p n x)
          pn-1 (legendre-p (dec n) x)]
      (* (/ (double n) (- (* x x) 1.0))
         (- (* x pn) pn-1)))))

(defn gauss-legendre-01
  "Return **{:nodes [x_i …], :weights [λ_i …]}** for **∫_0^1 f(x) dx ≈ Σ λ_i f(x_i)**."
  [n]
  (let [n (long n)]
    (when (< n 2)
      (throw (ex-info "gauss-legendre-01: n must be ≥ 2" {:n n})))
    (let [raw-nodes
          (vec
           (for [k (range 1 (inc n))]
             (let [guess (Math/cos (* Math/PI (/ (- (* 2.0 k) 1.0) (* 2.0 n))))]
               (loop [u (double guess) iter 0]
                 (let [pu (legendre-p n u)
                       dpu (legendre-dp n u)]
                   (if (or (>= iter 80) (< (Math/abs pu) 1e-15))
                     u
                     (recur (- u (/ pu dpu)) (inc iter))))))))
          nodes (mapv (fn [u] (* 0.5 (+ u 1.0))) raw-nodes)
          weights (mapv (fn [u]
                          (let [dpu (legendre-dp n u)]
                            (/ 2.0 (* (- 1.0 (* u u)) dpu dpu))))
                        raw-nodes)
          weights-01 (mapv #(/ % 2.0) weights)]
      {:nodes nodes :weights weights-01})))

;; --- Linear algebra ---

(defn- solve-linear [A b]
  (let [n (count A)
        Ab (reduce
            (fn [Ab col]
              (let [pivot-row (loop [r col best col best-val -1.0]
                                (if (>= r n) best
                                    (let [v (Math/abs (double (get (nth Ab r) col)))]
                                      (if (> v best-val) (recur (inc r) r v)
                                          (recur (inc r) best best-val)))))
                    pivot-val (double (get (nth Ab pivot-row) col))]
                (when (< (Math/abs pivot-val) 1e-18)
                  (throw (ex-info "solve-linear: singular matrix" {:col col})))
                (let [Ab (if (= pivot-row col) Ab
                           (let [tmp (nth Ab col)]
                             (assoc Ab col (nth Ab pivot-row) pivot-row tmp)))
                      scale (/ 1.0 pivot-val)
                      pivot (mapv #(* scale (double %)) (nth Ab col))
                      Ab (assoc Ab col pivot)]
                  (reduce
                   (fn [Ab r]
                     (if (= r col) Ab
                         (let [factor (double (get (nth Ab r) col))]
                           (assoc Ab r
                                  (mapv (fn [j]
                                          (- (double (get (nth Ab r) j))
                                             (* factor (double (get pivot j)))))
                                        (range (inc n)))))))
                   Ab (range n)))))
            (mapv (fn [i row bi]
                    (vec (concat (mapv double row) [(double bi)])))
                  (range n) A b)
            (range n))]
    (loop [x (vec (repeat n 0.0)) i (dec n)]
      (if (neg? i) x
          (let [row (nth Ab i)
                rhs (- (double (get row n))
                       (reduce + 0.0
                               (map (fn [j] (* (double (get x j)) (double (get row j))))
                                    (range (inc i) n))))]
            (recur (assoc x i (/ rhs (double (get row i)))) (dec i)))))))

(defn- mat-inv [A]
  (let [n (count A)]
    (apply mapv vector
           (mapv (fn [j] (solve-linear A (mapv #(if (= j %) 1.0 0.0) (range n))))
                 (range n)))))

;; --- Lagrange mesh (Appendix C) ---

(defn lagrange-kinetic-bloch-matrix
  "**⟨φ_i | T_0 + L(0) | φ_j⟩** in units of **ℏ²/(2μ)** (multiply by **1/mass-factor** for MeV)."
  [n a]
  (let [n (long n) a (double a) a2 (* a a)
        {:keys [nodes]} (gauss-legendre-01 n)
        nn (+ (* n n) n)]
    (mapv (fn [i]
            (mapv (fn [j]
                    (let [xi (nth nodes i) xj (nth nodes j)]
                      (if (= i j)
                        (/ (- (* (+ (* 4 n n) (* 4 n) 3) xi (- 1.0 xi))
                              (* 6.0 xi) 1.0)
                           (* 3.0 a2 xi xi (Math/pow (- 1.0 xi) 2)))
                        (let [dx (- xi xj)
                              denom (Math/sqrt (* xi xj (- 1.0 xi) (- 1.0 xj)))]
                          (* (/ (if (even? (+ i j)) 1.0 -1.0) (* a2 denom))
                             (+ (double nn)
                                (/ (- (+ xi xj) (* 2.0 xi xj)) (* dx dx))
                                (- (/ 1.0 (- 1.0 xi)))
                                (- (/ 1.0 (- 1.0 xj)))))))))
                  (range n)))
          (range n))))

(defn lagrange-phi-at-boundary
  "**φ_i(a)** for the Gauss–Legendre Lagrange basis (eq. 4.7 at **r = a**)."
  [n a nodes i]
  (let [xi (nth nodes i)
        sign (if (even? (+ n i 1)) 1.0 -1.0)]
    (* sign (/ 1.0 (Math/sqrt (* a xi (- 1.0 xi)))))))

(defn lagrange-C-matrix
  "**C_ij = ⟨φ_i | T + L(0) + V − E | φ_j⟩** on the Gauss–Legendre Lagrange mesh."
  [E l a n V]
  (let [n (long n) E (double E) l (long l) a (double a)
        {:keys [nodes]} (gauss-legendre-01 n)
        kin (lagrange-kinetic-bloch-matrix n a)
        scale (/ 1.0 f/mass-factor)]
    (mapv (fn [i]
            (mapv (fn [j]
                    (let [t (* scale (get-in kin [i j]))
                          r (* a (nth nodes i))]
                      (if (= i j)
                        (+ t (- (V r) E) (/ (* l (inc l)) (* f/mass-factor r r)))
                        t)))
                  (range n)))
          (range n))))

(defn- lagrange-Ra-from-C
  "**R a = u(a)/u'(a)** from eq. 3.15 (code convention; **B = 0**)."
  [n a nodes C]
  (let [scale (/ 1.0 f/mass-factor)
        phi-a (mapv #(lagrange-phi-at-boundary n a nodes %) (range n))
        Cinv (mat-inv C)]
    (* scale
       (reduce + 0.0
               (for [i (range n) j (range n)]
                 (* (phi-a i) (get-in Cinv [i j]) (phi-a j)))))))

(defn lagrange-Ra
  "Dimensionless **R a = u(a)/u'(a)** from the calculable R-matrix."
  [E l a n V]
  (let [{:keys [nodes]} (gauss-legendre-01 n)
        C (lagrange-C-matrix E l a n V)]
    (lagrange-Ra-from-C n a nodes C)))

(defn lagrange-r-matrix
  "**R = u(a)/(a u'(a))** (same as **`functions/r-matrix-from-numerov`**)."
  [E l a n V]
  (/ (lagrange-Ra E l a n V) (double a)))

(defn lagrange-s-matrix0
  "Neutral **S_L** from Lagrange **R a**, matching **`functions/s-matrix0`**."
  [E a L Ra]
  (let [k (Math/sqrt (* f/mass-factor (double E)))
        rho (* k (double a))
        h- (f/hankel0- L rho) h+ (f/hankel0+ L rho)
        h-- (f/deriv f/hankel0- L rho 1.0e-6)
        h++ (f/deriv f/hankel0+ L rho 1.0e-6)]
    (c/div (c/subt2 h- (c/mul Ra k h--))
           (c/subt2 h+ (c/mul Ra k h++)))))

(defn phase-shift-lagrange
  "Phase shift **δ_l(E)** from the Lagrange-mesh calculable R-matrix."
  ([E V a L] (phase-shift-lagrange E V a L 15))
  ([E V a L n]
   (let [V0 (first V) R0 (second V) a0 (last V)
         V-fn (fn [r] (f/WS r [V0 R0 a0]))
         Ra (lagrange-Ra E L a n V-fn)]
     (/ (c/arg (lagrange-s-matrix0 E a L Ra)) 2.0))))

(defn phase-shift-lagrange-numerov-Ra
  "Lagrange **S-matrix** pipeline with **R a** supplied by **`functions/r-matrix-a`** (reference)."
  [E V a L]
  (let [Ra (f/r-matrix-a E V (double a) L)]
    (/ (c/arg (lagrange-s-matrix0 E a L Ra)) 2.0)))

(defn lagrange-vs-numerov-phase-shift
  "Compare Lagrange-mesh and Numerov phase shifts at channel radius **a**.

  **Note:** the calculable **R a** from the discrete eq. 3.15 is still being tuned against
  **`functions/r-matrix-a`**; use **`:numerov-r-matrix-a`** in **`phase-shift-lagrange*`** for
  production parity with Numerov / **`phase-shift0`**."
  [E l V a n & {:keys [numerov-h] :or {numerov-h 0.001}}]
  (let [V0 (first V) R0 (second V) a0 (last V)
        delta-l (phase-shift-lagrange E V a l n)
        delta-ref (phase-shift-lagrange-numerov-Ra E V a l)
        u (f/solve-numerov E l V0 R0 a0 numerov-h a)
        delta-n (f/phase-shift-from-numerov u numerov-h a E l)]
    {:lagrange delta-l
     :numerov delta-n
     :numerov-r-matrix-a delta-ref
     :delta-diff (- delta-l delta-ref)}))

(defn phase-shift-lagrange*
  "Phase shift with selectable **R a** source: **`:lagrange`** (default) or **`:numerov-r-matrix-a`**."
  ([E V a L] (phase-shift-lagrange* E V a L 15 :lagrange))
  ([E V a L n method]
   (case method
     :numerov-r-matrix-a (phase-shift-lagrange-numerov-Ra E V a L)
     :lagrange (phase-shift-lagrange E V a L n)
     (throw (ex-info "phase-shift-lagrange*: unknown method" {:method method})))))
