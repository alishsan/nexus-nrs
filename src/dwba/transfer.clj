(ns dwba.transfer
  "DWBA calculations for single nucleon transfer reactions.
   
   This namespace implements bound state wavefunctions, transfer form factors,
   and cross-section calculations for reactions like (d,p), (p,d), etc.

   **Coulomb vs nuclear S (do not regress):** Elastic partial-wave **S_L^n** for the Woods–Saxon +
   R-matrix match lives only in **`functions/s-matrix`** — same quotient as **`functions/s-matrix0`**, with
   **`Hankel±`** when **η ≠ 0**. **Do not** define **S^n** by dividing some **“full”** amplitude by **e^{2iσ}**
   (that was a bad **M_L** convention). Coulomb **σ_L** for **(3.1.88)** belongs in **`functions/partial-wave-exp2sigma-Sn-minus-one`**
   / **`e^{2iσ}(S^n−1)`** only. Here, **`distorted-wave-optical`** **`:coulomb-tail`** only scales **|u|** to **|H_L^+|**;
   **`nuclear-phase-shifts-map`** uses **`functions/phase-shift`** (**½ arg(S^n)**). Austern **(5.6)** rows use
   **`coulomb-sigma-L`** + that **δ** map separately — keep those roles split.

   **DWBA distorted waves vs elastic:** Each partial wave **χ_L** here is the **same kind** of continuum solution as in
   elastic scattering: integrate with **U = V_C + V_opt(± iW)**, then **match** the interior to **Coulomb–Hankel**
   asymptotics via the **R-matrix** (log-derivative) at the matching radius — the elastic **`functions/s-matrix`**
   path exposes that match as **S_L^n** explicitly; **`distorted-wave-optical`** instead produces **χ(r)** on a grid
   and may align **|u|** to **|H_L^+|** at the outer edge (**`:coulomb-tail`**). Same matching idea; DWBA post
   form factors consume **χ** directly."
  (:require [fastmath.core :as m]
            [fastmath.special :as spec]
            [fastmath.polynomials :as poly]
            [functions :refer :all]
            [dwba.finite-well :refer :all]
            [dwba.form-factors :as ff]
            [dwba.angular-momentum :as jam]
            [complex :refer :all]))

;; ============================================================================
;; BOUND STATE WAVEFUNCTION SOLVER
;; ============================================================================

(defn bound-state-start [rho l]
  "Power series expansion for bound state near rho=0 (dimensionless).
   For bound states: u(rho) ≈ rho^(l+1) for small rho.
   This is the same as the naive start for scattering states."
  (m/pow rho (inc l)))

(defn woods-saxon-dimensionless [rho alpha]
  "Woods-Saxon potential in dimensionless form.
   
   Parameters:
   - rho: Dimensionless radius (r/R0)
   - alpha: Dimensionless diffuseness (a0/R0)
   
   Returns: Dimensionless potential v(rho) = V(r)/V0 = -1/(1 + exp((rho-1)/alpha))"
  (/ -1.0 (+ 1.0 (Math/exp (/ (- rho 1.0) alpha)))))

(defn f-rho-numerov-dimensionless [rho epsilon l lambda alpha]
  "Effective potential function for Numerov integration in dimensionless form.
   
   Parameters:
   - rho: Dimensionless radius (r/R0)
   - epsilon: Dimensionless energy (E/V0)
   - l: Orbital angular momentum quantum number
   - lambda: Dimensionless parameter λ = (2μ/ħ²) · V0 · R0²
   - alpha: Dimensionless diffuseness (a0/R0)
   
   Returns: f(rho) = λ · [v(rho) + l(l+1)/(rho²) - ε]
   
   The Schrödinger equation in dimensionless form is:
   -d²u/dρ² + f(ρ)u = 0
   where f(ρ) = λ[v(ρ) + l(l+1)/(ρ²) - ε]"
  (if (zero? rho)
    ;; At rho=0, centrifugal term dominates: l(l+1)/rho^2 -> infinity
    ;; But we never actually use rho=0 in Numerov (starts at rho=h_rho)
    Double/POSITIVE_INFINITY
    (let [v-rho (woods-saxon-dimensionless rho alpha)
          centrifugal (/ (* l (inc l)) (* rho rho))
          v-eff (+ v-rho centrifugal)]
      (* lambda (- v-eff epsilon)))))

(defn solve-bound-state-numerov 
([e l v0 rad diff m-f h r-max]
                                 "Solve the radial Schrödinger equation for a bound state using Numerov method
   in physical units.
   
   Parameters (all in physical units):
   - e: Energy in MeV (must be negative for bound states)
   - l: Orbital angular momentum quantum number
   - v0: Woods-Saxon potential depth in MeV
   - rad: R0 parameter (nuclear radius) in fm
   - diff: a0 parameter (surface diffuseness) in fm
   - h: Step size in fm
   - r-max: Maximum radius for integration in fm
   
   Returns: Vector of wavefunction values u(r) at each grid point (physical units).
   
   Note: For bound states, we expect u(r → ∞) → 0. This function
   just integrates; use find-bound-state-energy to find the correct energy."
                                 (let [steps (int (/ r-max h))
                                       ;; Initialize with bound state start: u(r) ≈ r^(l+1) for small r
                                       ;; For l=0: u(r) ≈ r, so u(h) ≈ h
                                       ;; For l=1: u(r) ≈ r^2, so u(h) ≈ h^2
                                       u0 0.0
                                       u1 (m/pow h (inc l))  ; u1 = h^(l+1) in physical units
                                       
                                       ;; Pre-calculate f(r) values for Numerov using f-r-numerov
                                       ;; f(r) = (2μ/ħ²) · [V_eff(r) - E]
                                       ;; For bound states, E < 0, so f(r) > 0 in classically allowed region
                                       fs (mapv (fn [r] 
                                                  (if (zero? r)
                                                    0.0  ; f(0) is infinite, but u(0)=0, so f(0)*u(0)=0
                                                    (f-r-numerov r e l v0 rad diff m-f)))
                                                (take (+ steps 2) (iterate #(+ % h) 0.0)))
                                       h2-12 (/ (* h h) 12.0)]
                                   
                                   (let [results (loop [n      (long 1)
                                                        u-prev (double u0)
                                                        u-curr (double u1)
                                                        acc    (transient [u0 u1])]
                                                   (if (>= n (dec steps))
                                                     (persistent! acc)
                                                     (let [;; Fast vector lookups with type hints
                                                           fn-1 (double (nth fs (dec n)))
                                                           fn   (double (nth fs n))
                                                           fn+1 (double (nth fs (inc n)))
                                                           
                                                           ;; Numerov step formula (physical units)
                                                           ;; Primitive arithmetic for performance
                                                           term1     (* 2.0 u-curr)
                                                           term2     (- u-prev)
                                                           inner-sum (+ (* 10.0 fn u-curr) 
                                                                        (* fn-1 u-prev))
                                                           term3     (* h2-12 inner-sum)
                                                           
                                                           numerator   (+ term1 term2 term3)
                                                           denominator (- 1.0 (* h2-12 fn+1))
                                                           
                                                           u-next      (/ numerator denominator)]
                                                       (recur (inc n)
                                                              u-curr
                                                              u-next
                                                              (conj! acc u-next)))))]
                                     results)))

  ([e l Vparams m-f]
    (solve-bound-state-numerov e l (first Vparams) (second Vparams) (last Vparams) m-f 0.01 20.)
    ))

(defn solve-bound-state-numerov-dimensionless 
  ([epsilon l lambda alpha h-rho rho-max]
   "Solve the radial Schrödinger equation for a bound state using Numerov method
   with DIMENSIONLESS variables.
   
   Parameters (all in physical units):
   - e: Energy in MeV (must be negative for bound states)
   - l: Orbital angular momentum quantum number
   - v0: Woods-Saxon potential depth in MeV
   - rad: R0 parameter (nuclear radius) in fm
   - diff: a0 parameter (surface diffuseness) in fm
   - h: Step size in fm
   

   
  Parameters: dimensionless variables:

   - ε = E/V0 (dimensionless energy)
  - l: Orbital angular momentum quantum number
   - λ = (2μ/ħ²) · V0 · R0² (dimensionless coupling)
   - α = a0/R0 (dimensionless diffuseness)
   - h_ρ = h/R0 (dimensionless step size)
   - rho-max: Maximum radius for integration (rho-max = r-max /rad)

   Returns: Vector of wavefunction values u(rho) at each grid point.

   Note:   - ρ = r/R0 (dimensionless radius)
   Note: For bound states, we expect u(r → ∞) → 0. This function
   just integrates; use find-bound-state-energy to find the correct energy."
   (let [               
         steps (int (/ rho-max h-rho))
         ;; Initialize with bound state start: u(rho) ≈ rho^(l+1)
         ;; Note: u(rho) is the dimensionless radial wavefunction
         ;; For l=0: u(rho) ≈ rho, so u(h_rho) ≈ h_rho
         ;; For l=1: u(rho) ≈ rho^2, so u(h_rho) ≈ h_rho^2
         u0 0.0
         u1 (bound-state-start h-rho l)  ; u1 = h_rho^(l+1) in dimensionless units
         
         ;; Pre-calculate f(rho) values for Numerov in dimensionless form
         ;; f(rho) = λ · [v(rho) + l(l+1)/(rho²) - ε]
         ;; For bound states, ε < 0, so f(rho) > 0 in classically allowed region
         fs (mapv (fn [rho] 
                    (if (zero? rho)
                      0.0  ; f(0) is infinite, but u(0)=0, so f(0)*u(0)=0
                      (f-rho-numerov-dimensionless rho epsilon l lambda alpha)))
                  (take (+ steps 2) (iterate #(+ % h-rho) 0.0)))
         h-rho2-12 (/ (* h-rho h-rho) 12.0)]
     
     (let [results (loop [n 1
                          results [u0 u1]]
                     (if (>= n (dec steps))
                       results
                       (let [un (get results n)
                             un-1 (get results (dec n))
                             fn-1 (get fs (dec n))
                             fn (get fs n)
                             fn+1 (get fs (inc n))
                             
                             ;; Numerov step formula (dimensionless)
                             numerator (+ (* 2.0 un) 
                                          (- un-1) 
                                          (* h-rho2-12 (+ (* 10.0 fn un) (* fn-1 un-1))))
                             denominator (- 1.0 (* h-rho2-12 fn+1))
                             un+1 (/ numerator denominator)]
                         (recur (inc n) (conj results un+1)))))]
       ;; The radial wavefunction u(r) should satisfy u(r) ≈ r^(l+1) for small r
       ;; In dimensionless: u(rho) ≈ rho^(l+1) where rho = r/R0
       ;; For l=0: u(rho) ≈ rho, so u(h_rho) = h_rho = h/R0
       ;; But we want u(h) = h in physical units
       ;; So we need: u(h) = R0 * u(h_rho) = R0 * (h/R0) = h ✓
       ;; Therefore, we scale by R0 to convert from dimensionless to physical
       results))))


(defn bound-state-boundary-value [u r-max h]
  "Check the boundary condition for a bound state.
   
   For a true bound state, u(r_max) should be approximately 0.
   Returns the value of u at r_max (should be close to 0 for bound state)."
  (let [idx (min (dec (count u)) (int (/ r-max h)))]
    (get u idx)))

(defn l-dot-s-nucleon
  "(l·s) for spin-½: (j(j+1)−l(l+1)−3/4)/2."
  [ell j]
  (/ (- (* j (+ j 1.0)) (* ell (+ ell 1.0)) 0.75)
     2.0))

(defn solve-bound-state-numerov-spin-orbit
  "Bound reduced radial **u(r)** on **r = n·h**: central Woods–Saxon **[V0 R a]** plus Thomas SO
   **(V_so, R_so, a_so)** for nucleon **j** (implicit **s = ½**). Same grid/Numerov loop as
   **`solve-bound-state-numerov`**."
  [e l v0 rad diff m-f h r-max v-so r-so a-so j]
  (let [lds (l-dot-s-nucleon (double l) (double j))
        steps (int (/ (double r-max) (double h)))
        u0 0.0
        u1 (m/pow (double h) (inc (long l)))
        fs (mapv (fn [r]
                   (if (zero? (double r))
                     0.0
                     (f-r-numerov-spin-orbit r e l v0 rad diff m-f v-so r-so a-so lds)))
                 (take (+ steps 2) (iterate #(+ % (double h)) 0.0)))
        h2-12 (/ (* (double h) (double h)) 12.0)]
    (let [results (loop [n      1
                          u-prev (double u0)
                          u-curr (double u1)
                          acc    (transient [u0 u1])]
                    (if (>= n (dec steps))
                      (persistent! acc)
                      (let [fn-1 (double (nth fs (dec n)))
                            fn   (double (nth fs n))
                            fn+1 (double (nth fs (inc n)))
                            term1      (* 2.0 u-curr)
                            term2      (- u-prev)
                            inner-sum  (+ (* 10.0 fn u-curr) (* fn-1 u-prev))
                            term3      (* h2-12 inner-sum)
                            numerator   (+ term1 term2 term3)
                            denominator (- 1.0 (* h2-12 fn+1))
                            u-next      (/ numerator denominator)]
                        (recur (inc n) u-curr u-next (conj! acc u-next)))))]
      results)))

(defn count-nodes [u]
  "Count the number of nodes (zeros) in the wavefunction.
   This helps identify the principal quantum number n.
   Number of radial nodes = n - l - 1
   (e.g., 1s: 0 nodes, 2s: 1 node, 2p: 0 nodes, 3s: 2 nodes, 3p: 1 node, 3d: 0 nodes)
   
   Note: 
   - We skip the initial region where u ≈ 0 (near r=0) to avoid
     counting the boundary condition as a node.
   - We also exclude the last few points near r_max to avoid counting
     spurious nodes from the shooting method boundary condition.
   A node is where the wavefunction crosses zero AFTER it has started
   and BEFORE it reaches the boundary region."
  (let [;; Find where wavefunction starts (becomes significantly non-zero)
        threshold 1e-6
        start-idx (loop [i 0]
                    (if (or (>= i (count u))
                            (> (Math/abs (get u i)) threshold))
                      i
                      (recur (inc i))))
        ;; Exclude last 5% of points near r_max to avoid boundary artifacts
        ;; This prevents counting spurious nodes from shooting method
        end-idx (max start-idx (- (count u) (max 10 (int (* 0.05 (count u))))))
        ;; Need at least 3 points between start and end to detect nodes reliably
        start-idx (max 2 (min start-idx (- end-idx 3)))]
    (if (>= start-idx end-idx)
      0  ; Can't count nodes if we don't have enough points
      (loop [n 0
             i (inc start-idx)
             prev-val (get u start-idx)
             prev-sign (m/signum prev-val)]
        (if (>= i end-idx)
          n
          (let [u-i (get u i)
                current-sign (m/signum u-i)
                ;; Count a node if:
                ;; 1. Sign changes (crossed zero)
                ;; 2. Previous value was non-zero (not starting from zero)
                ;; 3. Current value is non-zero (not exactly at zero, which we handle separately)
                ;; 4. The values on either side are significant (not noise)
                crossed-zero (and (not= prev-sign current-sign)
                                  (not (zero? prev-sign))
                                  (not (zero? current-sign))
                                  (> (Math/abs prev-val) threshold)
                                  (> (Math/abs u-i) threshold))]
            (recur (if crossed-zero (inc n) n)
                   (inc i)
                   u-i
                   (if (zero? u-i) prev-sign current-sign))))))))

(defn- boundary-value-spin-orbit-E
  [E l v0 rad diff m-f h r-max v-so r-so a-so j]
  (bound-state-boundary-value
   (solve-bound-state-numerov-spin-orbit E l v0 rad diff m-f h r-max v-so r-so a-so j)
   r-max h))

(defn find-eigenenergy-spin-orbit
  "**E** (MeV, negative) for fixed **V0** + SO by **E**-shooting (**u(r_max)→0**), optionally
   checking **radial node count** against **`target-nodes`** (pass **nil** to skip)."
  [v0 rad diff v-so r-so a-so l j m-f h r-max E-scan-lo E-scan-hi E-step target-nodes]
  (let [f #(boundary-value-spin-orbit-E % l v0 rad diff m-f h r-max v-so r-so a-so j)
        pts (mapv (fn [e] [e (f e)])
                  (take-while #(>= % (double E-scan-hi))
                              (iterate #(+ % (double E-step)) (double E-scan-lo))))
        br (first
            (for [i (range (dec (count pts)))
                  :let [[e0 b0] (nth pts i)
                        [e1 b1] (nth pts (inc i))]
                  :when (and (Double/isFinite b0) (Double/isFinite b1)
                             (not= (m/signum b0) (m/signum b1)))]
              [(min e0 e1) (max e0 e1)]))]
    (when br
      (let [[elo ehi] br
            ;; Interval-only bisection; |f(root)| for shooting is checked below (not |f|<tol).
            res (bisection f [elo ehi] 1.0e20 120)
            eroot (:root res)
            babs (when eroot (Math/abs (double (f eroot))))
            u (when (and eroot babs (< babs 5.0e-2))
                (solve-bound-state-numerov-spin-orbit eroot l v0 rad diff m-f h r-max v-so r-so a-so j))]
        (when (and u (or (nil? target-nodes) (= (long target-nodes) (long (count-nodes u)))))
          eroot)))))

(defn find-woods-saxon-v0-binding-spin-orbit
  "Bisect central **V0** so the **l,j** eigenenergy (Thomas SO + WS volume) matches **−E_bind**.
   Inner: **`find-eigenenergy-spin-orbit`** at each **V0** with **`target-radial-nodes`** (e.g. **0** for **1d_{5/2}**)."
  [E-bind l j v-so r-so a-so rad diff m-f h r-max v0-min v0-max v0-step target-radial-nodes]
  (let [e-t (- (Math/abs (double E-bind)))
        eigen (fn [v0]
                (find-eigenenergy-spin-orbit v0 rad diff v-so r-so a-so l j m-f h r-max
                                               (- (+ v0 5.0)) -0.05 0.25
                                               target-radial-nodes))
        pts (mapv (fn [v]
                    [v (when-let [ev (eigen v)] (- ev e-t))])
                  (take-while #(<= % (double v0-max))
                              (iterate #(+ % (double v0-step)) (double v0-min))))
        br (first
            (for [i (range (dec (count pts)))
                  :let [[v0a ga] (nth pts i)
                        [v0b gb] (nth pts (inc i))]
                  :when (and (some? ga) (some? gb)
                             (Double/isFinite ga) (Double/isFinite gb)
                             (not= (m/signum ga) (m/signum gb)))]
                  [(min v0a v0b) (max v0a v0b)]))]
    (if-not br
      {:converged? false :error "no V0 bracket for E-eigenvalue = -E_bind; widen scan or nodes"}
      (let [[a0 b0] br
            ;; If **eigen** fails at an interior **V0**, treat as unbound (**E→0**) so **g = −e_t > 0**.
            g (fn [v] (- (double (or (eigen v) 0.0)) e-t))
            res (bisection g [a0 b0] 0.02 80)
            v-root (:root res)
            e-found (when v-root (eigen v-root))]
        (-> res
            (assoc :E-bind E-bind :eigen-eV e-found :l l :j j
                   :target-radial-nodes target-radial-nodes)
            (assoc :v0-scan-range [v0-min v0-max]))))))

(defn scan-energy-range
  "Helper function to scan an energy range and compute wavefunctions.
   Returns a vector of candidate maps with :energy, :boundary-value, :nodes."
  [E-start E-end num-steps V-params l r-max h]
  (let [v0 (first V-params)
        rad (second V-params)
        diff (last V-params)
        E-step (/ (- E-end E-start) num-steps)]
    (for [i (range (inc num-steps))]
      (let [E (+ E-start (* i E-step))
            u (solve-bound-state-numerov E l v0 rad diff mass-factor h r-max)
            u-end (bound-state-boundary-value u r-max h)
            nodes (count-nodes u)]
        {:energy E
         :boundary-value u-end
         :nodes nodes}))))

(defn find-sign-changes
  "Find sign changes in boundary values across energy range.
   Returns candidates where boundary value changes sign (indicates bound state)."
  [candidates]
  (filter some?
    (for [i (range (dec (count candidates)))]
      (let [curr (nth candidates i)
            next (nth candidates (inc i))]
        (when (not= (m/signum (:boundary-value curr))
                   (m/signum (:boundary-value next)))
          (if (< (m/abs (:boundary-value curr))
                (m/abs (:boundary-value next)))
            curr
            next))))))


(defn find-sign-change-pairs
  "Find sign changes in boundary values across energy range.
   Returns a list of vectors, where each vector contains two energies:
   [E1 E2] where E1 and E2 are the energies on either side of a zero crossing.
   This indicates a bound state exists between E1 and E2."
  [candidates]
  (filter some?
    (for [i (range (dec (count candidates)))]
      (let [curr (nth candidates i)
            next (nth candidates (inc i))]
        (when (not= (m/signum (:boundary-value curr))
                   (m/signum (:boundary-value next)))
          [(:energy curr) (:energy next)])))))

(defn find-energy-with-nodes
  "Coarse scan to find bound state energies with specific number of nodes.
   
   Parameters:
   - E-start, E-end: Energy range to search
   - num-steps: Number of energy points to scan
   - target-nodes: Desired number of nodes
   - V-params: Woods-Saxon parameters [V0, R0, a0]
   - l: Orbital angular momentum
   - r-max: Maximum radius
   - h: Step size
   
   Returns: List of candidates with {:energy E, :wavefunction u, :boundary-value u-end, :nodes n}"
  [E-start E-end num-steps target-nodes V-params l r-max h]
  (let [all-candidates (scan-energy-range E-start E-end num-steps V-params l r-max h)
        negative-candidates (filter (fn [c] (< (:energy c) -0.01)) all-candidates)
        sign-changes (find-sign-changes negative-candidates)
        candidates-with-nodes (filter (fn [c] (= (:nodes c) target-nodes)) sign-changes)]
    (if (seq candidates-with-nodes)
      candidates-with-nodes
      (filter (fn [c] (= (:nodes c) target-nodes)) negative-candidates))))



(defn get-refinement-energy-range [E-guess v0 boundary-value]
  "Calculate energy range for refinement search.
   
   Parameters:
   - E-guess: Initial energy guess
   - v0: Potential depth
   - boundary-value: Current boundary value (if large, use wider range)
   
   Returns: [E-lo E-hi]
   
   Uses a narrower range to avoid jumping to different quantum states,
   but expands if boundary value is large (not a true bound state).
   For very large boundary values, search a much wider range to find sign changes."
  (let [;; Use much wider range if boundary value is very large (not a true bound state)
        ;; This helps find sign changes that might be far from the coarse candidate
        E-range (cond
                  (> (Math/abs boundary-value) 1e10)  ; Very large boundary value
                  (* v0 0.5)  ; Search ±50% of V0 to find sign changes
                  
                  (> (Math/abs boundary-value) 100.0)  ; Large boundary value
                  (min 50.0 (* v0 0.3))  ; Wider range: ±50 MeV or 30% of V0
                  
                  :else
                  (min 10.0 (* v0 0.2)))  ; Narrower range: ±10 MeV or 20% of V0
        E-lo (max (- E-guess E-range) (- v0))
        E-hi (min (+ E-guess E-range) -0.1)]
    [E-lo E-hi]))

(defn create-boundary-value-function [V-params l r-max h]
  "Create function f(E) = u(r_max) for root finding."
  (let [v0 (first V-params)
        rad (second V-params)
        diff (last V-params)]
    (fn [E]
      (bound-state-boundary-value 
       (solve-bound-state-numerov E l v0 rad diff mass-factor h r-max) r-max h))))

(defn validate-secant-root [root E-guess E-lo E-hi v0]
  "Validate and clamp secant root to valid range.
   
   Returns: validated root or nil if invalid"
  (let [max-deviation 20.0
        E-min-allowed (max E-lo (- E-guess max-deviation))
        E-max-allowed (min E-hi (+ E-guess max-deviation))
        clamped (max E-min-allowed (min E-max-allowed (min root -0.1)))]
    (when (and (< clamped (* -1.0 v0 0.05))  ; Must be at least 5% of well depth
               (< clamped -0.01)
               (> clamped (- v0)))
      clamped)))

(defn find-sign-change-range [f E-lo E-hi search-points]
  "Find energy range where boundary value changes sign.
   
   Returns: [E-lo E-hi] of sign change range, or nil if none found"
  (loop [i 0
         found-range nil]
    (if (or (>= i search-points) found-range)
      found-range
      (let [E-test (+ E-lo (* i (/ (- E-hi E-lo) (dec search-points))))
            u-test-val (f E-test)
            E-next (if (< i (dec search-points))
                     (+ E-lo (* (inc i) (/ (- E-hi E-lo) (dec search-points))))
                     E-hi)
            u-next-val (f E-next)]
        (if (not= (m/signum u-test-val) (m/signum u-next-val))
          (recur (inc i) [E-test E-next])
          (recur (inc i) found-range))))))

(defn create-refined-result [E-root target-nodes V-params l r-max h]
  "Create result map from refined energy."
  (let [v0 (first V-params)
        rad (second V-params)
        diff (last V-params)
        u-final (solve-bound-state-numerov E-root l v0 rad diff mass-factor h r-max)
        u-final-val (bound-state-boundary-value u-final r-max h)
        nodes (count-nodes u-final)]
    {:energy E-root
     :wavefunction u-final
     :boundary-value u-final-val
     :nodes nodes}))

(defn refine-with-secant [f E-guess E-lo E-hi u-lo-val u-mid-val u-hi-val v0 tolerance target-nodes V-params l r-max h]
  "Try to refine using secant method.
   
   Returns: result map or nil if secant fails.
   IGNORES node count - only checks if boundary value is small."
  (let [E-secant-0 (if (< (Math/abs u-mid-val) (Math/abs u-lo-val)) E-guess E-lo)
        E-secant-1 (if (< (Math/abs u-hi-val) (Math/abs u-mid-val)) E-hi E-guess)
        secant-result (secant f E-secant-0 E-secant-1 tolerance 50)
        ;; Check precision: |value| < tolerance
        secant-precise? (< (Math/abs (:value secant-result)) tolerance)
        E-secant-root (when secant-precise?
                       (validate-secant-root (:root secant-result) E-guess E-lo E-hi v0))
        result (when (and secant-precise? E-secant-root)
                 (create-refined-result E-secant-root target-nodes V-params l r-max h))]
    ;; Return result if boundary value is reasonable (ignore node count)
    (when (and result (< (Math/abs (:boundary-value result)) 1e6))
      result)))

(defn refine-with-bisection [f E-lo E-hi tolerance target-nodes V-params l r-max h]
  "Try to refine using bisection method.
   
   Returns: result map or nil if bisection fails.
   IGNORES node count - only checks if boundary value is small."
  (let [sign-change-range (find-sign-change-range f E-lo E-hi 20)
        [E-bisect-lo E-bisect-hi] (or sign-change-range [E-lo E-hi])
        bisection-result (bisection f [E-bisect-lo E-bisect-hi] tolerance 100)
        ;; Check precision: |value| < tolerance
        bisection-precise? (< (Math/abs (:value bisection-result)) tolerance)
        result (when bisection-precise?
                 (create-refined-result (:root bisection-result) target-nodes V-params l r-max h))]
    ;; Return result if boundary value is reasonable (ignore node count)
    (when (and result (< (Math/abs (:boundary-value result)) 1e6))
      result)))

(defn refine-with-grid-search [E-lo E-hi target-nodes V-params l r-max h]
  "Fallback: find minimum boundary value using grid search.
   
   Returns: result map with smallest boundary value, or nil if none found.
   IGNORES node count - only looks for minimum boundary value."
  (let [v0 (first V-params)
        rad (second V-params)
        diff (last V-params)
        search-points 100
        candidates (for [i (range (inc search-points))]
                     (let [E-test (+ E-lo (* i (/ (- E-hi E-lo) search-points)))
                           u-test (solve-bound-state-numerov E-test l v0 rad diff mass-factor h r-max)
                           u-test-val (bound-state-boundary-value u-test r-max h)
                           nodes-test (count-nodes u-test)]
                       {:energy E-test
                        :wavefunction u-test
                        :boundary-value u-test-val
                        :boundary-abs (Math/abs u-test-val)
                        :nodes nodes-test}))
        ;; Find candidate with smallest boundary value (ignore node count)
        best (when (seq candidates)
               (apply min-key :boundary-abs candidates))]
    (when (and best (< (:boundary-abs best) 1e6))
      {:energy (:energy best)
       :wavefunction (:wavefunction best)
       :boundary-value (:boundary-value best)
       :nodes (:nodes best)})))

(defn has-sign-change? [u-lo-val u-mid-val u-hi-val]
  "Check if boundary values have opposite signs."
  (or (not= (m/signum u-lo-val) (m/signum u-hi-val))
      (not= (m/signum u-lo-val) (m/signum u-mid-val))
      (not= (m/signum u-mid-val) (m/signum u-hi-val))))

(defn refine-bound-state-energy
  "Refine bound state energy using secant method (or bisection as fallback).
   
   Parameters:
   - E-guess: Initial energy guess
   - target-nodes: Desired number of nodes
   - V-params: Woods-Saxon parameters [V0, R0, a0]
   - l: Orbital angular momentum
   - r-max: Maximum radius
   - h: Step size
   - tolerance: Energy convergence tolerance
   
   Returns: {:energy E, :wavefunction u, :boundary-value u-end, :nodes n}
   Check precision using |boundary-value| < tolerance
   or nil if no valid refinement found with correct node count.
   
   Uses secant method first (faster for smooth functions), falls back to bisection
   if secant doesn't converge or if boundary values have opposite signs."
  [E-guess target-nodes V-params l r-max h tolerance]
  (let [v0 (first V-params)
        rad (second V-params)
        diff (last V-params)
        ;; Get initial boundary value to determine refinement range
        u-guess (solve-bound-state-numerov E-guess l v0 rad diff mass-factor h r-max)
        u-guess-val (bound-state-boundary-value u-guess r-max h)
        [E-lo E-hi] (get-refinement-energy-range E-guess v0 u-guess-val)
        f (create-boundary-value-function V-params l r-max h)
        u-lo (solve-bound-state-numerov E-lo l v0 rad diff mass-factor h r-max)
        u-hi (solve-bound-state-numerov E-hi l v0 rad diff mass-factor h r-max)
        u-mid (solve-bound-state-numerov E-guess l v0 rad diff mass-factor h r-max)
        u-lo-val (bound-state-boundary-value u-lo r-max h)
        u-hi-val (bound-state-boundary-value u-hi r-max h)
        u-mid-val (bound-state-boundary-value u-mid r-max h)
        secant-result (refine-with-secant f E-guess E-lo E-hi u-lo-val u-mid-val u-hi-val v0 tolerance target-nodes V-params l r-max h)
        bisection-result (when (and (not secant-result) (has-sign-change? u-lo-val u-mid-val u-hi-val))
                          (refine-with-bisection f E-lo E-hi tolerance target-nodes V-params l r-max h))
        grid-result (when (and (not secant-result) (not bisection-result))
                     (refine-with-grid-search E-lo E-hi target-nodes V-params l r-max h))
        result (or secant-result bisection-result grid-result)]
    ;; Return result if boundary value is reasonable (ignore node count)
    (when (and result (< (Math/abs (:boundary-value result)) 1e6))
      result)))

(defn get-energy-search-range [n l v0]
  "Calculate energy search range for principal quantum number n and orbital angular momentum l.
   
   For bound states, we need to search from deep in the well (near -V0) up to near zero.
   Higher l states can be either deeper or shallower depending on the potential.
   
   Returns: [E-min E-max] in MeV"
  (let [;; For bound states, search from deep in well to near zero
        ;; Start deeper for higher n (more nodes) and account for l
        ;; For l > 0, we need to search a wider range since centrifugal barrier affects energy
        base-E-min (cond
                     (= n 1) (- (* v0 0.8))  ; Ground: search from 80% of V0
                     (= n 2) (- (* v0 0.6))  ; 2s/2p: search from 60% of V0
                     (= n 3) (- (* v0 0.4))  ; 3s/3p/3d: search from 40% of V0
                     :else (- (* v0 (- 0.8 (* (- n 1) 0.15)))))  ; Higher n: progressively higher
        base-E-max (cond
                     (= n 1) (- (* v0 0.2))  ; Ground: up to 20% of V0
                     (= n 2) (- (* v0 0.15))  ; 2s/2p: up to 15% of V0
                     (= n 3) (- (* v0 0.1))  ; 3s/3p/3d: up to 10% of V0
                     :else (- (* v0 (- 0.2 (* (- n 1) 0.02)))))  ; Higher n: closer to zero
        ;; For l > 0, extend the search range deeper (l=1 states can be quite deep)
        ;; and also shallower (centrifugal barrier pushes some states up)
        l-deepening (* l v0 0.1)  ; Each l extends search deeper by 10% of V0
        E-min (max (+ base-E-min (- l-deepening)) (- v0))  ; Extend deeper (more negative), but not beyond well depth
        E-max (min base-E-max -0.1)]  ; Don't go above zero
    [E-min E-max]))

(defn try-wider-search [v0 expected-nodes V-params l r-max h tolerance]
  "Try a wider energy search range.
   
   Returns: best result found, or nil if nothing found"
  (let [;; Wider range: search from 80% to 10% of V0
        E-wide-min (- (* v0 0.8))
        ;; Adjust upper bound for l: higher l states are at higher energies
        base-E-wide-max (- (* v0 0.1))
        l-adjustment (* l v0 0.05)  ; Same adjustment as in get-energy-search-range
        E-wide-max (min (+ base-E-wide-max l-adjustment) -0.1)
        wide-candidates (find-energy-with-nodes E-wide-min E-wide-max 200 expected-nodes 
                                                V-params l r-max h)
        wide-result (first wide-candidates)]
    (if (and wide-result (< (:energy wide-result) -0.01))
      ;; Only refine if boundary value is reasonable (not huge) - IGNORE node count
      (let [refined (refine-bound-state-energy (:energy wide-result) expected-nodes V-params l r-max h tolerance)]
        ;; Only return if refined result has reasonable boundary value
        (when (and refined
                   (< (Math/abs (:boundary-value refined)) 1e6))  ; Reject huge boundary values
          refined))
      nil)))

(defn valid-energy? [E E-min E-max]
  "Check if energy is valid (negative and within range)."
  (and (< E -0.01)
       (>= E E-min)
       (<= E E-max)))

(defn refinement-improved? [refined coarse-result coarse-boundary expected-nodes]
  "Check if refinement significantly improved the result.
   
   IGNORES node count - only checks if boundary value improved."
  (let [refined-boundary (Math/abs (:boundary-value refined))]
    (< refined-boundary coarse-boundary)))    ; Must improve boundary value

(defn print-refinement-debug [coarse-boundary refined-boundary coarse-result refined-nodes expected-nodes]
  "Print debug information about refinement."
  (when (> coarse-boundary 10.0)
    (println (format "  Refinement: coarse boundary=%.2e, refined boundary=%.2e, improvement=%.1f%%"
                   coarse-boundary refined-boundary 
                   (* 100.0 (/ (- coarse-boundary refined-boundary) coarse-boundary)))
    (println (format "  Coarse nodes=%d, refined nodes=%d, expected=%d"
                   (:nodes coarse-result) refined-nodes expected-nodes)))))

(defn try-refinement-with-wide-search [E-guess expected-nodes coarse-result coarse-boundary
                                       E-search-min E-search-max v0 V-params l r-max h tolerance]
  "Try to refine the energy, with fallback to wider search if needed.
   
   Returns: refined result if boundary value improved (ignores node count)"
  (let [refined (refine-bound-state-energy E-guess expected-nodes V-params l r-max h tolerance)]
    (if refined
      (let [refined-boundary (Math/abs (:boundary-value refined))
            refined-energy (:energy refined)
            refined-nodes (:nodes refined)]
        (print-refinement-debug coarse-boundary refined-boundary coarse-result refined-nodes expected-nodes)
        ;; Accept refinement if it improves boundary value (ignore node count and energy range)
        ;; The energy might be outside the original search range but still valid
        (if (refinement-improved? refined coarse-result coarse-boundary expected-nodes)
          refined
          ;; If refinement doesn't improve, try wider search or return the better one
          (let [wide-refined (try-wider-search v0 expected-nodes V-params l r-max h tolerance)]
            (if (and wide-refined
                     (< (Math/abs (:boundary-value wide-refined)) coarse-boundary))
              wide-refined
              ;; Use the one with smaller boundary value (refined or coarse)
              (if (< refined-boundary coarse-boundary)
                refined
                coarse-result))))
      ;; If refinement returned nil, use coarse result
      coarse-result)))
)

(defn should-refine? [coarse-result expected-nodes]
  "Check if refinement should be attempted.
   IGNORES node count - only checks boundary value."
  (> (Math/abs (:boundary-value coarse-result)) 1.0))

(defn handle-invalid-energy [E-guess E-search-min E-search-max v0 expected-nodes V-params l r-max h tolerance coarse-result]
  "Handle case when coarse scan finds invalid energy."
  (println (format "Warning: Coarse scan found invalid energy %.6f MeV (expected range: [%.2f, %.2f] MeV)"
                 E-guess E-search-min E-search-max))
  (or (try-wider-search v0 expected-nodes V-params l r-max h tolerance)
      coarse-result))

(defn handle-wrong-nodes [v0 expected-nodes V-params l r-max h tolerance coarse-result]
  "Handle case when coarse scan finds wrong number of nodes."
  (or (try-wider-search v0 expected-nodes V-params l r-max h tolerance)
      (do
        (println (format "Warning: Could not find state with %d nodes. Found %d nodes at E=%.2f MeV"
                       expected-nodes (:nodes coarse-result) (:energy coarse-result)))
        coarse-result)))


(defn find-bound-state-energy
  "Find bound state energy using shooting method.
   
   Parameters:
   - V-params: Woods-Saxon parameters [V0, R0, a0]
   - l: Orbital angular momentum
   - n: (optional) Principal/radial quantum number; when provided (5-arg), used for search range
   - r-max: Maximum radius for integration
   - h: Step size
   - E-min: Minimum energy to search (default: -V0, the potential depth)
   - E-max: Maximum energy to search (default: -0.1 MeV, just below zero)
   - tolerance: Energy convergence tolerance (default: 0.01 MeV)
   
   Returns: {:energy E, :wavefunction u, :nodes n-nodes, :boundary-value u-end}
   Check precision using |boundary-value| < tolerance
   
   Algorithm:
   1. Coarse scan to find approximate energy with n radial nodes
   2. Refine using secant/bisection around best candidate
   3. For different n values, searches in different energy ranges
       because bound states are ordered: E(n=0) < E(n=1) < E(n=2) < ..."
  ([V-params l r-max h]
   (let [v0 (first V-params)]
     (find-bound-state-energy V-params l r-max h (- v0) -0.1 0.01)))
  ([V-params l n r-max h]
   ;; 5-arg: n = principal quantum number (1,2,...), target radial nodes = n - 1
   (let [v0 (first V-params)
         target-nodes (max 0 (dec n))
         tolerance 0.01
         [E-min E-max] (get-energy-search-range n l v0)
         candidates (find-energy-with-nodes E-min E-max 150 target-nodes V-params l r-max h)
         best (first (sort-by #(Math/abs (:boundary-value %)) (or candidates [])))
         E-guess (when best (:energy best))
         refined (when E-guess (refine-bound-state-energy E-guess target-nodes V-params l r-max h tolerance))]
     (or refined
         (try-wider-search v0 target-nodes V-params l r-max h tolerance)
         ;; fallback: 4-arg (first sign change, often ground state)
         (find-bound-state-energy V-params l r-max h))))
  ([V-params l r-max h E-min E-max tolerance]
   (let [v0 (first V-params)
         rad (second V-params)
         diff (last V-params)
         [E-search-min E-search-max] [(- v0) -0.1]
         coarse-candidates (scan-energy-range E-search-min E-search-max 100 V-params l r-max h)
         sign-change-pairs (find-sign-change-pairs coarse-candidates)]
     (when (seq sign-change-pairs)
       (let [first-pair (first sign-change-pairs)
             [E1 E2] first-pair
             root-result (bisection (fn [x] (last (solve-bound-state-numerov x l v0 rad diff mass-factor h r-max))) 
                                    [E1 E2] tolerance 100)
             E-refined (:root root-result)
             u-refined (solve-bound-state-numerov E-refined l v0 rad diff mass-factor h r-max)
             u-end (bound-state-boundary-value u-refined r-max h)
             nodes (count-nodes u-refined)]
         {:energy E-refined
          :wavefunction u-refined
          :boundary-value u-end
          :nodes nodes})))))


(defn normalize-bound-state [u h]
  "Normalize bound state wavefunction so that ∫₀^∞ |u(r)|² dr = 1.
   
   Uses Simpson's rule for integration: ∫ f(r) dr ≈ (h/3) * [f₀ + fₙ + 4∑f_odd + 2∑f_even]
   
   Parameters:
   - u: Wavefunction vector
   - h: Step size
   
   Returns: Normalized wavefunction vector"
  (when (or (nil? u) (empty? u))
    (throw (IllegalArgumentException. 
            (format "Cannot normalize empty or nil wavefunction. Wavefunction: %s" u))))
  (let [;; Calculate normalization integral: N² = ∫ u²(r) r² dr
        ;; For radial wavefunctions, normalization is ∫ u²(r) dr (not r²)
        ;; But in some conventions it's ∫ u²(r) r² dr - we'll use ∫ u²(r) dr
        integrand (mapv #(* % %) u)
        n (count integrand)]
    (when (< n 2)
      (throw (IllegalArgumentException. 
              (format "Wavefunction too short for normalization: %d points (need at least 2)" n))))
    (let [;; Simpson's rule
          simpson-sum (loop [i 1
                             sum 0.0]
                        (if (>= i (dec n))
                          sum
                          (let [coeff (if (odd? i) 4.0 2.0)
                                term (* coeff (get integrand i))]
                            (recur (inc i) (+ sum term)))))
          integral (* (/ h 3.0) 
                      (+ (first integrand) 
                         (last integrand) 
                         simpson-sum))
          norm-factor (Math/sqrt integral)
          ;; Avoid division by zero
          norm-factor (if (zero? norm-factor) 1.0 (/ 1.0 norm-factor))]
      (mapv #(* % norm-factor) u))))

(defn solve-bound-state
  "Main function to solve for a bound state wavefunction.
   
   Parameters:
   - V-params: Woods-Saxon parameters [V0, R0, a0]
   - n: Principal quantum number (1, 2, 3, ...)
   - l: Orbital angular momentum (0, 1, 2, ...)
   - j: Total angular momentum (l ± 1/2 for nucleons, but we'll ignore spin-orbit for now)
   - r-max: Maximum radius (default: 20.0 fm)
   - h: Step size (default: 0.01 fm)
   
   Returns: {:energy E, :wavefunction u, :normalized-wavefunction u-norm, 
             :nodes n-nodes, :quantum-numbers {:n n, :l l, :j j}}
   
   Example:
   (solve-bound-state [50.0 2.0 0.6] 1 0 nil)
   => Finds 1s bound state in Woods-Saxon well"
  ([V-params n l]
   (solve-bound-state V-params n l nil 20.0 0.01))
  ([V-params n l j]
   (solve-bound-state V-params n l j 20.0 0.01))
  ([V-params n l j r-max h]
   ;; For nuclear potentials, n represents the number of radial nodes (not principal quantum number)
   ;; No constraint like n > l applies - nuclear bound states are labeled by radial nodes
   (let [result (find-bound-state-energy V-params l n r-max h)
         wavefunction (:wavefunction result)]
     (when (or (nil? wavefunction) (empty? wavefunction))
       (throw (IllegalArgumentException. 
               (format "Failed to find bound state for n=%d, l=%d. Result: %s" 
                      n l (pr-str (select-keys result [:energy :nodes :boundary-value]))))))
     (let [u-norm (normalize-bound-state wavefunction h)]
       {:energy (:energy result)
        :wavefunction (:wavefunction result)
        :normalized-wavefunction u-norm
        :nodes (:nodes result)
        :boundary-value (:boundary-value result)
        :quantum-numbers {:n n, :l l, :j j}
        :r-max r-max
        :h h}))))

;; ============================================================================
;; ASYMPTOTIC NORMALIZATION COEFFICIENT (ANC)
;; ============================================================================

(defn calculate-anc
  "Calculate Asymptotic Normalization Coefficient (ANC) from a MODEL bound state wavefunction.
   
   The ANC C is the coefficient in the asymptotic form:
   u(r) → C · e^(-κr) / r^l  (for large r)
   
   This fits the tail of the wavefunction produced by the potential used in
   solve-bound-state (e.g. Woods-Saxon). Prefer using ANCs from experimental
   data (e.g. (d,p) analysis) when available; pass those directly to
   anc-normalized-overlap or other callers instead of calling this.
   
   Parameters:
   - u: Bound state wavefunction (vector)
   - E: Bound state energy (MeV, negative)
   - l: Orbital angular momentum
   - mass-factor: Mass factor (2μ/ħ²) in MeV⁻¹·fm⁻²
   - h: Step size (fm)
   - r-fit-min: Minimum radius for fitting (fm, should be outside nuclear potential)
   - r-fit-max: Maximum radius for fitting (fm)
   
   Returns: ANC value C (in fm^(-l-1/2) for normalized u)
   
   Example:
   (let [result (solve-bound-state [50.0 2.0 0.6] 1 0)
         u (:normalized-wavefunction result)
         E (:energy result)]
     (calculate-anc u E 0 mass-factor 0.01 5.0 15.0))"
  [u E l mass-factor h r-fit-min r-fit-max]
  (let [;; Decay constant: κ = sqrt(2μ|E|)/ħ = sqrt(mass-factor · |E|)
        kappa (Math/sqrt (* mass-factor (Math/abs E)))
        ;; Find indices for fitting region
        idx-min (int (/ r-fit-min h))
        idx-max (int (/ r-fit-max h))
        idx-max-safe (min idx-max (dec (count u)))
        ;; Extract wavefunction values in fitting region
        u-fit (mapv (fn [i]
                     (let [r (* i h)
                           u-val (get u i)]
                       {:r r :u u-val}))
                   (range idx-min (inc idx-max-safe)))
        ;; Fit to asymptotic form: u(r) = C · e^(-κr) / r^l
        ;; Taking logarithm: ln(u · r^l) = ln(C) - κr
        ;; This is a linear fit: y = a + b·r, where y = ln(u·r^l), a = ln(C), b = -κ
        fit-data (filter (fn [{:keys [r u]}]
                          (and (> (Math/abs u) 1e-10)  ; Avoid log of zero
                               (> r 0.1)))  ; Avoid r=0 issues
                        u-fit)
        ;; Perform linear fit: ln(u·r^l) = ln(C) - κr
        fit-values (mapv (fn [{:keys [r u]}]
                          (let [r-l (if (zero? l)
                                    1.0
                                    (Math/pow r l))
                                u-times-rl (* u r-l)
                                log-u (if (> u-times-rl 0)
                                       (Math/log u-times-rl)
                                       -100.0)]  ; Large negative for very small values
                            {:r r :y log-u}))
                        fit-data)]
    (if (and (seq fit-values) (> (count fit-values) 2))
      (let [;; Linear regression: y = a + b·r
            ;; a = ln(C), b = -κ
            n-fit (count fit-values)
            sum-r (reduce + (map :r fit-values))
            sum-y (reduce + (map :y fit-values))
            sum-r2 (reduce + (map #(* (:r %) (:r %)) fit-values))
            sum-ry (reduce + (map #(* (:r %) (:y %)) fit-values))
            ;; Calculate slope and intercept
            denominator (- (* n-fit sum-r2) (* sum-r sum-r))
            slope (if (> (Math/abs denominator) 1e-10)
                   (/ (- (* n-fit sum-ry) (* sum-r sum-y)) denominator)
                   (- kappa))  ; Fallback to expected value
            intercept (if (> (Math/abs denominator) 1e-10)
                       (/ (- (* sum-y sum-r2) (* sum-r sum-ry)) denominator)
                       (Math/log (Math/abs (first (map :u fit-data)))))
            ;; ANC: C = exp(intercept)
            anc (Math/exp intercept)
            ;; Verify: slope should be approximately -κ
            slope-check (Math/abs (+ slope kappa))]
        ;; Return ANC, but check if fit is reasonable
        (if (< slope-check (* 0.5 kappa))  ; Slope should be close to -κ
          anc
          ;; If fit is poor, estimate from single point at large r
          (let [last-point (last fit-values)
                r-last (:r last-point)
                u-last (get u (int (/ r-last h)))
                r-l-power (if (zero? l) 1.0 (Math/pow r-last l))
                exp-kappa-r (Math/exp (* kappa r-last))
                anc-estimate (* u-last r-l-power exp-kappa-r)]
            anc-estimate)))
      ;; If not enough points, estimate from last point
      (let [last-idx (min (dec (count u)) idx-max-safe)
            r-last (* last-idx h)
            u-last (get u last-idx)
            r-l-power (if (zero? l) 1.0 (Math/pow r-last l))
            exp-kappa-r (Math/exp (* kappa r-last))
            anc-estimate (* u-last r-l-power exp-kappa-r)]
        anc-estimate))))

(defn anc-normalized-overlap
  "Calculate overlap integral normalized by ANC product.
   
   O_ANC = ∫ φ*_f(r) φ_i(r) r² dr / (C_f · C_i)
   
   Prefer passing ANCs from experiment or literature (e.g. from (d,p) or
   transfer analyses) when available. If not, you can use calculate-anc
   from model wavefunctions, but that is less reliable than structure
   or experimental ANCs.
   
   Parameters:
   - phi-i: Initial bound state wavefunction
   - phi-f: Final bound state wavefunction
   - ANC-i: ANC of initial state (e.g. from experiment when available)
   - ANC-f: ANC of final state (e.g. from experiment when available)
   - r-max: Maximum radius (fm)
   - h: Step size (fm)
   
   Returns: Normalized overlap O_ANC
   
   Example (using experimental ANCs when available):
   (anc-normalized-overlap phi-i phi-f 0.78 0.14 r-max h)
   Example (fallback to model extraction):
   (let [ANC-i (calculate-anc phi-i E-i l-i mass-factor h 5.0 15.0)
         ANC-f (calculate-anc phi-f E-f l-f mass-factor h 5.0 15.0)]
     (anc-normalized-overlap phi-i phi-f ANC-i ANC-f r-max h))"
  [phi-i phi-f ANC-i ANC-f r-max h]
  (let [overlap (ff/overlap-integral phi-i phi-f r-max h)
        anc-product (* ANC-i ANC-f)]
    (if (> anc-product 1e-10)
      (/ overlap anc-product)
      0.0)))

;; ============================================================================
;; UTILITY FUNCTIONS
;; ============================================================================

(defn bound-state-energy-approx [V-params n l]
  "Rough estimate of bound state energy using infinite square well approximation.
   Useful for providing initial guess: E ≈ -V0 + (n²π²ħ²)/(2mR²)
   
   This is just a rough estimate; actual energy will be different due to
   finite well depth and Woods-Saxon shape."
  (let [V0 (first V-params)
        R0 (second V-params)
        ;; Infinite square well: E_n = n²π²ħ²/(2mR²) - V0
        ;; For bound state, E < 0, so: E ≈ -V0 + n²π²ħ²/(2mR²)
        ;; But this gives positive energy, so we need: E ≈ -V0 + ...
        ;; Actually, for a well of depth V0, ground state is roughly at -V0/2
        ;; For excited states: E_n ≈ -V0 + n² * (V0/4) for n=1,2,3,...
        ]
    (- V0 (* n n (/ V0 4.0)))))

(defn plot-bound-state-info [result]
  "Print information about a bound state solution.
   Useful for debugging and verification."
  (println "=== Bound State Information ===")
  (println (format "Quantum numbers: n=%d, l=%d, j=%s"
                   (get-in result [:quantum-numbers :n])
                   (get-in result [:quantum-numbers :l])
                   (get-in result [:quantum-numbers :j])))
  (println (format "Energy: %.6f MeV" (:energy result)))
  (println (format "Number of nodes: %d" (:nodes result)))
  (println (format "Boundary value at r_max: %.6e" (:boundary-value result)))
  (let [boundary-precise? (< (Math/abs (:boundary-value result)) 10.0)]
    (println (format "Precision: |boundary-value| < 10.0: %s" boundary-precise?)))
  (println (format "Wavefunction length: %d points" (count (:normalized-wavefunction result))))
  (println ""))

;; ============================================================================
;; Radial convention: reduced u(r)=r·R(r) → R(r) for volume overlaps
;; ============================================================================

(defn radial-R-from-reduced-u
  "Convert reduced radial amplitudes u(r)=r·R(r) (Numerov convention in this namespace)
  to the radial factor R(r)=u(r)/r used in the **three-dimensional** zero-range POST overlap

    ∫ χ_f^*(\\vec r) χ_i(\\vec r) φ_f^*(\\vec r) φ_i(\\vec r) d³r
    → 4π ∫ R_{χ_f} R_{χ_i} R_{φ_f} R_{φ_i} r² dr   (spherical, s-wave cluster / fixed L block).

  All χ and φ passed to `transfer-amplitude-post` are **u** on an r = i·h grid; using u in place
  of R mis-scales the integrand by O(r⁴) relative to the correct r² measure and blows the
  cross section vs. flux-normalized codes (e.g. DWUCK).

  At r=0 we return R=0: the Simpson integrand is weighted by r², so this point contributes
  nothing; avoids 0/0 for u(0)=0."
  [u-vec ^double h]
  (mapv (fn [^long i u-val]
          (let [r (* (double i) h)]
            (if (< r 1e-14)
              (if (number? u-val) 0.0 (complex-cartesian 0.0 0.0))
              (if (number? u-val)
                (/ (double u-val) r)
                (mul u-val (complex-cartesian (/ 1.0 r) 0.0))))))
        (range (count u-vec))
        u-vec))

(defn F-lsj-r-from-bound-reduced-u
  "Radial factor on **r = i·h** for **N. Austern**-style ZR **(5.5)** via **`austern-radial-integral-I-zr-eq-5-5-from-u`**.

  **This repo:** builds **F(r) = R_{φ_f}^*(r) R_{φ_i}(r)** with **R(r) = u(r)/r** (**`radial-R-from-reduced-u`**).
  For **(p,d)** pickup, **φ_i** is normally the **bound neutron** and **φ_f** the **bound deuteron** (cluster) —
  so **F** is a **two-factor overlap**, not the nucleon radial function alone.

  **Other references (e.g. *Handbook of direct nuclear reaction for retarded theorist* §5.4, pdf p. 61–62):**
  **F_{ℓsj}(r)** is defined as the **single-particle bound nucleon** radial wavefunction (Schrödinger at negative
  energy, e.g. WS + spin–orbit), normalized with **∫ F_{ℓsj}^2(r) r^2 dr = 1**; their radial integral then has
  **∫ χ_β F_{ℓsj} χ_α dr** with **F_{ℓsj}** only under the **distorted waves**. Same physics after folding
  cluster / **D₀** conventions; the **symbol** **F** is not identical to our **R_f^* R_i** product.

  Pass **φ_i**, **φ_f** as reduced **u = rR** from **`solve-bound-state-numerov`** and **`normalize-bound-state`**.
  Quantum numbers **(ℓ,s,j)** are **implicit** in those solutions.

  **Not included:** **D₀**, distorted **χ**, or **Y_{ℓ}^{m*}(\\hat r)** — use the returned vector as **`F-vec`**
  in **`austern-radial-integrand-zr-F-Ra-Rb-r2`** / **`austern-radial-integral-I-zr-eq-5-5-from-u`**.

  **Returns:** vector **F(r_i)** for **i = 0 … n−1**, **n = min** lengths of **φ_i**, **φ_f**; entries are
  complex-ready (real bound states → real **F**)."
  [phi-i-reduced-u phi-f-reduced-u ^double h]
  (let [Ri (radial-R-from-reduced-u phi-i-reduced-u h)
        Rf (radial-R-from-reduced-u phi-f-reduced-u h)
        n (min (count Ri) (count Rf))]
    (mapv (fn [^long i]
            (mul (complex-conjugate (get Rf i))
                 (get Ri i)))
          (range n))))

(defn handbook-F-lsj-radial-from-neutron-bound-u
  "**F_{ℓsj}(r_i) = R_n(r_i)** for the **transferred nucleon**: **R = u/r** from **`radial-R-from-reduced-u`** on the
  **normalized** reduced bound **u** (**`normalize-bound-state`** ⇒ **∫|u|² dr = 1 = ∫|R|² r² dr**, §5.4 style).

  **Not** multiplied by the residual cluster radial factor (deuteron **R_d**); that physics is folded into **D₀** in
  the handbook-style DWBA."
  [neutron-u-reduced ^double h]
  (radial-R-from-reduced-u neutron-u-reduced h))

(defn- transfer-radial-R-at-r-linear
  "Linear interpolation of radial factor **R(r)** on uniform grid **r = i·h**.
  **R-vec** entries are real or complex. **r** is clamped to **[0, (n−1)h]**."
  [R-vec ^double h ^double r]
  (let [n (count R-vec)]
    (if (zero? n)
      (complex-cartesian 0.0 0.0)
      (let [r-max (* (double (dec n)) h)
            r (max 0.0 (min r r-max))
            f (/ r h)
            i0 (long (Math/floor f))
            i0 (max 0 (min i0 (dec n)))
            i1 (min (dec n) (inc i0))
            t (- f (double i0))
            v0 (get R-vec i0)
            v1 (get R-vec i1)
            c0 (if (number? v0) (complex-cartesian (double v0) 0.0) v0)
            c1 (if (number? v1) (complex-cartesian (double v1) 0.0) v1)]
        (if (< (Math/abs t) 1e-15)
          c0
          (add (mul c0 (complex-cartesian (- 1.0 t) 0.0))
               (mul c1 (complex-cartesian t 0.0))))))))

;; ============================================================================
;; ZERO-RANGE AND FINITE-RANGE INTERACTIONS
;; ============================================================================

(defn zero-range-constant
  "Get zero-range constant D₀ for specific reaction type.
   
   Parameters:
   - reaction-type: Keyword indicating reaction type
     :d-p  - (d,p) reaction (neutron transfer)
     :d-n  - (d,n) reaction (proton transfer)
     :p-d  - (p,d) reaction (neutron pickup)
     :n-d  - (n,d) reaction (proton pickup)
     :alpha-t - (α,t) reaction (proton transfer)
     :alpha-he3 - (α,³He) reaction (neutron transfer)
   
   Returns: D₀ in MeV·fm^(3/2) (typical units)
   
   Typical values:
   - (d,p): D₀ ≈ -122.4 MeV·fm^(3/2) (for neutron transfer)
   - (p,d): D₀ ≈ -122.4 MeV·fm^(3/2) (same magnitude, sign depends on convention)
   - (α,t): D₀ ≈ -400 MeV·fm^(3/2) (approximate, depends on α-t interaction)
   
   Note: The sign convention varies in literature. We use the standard convention
   where D₀ is negative for attractive interactions."
  [reaction-type]
  (case reaction-type
    :d-p -122.4    ; (d,p) neutron transfer
    :d-n -122.4    ; (d,n) proton transfer
    :p-d -122.4    ; (p,d) neutron pickup
    :n-d -122.4    ; (n,d) proton pickup
    :alpha-t -400.0  ; (α,t) proton transfer (approximate)
    :alpha-he3 -400.0 ; (α,³He) neutron transfer (approximate)
    (throw (IllegalArgumentException. 
            (format "Unknown reaction type: %s. Supported types: :d-p, :d-n, :p-d, :n-d, :alpha-t, :alpha-he3" 
                   reaction-type)))))

(defn transfer-amplitude-zero-range
  "Calculate transfer amplitude in zero-range approximation.
   
   In zero-range approximation, the transfer interaction is:
   V_transfer(r) = D₀ δ(r)
   
   This simplifies the transfer amplitude to:
   T = D₀ · ∫ χ*_f(r) φ*_f(r) φ_i(r) χ_i(r) d³r
   
   For zero-range, this reduces to evaluating the form factor at r=0:
   T = D₀ · φ*_f(0) · φ_i(0)
   
   However, for l > 0 states, φ(0) = 0, so we need the full overlap integral.
   The correct zero-range amplitude is:
   T = D₀ · ∫ φ*_f(r) φ_i(r) r² dr
   
   Parameters:
   - overlap-integral: Overlap integral O = ∫ φ*_f(r) φ_i(r) r² dr
   - D0: Zero-range constant (from zero-range-constant function)
   
   Returns: Transfer amplitude T (complex number, in general)
   
   Example:
   (require '[dwba.form-factors :as ff])
   (let [overlap (ff/overlap-integral phi-i phi-f r-max h)
         D0 (zero-range-constant :d-p)]
     (transfer-amplitude-zero-range overlap D0))"
  [overlap-integral D0]
  (* D0 overlap-integral))

(defn yukawa-form-factor
  "Yukawa form factor for finite-range interaction.
   
   The Yukawa form factor is commonly used for finite-range interactions:
   F_Y(r) = exp(-μr) / r
   
   This corresponds to the standard Yukawa potential form.
   
   Parameters:
   - r: Radial distance (fm)
   - mu: Range parameter (fm⁻¹), typically μ ≈ 0.7 fm⁻¹ for nucleon-nucleon interaction
   
   Returns: F_Y(r)
   
   Note: For r = 0, we use a small value to avoid division by zero.
   The limit as r → 0 is handled by using r = h (step size) as minimum."
  [r mu]
  (let [r-safe (if (zero? r) 1e-10 r)]  ; Avoid division by zero
    (/ (Math/exp (* (- mu) r-safe)) r-safe)))

(defn gaussian-form-factor
  "Gaussian form factor for finite-range interaction.
   
   Alternative to Yukawa form factor:
   F_G(r) = exp(-(r/β)²)
   
   Parameters:
   - r: Radial distance (fm)
   - beta: Range parameter (fm), typically β ≈ 1.0-1.5 fm
   
   Returns: F_G(r)
   
   This form is sometimes used for finite-range corrections."
  [r beta]
  (Math/exp (- (/ (* r r) (* beta beta)))))

(defn finite-range-interaction
  "Calculate finite-range interaction potential.
   
   The finite-range interaction is:
   V_transfer(r) = V₀ · F(r)
   
   where F(r) is a form factor (Yukawa or Gaussian).
   
   Parameters:
   - r: Radial distance (fm)
   - V0: Interaction strength (MeV)
   - form-factor-type: :yukawa or :gaussian
   - range-param: Range parameter
     - For Yukawa: μ (fm⁻¹), typically 0.7 fm⁻¹
     - For Gaussian: β (fm), typically 1.0-1.5 fm
   
   Returns: V_transfer(r) in MeV
   
   Example:
   (finite-range-interaction 2.0 50.0 :yukawa 0.7)  ; Yukawa at r=2 fm
   (finite-range-interaction 2.0 50.0 :gaussian 1.2)  ; Gaussian at r=2 fm"
  [r V0 form-factor-type range-param]
  (let [F-r (case form-factor-type
              :yukawa (yukawa-form-factor r range-param)
              :gaussian (gaussian-form-factor r range-param)
              (throw (IllegalArgumentException. 
                      (format "Unknown form factor type: %s. Use :yukawa or :gaussian" 
                             form-factor-type))))]
    (* V0 F-r)))

(defn finite-range-overlap-integral
  "Calculate overlap integral with finite-range form factor.
   
   This is the overlap integral weighted by the finite-range interaction:
   O_FR = ∫₀^∞ φ*_f(r) φ_i(r) F(r) r² dr
   
   where F(r) is the form factor (Yukawa or Gaussian).
   
   Parameters:
   - phi-i: Initial bound state wavefunction (vector)
   - phi-f: Final bound state wavefunction (vector)
   - r-max: Maximum radius (fm)
   - h: Step size (fm)
   - form-factor-type: :yukawa or :gaussian
   - range-param: Range parameter (μ for Yukawa, β for Gaussian)
   
   Returns: O_FR = ∫ φ*_f(r) φ_i(r) F(r) r² dr
   
   Example:
   (finite-range-overlap-integral phi-i phi-f 20.0 0.01 :yukawa 0.7)"
  [phi-i phi-f _r-max h form-factor-type range-param]
  (let [n (min (count phi-i) (count phi-f))
        integrand (mapv (fn [i]
                         (let [r (* i h)
                               phi-i-val (get phi-i i)
                               phi-f-val (get phi-f i)
                               phi-f-conj (if (number? phi-f-val)
                                           phi-f-val
                                           (complex-cartesian (re phi-f-val) (- (im phi-f-val))))
                               F-r (case form-factor-type
                                     :yukawa (yukawa-form-factor r range-param)
                                     :gaussian (gaussian-form-factor r range-param)
                                     (throw (IllegalArgumentException. 
                                             (format "Unknown form factor type: %s" form-factor-type))))]
                           (* phi-f-conj phi-i-val F-r r r)))
                       (range n))
        ;; Simpson's rule integration
        simpson-sum (loop [i 1 sum 0.0]
                     (if (>= i (dec n))
                       sum
                       (let [coeff (if (odd? i) 4.0 2.0)
                             term (* coeff (get integrand i))]
                         (recur (inc i) (+ sum term)))))
        integral (* (/ h 3.0)
                   (+ (first integrand)
                      (last integrand)
                      simpson-sum))]
    integral))

(defn transfer-amplitude-finite-range
  "Calculate transfer amplitude with finite-range interaction.
   
   Parameters:
   - finite-range-overlap: Overlap integral with form factor (from finite-range-overlap-integral)
   - V0: Interaction strength (MeV)
   
   Returns: Transfer amplitude T = V₀ · O_FR
   
   Example:
   (let [overlap-fr (finite-range-overlap-integral phi-i phi-f 20.0 0.01 :yukawa 0.7)]
     (transfer-amplitude-finite-range overlap-fr 50.0))"
  [finite-range-overlap V0]
  (* V0 finite-range-overlap))

(defn transfer-amplitude-post
  "Calculate transfer amplitude in POST formulation.
   
   In the post formulation, the interaction is in the exit channel:
   T_post = <χ_f|V_transfer|χ_i>
   
   where:
   - χ_i: Distorted wave in entrance channel
   - χ_f: Distorted wave in exit channel
   - V_transfer: Transfer interaction
   
   For zero-range: T_post = D₀ · ∫ χ*_f(r) φ*_f(r) φ_i(r) χ_i(r) d³r
   
   For finite-range: T_post = ∫ χ*_f(r) V_transfer(r) φ*_f(r) φ_i(r) χ_i(r) d³r

   **Austern (1970), Eq. (5.3) — ZR reduced amplitude volume element:** in zero-range DWBA the
   double coordinate integral **(4.60)** collapses to **∫ d³r** in the **entrance** relative coordinate
   **r**, with the **exit** distorted wave evaluated at **(M_A/M_B) r** (Austern’s mass ratio;
   **M_B = M_A + M_x**). Pass **`:zr-chi-exit-mass-ratio`** `(M_A/M_B)` in the optional **opts** map
   (9th argument). Default **1.0** recovers the legacy single-radius sampling **χ_f(r)**. See
   `austern-zr-chi-exit-mass-ratio`.
   
   Parameters:
   - chi-i, chi-f, phi-i, phi-f: **Reduced** radial functions u(r)=r·R(r) on grid r=i·h, as
     produced by `distorted-wave-optical`, `solve-bound-state-numerov` + `normalize-bound-state`.
     **φ_i**, **φ_f** are **bound** states (**negative** eigenenergy). The ZR integrand uses **handbook**
     single-particle **F_{ℓsj}(r)=R_n(r)=u_n/r** (**`handbook-F-lsj-radial-from-neutron-bound-u`**), i.e.
     **∫|R_n|² r² dr = 1** when **u_n** is **`normalize-bound-state`** — not Austern’s **R_{φ_f}^* R_{φ_i}**
     product (**`F-lsj-r-from-bound-reduced-u`**, kept for low-level benchmarks only).
   - r-max: Maximum radius (fm)
   - h: Step size (fm)
   - interaction-type: :zero-range or :finite-range
   - interaction-params: Parameters for interaction
     - For zero-range: {:D0 value} or just D0 value
     - For finite-range: {:V0 value, :form-factor :yukawa/:gaussian, :range-param value}
   - opts (optional): **`:zr-chi-exit-mass-ratio`** for **:zero-range**; **`:handbook-F-from`** **`:phi-i`**
     (default, e.g. nucleus nucleon in **(p,d)**) or **`:phi-f`** (e.g. captured nucleon in **(d,p)**).
   
   Returns: Transfer amplitude T_post
   
   Note: Default distorted waves are max-normalized (`distorted-wave-optical` **`:normalize-mode :max`**); optional
   **`:coulomb-tail`** (see **`distorted-wave-optical`** doc) gives better partial-wave ratios for Austern **(5.6)**; absolute
   scale vs. unit-flux DWUCK may need asymptotic matching — but u→R fixes the dominant r-measure bug."
  ([chi-i chi-f phi-i phi-f r-max h interaction-type interaction-params]
   (transfer-amplitude-post chi-i chi-f phi-i phi-f r-max h interaction-type interaction-params {}))
  ([chi-i chi-f phi-i phi-f _r-max h interaction-type interaction-params opts]
  (let [ratio (double (or (:zr-chi-exit-mass-ratio opts) 1.0))
        ratio-1? (< (Math/abs (- ratio 1.0)) 1e-14)
        chi-i-R (radial-R-from-reduced-u chi-i h)
        chi-f-R (radial-R-from-reduced-u chi-f h)
        nucleon-u (case (or (:handbook-F-from opts) :phi-i)
                     :phi-i phi-i
                     :phi-f phi-f
                     (throw (IllegalArgumentException.
                             (format "transfer-amplitude-post: :handbook-F-from must be :phi-i or :phi-f, got %s"
                                     (pr-str (:handbook-F-from opts))))))
        F-lsj (handbook-F-lsj-radial-from-neutron-bound-u nucleon-u h)
        n (min (count chi-i-R) (count chi-f-R) (count F-lsj))
        integrand (mapv (fn [i]
                         (let [r (* i h)
                               chi-i-val (get chi-i-R i)
                               chi-f-val (if (and ratio-1? (= :zero-range interaction-type))
                                           (get chi-f-R i)
                                           (if (= :zero-range interaction-type)
                                             (transfer-radial-R-at-r-linear chi-f-R h (* ratio r))
                                             (get chi-f-R i)))
                               F-val (get F-lsj i)
                               ;; Compute conjugates (complex-conjugate works for both real and complex)
                               chi-f-conj (complex-conjugate chi-f-val)
                               ;; Transfer interaction
                               V-transfer (case interaction-type
                                            :zero-range (let [D0 (if (map? interaction-params)
                                                                  (:D0 interaction-params)
                                                                  interaction-params)]
                                                          ;; Zero-range: V = D₀ δ(r), so at r we use D₀
                                                          ;; In practice, we integrate the overlap separately
                                                          D0)
                                            :finite-range (let [params (if (map? interaction-params)
                                                                        interaction-params
                                                                        {:V0 interaction-params})]
                                                           (finite-range-interaction r 
                                                                                     (:V0 params)
                                                                                     (:form-factor params :yukawa)
                                                                                     (:range-param params 0.7)))
                                            (throw (IllegalArgumentException. 
                                                    (format "Unknown interaction type: %s" interaction-type))))
                               r-squared (* r r)]
                           ;; Use complex multiplication: χ*_f · V · F_{ℓsj} · χ_i · r² (**F** = handbook **R_n**)
                           ;; Convert real numbers to complex for safe multiplication
                           (let [V-complex (if (number? V-transfer)
                                           (complex-cartesian V-transfer 0.0)
                                           V-transfer)
                                 r-sq-complex (complex-cartesian r-squared 0.0)]
                             (mul chi-f-conj V-complex F-val chi-i-val r-sq-complex))))
                       (range n)) ;integrand
        ;; Simpson's rule integration with complex numbers
        simpson-sum (loop [i 1 sum (complex-cartesian 0.0 0.0)]
                     (if (>= i (dec n))
                       sum
                       (let [coeff (if (odd? i) 4.0 2.0)
                             term-val (get integrand i)
                             coeff-complex (complex-cartesian coeff 0.0)
                             term (mul coeff-complex term-val)]
                         (recur (inc i) (add sum term)))))
        h-over-3 (/ h 3.0)
        first-term (get integrand 0)
        last-term (get integrand (dec n))
        h-over-3-complex (complex-cartesian h-over-3 0.0)
        integral (mul h-over-3-complex
                     (add first-term last-term simpson-sum))]
    ;; Distorted u(r) are max-normalized, not unit asymptotic flux — remaining scale offset vs
    ;; DWUCK is typically O(1–10²) unless we re-match χ to Coulomb/nuclear asymptotics.
    integral)))

(defn transfer-amplitude-prior
  "Calculate transfer amplitude in PRIOR formulation.
   
   In the prior formulation, the interaction is in the entrance channel:
   T_prior = <χ_f|V_transfer|χ_i>
   
   The mathematical form is the same as post, but the physical interpretation
   is different. For exact calculations, T_post = T_prior (post-prior equivalence).
   
   Parameters: Same as `transfer-amplitude-post` (optional **opts** for **Austern (5.3)**).
   
   Returns: Transfer amplitude T_prior
   
   Note: For testing, you can verify that T_post ≈ T_prior numerically."
  ([chi-i chi-f phi-i phi-f r-max h interaction-type interaction-params]
   (transfer-amplitude-prior chi-i chi-f phi-i phi-f r-max h interaction-type interaction-params {}))
  ([chi-i chi-f phi-i phi-f r-max h interaction-type interaction-params opts]
   ;; Prior form is mathematically identical to post form
   (transfer-amplitude-post chi-i chi-f phi-i phi-f r-max h interaction-type interaction-params opts)))

;; ============================================================================
;; LOCAL ENERGY APPROXIMATION (LEA) WITH HULTHEN POTENTIAL
;; ============================================================================

(defn hulthen-potential
  "Hulthen potential for deuteron model.
   
   The Hulthen potential is:
   V(r) = -V₀ · (e^(-αr) / (1 - e^(-βr)))
   
   For the standard Hulthen form (β = α):
   V(r) = -V₀ · (e^(-αr) / (1 - e^(-αr)))
   
   Parameters:
   - r: Radial distance (fm)
   - V0: Potential depth (MeV), typically 50-70 MeV for deuteron
   - alpha: Range parameter (fm⁻¹), typically α ≈ 0.23 fm⁻¹ for deuteron
   - beta: Optional second range parameter (fm⁻¹). If nil, uses β = α (standard form)
   
   Returns: V(r) in MeV
   
   Typical deuteron parameters:
   - V₀ ≈ 50-70 MeV
   - α ≈ 0.23 fm⁻¹
   - Binding energy: E_d ≈ -2.225 MeV
   
   Example:
   (hulthen-potential 1.0 60.0 0.23)  ; Standard form at r=1 fm"
  ([r V0 alpha]
   (hulthen-potential r V0 alpha alpha))
  ([r V0 alpha beta]
   (let [r-safe (if (zero? r) 1e-10 r)
         exp-alpha-r (Math/exp (* (- alpha) r-safe))
         exp-beta-r (Math/exp (* (- beta) r-safe))
         denominator (- 1.0 exp-beta-r)]
     (if (< (Math/abs denominator) 1e-10)
       (* -1.0 V0)  ; Limit as r → 0: V(0) = -V₀
       (* -1.0 V0 (/ exp-alpha-r denominator))))))

(defn hulthen-wavefunction
  "Hulthen wavefunction for deuteron (analytical form).
   
   The Hulthen wavefunction for the deuteron (l=0, S-state) is:
   u(r) = N · (e^(-αr) - e^(-βr)) / (1 - e^(-βr))
   
   For the standard form (β = α), this simplifies, but we use the general form.
   
   Parameters:
   - r: Radial distance (fm)
   - alpha: Range parameter (fm⁻¹), typically α ≈ 0.23 fm⁻¹
   - beta: Second range parameter (fm⁻¹), typically β ≈ 1.0-1.5 fm⁻¹
   - normalization: Normalization constant N (optional, calculated if nil)
   
   Returns: u(r) - radial wavefunction
   
   Note: This is the analytical form. For l > 0, numerical solution is needed."
  ([r alpha beta]
   (hulthen-wavefunction r alpha beta nil))
  ([r alpha beta normalization]
   (let [r-safe (if (zero? r) 1e-10 r)
         exp-alpha-r (Math/exp (* (- alpha) r-safe))
         exp-beta-r (Math/exp (* (- beta) r-safe))
         numerator (- exp-alpha-r exp-beta-r)
         denominator (- 1.0 exp-beta-r)]
     (if (< (Math/abs denominator) 1e-10)
       0.0  ; u(0) = 0
       (let [u-unnormalized (/ numerator denominator)
             N (or normalization 1.0)]  ; Default normalization, should be calculated
         (* N u-unnormalized))))))

(defn hulthen-wavefunction-normalized
  "Calculate normalized Hulthen wavefunction.
   
   First calculates the normalization constant, then returns normalized wavefunction.
   
   Parameters:
   - r: Radial distance (fm)
   - alpha: Range parameter (fm⁻¹)
   - beta: Second range parameter (fm⁻¹)
   - r-max: Maximum radius for normalization (fm)
   - h: Step size (fm)
   
   Returns: Normalized wavefunction value at r
   
   The normalization is: ∫₀^∞ |u(r)|² dr = 1"
  [r alpha beta r-max h]
  (let [n (int (/ r-max h))
        ;; Calculate unnormalized wavefunction at all points
        u-unnorm (mapv (fn [i]
                        (let [r-val (* i h)]
                          (hulthen-wavefunction r-val alpha beta 1.0)))
                      (range n))
        ;; Calculate norm squared
        integrand (mapv #(* % %) u-unnorm)
        simpson-sum (loop [i 1 sum 0.0]
                     (if (>= i (dec n))
                       sum
                       (let [coeff (if (odd? i) 4.0 2.0)
                             term (* coeff (get integrand i))]
                         (recur (inc i) (+ sum term)))))
        norm-squared (* (/ h 3.0)
                       (+ (first integrand)
                          (last integrand)
                          simpson-sum))
        normalization (Math/sqrt norm-squared)
        idx (int (/ r h))
        idx-safe (min idx (dec n))]
    (if (and (>= idx-safe 0) (< idx-safe n))
      (/ (get u-unnorm idx-safe) normalization)
      0.0)))

(defn lea-transfer-amplitude
  "Calculate transfer amplitude using Local Energy Approximation (LEA) with Hulthen potential.
   
   LEA approximates the transfer amplitude by evaluating the overlap at the local energy
   of the transferred particle. For (d,p) reactions, the deuteron is described by a
   Hulthen wavefunction.
   
   The LEA amplitude is:
   T_LEA = D₀ · ∫ φ*_f(r) φ_Hulthen(r) r² dr
   
   where φ_Hulthen is the Hulthen wavefunction for the deuteron.
   
   Parameters:
   - phi-f: Final bound state wavefunction (vector) - nucleon bound in final nucleus
   - alpha: Hulthen range parameter α (fm⁻¹), typically 0.23 fm⁻¹ for deuteron
   - beta: Hulthen range parameter β (fm⁻¹), typically 1.0-1.5 fm⁻¹
   - r-max: Maximum radius (fm)
   - h: Step size (fm)
   - D0: Zero-range constant (MeV·fm^(3/2)), typically from zero-range-constant function
   
   Returns: Transfer amplitude T_LEA
   
   Example:
   (let [phi-f (:normalized-wavefunction (solve-bound-state ...))
         D0 (zero-range-constant :d-p)]
     (lea-transfer-amplitude phi-f 0.23 1.4 20.0 0.01 D0))"
  [phi-f alpha beta r-max h D0]
  (let [n (min (count phi-f) (int (/ r-max h)))
        integrand (mapv (fn [i]
                         (let [r (* i h)
                               phi-f-val (get phi-f i)
                               phi-f-conj (if (number? phi-f-val)
                                           phi-f-val
                                           (complex-cartesian (re phi-f-val) (- (im phi-f-val))))
                               phi-hulthen (hulthen-wavefunction-normalized r alpha beta r-max h)]
                           (* phi-f-conj phi-hulthen r r)))
                       (range n))
        simpson-sum (loop [i 1 sum 0.0]
                     (if (>= i (dec n))
                       sum
                       (let [coeff (if (odd? i) 4.0 2.0)
                             term (* coeff (get integrand i))]
                         (recur (inc i) (+ sum term)))))
        overlap (* (/ h 3.0)
                  (+ (first integrand)
                     (last integrand)
                     simpson-sum))]
    (* D0 overlap)))

(defn lea-transfer-amplitude-simplified
  "Simplified LEA transfer amplitude using standard Hulthen parameters for deuteron.
   
   Uses typical deuteron parameters:
   - α = 0.23 fm⁻¹
   - β = 1.4 fm⁻¹
   
   Parameters:
   - phi-f: Final bound state wavefunction (vector)
   - r-max: Maximum radius (fm)
   - h: Step size (fm)
   - reaction-type: Reaction type keyword (:d-p, :p-d, etc.) for D₀
   
   Returns: Transfer amplitude T_LEA
   
   Example:
   (let [phi-f (:normalized-wavefunction (solve-bound-state ...))]
     (lea-transfer-amplitude-simplified phi-f 20.0 0.01 :d-p))"
  [phi-f r-max h reaction-type]
  (let [alpha 0.23  ; Standard deuteron parameter
        beta 1.4    ; Standard deuteron parameter
        D0 (zero-range-constant reaction-type)]
    (lea-transfer-amplitude phi-f alpha beta r-max h D0)))

;; ============================================================================
;; PHASE 5: ANGULAR MOMENTUM COUPLING
;; ============================================================================

(defn transfer-L-selection-weight
  "Selection weight for angular momentum transfer L in single-nucleon transfer.
   
   For transfer from initial orbital l_i to final orbital l_f:
   - Triangle rule: |l_i - l_f| ≤ L ≤ l_i + l_f
   - Parity: (-1)^L = (-1)^(l_i + l_f)
   
   Returns 1.0 if L is allowed, 0.0 if forbidden. Use when building dσ/dΩ
   so that only allowed L contribute (some L dominate, others suppressed).
   
   Example: l_i=1, l_f=0 → only L=1 allowed (L=0 parity-forbidden, L>1 triangle-forbidden)."
  [L l-i l-f]
  (let [l-i (long l-i)
        l-f (long l-f)
        L (long L)
        triangle-ok (and (>= L (Math/abs (- l-i l-f)))
                        (<= L (+ l-i l-f)))
        parity-l (+ l-i l-f)  ; (-1)^(l_i + l_f)
        parity-L L            ; (-1)^L
        parity-ok (= (mod parity-l 2) (mod parity-L 2))]
    (if (and triangle-ok parity-ok) 1.0 0.0)))

(defn clebsch-gordan
  "Calculate Clebsch-Gordan coefficient <j1 m1 j2 m2 | J M>.
   
   Clebsch-Gordan coefficients couple two angular momenta:
   |J M> = Σ_{m1,m2} <j1 m1 j2 m2 | J M> |j1 m1> |j2 m2>
   
   Parameters:
   - j1: First angular momentum (half-integer, e.g., 1/2, 1, 3/2, ...)
   - m1: First magnetic quantum number (-j1 ≤ m1 ≤ j1)
   - j2: Second angular momentum
   - m2: Second magnetic quantum number (-j2 ≤ m2 ≤ j2)
   - J: Total angular momentum
   - M: Total magnetic quantum number (M = m1 + m2)
   
   Returns: Clebsch-Gordan coefficient (real number)
   
   Note: This implementation uses selection rules and simplified formulas.
   For production use, consider using a specialized library for more accurate calculations.
   
   Example:
   (clebsch-gordan 1 0 1 0 2 0)  ; <1 0 1 0 | 2 0>"
  [j1 m1 j2 m2 J M]
  ;; Check selection rules
  (if (not= M (+ m1 m2))
    0.0  ; M must equal m1 + m2
    (if (or (< J (Math/abs (- j1 j2)))
            (> J (+ j1 j2)))
      0.0  ; Triangle inequality: |j1 - j2| ≤ J ≤ j1 + j2
      ;; Simplified calculation using Wigner 3-j symbol relation
      ;; CG = (-1)^(j1-j2+M) * sqrt(2J+1) * Wigner3j(j1 j2 J; m1 m2 -M)
      ;; For now, use a simplified approximation
      (let [sign (if (even? (int (- j1 j2 M))) 1.0 -1.0)
            sqrt-factor (Math/sqrt (inc (* 2 J)))
            ;; Simplified: use approximate formula
            ;; For specific cases, we can calculate exactly
            approx-value (cond
                          ;; Special case: J = j1 + j2, M = m1 + m2
                          (= J (+ j1 j2))
                          (* sign sqrt-factor 1.0)
                          ;; Special case: J = |j1 - j2|, M = m1 + m2
                          (= J (Math/abs (- j1 j2)))
                          (* sign sqrt-factor 1.0)
                          ;; General case: use approximate formula
                          :else
                          (let [j1-int (int j1)
                                j2-int (int j2)
                                J-int (int J)
                                fact-sum (m/factorial (+ j1-int j2-int J-int))
                                fact-diff1 (m/factorial (Math/abs (- j1-int j2-int J-int)))
                                fact-diff2 (m/factorial (Math/abs (- j2-int j1-int J-int)))]
                            (* sign sqrt-factor 
                               (Math/sqrt (/ (* (inc (* 2 j1)) (inc (* 2 j2)))
                                            (* (inc (* 2 J)) fact-sum fact-diff1 fact-diff2))))))]
        approx-value))))

(defn wigner-3j
  "Calculate Wigner 3-j symbol (j1 j2 j3; m1 m2 m3).
   
   The Wigner 3-j symbol is related to Clebsch-Gordan coefficients:
   (j1 j2 j3; m1 m2 m3) = (-1)^(j1-j2-m3) / sqrt(2j3+1) * <j1 m1 j2 m2 | j3 -m3>
   
   Parameters:
   - j1, j2, j3: Angular momenta (half-integers)
   - m1, m2, m3: Magnetic quantum numbers
   
   Returns: Wigner 3-j symbol value
   
   Note: This is a simplified implementation. For accurate calculations,
   use specialized libraries or more sophisticated algorithms."
  [j1 j2 j3 m1 m2 m3]
  (if (not= (+ m1 m2 m3) 0)
    0.0  ; Selection rule: m1 + m2 + m3 = 0
    (if (or (< j3 (Math/abs (- j1 j2)))
            (> j3 (+ j1 j2)))
      0.0  ; Triangle inequality
      (let [sign (if (even? (int (- j1 j2 m3))) 1.0 -1.0)
            sqrt-factor (/ 1.0 (Math/sqrt (inc (* 2 j3))))
            cg (clebsch-gordan j1 m1 j2 m2 j3 (- m3))]
        (* sign sqrt-factor cg)))))

(defn racah-coefficient
  "Calculate Racah coefficient W(j1 j2 J j3; J12 J23).
   
   Racah coefficients are used for recoupling three angular momenta.
   They are related to 6-j symbols:
   W(j1 j2 J j3; J12 J23) = (-1)^(j1+j2+j3+J) * {j1 j2 J12; j3 J J23}
   
   Parameters:
   - j1, j2, j3, J: Angular momenta
   - J12: Coupled angular momentum of j1 and j2
   - J23: Coupled angular momentum of j2 and j3
   
   Returns: Racah coefficient value
   
   Note: This is a simplified implementation. For accurate calculations,
   use specialized libraries or more sophisticated algorithms."
  [j1 j2 j3 J J12 J23]
  ;; Check triangle inequalities
  (if (or (< J12 (Math/abs (- j1 j2)))
          (> J12 (+ j1 j2))
          (< J23 (Math/abs (- j2 j3)))
          (> J23 (+ j2 j3))
          (< J (Math/abs (- J12 j3)))
          (> J (+ J12 j3))
          (< J (Math/abs (- j1 J23)))
          (> J (+ j1 J23)))
    0.0
    ;; Simplified calculation using sum over magnetic quantum numbers
    ;; W = Σ_m1 m2 m3 M12 M23 (-1)^(j1+j2+j3+J) * CG(j1 m1 j2 m2 | J12 M12)
    ;;     * CG(J12 M12 j3 m3 | J M) * CG(j1 m1 J23 M23 | J M) * CG(j2 m2 j3 m3 | J23 M23)
    ;; This is computationally expensive, so we use an approximation
    (let [sign (if (even? (int (+ j1 j2 j3 J))) 1.0 -1.0)
          ;; Simplified approximation based on symmetry properties
          approx-value (* sign
                         (Math/sqrt (/ (* (inc (* 2 J12)) (inc (* 2 J23)))
                                      (* (inc (* 2 j1)) (inc (* 2 j2)) (inc (* 2 j3)) (inc (* 2 J))))))]
      approx-value)))

(defn- assoc-legendre-double-factorial-odd
  ^double [^long n]
  (if (<= n 1) 1.0 (* (double n) (assoc-legendre-double-factorial-odd (- n 2)))))

(defn- assoc-legendre-P-mm
  "**P_m^m(x)**, Ferrers / Condon–Shortley (**m ≥ 0**): **(−1)^m (2m−1)!! (1−x²)^{m/2}**."
  ^double [^long m ^double x]
  (if (zero? m)
    1.0
    (let [omx2 (max 0.0 (- 1.0 (* x x)))
          sm (Math/sqrt omx2)
          df (assoc-legendre-double-factorial-odd (dec (* 2 m)))
          phase (if (even? m) 1.0 -1.0)]
      (* phase df (Math/pow sm m)))))

(defn- assoc-legendre-P-lm
  "Associated Legendre **P_l^m(x)**, Condon–Shortley, integers **0 ≤ m ≤ l**, **x ∈ [−1,1]**.
  Stable upward recurrence in **l** (Abramowitz & Stegun style)."
  ^double [^long l ^long m ^double x]
  (cond
    (or (< m 0) (> m l)) 0.0
    (= l m) (assoc-legendre-P-mm m x)
    (= l (inc m)) (* x (double (+ (* 2 m) 1)) (assoc-legendre-P-mm m x))
    :else
    (loop [k (+ m 2)
           p0 (assoc-legendre-P-mm m x)
           p1 (* x (double (+ (* 2 m) 1)) (assoc-legendre-P-mm m x))]
      (let [p2 (/ (- (* (double (dec (* 2 k))) x p1)
                   (* (double (+ (+ k m) -1)) p0))  ;; (k+m-1)
                (double (- k m)))]
        (if (= k l)
          p2
          (recur (inc k) p1 p2))))))

(defn spherical-harmonic
  "Spherical harmonic **Y_l^m(θ, φ)** (Condon–Shortley):

  **Y_l^m = (−1)^m √[(2l+1)/(4π) · (l−m)!/(l+m)!] · P_l^m(cos θ) · e^{imφ}**,

  with **P_l^m** from **`assoc-legendre-P-lm`** (replaces the old incorrect **(sin θ)^|m| P_l** patch).

  Parameters: **l ≥ 0**, **|m| ≤ l**, **θ**, **φ** (radians). Returns a complex number.

  Example: `(spherical-harmonic 1 0 (/ Math/PI 2) 0)`."
  [l m theta phi]
  (let [cos-theta (Math/cos theta)
        m-i (int m)
        ma (long (Math/abs m-i))
        ;; Normalization √( (2l+1)/(4π) · (l−|m|)! / (l+|m|)! )
        norm-factor (if (and (zero? l) (zero? ma))
                     (/ 1.0 (Math/sqrt (* 4.0 Math/PI)))
                     (Math/sqrt (/
                                  (* (inc (* 2 (long l)))
                                     (m/factorial (- (long l) ma)))
                                  (* 4.0 Math/PI (m/factorial (+ (long l) ma))))))
        leg-value (if (and (zero? l) (zero? ma))
                   1.0
                   (assoc-legendre-P-lm (long l) ma cos-theta))
        ;; Condon–Shortley (−1)^m on the full harmonic (signed m)
        condon-phase (if (even? m-i) 1.0 -1.0)
        exp-factor (if (or (zero? m-i) (zero? phi))
                    (complex-cartesian 1.0 0.0)
                    (complex-polar (* m-i phi) 1.0))
        pref (* condon-phase norm-factor (double leg-value))
        result (mul (complex-cartesian pref 0.0) exp-factor)]
    result))

(defn transfer-angular-distribution
  "Calculate angular distribution for transfer reactions.
   
   The angular distribution for transfer reactions depends on the angular
   momentum coupling between the initial and final states.
   
   For a transfer reaction with angular momentum transfer L:
   dσ/dΩ(θ) ∝ Σ_L |w_L·T_L|² · |Y_L0(θ, 0)|² (**incoherent** multipoles: no **L–L′** cross terms).
   This is **not** isotropic when **several** **L** contribute — each term varies with **θ** via **|Y_L0|² ∝ P_L(cos θ)²**;
   only the **interference** pattern differs from **`transfer-angular-distribution-coherent`**
   (**|Σ_L w_L T_L Y_L0|²**). When l-i and l-f are given, w_L is the selection weight (1 if L allowed, 0 if forbidden
   by triangle/parity), so some L dominate, others suppressed.
   
   Parameters:
   - T-amplitudes: Map of {L → T_L} transfer amplitudes for each angular momentum
   - theta: Scattering angle (radians)
   - phi: Azimuthal angle (radians, default 0 for coplanar scattering)
   - l-i: Optional. Initial bound-state orbital (e.g. 1 for p-wave). Enables L selection.
   - l-f: Optional. Final bound-state orbital (e.g. 0 for s-wave). Enables L selection.
   
   Returns: Angular distribution value
   
   Example:
   (let [T-map {0 1.0, 1 0.5, 2 0.2}]
     (transfer-angular-distribution T-map (/ Math/PI 2) 0))
   ;; With L selection (l_i=1, l_f=0 → only L=1 allowed):
   (transfer-angular-distribution T-map (/ Math/PI 2) 0.0 1 0)"
  ([T-amplitudes theta]
   (transfer-angular-distribution T-amplitudes theta 0.0 nil nil))
  ([T-amplitudes theta phi]
   (transfer-angular-distribution T-amplitudes theta phi nil nil))
  ([T-amplitudes theta phi l-i l-f]
   (let [weight-fn (if (and (number? l-i) (number? l-f))
                     (fn [L] (transfer-L-selection-weight L l-i l-f))
                     (fn [_L] 1.0))
         sum-over-L (reduce + (map (fn [[L T-L]]
                                   (let [w (weight-fn L)
                                         Y-L0 (spherical-harmonic L 0 theta phi)
                                         Y-L0-mag-squared (+ (* (re Y-L0) (re Y-L0))
                                                            (* (im Y-L0) (im Y-L0)))
                                         T-L-mag-squared (if (number? T-L)
                                                          (* T-L T-L)
                                                          (+ (* (re T-L) (re T-L))
                                                            (* (im T-L) (im T-L))))]
                                     (* w w T-L-mag-squared Y-L0-mag-squared)))
                                 T-amplitudes))]
     sum-over-L)))

(defn- transfer-cg-orbitals-exact
  "Orbital Clebsch–Gordan **<l_i m_i L M | l_f m_f>** (integer **l**), exact via **`wigner-3j`**:

  **<j1 m1 j2 m2 | j m>** = (−1)^{j1−j2+m} √(2j+1) **(j1 j2 j; m1 m2 −m)**."
  [l-i m-i L M l-f m-f]
  (let [l-i (long l-i) m-i (long m-i) L (long L) M (long M) l-f (long l-f) m-f (long m-f)]
    (if (not= (+ m-i M) m-f)
      0.0
      (let [phase (if (even? (+ l-i (- L) m-f)) 1.0 -1.0)
            pref (* phase (Math/sqrt (double (inc (* 2 l-f)))))
            w3 (jam/wigner-3j (double l-i) (double L) (double l-f)
                 (double m-i) (double M) (double (- m-f)))]
        (* pref w3)))))

(defn transfer-orbital-cg-bilinear-sum-mi
  "**Scalar recoupling factor** **S(L,L′) = Σ_{m_i} <l_i m_i L M | l_f m_i+M> <l_i m_i L′ M | l_f m_i+M>**
  with integer orbitals **l_i, l_f, L, L′** and **exact** orbital CGs (`transfer-cg-orbitals-exact` / **3j**).

  For allowed triangles this sum is **independent of M** (standard angular-momentum recoupling); it is
  evaluated at **M = 0** for efficiency. Equivalent to a **6j**-type contraction of four **3j** symbols
  (Varshalovich *et al.*; DLMF §34.3 orthogonality sums).

  Returns a **real double** (0 when no **m_i** contributes)."
  [l-i L L-prime l-f]
  (let [l-i (long l-i) L (long L) Lp (long L-prime) l-f (long l-f)
        M 0]
    (double
     (reduce
      (fn [^double acc ^long m-i]
        (let [m-f (+ m-i M)]
          (if (> (Math/abs m-f) l-f)
            acc
            (let [cg1 (transfer-cg-orbitals-exact l-i m-i L M l-f m-f)
                  cg2 (transfer-cg-orbitals-exact l-i m-i Lp M l-f m-f)]
              (+ acc (* (double cg1) (double cg2)))))))
      0.0
      (range (- l-i) (inc l-i))))))

(defn- transfer-Ylm-bilinear-sum-over-M
  "**χ_{L,L′}(θ,φ) = Σ_M Y_{LM}(θ,φ) Y*_{L′M}(θ,φ)** over **M** with **|M| ≤ min(L,L′)** (complex)."
  [^long L ^long L-prime theta phi]
  (let [mmax (min L L-prime)]
    (reduce
     (fn [acc ^long M]
       (let [y1 (spherical-harmonic L M theta phi)
             y2 (spherical-harmonic L-prime M theta phi)]
         (add acc (mul y1 (complex-conjugate y2)))))
     (complex-cartesian 0.0 0.0)
     (range (- mmax) (inc mmax)))))

(defn- transfer-T->complex [T-L]
  (if (number? T-L)
    (complex-cartesian (double T-L) 0.0)
    T-L))

(defn- transfer-angular-distribution-m-sum-unpolarized-legacy-mi-mf
  "Reference implementation: double loop over **m_i**, **m_f** (for regression tests only)."
  [T-amplitudes theta phi l-i l-f]
  (let [theta (double theta)
        phi (double phi)
        weight-fn (fn [^long L] (transfer-L-selection-weight L l-i l-f))
        t-entries (sort-by first (or (seq T-amplitudes) []))
        norm-i (double (inc (* 2 l-i)))
        mis (range (- l-i) (inc l-i))
        mfs (range (- l-f) (inc l-f))
        amp-for-pair
        (fn [^long m-i ^long m-f]
          (let [M (- m-f m-i)]
            (reduce
             (fn [acc [L T-L]]
               (let [L (long L)
                     w (double (weight-fn L))]
                 (if (< (Math/abs w) 1e-15)
                   acc
                   (if (> (Math/abs M) L)
                     acc
                     (let [cg (transfer-cg-orbitals-exact l-i m-i L M l-f m-f)]
                       (if (< (Math/abs cg) 1e-15)
                         acc
                         (let [Ylm (spherical-harmonic L M theta phi)
                               T-c (transfer-T->complex T-L)
                               wc (complex-cartesian (* w cg) 0.0)
                               term (mul T-c (mul wc Ylm))]
                           (add acc term))))))))
             (complex-cartesian 0.0 0.0)
             t-entries)))
        mag2 (fn [z] (+ (* (re z) (re z)) (* (im z) (im z))))
        sum-m (reduce
               (fn [^double s [^long mi ^long mf]]
                 (+ s (mag2 (amp-for-pair mi mf))))
               0.0
               (for [mi mis mf mfs] [mi mf]))]
    (max 0.0 (/ sum-m norm-i))))

(defn transfer-angular-distribution-m-sum-unpolarized
  "Unpolarized **orbital** angular factor: **(1/(2l_i+1)) Σ_{m_i,m_f} |A(m_i,m_f)|²** with

  **A(m_i,m_f) = Σ_L w_L · T_L · <l_i m_i L (m_f−m_i) | l_f m_f> · Y_{L,(m_f−m_i)}(θ,φ)**.

  **w_L** = `transfer-L-selection-weight`. CGs use **`dwba.angular-momentum` / Wigner 3j**.

  **Legendre vs m-sum:** For each **l**, **Y_{l0}(θ) ∝ P_l(cos θ)** (φ=0, real). The **unpolarized**
  sum that enters here is **(1/(2l_i+1)) Σ_{m_i,m_f} |A|²**; for **l_i=0** it collapses to a factor times
  **Σ_m |Y_{l_f m}(θ)|²** at a **single** θ. By the **addition theorem**, **Σ_m |Y_{l_f m}(θ,φ)|² =
  (2l_f+1)/(4π)** for **all** θ — that is **not** **P_{l_f}(cos θ)**; the **P_l** dependence cancels
  when you add **all** **m** at one point. (Contrast: **Σ_m Y*_{lm}(θ')Y_{lm}(θ) = ((2l+1)/4π) P_l(cos γ)**,
  two directions **γ** apart — that **is** **P_l**, but a different object than **Σ_m |Y_lm(θ)|²**.)
  So for **l_i=0** this path is **|w T_{l_f}|² (2l_f+1)/(4π)**, **θ-flat**. For **θ** shape from **P_L**
  in a reduced model use **`transfer-angular-distribution-coherent`** (**Y_L0** only, i.e. **|Σ w T Y_L0|² ∝ |P_L|²**).

  **l_i ≥ 1 — recoupled (same algebra as DWUCK-style m-averages):** factor
  **Σ_{m_i} |A|² = Σ_{L,L′} w_L w_{L′} T_L T_{L′}* S(L,L′) χ_{L,L′}(θ,φ)** for fixed **M = m_f−m_i**, with
  **S(L,L′) = `transfer-orbital-cg-bilinear-sum-mi`** and **χ_{L,L′} = Σ_M Y_{LM} Y*_{L′M}**
  (internal **χ** sum over **M**). No double loop over **(m_i,m_f)**; unit tests match the legacy
  **m_i**, **m_f** implementation to round-off.

  **Parameters:** `T-amplitudes` map **{L → T_L}**; **l-i**, **l-f** integer bound orbitals (required).
  **theta**, **phi** unused when **l_i=0**."
  [T-amplitudes theta phi l-i l-f]
  (let [l-i (long l-i)
        l-f (long l-f)]
    (if (zero? l-i)
      (let [w (double (transfer-L-selection-weight l-f l-i l-f))
            T-L (get T-amplitudes l-f)]
        (if (or (< (Math/abs w) 1e-15) (nil? T-L))
          0.0
          (let [T-mag-sq (if (number? T-L)
                          (let [x (double T-L)] (* x x))
                          (+ (* (re T-L) (re T-L)) (* (im T-L) (im T-L))))]
            (max 0.0 (* w w T-mag-sq (/ (double (inc (* 2 l-f))) (* 4.0 Math/PI)))))))
      (let [theta (double theta)
            phi (double phi)
            weight-fn (fn [^long L] (transfer-L-selection-weight L l-i l-f))
            t-entries (sort-by first (or (seq T-amplitudes) []))
            norm-i (double (inc (* 2 l-i)))
            zsum
            (reduce
             (fn [acc [[L T-L] [Lp T-Lp]]]
               (let [L (long L) Lp (long Lp)
                     wL (double (weight-fn L))
                     wLp (double (weight-fn Lp))]
                 (if (or (< (Math/abs wL) 1e-15) (< (Math/abs wLp) 1e-15))
                   acc
                   (let [Sbf (transfer-orbital-cg-bilinear-sum-mi l-i L Lp l-f)
                         chi (transfer-Ylm-bilinear-sum-over-M L Lp theta phi)
                         TLC (transfer-T->complex T-L)
                         TLpC (transfer-T->complex T-Lp)
                         ;; wL wLp S T_L T_L'* χ
                         pref (complex-cartesian (* wL wLp Sbf) 0.0)
                         tprod (mul TLC (complex-conjugate TLpC))
                         term (mul pref (mul tprod chi))]
                     (add acc term)))))
             (complex-cartesian 0.0 0.0)
             (for [a t-entries b t-entries] [a b]))
            q (re zsum)]
        (max 0.0 (/ (double q) norm-i))))))

(defn transfer-angular-distribution-coherent
  "Coherent DWBA-style multipole sum (coplanar, m = 0 only), then square:

  dσ/dΩ ∝ | Σ_L w_L · T_L · Y_L0(θ,φ) |²

  where w_L = `transfer-L-selection-weight` when `l-i` and `l-f` are supplied, else 1.
  For **l_i=0** the strict **`transfer-angular-distribution-m-sum-unpolarized`** factor is **θ-flat**;
  Ca40 benchmark uses **this** coherent path for a **shaped** σ(θ).

  **CM symmetry (θ ↔ π−θ):** With φ=0, **Y_{L0} ∝ P_L(cos θ)**. Parity-conserving transfer forces every
  allowed **L** to share the same parity **(−1)^L = (−1)^{l_f−l_i}**, so
  **f(θ)=Σ_L c_L P_L(cos θ)** obeys **f(π−θ)=(−1)^{l_f−l_i} f(θ)** and **|f(θ)|² = |f(π−θ)|²** even when
  **c_L** are **complex**. So this **Y_L0-only** angular factor is always **symmetric** about **θ=90°** CM for
  one-step transfer. **Incoherent** **Σ_L |w T Y_L0|²** is symmetric term-by-term. **Asymmetry** (DWUCK-style)
  needs **m≠0** physics via **`transfer-angular-distribution-m-sum-unpolarized`** (and **l_i≥1** with several
  multipoles — see its docstring); that path is **θ-flat** when **l_i=0**."
  ([T-amplitudes theta]
   (transfer-angular-distribution-coherent T-amplitudes theta 0.0 nil nil))
  ([T-amplitudes theta phi]
   (transfer-angular-distribution-coherent T-amplitudes theta phi nil nil))
  ([T-amplitudes theta phi l-i l-f]
   (let [weight-fn (if (and (number? l-i) (number? l-f))
                     (fn [L] (transfer-L-selection-weight L (long l-i) (long l-f)))
                     (fn [_L] 1.0))
         sum-amp
         (reduce
          (fn [acc [L T-L]]
            (let [w (double (weight-fn L))]
              (if (< (Math/abs w) 1e-15)
                acc
                (let [Y-L0 (spherical-harmonic (long L) 0 theta phi)
                      T-c (if (number? T-L)
                            (complex-cartesian (double T-L) 0.0)
                            T-L)
                      w-c (complex-cartesian w 0.0)
                      term (mul T-c (mul w-c Y-L0))]
                  (add acc term)))))
          (complex-cartesian 0.0 0.0)
          (sort-by first (or (seq T-amplitudes) [])))]
     (+ (* (re sum-amp) (re sum-amp))
        (* (im sum-amp) (im sum-amp))))))

(defn transfer-angular-distribution-function
  "Calculate angular distribution as a function of angle.
   
   Returns a vector of [theta, dσ/dΩ(theta)] pairs.
   
   Parameters:
   - T-amplitudes: Map of {L → T_L} transfer amplitudes
   - theta-min: Minimum angle (radians)
   - theta-max: Maximum angle (radians)
   - n-points: Number of points to calculate
   - phi: Azimuthal angle (radians, default 0)
   
   Returns: Vector of [theta, dσ/dΩ(theta)] pairs
   
   Example:
   (let [T-map {0 1.0, 1 0.5}]
     (transfer-angular-distribution-function T-map 0 Math/PI 100))"
  ([T-amplitudes theta-min theta-max n-points]
   (transfer-angular-distribution-function T-amplitudes theta-min theta-max n-points 0.0 nil nil))
  ([T-amplitudes theta-min theta-max n-points phi]
   (transfer-angular-distribution-function T-amplitudes theta-min theta-max n-points phi nil nil))
  ([T-amplitudes theta-min theta-max n-points phi l-i l-f]
   (let [d-theta (/ (- theta-max theta-min) (dec n-points))
         thetas (map #(+ theta-min (* % d-theta)) (range n-points))]
     (mapv (fn [theta]
             [theta (transfer-angular-distribution T-amplitudes theta phi l-i l-f)])
           thetas))))

(defn sum-over-magnetic-substates
  "Sum transfer amplitude over all magnetic substates.
   
   For transfer reactions, we need to sum over all possible magnetic quantum
   number projections. This function performs the sum:
   T_total = Σ_{m1,m2,m3,m4} |<j1 m1 j2 m2 | J M>|² · |T(m1,m2,m3,m4)|²
   
   Parameters:
   - T-function: Function that takes magnetic quantum numbers and returns amplitude
                 T(m1, m2, m3, m4)
   - j1, j2, j3, j4: Angular momenta for the four particles
   - J: Total angular momentum
   - M: Total magnetic quantum number
   
   Returns: Sum over all magnetic substates
   
   Note: This is a simplified version. Full implementation would require
   proper angular momentum coupling for all particles involved."
  [T-function j1 j2 j3 j4 J M]
  (let [;; Generate all possible magnetic quantum number combinations
        m1-values (range (- j1) (+ j1 1))
        m2-values (range (- j2) (+ j2 1))
        m3-values (range (- j3) (+ j3 1))
        m4-values (range (- j4) (+ j4 1))
        ;; Sum over all combinations
        total-sum (reduce + (for [m1 m1-values
                                 m2 m2-values
                                 m3 m3-values
                                 m4 m4-values
                                 :when (= M (+ m1 m2 m3 m4))]
                             (let [cg (clebsch-gordan j1 m1 j2 m2 J M)
                                   T-val (T-function m1 m2 m3 m4)
                                   T-mag-squared (if (number? T-val)
                                                  (* T-val T-val)
                                                  (+ (* (re T-val) (re T-val))
                                                    (* (im T-val) (im T-val))))]
                               (* cg cg T-mag-squared))))]
    total-sum))

;; ============================================================================
;; PHASE 6: DIFFERENTIAL CROSS-SECTION
;; ============================================================================

(defn transfer-differential-cross-section
  "Calculate differential cross-section for transfer reactions.
   
   Using the standard DWBA prefactor in terms of the mass-factors already used
   throughout this code:
   
   - mass-factor ≡ 2μ/ħ²  (in MeV⁻¹·fm⁻²)
   - Textbook DWBA form:
       dσ/dΩ = (μ_i μ_f / (2πħ²)²) · (k_f/k_i) · |T_transfer|² · S
   
   Expressed in terms of mass-factors this becomes:
   
       dσ/dΩ = (mass-factor-i · mass-factor-f)/(16π²)
               · (k_f/k_i) · |T_transfer|² · S
   
   Parameters:
   - T-amplitude: Transfer amplitude (can be complex or real number)
   - S-factor: Spectroscopic factor (dimensionless, typically 0 < S < 1)
   - k-i: Wavenumber in entrance channel (fm⁻¹)
   - k-f: Wavenumber in exit channel (fm⁻¹)
   - mass-factor-i: Entrance-channel mass factor (2μ_i/ħ²)
   - mass-factor-f: Exit-channel mass factor (2μ_f/ħ²)
   
   Optional parameters:
   - E-i: Incident energy (MeV) - for documentation/validation
   - E-f: Final energy (MeV) - for documentation/validation
   
   Returns: dσ/dΩ in **mb/sr** (1 fm² = 10 mb; internal kinematic result × 10)."
  ([T-amplitude S-factor k-i k-f mass-factor-i mass-factor-f]
   (transfer-differential-cross-section T-amplitude S-factor k-i k-f mass-factor-i mass-factor-f nil nil))
  ([T-amplitude S-factor k-i k-f mass-factor-i mass-factor-f _E-i _E-f]
   (let [;; Combined prefactor (mass-factor-i * mass-factor-f) / (16π²)
         prefactor (/ (* mass-factor-i mass-factor-f)
                      (* 16.0 Math/PI Math/PI))
         ;; Wavenumber ratio: k_f/k_i
         k-ratio (/ k-f k-i)
         ;; Amplitude squared: |T|²
         T-squared (if (number? T-amplitude)
                    (* T-amplitude T-amplitude)
                    (let [T-mag-val (mag T-amplitude)]
                      (* T-mag-val T-mag-val)))]
     (* 10.0 prefactor k-ratio T-squared S-factor))))

(defn transfer-differential-cross-section-angular
  "Calculate differential cross-section as a function of angle.
   
   This combines the transfer amplitude with angular distribution:
   dσ/dΩ(θ) = (μ_i μ_f/(2πħ²)²) · (k_f/k_i) · |T(θ)|² · S
   
   where T(θ) includes the angular momentum coupling and spherical harmonics.
   
   Parameters:
   - T-amplitudes: Map of {L → T_L} transfer amplitudes for each angular momentum
   - S-factor: Spectroscopic factor
   - k-i: Wavenumber in entrance channel (fm⁻¹)
   - k-f: Wavenumber in exit channel (fm⁻¹)
   - theta: Scattering angle in center-of-mass frame (radians)
   - mass-factor-i: Entrance channel mass factor (2μ_i/ħ²)
   - mass-factor-f: Exit channel mass factor (2μ_f/ħ²)
   - phi: Azimuthal angle (radians, default 0)
   
   Returns: dσ/dΩ(θ) in **mb/sr** (CM frame). Use transfer-cm-to-lab to convert to lab.
   
   Example:
   (let [T-map {0 1.0, 1 0.5, 2 0.2}
         k-i (Math/sqrt (* mass-factor-i 10.0))
         k-f (Math/sqrt (* mass-factor-f 8.0))
         S 0.5]
     (transfer-differential-cross-section-angular T-map S k-i k-f (/ Math/PI 2) mass-factor-i mass-factor-f))"
  ([T-amplitudes S-factor k-i k-f theta mass-factor]
   (transfer-differential-cross-section-angular T-amplitudes S-factor k-i k-f theta mass-factor mass-factor))
  ([T-amplitudes S-factor k-i k-f theta mass-factor-i mass-factor-f]
   (transfer-differential-cross-section-angular T-amplitudes S-factor k-i k-f theta mass-factor-i mass-factor-f 0.0 nil nil))
  ([T-amplitudes S-factor k-i k-f theta mass-factor-i mass-factor-f phi]
   (transfer-differential-cross-section-angular T-amplitudes S-factor k-i k-f theta mass-factor-i mass-factor-f phi nil nil))
  ([T-amplitudes S-factor k-i k-f theta mass-factor-i mass-factor-f phi l-i l-f]
   (let [;; Get angular distribution (sum over L: |w_L·T_L|² · |Y_L0|²; w_L from l-i, l-f when provided)
         angular-dist (transfer-angular-distribution T-amplitudes theta phi l-i l-f)
         ;; Combined prefactor (mass-factor-i * mass-factor-f) / (16π²)
         prefactor (/ (* mass-factor-i mass-factor-f)
                     (* 16.0 Math/PI Math/PI))
         ;; Wavenumber ratio: k_f/k_i
         k-ratio (/ k-f k-i)
         ;; Multiply all factors: (μ_i μ_f/(2πħ²)²) · (k_f/k_i) · angular_dist · S → mb/sr
         dsigma (* 10.0 prefactor k-ratio angular-dist S-factor)]
     dsigma)))

(defn transfer-differential-cross-section-angular-coherent
  "Like `transfer-differential-cross-section-angular` but uses
  `transfer-angular-distribution-coherent` (interference between multipoles)."
  ([T-amplitudes S-factor k-i k-f theta mass-factor]
   (transfer-differential-cross-section-angular-coherent T-amplitudes S-factor k-i k-f theta mass-factor mass-factor))
  ([T-amplitudes S-factor k-i k-f theta mass-factor-i mass-factor-f]
   (transfer-differential-cross-section-angular-coherent T-amplitudes S-factor k-i k-f theta mass-factor-i mass-factor-f 0.0 nil nil))
  ([T-amplitudes S-factor k-i k-f theta mass-factor-i mass-factor-f phi]
   (transfer-differential-cross-section-angular-coherent T-amplitudes S-factor k-i k-f theta mass-factor-i mass-factor-f phi nil nil))
  ([T-amplitudes S-factor k-i k-f theta mass-factor-i mass-factor-f phi l-i l-f]
   (let [angular-dist (transfer-angular-distribution-coherent T-amplitudes theta phi l-i l-f)
         prefactor (/ (* mass-factor-i mass-factor-f)
                      (* 16.0 Math/PI Math/PI))
         k-ratio (/ k-f k-i)
         dsigma (* 10.0 prefactor k-ratio angular-dist S-factor)]
     dsigma)))

(defn transfer-differential-cross-section-angular-m-sum
  "Like `transfer-differential-cross-section-angular-coherent` but uses
  `transfer-angular-distribution-m-sum-unpolarized` (exact orbital CGs; **l_i≥1** via **S(L,L′)** recoupling
  and **χ_{L,L′} = Σ_M Y_{LM} Y*_{L′M}**, same algebra as explicit **m_i**, **m_f** average)."
  ([T-amplitudes S-factor k-i k-f theta mass-factor]
   (transfer-differential-cross-section-angular-m-sum T-amplitudes S-factor k-i k-f theta mass-factor mass-factor))
  ([T-amplitudes S-factor k-i k-f theta mass-factor-i mass-factor-f]
   (transfer-differential-cross-section-angular-m-sum T-amplitudes S-factor k-i k-f theta mass-factor-i mass-factor-f 0.0 nil nil))
  ([T-amplitudes S-factor k-i k-f theta mass-factor-i mass-factor-f phi]
   (transfer-differential-cross-section-angular-m-sum T-amplitudes S-factor k-i k-f theta mass-factor-i mass-factor-f phi nil nil))
  ([T-amplitudes S-factor k-i k-f theta mass-factor-i mass-factor-f phi l-i l-f]
   (let [angular-dist (transfer-angular-distribution-m-sum-unpolarized T-amplitudes theta phi l-i l-f)
         prefactor (/ (* mass-factor-i mass-factor-f)
                      (* 16.0 Math/PI Math/PI))
         k-ratio (/ k-f k-i)
         dsigma (* 10.0 prefactor k-ratio angular-dist S-factor)]
     dsigma)))

;; -----------------------------------------------------------------------------
;; Point-charge **elastic** Rutherford (optional reference / ratios only — not used in transfer DWBA)
;;
;; Single-nucleon **transfer** dσ/dΩ uses `distorted-wave-optical` (Coulomb + nuclear in χ), radial
;; POST integrals, and an angular factor (coherent multipole or orbital m-sum — see `dwba.benchmark.ca40-dwuck`) — no
;; separate “Coulomb tail” on σ.
;; -----------------------------------------------------------------------------

(def ^:private transfer-coulomb-e2-mev-fm
  "Coulomb constant **e²** in **MeV·fm** (same convention as optical / examples)."
  1.44)

(def ^:private transfer-rutherford-theta-floor-rad
  "Minimum CM **θ** (rad) for analytic **1/sin⁴(θ/2)** so θ→0 stays finite (~2°)."
  (* 2.0 (/ Math/PI 180.0)))

(defn- transfer-rutherford-effective-theta-rad
  ^double [^double theta-rad]
  (let [t (Math/abs (double theta-rad))
        lo transfer-rutherford-theta-floor-rad]
    (if (< t lo) lo t)))

(defn transfer-rutherford-dsigma-mb-sr
  "Point-charge **elastic** Rutherford **dσ/dΩ** in **mb/sr**, non-relativistic CM:
  **(Z₁ Z₂ e² / (4 E_cm))² / sin⁴(θ/2)**, then ×10 (**fm²/sr → mb/sr**). For **elastic** comparisons
  only — **transfer** DWBA does not add this to stripping/pickup σ.

  **θ** in radians. Below **`transfer-rutherford-theta-floor-rad`** (~2°), **θ** is floored so the
  pole at θ→0 is finite."
  [Z1 Z2 ^double e-cm-mev ^double theta-rad]
  (let [e2 transfer-coulomb-e2-mev-fm
        zze (* (double Z1) (double Z2) e2)
        th (transfer-rutherford-effective-theta-rad theta-rad)
        s (max 1e-15 (Math/abs (Math/sin (* 0.5 th))))
        p (/ zze (* 4.0 e-cm-mev))]
    (* 10.0 (/ (* p p) (* s s s s)))))

;; -----------------------------------------------------------------------------
;; Nuclear spin / j-coupling for single-nucleon transfer (optional on top of |T|² × S)
;; -----------------------------------------------------------------------------

(defn transfer-nuclear-spin-statistical-factor
  "Unpolarized target + unpolarized residual m-average/sum: (2J_f+1)/(2J_i+1).
  J_i, J_f are total angular momenta of initial and final nuclei (integer or half-integer)."
  [J-i J-f]
  (/ (+ (* 2.0 (double J-f)) 1.0) (+ (* 2.0 (double J-i)) 1.0)))

(defn transfer-unpolarized-deuteron-spin-factor
  "Spin average for an unpolarized deuteron beam: 1/3."
  []
  (/ 1.0 3.0))

(defn transfer-one-nucleon-recoupling-6j
  "Wigner 6j factor {J_i J_f j; l ½ j} for coupling core spins to the transferred
  nucleon's (l,½)j single-particle state (stripping/pickup). Returns 0 if triads fail.
  Uses exact 6j from dwba.angular-momentum."
  [J-i J-f l j]
  (jam/wigner-6j J-i J-f j l 0.5 j))

(defn transfer-one-nucleon-spin-prefactor
  "Dimensionless multiplier for dσ/dΩ when folding in nuclear spins (multiply after
  the usual kinematic × |T(θ)|² × S from transfer-differential-cross-section-angular).

  Combines:
    • (2J_f+1)/(2J_i+1) — initial/final spin statistical factor
    • 1/3 if :deuteron-unpolarized? true (default) — unpolarized deuteron entrance
    • (2j+1) × {J_i J_f j; l ½ j}² — lsj recoupling (nucleon spin ½)

  Convention note: some works absorb (2j+1) or pieces of the 6j² into the definition of
  the spectroscopic factor C²S; compare your S before interpreting absolute scale.

  opts: :deuteron-unpolarized? (default true)."
  [J-i J-f l j & {:keys [deuteron-unpolarized?] :or {deuteron-unpolarized? true}}]
  (let [stat (transfer-nuclear-spin-statistical-factor J-i J-f)
        d-fac (if deuteron-unpolarized? (transfer-unpolarized-deuteron-spin-factor) 1.0)
        sixj (transfer-one-nucleon-recoupling-6j J-i J-f l j)
        rec (* (inc (* 2.0 (double j))) (* sixj sixj))]
    (* stat d-fac rec)))

(defn- transfer-imaginary-unit-to-integer-power
  "i^ℓ for integer ℓ (Edmonds / Austern phase)."
  [^long ell]
  (case (Math/floorMod ell 4)
    0 (complex-cartesian 1.0 0.0)
    1 (complex-cartesian 0.0 1.0)
    2 (complex-cartesian -1.0 0.0)
    3 (complex-cartesian 0.0 -1.0)))

(defn austern-zr-chi-exit-mass-ratio
  "Mass ratio **M_A/M_B** for **N. Austern**, *Direct Nuclear Reaction Theories*, **Eq. (5.3)**.
  After the zero-range δ-contraction, the volume integral is **∫ d³ r** in the **entrance**
  relative coordinate **r**, with the **exit** distorted wave **χ_β^{(-)}** sampled at
  **(M_A/M_B) r** (often **M_B = M_A + M_x**, residual = core + transferred cluster). Masses in
  any consistent units (e.g. MeV/c²).

  Pass as **`:zr-chi-exit-mass-ratio`** to `transfer-amplitude-post` (**`:zero-range` only**)."
  [M-A M-B]
  (/ (double M-A) (double M-B)))

(defn handbook-zr-chi-exit-mass-ratio
  "Same as **`austern-zr-chi-exit-mass-ratio`** — **M_A/M_B** for ZR **χ_β** radius scaling in handbook DWBA."
  [M-A M-B]
  (austern-zr-chi-exit-mass-ratio M-A M-B))

(defn- austern-real-scalar
  ^double [x]
  (if (number? x) (double x) (double (re x))))

(defn- austern-radial-simpson-integrate-real
  "Composite Simpson's rule, uniform **h**, real samples **f(i·h)**, **i = 0 … n−1**."
  ^double [^double h integrand-vec]
  (let [n (count integrand-vec)]
    (cond
      (zero? n) 0.0
      (= n 1) (* h (austern-real-scalar (first integrand-vec)))
      :else
      (let [simpson-sum (loop [i 1 sum 0.0]
                          (if (>= i (dec n))
                            sum
                            (let [coeff (if (odd? i) 4.0 2.0)
                                  term (* coeff (austern-real-scalar (get integrand-vec i)))]
                              (recur (inc i) (+ sum term)))))
            h3 (/ h 3.0)
            f0 (austern-real-scalar (first integrand-vec))
            fl (austern-real-scalar (nth integrand-vec (dec n)))]
        (* h3 (+ f0 fl simpson-sum))))))

(defn- austern-radial-simpson-integrate-complex
  "Composite Simpson's rule, uniform **h**, complex (or real) samples **f(i·h)**."
  [^double h integrand-vec]
  (let [n (count integrand-vec)]
    (cond
      (zero? n) (complex-cartesian 0.0 0.0)
      (= n 1) (mul (complex-cartesian h 0.0)
                   (if (number? (first integrand-vec))
                     (complex-cartesian (double (first integrand-vec)) 0.0)
                     (first integrand-vec)))
      :else
      (let [simpson-sum (loop [i 1 sum (complex-cartesian 0.0 0.0)]
                          (if (>= i (dec n))
                            sum
                            (let [coeff (if (odd? i) 4.0 2.0)
                                  term-val (get integrand-vec i)
                                  term (if (number? term-val)
                                         (complex-cartesian (* coeff (double term-val)) 0.0)
                                         (mul (complex-cartesian coeff 0.0) term-val))]
                              (recur (inc i) (add sum term)))))
            h3 (/ h 3.0)
            f0 (let [v (first integrand-vec)]
                 (if (number? v) (complex-cartesian (double v) 0.0) v))
            fl (let [v (nth integrand-vec (dec n))]
                 (if (number? v) (complex-cartesian (double v) 0.0) v))]
        (mul (complex-cartesian h3 0.0) (add f0 fl simpson-sum))))))

(defn austern-radial-integral-prefactor-eq-5-5
  "Prefactor **(M_B/M_A)(4π/(k_α k_β))** from **N. Austern**, **Eq. (5.5)**.
  **M_B = M_A + M_x** in the book; **k_α**, **k_β** in **fm⁻¹** (same convention as your radial waves)."
  ^double [^double M-A ^double M-B ^double k-alpha ^double k-beta]
  (* (/ M-B M-A) (* 4.0 Math/PI) (/ 1.0 (* k-alpha k-beta))))

(defn f-alphaL
  "Austern **Eq. (5.5)** entrance factor **f_{α L_α}(k_α, r)** on **r = i·h**.

  **u_α** is the **reduced** distorted wave (**u = rR**) from **`distorted-wave-optical`** for partial wave
  **L_α**. Returns a vector **f(r_i)** with **(count u_α)** entries:

  **f_{α L_α}(k_α, r_i) := R_α(r_i) = u_α(r_i)/r_i** via the same **u → R** rule as **`transfer-amplitude-post`**.

  **k_α** is implicit in which solution you solved. If the book’s **f** differs from **R** (e.g. **u/(ikr)** or Coulomb
  normalization), multiply this vector or replace with your convention before **`austern-radial-integral-I-Lb-La-eq-5-5`**."
  [u-alpha ^double h]
  (radial-R-from-reduced-u u-alpha h))

(defn f-betaL
  "Austern **Eq. (5.5)** exit factor **f_{β L_β}(k_β, (M_A/M_B)\\, r)** on the **entrance** grid **r_i = i·h**.

  **u_β** is the reduced exit distorted wave for partial wave **L_β**. For each **i = 0 … n−1**, samples
  **R_β** at argument **ρ \\, r_i** with **ρ = zr-mass-ratio** (**M_A/M_B**, same as **`austern-zr-chi-exit-mass-ratio`** /
  **`:zr-chi-exit-mass-ratio`**: exit χ evaluated at scaled radius in ZR). **ρ = 1** gives **f_{β L_β}(k_β, r_i)**.

  Returns a vector of length **n-points** (use **n = min(count F, count u_α, count u_β)** for ZR integrands).
  Linear interpolation: **`transfer-radial-R-at-r-linear`** on **R_β = u_β/r**.

  **Convention:** same **f ≡ R** note as **`f-alphaL`**."
  [u-beta ^double h zr-mass-ratio n-points]
  (let [R-beta (radial-R-from-reduced-u u-beta h)
        rho (double zr-mass-ratio)
        n (long n-points)]
    (mapv (fn [^long i]
            (transfer-radial-R-at-r-linear R-beta h (* rho (* (double i) h))))
          (range n))))

(defn austern-radial-integrand-zr-F-Ra-Rb-r2-complex
  "Like **`austern-radial-integrand-zr-F-Ra-Rb-r2`**, but **f_α**, **f_β** may be **complex** (**R = u/r**).
  Entries: **F(r) f_α(r) f_β((M_A/M_B)r) r²** (same indexing as real ZR DWBA)."
  [F-vec u-alpha u-beta h zr-mass-ratio]
  (let [h* (double h)
        n (min (count F-vec) (count u-alpha) (count u-beta))
        fa (radial-R-from-reduced-u u-alpha h*)
        fb (f-betaL u-beta h* (double zr-mass-ratio) n)]
    (mapv (fn [^long i]
            (let [r (* (double i) h*)
                  fa-i (get fa i)
                  fb-i (get fb i)
                  Fi (get F-vec i)
                  Fi-c (if (number? Fi) (complex-cartesian (double Fi) 0.0) Fi)
                  fa-c (if (number? fa-i) (complex-cartesian (double fa-i) 0.0) fa-i)
                  fb-c (if (number? fb-i) (complex-cartesian (double fb-i) 0.0) fb-i)
                  rsq (* r r)]
              (mul Fi-c fa-c fb-c (complex-cartesian rsq 0.0))))
          (range n))))

(defn austern-radial-integral-I-zr-eq-5-5-from-u-complex
  "Same as **`austern-radial-integral-I-zr-eq-5-5-from-u`** (**Austern (5.5)** prefactor × Simpson), but
  allows **complex** **f_{αL}**, **f_{βL}** from **complex** optical **`distorted-wave-optical`** **u** (absorption)."
  [F-vec u-alpha u-beta h M-A M-B k-alpha k-beta zr-mass-ratio]
  (let [integrand (austern-radial-integrand-zr-F-Ra-Rb-r2-complex
                    F-vec u-alpha u-beta (double h) (double zr-mass-ratio))
        sum (austern-radial-simpson-integrate-complex (double h) integrand)
        pref (austern-radial-integral-prefactor-eq-5-5 (double M-A) (double M-B)
              (double k-alpha) (double k-beta))]
    (mul (complex-cartesian pref 0.0) sum)))

(defn austern-radial-integral-I-Lb-La-eq-5-5
  "Radial integral **I_{L_β L_α}^{ℓsj}** from **N. Austern**, *Direct Nuclear Reaction Theories*,
  **Eq. (5.5)**:

  \\[
  I_{L_\\beta L_\\alpha}^{\\ell s j}
  = \\frac{M_B}{M_A}\\,\\frac{4\\pi}{k_\\alpha k_\\beta}
  \\int_0^\\infty \\mathrm{d}r\\,
  F_{\\ell s j}(r)\\,
  f_{\\beta L_\\beta}\\!\\left(k_\\beta,\\frac{M_A}{M_B}r\\right)\\,
  f_{\\alpha L_\\alpha}(k_\\alpha,r) \\, .
  \\]

  Pass **integrand-real-vec** with entry **i** = **F(r_i) f_β((M_A/M_B)r_i) f_α(r_i)** at **r_i = i·h**
  (the **integrand** of the book's **dr** integral — **not** necessarily including **r²**; see below).

  **Convention bridge:** partial-wave codes often reduce **∫ d³r** to **4π ∫ R_β R_α r² dr** for fixed
  **(L_α,L_β)**. Build **f_α**, **f_β** with **`f-alphaL`** and **`f-betaL`**, or use **`austern-radial-integrand-zr-F-Ra-Rb-r2`**
  and then this function. If your optical code defines **f** differently (e.g. **u/(kr)**), fold that
  into **integrand-real-vec** before calling here.

  **Returns:** real **I** for **`austern-reduced-amplitude-beta-sum-eq-5-6`** (Austern **(5.5)** prefactor). Handbook §5.5.2 **I** on the same integrand: **`handbook-radial-integral-I-zr`** + **`handbook-zr-multipole-amplitude-sum`**."
  [integrand-real-vec h M-A M-B k-alpha k-beta]
  {:pre [(sequential? integrand-real-vec)
         (every? number? [h M-A M-B k-alpha k-beta])
         (pos? (double k-alpha))
         (pos? (double k-beta))]}
  (* (austern-radial-integral-prefactor-eq-5-5 (double M-A) (double M-B)
      (double k-alpha) (double k-beta))
     (austern-radial-simpson-integrate-real (double h) integrand-real-vec)))

(defn austern-radial-integrand-zr-F-Ra-Rb-r2
  "Per-point samples **F(r) f_{αL}(r) f_{βL}((M_A/M_B)r) r²** on **r = i·h**, with **f** from **`f-alphaL`**
  / **`f-betaL`** (**R** convention). **u_α**, **u_β** reduced (**u=rR**). **zr-mass-ratio** = **M_A/M_B**."
  [F-vec u-alpha u-beta h zr-mass-ratio]
  (let [h* (double h)
        n (min (count F-vec) (count u-alpha) (count u-beta))
        fa (f-alphaL u-alpha h*)
        fb (f-betaL u-beta h* zr-mass-ratio n)]
    (mapv (fn [^long i]
            (let [r (* (double i) h*)
                  fa-i (austern-real-scalar (get fa i))
                  fb-i (austern-real-scalar (get fb i))
                  Fi (austern-real-scalar (get F-vec i))
                  rsq (* r r)]
              (* Fi fa-i fb-i rsq)))
          (range n))))

(defn austern-radial-integral-I-zr-eq-5-5-from-u
  "Convenience: **(5.5)** with ZR integrand **F f_{αL} f_{βL} r²** (`**f-alphaL**`, **`f-betaL`**) when
  **F_{ℓsj}** is already on a vector (**one** **L_α**, **one** **L_β** partial wave per **u_α**, **u_β**).
  To build **F** from bound states automatically, use **`austern-radial-integral-I-eq-5-5-from-F-lsj`**."
  [F-vec u-alpha u-beta h M-A M-B k-alpha k-beta zr-mass-ratio]
  (let [integrand (austern-radial-integrand-zr-F-Ra-Rb-r2
                    F-vec u-alpha u-beta (double h) (double zr-mass-ratio))]
    (austern-radial-integral-I-Lb-La-eq-5-5 integrand h M-A M-B k-alpha k-beta)))

(defn austern-radial-integral-I-eq-5-5-from-F-lsj
  "Full **N. Austern**, **Eq. (5.5)** radial integral **I_{L_β L_α}^{ℓsj}** from bound + distorted **u** grids.

  1. **F_{ℓsj}(r_i) = R_{φ_f}^* R_{φ_i}** — **`F-lsj-r-from-bound-reduced-u`** on **φ_i**, **φ_f** (reduced bound **u**, **E < 0**).
  2. **f_{α L_α}(k_α, r)**, **f_{β L_β}(k_β, (M_A/M_B)r)** — **`f-alphaL`**, **`f-betaL`** inside **`austern-radial-integrand-zr-F-Ra-Rb-r2`**.
  3. **(M_B/M_A)(4π/(k_α k_β)) ∫ dr F f_α f_β** (Simpson) — **`austern-radial-integral-I-Lb-La-eq-5-5`**.
  
  Uses **n = min(count φ_i, count φ_f, count u_α, count u_β)** samples via the shared integrand; align grids.
  **k_α**, **k_β** (**fm⁻¹**). **zr-mass-ratio** = **M_A/M_B** (`**austern-zr-chi-exit-mass-ratio**`).
  
  Equivalent to **`(austern-radial-integral-I-zr-eq-5-5-from-u (F-lsj-r-from-bound-reduced-u φ_i φ_f h) u_α u_β …)`**."
  [phi-i-reduced-u phi-f-reduced-u u-alpha u-beta h M-A M-B k-alpha k-beta zr-mass-ratio]
  (let [F-vec (F-lsj-r-from-bound-reduced-u phi-i-reduced-u phi-f-reduced-u h)]
    (austern-radial-integral-I-zr-eq-5-5-from-u
      F-vec u-alpha u-beta h M-A M-B k-alpha k-beta zr-mass-ratio)))

;; -----------------------------------------------------------------------------
;; **Handbook** single-nucleon **F** (§5.4): **F_{ℓsj} = R_n = u/r** on **`austern-radial-integrand-zr-F-Ra-Rb-r2`**.
;;
;; **I** for **`handbook-zr-multipole-amplitude-sum`** + **`T = D₀ √(2ℓ+1) β`** + **`transfer-differential-cross-section`**
;; must use the same **(5.5)** prefactor **(M_B/M_A)(4π/(k_α k_β))** as **`austern-radial-integral-I-zr-eq-5-5-from-u`**
;; (LNPS / Austern angular reduction is derived for that **I**).  A printed §5.5.2 line often shows **√(4π)(M_B/M_A)/(k_α k_β)** —
;; that matches a **reduced** amplitude in the book, not this **β → T → dσ/dΩ** chain; using **√(4π)** here undershoots **|T|²** by **4π**.
;;
;; **Display only:** **`handbook-radial-integral-prefactor-handbook-display-zr`** (ratio **√(4π)/(4π)** to **(5.5)**).
;; -----------------------------------------------------------------------------

(defn handbook-radial-integral-prefactor-handbook-display-zr
  "Handbook §5.5.2 **printed** factor **√(4π)\\,(M_B/M_A)/(k_α k_β)** — for comparison to tables only.

  **Do not** use for **`handbook-radial-integral-I-zr`**; that fn uses **(5.5)** **`austern-radial-integral-prefactor-eq-5-5`** so **dσ/dΩ** matches the **β** multipole sum.

  Ratio **display / (5.5)** = **√(4π)/(4π)**."
  ^double [^double M-A ^double M-B ^double k-alpha ^double k-beta]
  (* (Math/sqrt (* 4.0 Math/PI)) (/ M-B M-A) (/ 1.0 (* k-alpha k-beta))))

(defn handbook-radial-integral-prefactor-zr
  "Deprecated alias for **`handbook-radial-integral-prefactor-handbook-display-zr`** (§5.5.2 display)."
  ^double [^double M-A ^double M-B ^double k-alpha ^double k-beta]
  (handbook-radial-integral-prefactor-handbook-display-zr M-A M-B k-alpha k-beta))

(defn handbook-radial-integral-I-zr
  "ZR radial **I_{L_β L_α}** with handbook **F = R_n** and **N. Austern (5.5)** prefactor **(M_B/M_A)(4π/(k_α k_β))** × Simpson on **`austern-radial-integrand-zr-F-Ra-Rb-r2`**.

  Equals **`austern-radial-integral-I-zr-eq-5-5-from-u`** for the same **F**-vector (only **F** construction differs: **R_n** vs overlaps). Use with **`handbook-zr-multipole-amplitude-sum`**."
  [F-vec u-alpha u-beta h M-A M-B k-alpha k-beta zr-mass-ratio]
  (let [h* (double h)
        integrand (austern-radial-integrand-zr-F-Ra-Rb-r2
                   F-vec u-alpha u-beta h* (double zr-mass-ratio))]
    (* (austern-radial-integral-prefactor-eq-5-5 (double M-A) (double M-B)
         (double k-alpha) (double k-beta))
       (austern-radial-simpson-integrate-real h* integrand))))

(defn handbook-radial-integral-I-zr-from-neutron-bound
  "Convenience: build **F** from **neutron** bound **u** only, then **`handbook-radial-integral-I-zr`**."
  [neutron-u-reduced u-alpha u-beta h M-A M-B k-alpha k-beta zr-mass-ratio]
  (handbook-radial-integral-I-zr
    (handbook-F-lsj-radial-from-neutron-bound-u neutron-u-reduced h)
    u-alpha u-beta h M-A M-B k-alpha k-beta zr-mass-ratio))

(defn austern-reduced-amplitude-beta-sj-ellm
  "Reduced amplitude **β_{sj}^{ℓm}** exactly as **N. Austern**, *Direct Nuclear Reaction
  Theories* (Wiley, 1970), **Eq. (4.60), p. 84**:

  \\[
  (2\\ell+1)^{1/2}\\, i^{\\ell}\\, \\beta_{sj}^{\\ell m}
  = \\int \\mathrm{d}^3 r_\\alpha \\int \\mathrm{d}^3 r_\\beta\\,
  \\chi_\\beta^{(-)*}(\\mathbf{k}_\\beta, \\mathbf{r}_\\beta)\\,
  f_{\\ell s j,\\,m}(\\mathbf{r}_\\beta, \\mathbf{r}_\\alpha)\\,
  \\chi_\\alpha^{(+)}(\\mathbf{k}_\\alpha, \\mathbf{r}_\\alpha) \\, .
  \\]

  **Zero-range (Eq. (5.3)):** the same left-hand side equals **∫ d³ r χ_β^*((M_A/M_B)\\mathbf r)
  F_{\\ell sj}(r) Y_{\\ell}^{m*}(\\hat r) χ_α(\\mathbf r)** (Austern’s notation). The **bound** radial
  content in **F_{ℓsj}** is the **handbook** single-particle **R_n = u_n/r** (**`handbook-F-lsj-radial-from-neutron-bound-u`**):
  **`transfer-amplitude-post`** picks **u_n** from **`:handbook-F-from`** **`:phi-i`** (default) or **`:phi-f`**.
  Angular **Y_{ℓm}** and full **f_{ℓsj,m}** are not expanded in that shortcut — partial wave and **D₀**
  carry the rest. ( **`F-lsj-r-from-bound-reduced-u`** keeps Austern’s **R_{φ_f}^* R_{φ_i}** for benchmarks.)

  **Partial-wave / angle:** after ZR radial integrals **I_{L_β L_α}** are known,
  **`austern-reduced-amplitude-beta-sum-eq-5-6`** (Austern **(5.6)**) or **`handbook-zr-multipole-amplitude-sum`**
  (same sum; use with **`handbook-radial-integral-I-zr`**) in the coplanar frame (**ẑ ‖ k_α**). Each **(L_α, L_β)** term carries **Y_{L_β}^{-m}(Θ, 0)** (Austern’s **Θ** between
  **k_α** and **k_β**). It is simpler than **(5.4)** (no **M_α**, **M_β** sums) but still requires
  that radial table and consistent **σ_{αL}**, **σ_{βL}**.

  **No Clebsch–Gordan factors** appear in (4.60); those enter the DW amplitude **T** in
  **Eq. (4.59)** — use `austern-dw-transition-amplitude-T-term-4-59`.

  **Parameters:** **ℓ** = orbital **l** in the book; **m** = projection on **f_{ℓsj,m}**;
  **s, j** label the multipole (book subscript on **β**); **I** = value of the double
  integral (4.60), or the single volume integral (5.3) in ZR — typically from your DWBA radial code.

  **Returns:** complex **β_{sj}^{ℓm} = I / (√(2ℓ+1) · i^ℓ)**."
  [ell m-ell s j I]
  {:pre [(every? number? [ell m-ell s j])]}
  (let [Ic (if (number? I) (complex-cartesian (double I) 0.0) I)
        sqrt-2l1 (Math/sqrt (inc (* 2.0 (double ell))))
        denom (mul (complex-cartesian sqrt-2l1 0.0)
                   (transfer-imaginary-unit-to-integer-power (long ell)))]
    (div Ic denom)))

(defn- zr-dwba-cg-product-Lbeta-l-Lalpha
  "Shared **⟨L_β ℓ;0,0|L_α 0⟩ ⟨L_β ℓ;−m,m|L_α 0⟩** (Edmonds order via **`jam/clebsch-gordan-exact`**)."
  ^double [L-alpha L-beta ell m]
  (let [L-a (double L-alpha)
        L-b (double L-beta)
        l (double ell)
        mproj (double m)]
    (* (double (jam/clebsch-gordan-exact L-b 0.0 l 0.0 L-a 0.0))
       (double (jam/clebsch-gordan-exact L-b (- mproj) l mproj L-a 0.0)))))

(defn austern-eq-5-6-cg-product-Lbeta-l-Lalpha
  "Clebsch–Gordan **product** from **N. Austern**, **Eq. (5.6)**:

  \\[
  \\langle L_\\beta \\ell; 0,\\, 0 \\,|\\, L_\\alpha\\, 0 \\rangle\\,
  \\langle L_\\beta \\ell; -m,\\, m \\,|\\, L_\\alpha\\, 0 \\rangle \\, .
  \\]

  **Notation:** Austern writes **⟨L_β ℓ; m_β, m_ℓ | L_α M⟩** (two **j**’s, then two **m**’s after the semicolon). That is
  the same coefficient as **⟨L_β m_β; ℓ m_ℓ | L_α M⟩** or Racah/Edmonds **⟨L_β m_β ℓ m_ℓ | L_α M⟩** — not a different coupling.
  Only **swapping** **L_β ↔ ℓ** would introduce **(−1)^{L_β+ℓ−L_α}**; we keep **L_β** first to match the book.

  **Convention:** **`jam/clebsch-gordan-exact`** implements **⟨j₁ m₁ j₂ m₂ | j₃ m₃⟩** with **j₁ = L_β**, **j₂ = ℓ**,
  **j₃ = L_α** (matches **`austern-reduced-amplitude-beta-sum-eq-5-6`**).

  **Parameters:** **L_α**, **L_β**, **ℓ** (integers or half-integers supported by **3j**), **m** (projection on **f_{ℓsj,m}**).

  **Returns:** real **double**."
  ^double [L-alpha L-beta ell m]
  (zr-dwba-cg-product-Lbeta-l-Lalpha L-alpha L-beta ell m))

(defn handbook-zr-cg-product-Lbeta-l-Lalpha
  "Same CG product as **`austern-eq-5-6-cg-product-Lbeta-l-Lalpha`** — ZR DWBA coplanar multipole geometry."
  ^double [L-alpha L-beta ell m]
  (zr-dwba-cg-product-Lbeta-l-Lalpha L-alpha L-beta ell m))

(defn- zr-dwba-admissible-L-beta-values
  [^long L-alpha ^long ell ^long L-max]
  (let [la L-alpha
        l ell
        lb-min (max 0 (Math/abs (- la l)))
        lb-max (min L-max (+ la l))]
    (if (> lb-min lb-max)
      ()
      (filter (fn [^long lb] (even? (+ la lb l)))
              (range lb-min (inc lb-max))))))

(defn austern-eq-5-6-admissible-L-beta-values
  "**L_β** that can contribute to **Austern (5.6)** for fixed entrance partial wave **L_α** and transferred orbital **ℓ**.

  The book writes **∑_{L_α,L_β}**, but **L⃗_α = L⃗_β + ℓ⃗** in angular-momentum coupling fixes **L_β** once **L_α**, **ℓ** are chosen: triangle
  **|L_α − ℓ| ≤ L_β ≤ L_α + ℓ** (with **L_β ≥ 0**) and **L_α + L_β + ℓ** even so **⟨L_β ℓ; 0,0 | L_α 0⟩** need not vanish.
  The second **CG** in **(5.6)** uses projections **(−m, m)** on the **ℓ** leg but does not enlarge this **L_β** span for integer **L**.

  **L-max** — keep **L_β ≤ L-max**. Returns a seq of long **L_β** (possibly empty)."
  [^long L-alpha ^long ell ^long L-max]
  (zr-dwba-admissible-L-beta-values L-alpha ell L-max))

(defn handbook-zr-partial-wave-L-beta-values
  "Handbook / standard ZR DWBA: admissible **L_β** for fixed **L_α**, transferred **ℓ**, cap **L_max** (triangle + parity)."
  [^long L-alpha ^long ell ^long L-max]
  (zr-dwba-admissible-L-beta-values L-alpha ell L-max))

(defn nuclear-phase-shifts-map
  "Partial-wave **nuclear** (short-range) phase shifts **δ^n_L** (**radians**) for **L = 0 … L-max**
  (**S_L^n = e^{2iδ^n_L}** in T&N **(3.1.84)**).

  **`functions/phase-shift`** = **½ arg(S_L^n)** with **S_L^n = `s-matrix`** (same as neutral **`phase-shift0`**, no added **σ**).
  This map is **`{L → δ^n_L}`** = **`phase-shift`** per **L**.

  **`austern-radial-rows-with-sigma`** forms **σ_L^total = σ_L^Coulomb + δ^n_L** for **(5.6)**.

  Bind **`mass-factor`**, **`Z1Z2ee`**, and **`functions/*elastic-imag-ws-params*`** like **`s-matrix`**.

  **Returns:** map **`{L → δ^n_L}`** (long keys)."
  [^double e-cm v-params ^long l-max]
  (into {}
        (map (fn [^long L]
               [L (double (phase-shift e-cm v-params L))])
             (range 0 (inc l-max)))))

(defn- austern-delta-branch
  ^double [delta-spec ^long L]
  (cond
    (nil? delta-spec) 0.0
    (number? delta-spec) (double delta-spec)
    (map? delta-spec) (double (or (get delta-spec L) (get delta-spec (int L)) 0.0))
    (fn? delta-spec) (double (delta-spec L))
    :else 0.0))

(defn austern-radial-rows-with-sigma
  "Set **:sigma-alpha** and **:sigma-beta** for Austern **Eq. (5.6)** as **Coulomb + nuclear**:

  **σ_{α L_α} = σ_L^Coulomb(η_α) + δ_α(L_α)**, **σ_{β L_β} = σ_L^Coulomb(η_β) + δ_β(L_β)** (**radians**).

  **δ_α** / **δ_β**: **`nil`** → **0**; a **number** → same for all **L**; a **map** **`{L → δ}`**; or **`(fn [L] δ)`**.
  Build maps with **`nuclear-phase-shifts-map`**.

  **η_α**, **η_β**: **`channel-sommerfeld-eta`** with **`binding`** **`mass-factor`**, **`Z1Z2ee`** per channel.

  **Coulomb only:** **`austern-radial-rows-with-coulomb-sigma`** (**`nil`** nuclear parts)."
  [radial-rows eta-alpha eta-beta delta-alpha delta-beta]
  (let [dα #(austern-delta-branch delta-alpha %)
        dβ #(austern-delta-branch delta-beta %)]
    (mapv (fn [row]
            (let [La (long (:L-alpha row))
                  Lb (long (:L-beta row))]
              (-> row
                  (assoc :sigma-alpha (+ (coulomb-sigma-L La (double eta-alpha)) (dα La)))
                  (assoc :sigma-beta (+ (coulomb-sigma-L Lb (double eta-beta)) (dβ Lb))))))
          radial-rows)))

(defn austern-radial-rows-with-coulomb-sigma
  "Coulomb **σ_L(η)** only (**no** **`phase-shift`** nuclear part). See **`austern-radial-rows-with-sigma`** for **σ_C + δ**.

  **σ_{α L_α} = arg Γ(L_α+1+iη_α)**, same **`coulomb-sigma-L`** as **`Hankel±`**. **η** from **`channel-sommerfeld-eta`**."
  [radial-rows eta-alpha eta-beta]
  (austern-radial-rows-with-sigma radial-rows eta-alpha eta-beta nil nil))

(defn handbook-zr-rows-with-coulomb-sigma
  "Attach Coulomb **σ_L(η)** only to **`{:L-alpha :L-beta :I}`** rows (*Handbook* / ZR DWBA partial-wave sums)."
  [radial-rows eta-alpha eta-beta]
  (austern-radial-rows-with-sigma radial-rows eta-alpha eta-beta nil nil))

(defn- zr-multipole-amplitude-sum-rows
  "Coplanar ZR DWBA multipole sum over **`radial-rows`** (**I**, **σ**, **L_α**, **L_β**). Shared by handbook and Austern APIs."
  [ell m-ell theta-rad radial-rows]
  {:pre [(every? number? [ell m-ell theta-rad]) (sequential? radial-rows)]}
  (let [l (double ell)
        mproj (double m-ell)
        th (double theta-rad)]
    (reduce
     (fn [acc row]
       (let [L-a (double (:L-alpha row))
             L-b (double (:L-beta row))
             Ival (:I row)
             sa (double (or (:sigma-alpha row) 0.0))
             sb (double (or (:sigma-beta row) 0.0))
             cg-prod (zr-dwba-cg-product-Lbeta-l-Lalpha L-a L-b l mproj)
             La (long (Math/round L-a))
             Lb (long (Math/round L-b))
             ll (long (Math/round l))
             mm (long (Math/round mproj))]
         (if (< (Math/abs cg-prod) 1e-20)
           acc
           (let [pow (Math/floorMod (- La Lb ll) 4)
                 iphase (transfer-imaginary-unit-to-integer-power pow)
                 eph (complex-polar (+ sa sb) 1.0)
                 sqrt2lb (Math/sqrt (inc (* 2.0 L-b)))
                 ybm (spherical-harmonic Lb (- mm) th 0.0)
                 Ic (if (number? Ival) (complex-cartesian (double Ival) 0.0) Ival)
                 pref (complex-cartesian (* sqrt2lb cg-prod) 0.0)
                 term (mul iphase eph pref Ic ybm)]
             (add acc term)))))
     (complex-cartesian 0.0 0.0)
     radial-rows)))

(defn handbook-zr-multipole-amplitude-sum
  "ZR DWBA angular weight **∑_{L_α,L_β} …** for coplanar kinematics (handbook multipole reduction).

  **radial-rows** — maps **`{:L-alpha :L-beta :I :sigma-alpha :sigma-beta}`** with **I** from **`handbook-radial-integral-I-zr-from-neutron-bound`** (**(5.5)** prefactor **(M_B/M_A)(4π/(k_α k_β))** on **∫ F R_α R_β r² dr** — not the §5.5.2 **√(4π)** display line). Returns complex **∝ β^{ℓm}(θ)** for **T = D₀ √(2ℓ+1) β** and **`transfer-differential-cross-section`**.

  **Same algebra** as **`austern-reduced-amplitude-beta-sum-eq-5-6`**."
  [ell m-ell theta-rad radial-rows]
  (zr-multipole-amplitude-sum-rows ell m-ell theta-rad radial-rows))

(defn austern-reduced-amplitude-beta-sum-eq-5-6
  "Partial-wave form of **β_{sj}^{ℓm}** from **N. Austern**, *Direct Nuclear Reaction Theories*,
  **Eq. (5.6)** — simpler than **(5.4)** because the **z**-axis is along **k_α** and **y** along
  **k_α × k_β**, so **k_β** lies in the **x–z** plane at polar angle **Θ** (CM angle between
  **k_α** and **k_β**).

  **Partial-wave bookkeeping:** for fixed **ℓ**, the **CG** **⟨L_β ℓ; 0,0 | L_α 0⟩** leaves only one independent orbital index — e.g. sum **L_α** and **L_β** in **`austern-eq-5-6-admissible-L-beta-values`** — not an unconstrained **(L_α,L_β)** grid.

  \\[
  \\beta_{sj}^{\\ell m}
  = \\sum_{L_\\alpha,\\,L_\\beta}
  i^{L_\\alpha - L_\\beta - \\ell}\\,
  e^{i(\\sigma_{\\alpha L_\\alpha} + \\sigma_{\\beta L_\\beta})}\\,
  \\sqrt{2L_\\beta+1}\\,
  I_{L_\\beta L_\\alpha}^{\\ell s j}\\,
  \\langle L_\\beta \\ell; 0, 0 \\,|\\, L_\\alpha\\, 0 \\rangle\\,
  \\langle L_\\beta \\ell; -m,\\, m \\,|\\, L_\\alpha\\, 0 \\rangle\\,
  Y_{L_\\beta}^{-m}(\\Theta,\\, 0) \\, .
  \\]

  **Angular piece (your Eq. (5.6) reading):** For each **(L_α, L_β)** pair, the contribution to
  **β_{sj}^{ℓm}** is proportional to **Y_{L_β}^{-m}(Θ, 0)** — same **Θ** as above, **φ = 0** (coplanar;
  **ẑ ‖ k_α**, **k_β** in the **x–z** plane). Implementation: **`spherical-harmonic` L_β (−m) Θ 0**. The full
  **β^{ℓm}(Θ)** is the **sum** over **L_α, L_β** of that factor times **I**, **√(2L_β+1)**, the two CGs,
  and **i^{L_α−L_β−ℓ} e^{i(σ_{αL_α}+σ_{βL_β})}**.

  **When to use:** Yes — **(5.6)** is the right angular reduction **once** you supply **I_{L_β L_α}**
  from **Eq. (5.5)** (radial factors **f_{α L_α}**, **f_{β L_β}**, **F_{ℓsj}**, and the **M_B/M_A**
  prefactor there). It is **not** a substitute for that radial work: it removes the **M_α**, **M_β**
  sums of **(5.4)** only. You still need **σ_{α L}**, **σ_{β L}** consistent with your distorted waves
  (Coulomb ± nuclear, same convention as the code that produces **f**).

  **radial-rows** — collection of maps **`{:L-alpha :L-beta :I :sigma-alpha :sigma-beta}`**:
  **:I** = **I_{L_β L_α}^{ℓsj}** with Austern **(5.5)** prefactor; **σ** (**radians**). Use **`austern-radial-rows-with-sigma`** (Coulomb + nuclear via **`phase-shift`**) or **`austern-radial-rows-with-coulomb-sigma`** (Coulomb only). Handbook §5.5.2 **I**: **`handbook-zr-multipole-amplitude-sum`** (same sum).

  **Clebsch–Gordan product:** **`austern-eq-5-6-cg-product-Lbeta-l-Lalpha`**. **Y:** `spherical-harmonic`.

  **Returns:** complex **β_{sj}^{ℓm}** (same object as in **(4.59)** — **s**, **j** are implicit in each **I**)."
  [ell m-ell theta-rad radial-rows]
  (zr-multipole-amplitude-sum-rows ell m-ell theta-rad radial-rows))

(defn austern-dw-transition-amplitude-T-term-4-59
  "One **(ℓ,s,j)** contribution to the distorted-wave amplitude **T_{αβ}^{DW}** from
  **Austern**, **Eq. (4.59), p. 84**:

  \\[
  T = \\sum_{\\ell,s,j} (2\\ell+1)^{1/2} A_{\\ell s j}\\, (-)^{s_b - m_b}\\,
  \\langle J_A\\, j;\\, M_A,\\, M_B - M_A \\,|\\, J_B\\, M_B \\rangle\\,
  \\langle \\ell\\, s;\\, m,\\, m_a - m_b \\,|\\, j,\\, m - m_b + m_a \\rangle\\,
  \\langle s_a\\, s_b;\\, m_a,\\, -m_b \\,|\\, s,\\, m_a - m_b \\rangle\\,
  \\beta_{sj}^{\\ell m} \\, .
  \\]

  Pass **beta** = **β_{sj}^{ℓm}** from `austern-reduced-amplitude-beta-sj-ellm`. **A-lsj** is
  the dynamical factor **A_{ℓsj}** in the book. **Clebsch–Gordan** via
  `clebsch-gordan-exact` (Edmonds **⟨j₁m₁ j₂m₂|j₃m₃⟩** order).

  **Note:** **(2ℓ+1)^{1/2} β** = **I/i^ℓ** from (4.60), so this term could also be written as
  **A (-)^{s_b-m_b} (CG_J)(CG_{ℓs})(CG_{ss}) I / i^ℓ** (the **√(2ℓ+1)** cancels).

  **Channel vs nucleon spin:** In **⟨s_a s_b | s⟩**, **s** is often **1** (deuteron triplet). The factor
  **⟨ℓ s; … | j …⟩** couples the transferred nucleon’s **ℓ** to **intrinsic spin ½** to form **j = ℓ ± ½**,
  so that leg’s second angular momentum is **½**, not **1**. Pass **`:s-for-ls-coupling ½`** when **s = 1**
  for the deuteron leg; default **`s-for-ls-coupling`** = **s** (backward compatible)."
  ([A-lsj ell m-ell s j J-A M-A J-B M-B s-a m-a s-b m-b beta]
   (austern-dw-transition-amplitude-T-term-4-59
     A-lsj ell m-ell s j J-A M-A J-B M-B s-a m-a s-b m-b beta nil))
  ([A-lsj ell m-ell s j J-A M-A J-B M-B s-a m-a s-b m-b beta opts]
   (let [s-ls (double (or (:s-for-ls-coupling opts) s))
         beta-c (if (number? beta) (complex-cartesian (double beta) 0.0) beta)
         phase (Math/pow -1.0 (- (double s-b) (double m-b)))
         sqrt-2l1 (Math/sqrt (inc (* 2.0 (double ell))))
         ;; ⟨J_A j; M_A, M_B - M_A | J_B M_B⟩
         cg-J (jam/clebsch-gordan-exact (double J-A) (double M-A) (double j)
                (- (double M-B) (double M-A)) (double J-B) (double M-B))
         ;; ⟨ℓ s_ls; m, m_a - m_b | j, m - m_b + m_a⟩
         cg-ls (jam/clebsch-gordan-exact (double ell) (double m-ell) s-ls
                 (- (double m-a) (double m-b)) (double j)
                 (+ (double m-ell) (double m-a) (- (double m-b))))
         ;; ⟨s_a s_b; m_a, -m_b | s, m_a - m_b⟩
         cg-ss (jam/clebsch-gordan-exact (double s-a) (double m-a) (double s-b)
                 (- (double m-b)) (double s) (- (double m-a) (double m-b)))
         pref (* (double A-lsj) phase sqrt-2l1 (double cg-J) (double cg-ls) (double cg-ss))]
     (mul (complex-cartesian pref 0.0) beta-c))))

(defn satchler-reduced-amplitude-eq13-diagonal-spin
  "**Not** Austern (4.60). **G. R. Satchler**, *Nucl. Phys.* **55** (1964) **Eq. (13)** in the
  **diagonal-spin** limit (no spin-flip in χ: only **m_a, m_b** term in the sum). Relates
  the double integral **\\mathcal{I}** to a **different** reduced quantity with **√(2j+1)**
  and two Clebsch factors (not **√(2ℓ+1)** as in Austern):

  \\[
  \\sqrt{2j+1}\\; i^{\\ell}\\, \\tilde\\beta
  = (-1)^{s_b-m_b}\\,
  \\langle \\ell\\, m_\\ell;\\, s,\\, M_J-m_\\ell \\,|\\, j\\, M_J \\rangle\\,
  \\langle s_a\\, m_a;\\, s_b,\\, -m_b \\,|\\, s,\\, m_a-m_b \\rangle\\;
  \\mathcal{I} \\,,
  \\quad M_J \\equiv M_B-M_A \\, .
  \\]

  Use **`austern-reduced-amplitude-beta-sj-ellm`** when you need **β** as in the textbook
  **Eq. (4.60)**."
  [ell m-ell s j s-a m-a s-b m-b M-J I]
  (let [Ic (if (number? I) (complex-cartesian (double I) 0.0) I)
        cg-ls (jam/clebsch-gordan-exact (double ell) (double m-ell) (double s)
                (- (double M-J) (double m-ell)) (double j) (double M-J))
        cg-sab (jam/clebsch-gordan-exact (double s-a) (double m-a) (double s-b)
                 (- (double m-b)) (double s) (- (double m-a) (double m-b)))
        phase (Math/pow -1.0 (- (double s-b) (double m-b)))
        numer-sc (* phase cg-ls cg-sab)
        numer (mul Ic (complex-cartesian numer-sc 0.0))
        sqrt-2j1 (Math/sqrt (inc (* 2.0 (double j))))
        denom (mul (complex-cartesian sqrt-2j1 0.0)
                   (transfer-imaginary-unit-to-integer-power (long ell)))]
    (div numer denom)))

(defn transfer-differential-cross-section-angular-with-spin
  "Same as transfer-differential-cross-section-angular, multiplied by
  transfer-one-nucleon-spin-prefactor for target/residual spins J_i, J_f, bound orbital
  l and total single-particle j of the transferred nucleon. Pass optional kwargs for
  transfer-one-nucleon-spin-prefactor (e.g. :deuteron-unpolarized? false)."
  [T-amplitudes S-factor k-i k-f theta mass-factor-i mass-factor-f phi l-i l-f J-i J-f l j & opts]
  (* (transfer-differential-cross-section-angular T-amplitudes S-factor k-i k-f theta
                                                 mass-factor-i mass-factor-f phi l-i l-f)
     (double (apply transfer-one-nucleon-spin-prefactor J-i J-f l j opts))))

(defn transfer-differential-cross-section-angular-coherent-with-spin
  "Coherent multipole angular distribution (`transfer-differential-cross-section-angular-coherent`)
  times `transfer-one-nucleon-spin-prefactor` — full single-nucleon spin / 6j folding on top of
  interfering L amplitudes."
  [T-amplitudes S-factor k-i k-f theta mass-factor-i mass-factor-f phi l-i l-f J-i J-f l j & opts]
  (* (transfer-differential-cross-section-angular-coherent T-amplitudes S-factor k-i k-f theta
                                                         mass-factor-i mass-factor-f phi l-i l-f)
     (double (apply transfer-one-nucleon-spin-prefactor J-i J-f l j opts))))

(defn transfer-total-cross-section
  "Calculate total cross-section by integrating differential cross-section.
   
   The total cross-section is:
   σ_total = ∫ dσ/dΩ dΩ = ∫₀^π ∫₀^2π dσ/dΩ(θ,φ) sin(θ) dθ dφ
   
   For coplanar scattering (φ=0), this simplifies to:
   σ_total = 2π ∫₀^π dσ/dΩ(θ) sin(θ) dθ
   
   Uses Simpson's rule for numerical integration.
   
   Parameters:
   - T-amplitudes: Map of {L → T_L} transfer amplitudes
   - S-factor: Spectroscopic factor
   - k-i: Wavenumber in entrance channel (fm⁻¹)
   - k-f: Wavenumber in exit channel (fm⁻¹)
   - mass-factor: Mass factor (2μ/ħ²)
   - n-points: Number of integration points (default: 100)
   - theta-min: Minimum angle (default: 0)
   - theta-max: Maximum angle (default: π)
   
   Returns: σ_total in **mb** (1 fm² = 10 mb; integral result × 10)
   
   Example:
   (let [T-map {0 1.0, 1 0.5}
         k-i (Math/sqrt (* mass-factor 10.0))
         k-f (Math/sqrt (* mass-factor 8.0))
         S 0.5]
     (transfer-total-cross-section T-map S k-i k-f mass-factor))"
  ([T-amplitudes S-factor k-i k-f mass-factor]
   (transfer-total-cross-section T-amplitudes S-factor k-i k-f mass-factor 100 0.0 Math/PI))
  ([T-amplitudes S-factor k-i k-f mass-factor n-points]
   (transfer-total-cross-section T-amplitudes S-factor k-i k-f mass-factor n-points 0.0 Math/PI))
  ([T-amplitudes S-factor k-i k-f mass-factor n-points theta-min theta-max]
   (let [;; Get angular distribution function
         angular-dist-fn (transfer-angular-distribution-function T-amplitudes theta-min theta-max n-points)
         ;; Reduced mass factor: μ/(2πħ²) = mass-factor/(4π)
         mu-factor (/ mass-factor (* 4.0 Math/PI))
         ;; Wavenumber ratio: k_f/k_i
         k-ratio (/ k-f k-i)
         ;; Integrate: ∫ dσ/dΩ(θ) sin(θ) dθ using Simpson's rule
         integral (loop [i 1 sum 0.0]
                   (if (>= i (dec (count angular-dist-fn)))
                     sum
                     (let [[theta dsigma-angular] (nth angular-dist-fn i)
                           coeff (if (odd? i) 4.0 2.0)
                           sin-theta (Math/sin theta)
                           term (* coeff dsigma-angular sin-theta)]
                       (recur (inc i) (+ sum term)))))
         ;; Final integral value: (h/3) * [f₀ + fₙ + sum]
         h-theta (/ (- theta-max theta-min) (dec n-points))
         [theta-first dsigma-first] (first angular-dist-fn)
         [theta-last dsigma-last] (last angular-dist-fn)
         final-integral (* (/ h-theta 3.0)
                          (+ (* dsigma-first (Math/sin theta-first))
                             (* dsigma-last (Math/sin theta-last))
                             integral))
         ;; Total cross-section: 2π · (μ/(2πħ²))² · (k_f/k_i) · integral · S → mb
         sigma-total (* 10.0 2.0 Math/PI mu-factor mu-factor k-ratio final-integral S-factor)]
     sigma-total)))

(defn transfer-kinematic-factors
  "Calculate kinematic factors for transfer reactions.
   
   Returns wavenumbers and energy factors needed for cross-section calculations.
   
   Parameters:
   - E-i: Incident energy in CM frame (MeV)
   - E-f: Final energy in CM frame (MeV)
   - mass-factor: Mass factor (2μ/ħ²) in MeV⁻¹·fm⁻²
   
   Returns: Map with {:k-i, :k-f, :k-ratio, :E-i, :E-f}
   
   Example:
   (transfer-kinematic-factors 10.0 8.0 mass-factor)
   => {:k-i 0.5, :k-f 0.45, :k-ratio 0.9, :E-i 10.0, :E-f 8.0}"
  [E-i E-f mass-factor]
  (let [k-i (Math/sqrt (* mass-factor E-i))
        k-f (Math/sqrt (* mass-factor E-f))
        k-ratio (/ k-f k-i)]
    {:k-i k-i
     :k-f k-f
     :k-ratio k-ratio
     :E-i E-i
     :E-f E-f}))

(defn transfer-lab-to-cm
  "Convert transfer cross-section from lab frame to CM frame.
   
   For transfer reactions, the transformation is:
   dσ/dΩ_CM = dσ/dΩ_lab · (dΩ_lab/dΩ_CM)
   
   where the Jacobian depends on the masses and angles.
   
   Parameters:
   - dsigma-lab: Differential cross-section in lab frame (mb/sr)
   - theta-lab: Scattering angle in lab frame (radians)
   - theta-cm: Scattering angle in CM frame (radians)
   - m-a: Mass of projectile a (amu)
   - m-A: Mass of target A (amu)
   - m-b: Mass of outgoing particle b (amu)
   - m-B: Mass of residual nucleus B (amu)
   
   Returns: dσ/dΩ_CM in mb/sr
   
   Note: This is a simplified version. Full transformation requires
   solving the kinematic equations for the specific reaction.
   For now, we use a simple sin ratio approximation."
  [dsigma-lab theta-lab theta-cm _m-a _m-A _m-b _m-B]
  (let [;; Simplified Jacobian: dΩ_lab/dΩ_CM ≈ sin(θ_lab)/sin(θ_cm)
        ;; This is an approximation - full calculation requires solving kinematic equations
        sin-lab (Math/sin theta-lab)
        sin-cm (Math/sin theta-cm)
        jacobian (if (< sin-cm 1e-10)
                  1.0  ; Avoid division by zero
                  (/ sin-lab sin-cm))]
    (* dsigma-lab jacobian)))

(defn transfer-cm-to-lab
  "Convert transfer cross-section from CM frame to lab frame.
   
   Inverse of transfer-lab-to-cm.
   
   Parameters:
   - dsigma-cm: Differential cross-section in CM frame (mb/sr)
   - theta-lab: Scattering angle in lab frame (radians)
   - theta-cm: Scattering angle in CM frame (radians)
   - m-a: Mass of projectile a (amu)
   - m-A: Mass of target A (amu)
   - m-b: Mass of outgoing particle b (amu)
   - m-B: Mass of residual nucleus B (amu)
   
   Returns: dσ/dΩ_lab in mb/sr"
  [dsigma-cm theta-lab theta-cm _m-a _m-A _m-b _m-B]
  (let [;; Inverse Jacobian: dΩ_CM/dΩ_lab ≈ sin(θ_cm)/sin(θ_lab)
        ;; This is an approximation - full calculation requires solving kinematic equations
        sin-lab (Math/sin theta-lab)
        sin-cm (Math/sin theta-cm)
        jacobian (if (< sin-lab 1e-10)
                  1.0  ; Avoid division by zero
                  (/ sin-cm sin-lab))]
    (* dsigma-cm jacobian)))

;; ============================================================================
;; OPTICAL POTENTIALS FOR TRANSFER REACTIONS
;; ============================================================================

(defn optical-potential-woods-saxon
  "Calculate optical potential with Woods-Saxon form factor.
   
   The optical potential has the form:
   U(r) = -V · f_V(r) - iW · f_W(r) + V_so · (l·s) · f_so(r) + V_C(r)
   
   where:
   - V is the real potential depth (MeV)
   - W is the imaginary potential depth (MeV) - accounts for absorption
   - V_so is the spin-orbit coupling strength (MeV)
   - f_V, f_W, f_so are Woods-Saxon form factors
   - V_C is the Coulomb potential
   
   Parameters:
   - r: Radial distance (fm)
   - V-params: Real potential parameters [V0, R_V, a_V]
   - W-params: Imaginary potential parameters [W0, R_W, a_W] (optional)
   - V-so: Spin-orbit potential depth (MeV, optional)
   - R-so: Spin-orbit radius (fm, optional)
   - a-so: Spin-orbit diffuseness (fm, optional)
   - l: Orbital angular momentum
   - s: Spin (1/2 for nucleons)
   - j: Total angular momentum
   - Z1, Z2: Charges of projectile and target (for Coulomb, optional)
   - R-C: Coulomb radius (fm, optional)
   
   Returns: Complex potential U(r) in MeV
   
   Example:
   (optical-potential-woods-saxon 2.0 
                                   [50.0 2.0 0.6]  ; Real potential
                                   [10.0 2.0 0.6]  ; Imaginary potential
                                   7.0 2.0 0.6     ; Spin-orbit
                                   1 0.5 1.5       ; l, s, j
                                   1 8 2.0)        ; Z1, Z2, R_C)"
  ([r V-params]
   (optical-potential-woods-saxon r V-params nil nil nil nil nil nil nil nil nil nil))
  ([r V-params W-params]
   (optical-potential-woods-saxon r V-params W-params nil nil nil nil nil nil nil nil nil))
  ([r V-params W-params V-so R-so a-so l s j Z1 Z2 R-C]
   (let [;; Real Woods-Saxon potential: -V0/(1+exp((r-R_V)/a_V))
         [V0 R-V a-V] V-params
         f-V (/ 1.0 (+ 1.0 (Math/exp (/ (- r R-V) a-V))))
         V-real (* -1.0 V0 f-V)
         
         ;; Imaginary Woods-Saxon potential: -iW0/(1+exp((r-R_W)/a_W))
         W-imag (if W-params
                 (let [[W0 R-W a-W] W-params
                       f-W (/ 1.0 (+ 1.0 (Math/exp (/ (- r R-W) a-W))))]
                   (* -1.0 W0 f-W))
                 0.0)
         
         ;; Spin-orbit coupling: V_so · (l·s) · f_so(r)
         ;; (l·s) = (j(j+1) - l(l+1) - s(s+1))/2
         V-so-term (if (and V-so R-so a-so l s j)
                    (if (zero? r)
                      0.0  ; Spin-orbit term is zero at r=0
                      (let [f-so (/ 1.0 (+ 1.0 (Math/exp (/ (- r R-so) a-so))))
                            ;; Derivative of f-so for spin-orbit form factor
                            df-so-dr (/ (* f-so (- 1.0 f-so)) a-so)
                            l-dot-s (/ (- (* j (+ j 1.0))
                                         (* l (+ l 1.0))
                                         (* s (+ s 1.0)))
                                      2.0)
                            V-so-radial (* V-so l-dot-s df-so-dr (/ 1.0 r))]
                        V-so-radial))
                    0.0)
         
         ;; Coulomb potential: Z1*Z2*e²/r for r > R_C, else Z1*Z2*e²*(3 - r²/R_C²)/(2*R_C)
         V-coulomb (if (and Z1 Z2 R-C)
                    (let [Z1Z2e2 (* Z1 Z2 1.44)]  ; e² = 1.44 MeV·fm
                      (if (> r R-C)
                        (/ Z1Z2e2 r)
                        (* Z1Z2e2 (/ (- 3.0 (/ (* r r) (* R-C R-C))) (* 2.0 R-C)))))
                   0.0)
         
         ;; Total potential: V_real + i*W_imag + V_so + V_coulomb
         total-real (+ V-real V-so-term V-coulomb)]
     (complex-cartesian total-real W-imag))))

(defn optical-potential-parameters
  "Get standard optical potential parameters for common reactions.
   
   Returns parameter sets for:
   - Real potential: [V0, R_V, a_V]
   - Imaginary potential: [W0, R_W, a_W]
   - Spin-orbit: [V_so, R_so, a_so]
   
   Parameters:
   - projectile: Projectile type (:p, :n, :d, :alpha, etc.)
   - target-A: Target mass number
   - E-lab: Lab energy (MeV)
   
   Returns: Map with {:V-params, :W-params, :V-so, :R-so, :a-so}
   
   Note: These are typical values. Actual parameters should be fitted
   to experimental data for specific reactions.
   
   Example:
   (optical-potential-parameters :p 16 10.0)
   => {:V-params [50.0 2.0 0.6], :W-params [10.0 2.0 0.6], ...}"
  [projectile target-A E-lab]
  (let [;; Radius parameter: R = r0 * A^(1/3)
        r0 1.25  ; fm
        R-base (* r0 (Math/pow target-A (/ 1.0 3.0)))
        a-diff 0.6  ; Typical diffuseness (fm)
        
        ;; Real potential depth (MeV) - energy dependent
        V0-base (case projectile
                 :p 50.0
                 :n 50.0
                 :d 100.0
                 :alpha 200.0
                 50.0)  ; Default
        V0-energy-dep (* V0-base (- 1.0 (* 0.003 E-lab)))  ; Rough energy dependence
        
        ;; Imaginary potential depth (MeV)
        W0-base (case projectile
                 :p 10.0
                 :n 10.0
                 :d 20.0
                 :alpha 30.0
                 10.0)  ; Default
        W0-energy-dep (* W0-base (Math/sqrt E-lab))  ; Rough energy dependence
        
        ;; Spin-orbit parameters
        V-so (case projectile
              :p 7.0
              :n 7.0
              :d 0.0  ; No spin-orbit for deuterons typically
              :alpha 0.0
              0.0)  ; Default: no spin-orbit
        
        R-so R-base
        a-so a-diff]
    {:V-params [V0-energy-dep R-base a-diff]
     :W-params [W0-energy-dep R-base a-diff]
     :V-so V-so
     :R-so R-so
     :a-so a-so}))

(defn optical-potential-energy-dependent
  "Calculate energy-dependent optical potential.
   
   Optical potential parameters often depend on energy. This function
   implements common energy-dependent parameterizations. When :global-set
   is provided (e.g. :ch89 for Chapel Hill 89), uses that global potential
   instead of the built-in generic one.
   
   Parameters:
   - r: Radial distance (fm)
   - projectile: Projectile type (:p, :n, :d, :alpha)
   - target-A: Target mass number
   - E-lab: Lab energy (MeV)
   - l: Orbital angular momentum
   - s: Spin
   - j: Total angular momentum
   - Optional: :Z1, :Z2 charges, :R-C Coulomb radius, :global-set (e.g. :ch89)
   
   Returns: Complex potential U(r) in MeV
   
   Example:
   (optical-potential-energy-dependent 2.0 :p 16 10.0 1 0.5 1.5 :Z1 1 :Z2 8 :global-set :ch89)"
  [r projectile target-A E-lab l s j & {:keys [Z1 Z2 R-C global-set]}]
  (let [;; Auto-select global potential if not explicitly provided
        auto-global-set (cond
                         (and (#{:p :n} projectile) (nil? global-set)) :ch89
                         (and (= projectile :d) (nil? global-set)) :daehnick80
                         :else global-set)
        effective-global-set (or global-set auto-global-set)]
    (cond
      ;; Use CH89 for protons and neutrons
      (and (#{:ch89} effective-global-set) (#{:p :n} projectile))
      (do (require 'dwba.global-potentials)
          ((resolve 'dwba.global-potentials/optical-potential-ch89)
           r projectile target-A (long (or Z2 0)) E-lab l s j))
      
      ;; Use Daehnick80 for deuterons
      (and (#{:daehnick80} effective-global-set) (= projectile :d))
      (do (require 'dwba.global-potentials)
          ((resolve 'dwba.global-potentials/optical-potential-daehnick80)
           r target-A (long (or Z2 0)) E-lab l s j))
      
      ;; Fall back to generic parameters
      :else
      (let [params (optical-potential-parameters projectile target-A E-lab)
            V-params (:V-params params)
            W-params (:W-params params)
            V-so (:V-so params)
            R-so (:R-so params)
            a-so (:a-so params)]
        (optical-potential-woods-saxon r V-params W-params V-so R-so a-so l s j Z1 Z2 R-C)))))

(defn optical-potential-entrance-channel
  "Calculate optical potential for entrance channel of transfer reaction.
   
   This is a convenience function that sets up the optical potential
   for the entrance channel (projectile + target). 
   
   Automatically uses:
   - CH89 (Chapel Hill 89) for protons (:p) and neutrons (:n)
   - Daehnick 1980 for deuterons (:d)
   
   You can override by passing :global-set explicitly.
   
   Parameters:
   - r: Radial distance (fm)
   - projectile-type: Type of projectile (:p, :n, :d, :alpha)
   - target-A: Target mass number
   - target-Z: Target charge number
   - E-lab: Lab energy (MeV) or energy per nucleon for deuterons
   - l: Orbital angular momentum
   - s: Spin
   - j: Total angular momentum
   - Optional: :global-set to override auto-selection (:ch89, :daehnick80, or nil for generic)
   
   Returns: Complex potential U(r) in MeV
   
   Example:
   (optical-potential-entrance-channel 2.0 :p 16 8 20.0 1 0.5 1.5)
   => Automatically uses CH89 for proton
   
   (optical-potential-entrance-channel 2.0 :d 16 8 10.0 1 1.0 2.0)
   => Automatically uses Daehnick80 for deuteron"
  [r projectile-type target-A target-Z E-lab l s j & {:keys [global-set]}]
  (let [projectile-Z (case projectile-type
                      :p 1
                      :n 0
                      :d 1
                      :alpha 2
                      0)]
    (apply optical-potential-energy-dependent
           r projectile-type target-A E-lab l s j
           (concat [:Z1 projectile-Z
                    :Z2 target-Z
                    :R-C (* 1.25 (Math/pow target-A (/ 1.0 3.0)))]
                   (when global-set [:global-set global-set])))))

(defn optical-potential-exit-channel
  "Calculate optical potential for exit channel of transfer reaction.
   
   Sets up the optical potential for the exit channel.
   
   Automatically uses:
   - CH89 (Chapel Hill 89) for protons (:p) and neutrons (:n)
   - Daehnick 1980 for deuterons (:d)
   
   You can override by passing :global-set explicitly.
   
   Parameters:
   - r: Radial distance (fm)
   - outgoing-type: Type of outgoing particle (:p, :n, :d, :alpha)
   - residual-A: Residual nucleus mass number
   - residual-Z: Residual nucleus charge number
   - E-lab: Lab energy in exit channel (MeV) or energy per nucleon for deuterons
   - l: Orbital angular momentum
   - s: Spin
   - j: Total angular momentum
   - Optional: :global-set to override auto-selection (:ch89, :daehnick80, or nil for generic)
   
   Returns: Complex potential U(r) in MeV
   
   Note: For transfer reactions, typically use CH89 for proton channel and Daehnick80 for deuteron channel."
  [r outgoing-type residual-A residual-Z E-lab l s j & {:keys [global-set]}]
  (let [outgoing-Z (case outgoing-type
                    :p 1
                    :n 0
                    :d 1
                    :alpha 2
                    0)]
    (apply optical-potential-energy-dependent
           r outgoing-type residual-A E-lab l s j
           (concat [:Z1 outgoing-Z
                    :Z2 residual-Z
                    :R-C (* 1.25 (Math/pow residual-A (/ 1.0 3.0)))]
                   (when global-set [:global-set global-set])))))

(defn f-r-numerov-complex
  "Calculate f(r) for Numerov method with complex optical potential.
   
   f(r) = (2μ/ħ²) * [E - U(r)] - l(l+1)/r²
   
   Parameters:
   - r: Radial distance (fm)
   - E: Energy (MeV)
   - l: Orbital angular momentum
   - U: Complex potential U(r) (MeV)
   - mass-factor: Mass factor (2μ/ħ²)
   
   Returns: Complex f(r) value"
  [r E l U mass-factor]
  (let [;; Check if U is NaN or Inf
        U-valid (if (number? U)
                  (and (not (Double/isNaN U)) (not (Double/isInfinite U)))
                  (and (not (Double/isNaN (re U))) (not (Double/isInfinite (re U)))
                       (not (Double/isNaN (im U))) (not (Double/isInfinite (im U)))))
        ;; Centrifugal term: l(l+1)/r²
        centrifugal (if (zero? r)
                     0.0  ; Avoid division by zero
                     (/ (* l (+ l 1.0)) (* r r)))
        ;; Effective potential term: (2μ/ħ²) * [E - U(r)]
        U-effective (if (not U-valid)
                     (complex-cartesian 0.0 0.0)  ; Return zero if U is invalid
                     (if (number? U)
                       (- E U)
                       (subt E U)))
        potential-term (mul mass-factor U-effective)
        ;; Total: f(r) = (2μ/ħ²)[E-U] - l(l+1)/r²
        centrifugal-complex (complex-cartesian centrifugal 0.0)
        f-total (subt potential-term centrifugal-complex)]
    f-total))

(defn distorted-wave-optical
  "Calculate distorted wave using optical potential with complex Numerov.

   **Context (DWBA):** This **χ_L** is the **same physical object** as the elastic distorted wave in that partial
   wave: solution of the radial equation with **U = V_C + V_opt**, matched at the nuclear surface to **outgoing**
   Coulomb behavior (**R-matrix / log-derivative** ↔ **Hankel±**), as **`functions/s-matrix`** does for the
   elastic **S_L^n** quotient. Here you get **χ(r)** on a grid for overlap integrals; optional **`:coulomb-tail`**
   enforces Coulomb normalization at **r_max** only (see **ns** doc).

   Solves the Schrödinger equation with optical potential:
   -∇²/2μ · χ + U(r) · χ = E · χ
   
   Parameters:
   - E: Energy (MeV)
   - l: Orbital angular momentum
   - s: Spin
   - j: Total angular momentum
   - optical-potential-fn: Function r → U(r) (complex potential)
   - r-max: Maximum radius (fm)
   - h: Step size (fm)
   - mass-factor: Mass factor (2μ/ħ²)
   - Optional keywords:
     - **`:normalize-mode`** — **`:max`** (default): scale χ so **max |χ| = 1** (legacy).
       **`:coulomb-tail`**: scale χ so **|u_last|** matches **|H_L^+(η, ρ)|** at **ρ = k r_max** with
       **u = r χ** reduced wave on the last grid point (same **η**, **k** as **`channel-sommerfeld-eta`** /
       **`Hankel+`**). This is **magnitude-only** tail matching for **χ** — **not** elastic **S_L^n** and **not**
       a substitute for **`functions/s-matrix`** (R-matrix / **Hankel** quotient). Improves **relative** partial-wave
       weights vs interior **`:max`** norm and tames **θ → π** blow-up in multi-**L** Austern **(5.6)** sums when
       high **L** radial integrals stay large. Pass **`:tail-eta`** (Sommerfeld **η**) and **`:tail-rho`** (**k r_max**, dimensionless).
   
   Returns: Vector of complex distorted wave values χ(r)
   
   Note: Uses complex Numerov integration. The wavefunction will be complex
   due to the imaginary part of the optical potential (absorption).
   
   Example:
   (let [U-fn (fn [r] (optical-potential-entrance-channel r :d 16 8 10.0 1 0.5 1.5))]
     (distorted-wave-optical 10.0 1 0.5 1.5 U-fn 20.0 0.01 mass-factor))"
  [E l s j optical-potential-fn r-max h mass-factor
   & {:keys [normalize-mode tail-eta tail-rho]
      :or {normalize-mode :max}}]
  (let [steps (int (/ r-max h))
        h2-12 (/ (* h h) 12.0)
        ;; Initial conditions: u(0) = 0, u(h) ≈ h^(l+1)
        u0 (complex-cartesian 0.0 0.0)
        u1-init (Math/pow h (inc l))
        u1 (complex-cartesian u1-init 0.0)
        
        ;; Pre-calculate f(r) values at all radial points
        rs (take (+ steps 2) (iterate #(+ % h) 0.0))
        fs (mapv (fn [r]
                  (let [U (optical-potential-fn r)]
                    (f-r-numerov-complex r E l U mass-factor)))
                rs)]
    ;; Complex Numerov integration
    (let [results (loop [n 1
                         results-vec [u0 u1]]
                    (if (>= n (dec steps))
                      results-vec
                      (let [un (get results-vec n)
                            un-1 (get results-vec (dec n))
                            fn-1 (get fs (dec n))
                            fn (get fs n)
                            fn+1 (get fs (inc n))
                            ;; Use complex Numerov step
                            un+1 (numerov-step-complex un un-1 fn-1 fn fn+1 h)]
                        (recur (inc n) (conj results-vec un+1)))))
          u-last (last results)
          ul (if (number? u-last) (Math/abs u-last) (mag u-last))
          tail-norm
          (when (= normalize-mode :coulomb-tail)
            (let [eta (double tail-eta)
                  rho (double tail-rho)
                  Hc (Hankel+ l eta rho)
                  Hmag (if (number? Hc) (Math/abs Hc) (mag Hc))]
              (when (and (pos? ul) (pos? Hmag) (Double/isFinite Hmag) (Double/isFinite ul))
                (/ Hmag ul))))
          max-mag (transduce (map #(if (number? %)
                                    (Math/abs %)
                                    (mag %)))
                            (completing max)
                            Double/NEGATIVE_INFINITY
                            results)
          norm-factor (or tail-norm
                          (if (and (> max-mag 1e-10) (< max-mag 1e20))
                            (/ 1.0 max-mag)
                            1.0))]
      (mapv (fn [u]
              (if (number? u)
                (* norm-factor u)
                (mul norm-factor u)))
            results))))

(defn optical-potential-summary
  "Get summary of optical potential parameters.
   
   Returns a formatted string with all optical potential parameters.
   
   Parameters:
   - projectile: Projectile type
   - target-A: Target mass number
   - E-lab: Lab energy (MeV)
   
   Returns: String summary
   
   Example:
   (optical-potential-summary :p 16 10.0)"
  [projectile target-A E-lab]
  (let [params (optical-potential-parameters projectile target-A E-lab)
        V-params (:V-params params)
        W-params (:W-params params)
        V-so (:V-so params)
        R-so (:R-so params)
        a-so (:a-so params)]
    (format "Optical Potential Parameters:
Projectile: %s
Target A: %d
Lab Energy: %.2f MeV
Real Potential: V0=%.2f MeV, R=%.2f fm, a=%.2f fm
Imaginary Potential: W0=%.2f MeV, R=%.2f fm, a=%.2f fm
Spin-Orbit: V_so=%.2f MeV, R=%.2f fm, a=%.2f fm"
            (name projectile) target-A E-lab
            (first V-params) (second V-params) (last V-params)
            (first W-params) (second W-params) (last W-params)
            V-so R-so a-so)))
