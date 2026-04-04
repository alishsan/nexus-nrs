;; **⁹Li + n** two-body system: search **negative-energy bound states** and **positive-energy
;; resonance candidates** from **partial-wave phase shifts** (neutral Woods–Saxon + Numerov).
;;
;; **Pauli / valence channels:** **⁹Li** has **N = 6** neutrons (Z = 3). In a spherical **j-j** cartoon,
;; **2** fill **1s₁/₂** and **4** fill **1p₃/₂** (spin–orbit lowers **p₃/₂** first), so those subshells are
;; **full**. The next neutron (valence in **¹⁰Li** relative to a **⁹Li** core) then occupies **1p₁/₂**—the
;; natural ℓ = 1 valence—and that matches transfer work assigning **p₁/₂** removal strength to **¹⁰Li**.
;; By default the **ℓ = 0** channel is omitted from the **continuum** scan (shell note; optional block probes **s**).
;; A separate block still asks: **if ℓ = 0 were allowed, would this Woods–Saxon bind an
;; s-wave?** With the default **V₀ = 50 MeV** geometry and **⁹Li+n** reduced mass, the shooting
;; solver finds **no** small-|u(R)| bound state (same story as **¹⁰Li** being unbound—deeper
;; **V₀** would eventually bind **s**). **ℓ = 1** (**p** waves) **is included** in main scans.
;;
;; **Literature:** **¹⁰Li** is unbound; **¹¹Li(p,d)¹⁰Li** (IRIS, TRIUMF) sees a strong **p₁/₂**-like
;; resonance at **E_r = 0.62 ± 0.04 MeV**, **Γ = 0.33 ± 0.07 MeV** (Sanetullaev *et al.*,
;; *Phys. Lett. B* **755** (2016) 481–485). For a fuller experimental summary and kinematics see
;; **`examples/example_11Li_pd_10Li.clj`**. A simple real potential cannot replace three-body /
;; continuum-Gamow treatments—this is a **schematic** partial-wave survey only.
;;
;; **WS parameter scan:** The script **coarse-scans** **V₀, R** (and fixed **a** from **`-main`**) over a grid,
;; searches **0.05–5.5 MeV** for the **lowest ℓ = 1** **δ = π/2** crossing when it exists, **bisects** **E_res**,
;; and sets **Γ_est ≈ 2/|dδ/dE|** (Wigner). If **no π/2 crossing** appears (common for shallow wells), it falls
;; back to the **interior** local maximum of **(2ℓ+1)sin²δ** closest to **E_r** (peak centroid, no **Γ_est**).
;; This **does not** replace DWBA / Gamow fits—edit **`:v0-list`**, **`:r-list`**, **`:a-list`**, **`:e-max`** in
;; **`-main`** to refine the search around a hit (e.g. **V₀ ≈ 36 MeV**, **R ≈ 2.35 fm** gives **E ≈ 0.65 MeV**
;; in the default grid).
;;
;; **Run (repo root):** `lein run -m clojure.main examples/example_9Li_n_bound_resonances.clj`
;;
(ns examples.example-9Li-n-bound-resonances
  (:require [functions :refer [mass-factor mass-factor-from-mu
                              solve-numerov phase-shift-from-numerov]]
            [dwba.transfer :as tr]))

(def ^:private m-n 939.5654205)
;; M = A·u + Δ (MeV/c²), u = 931.494 MeV, Δ from mass tables (rounded, same spirit as `example_11Li_pd_10Li.clj`).
(def ^:private m-9Li (+ (* 9.0 931.494) 12.415))

(defn reduced-mass-9Li-n
  ^double []
  (/ (* m-9Li m-n) (+ m-9Li m-n)))

(defn local-max-indices
  "Indices i in (1..n-2) where ys[i] is strictly greater than neighbors (interior local maxima)."
  [ys]
  (let [n (count ys)]
    (if (< n 3)
      []
      (filter (fn [i]
                (let [y (nth ys i)
                      y- (nth ys (dec i))
                      y+ (nth ys (inc i))]
                  (and (> y y-) (> y y+)
                       (every? Double/isFinite [y y- y+]))))
              (range 1 (dec n))))))

(defn resonance-candidates-from-sin2
  "\"Resonance\" proxy: local maxima of (2ℓ+1)·sin²δ(E) on the energy grid (narrow peaks)."
  [l curve]
  (let [weights (map (fn [{:keys [delta]}]
                       (* (inc (* 2 l))
                          (let [s (Math/sin delta)]
                            (* s s))))
                     curve)
        idxs (local-max-indices (vec weights))]
    (map (fn [i]
           (let [{:keys [E delta]} (nth curve i)]
             {:E-MeV E :l l :delta-rad delta :sin2-weighted (nth weights i)}))
         idxs)))

(defn unwrap-delta-sequence
  "Cumulative unwrap of principal δ(E): step increments constrained to (−π, π]."
  [raw-deltas]
  (if (empty? raw-deltas)
    []
    (loop [i 1
           out [(double (first raw-deltas))]]
      (if (>= i (count raw-deltas))
        out
        (let [yim1 (double (nth raw-deltas (dec i)))
              yi (double (nth raw-deltas i))
              d (loop [x (- yi yim1)]
                  (cond (> x Math/PI) (recur (- x (* 2.0 Math/PI)))
                        (< x (- Math/PI)) (recur (+ x (* 2.0 Math/PI)))
                        :else x))
              ui (+ (double (peek out)) d)]
          (recur (inc i) (conj out ui)))))))

(defn curve-with-unwrapped
  [curve]
  (let [raw (mapv :delta curve)
        unw (unwrap-delta-sequence raw)]
    (mapv #(assoc %1 :delta-unwrapped %2) curve unw)))

(defn wigner-slope-peaks
  "Largest |dδ_unwrap/dE| between adjacent points (avoids false 'peaks' from −π/π glitches)."
  [curve top-k]
  (let [c (curve-with-unwrapped curve)
        slopes
        (for [i (range (dec (count c)))
              :let [a (nth c i)
                    b (nth c (inc i))
                    dE (- (:E b) (:E a))]
              :when (and (pos? dE)
                         (every? #(Double/isFinite (:delta-unwrapped %)) [a b]))]
            {:E-mid (* 0.5 (+ (:E a) (:E b)))
             :dE dE
             :slope (/ (- (:delta-unwrapped b) (:delta-unwrapped a)) dE)
             :delta-a (:delta a)
             :delta-b (:delta b)})]
    (->> slopes
         (sort-by (comp - #(Math/abs (double %)) :slope))
         (take top-k))))

(defn global-sin2-peak
  "Global maximum of (2ℓ+1)·sin²δ on the grid. Includes :at-edge? when max sits on first/last point."
  [l curve]
  (when (seq curve)
    (let [ws (map-indexed (fn [idx {:keys [E delta]}]
                            (let [s (Math/sin delta)
                                  w (* (inc (* 2 l)) (* s s))]
                              {:idx idx :E E :weight w :delta delta}))
                          curve)
          best (apply max-key (fn [r] (:weight r)) ws)
          n (count curve)
          edge? (or (zero? (:idx best)) (= (:idx best) (dec n)))]
      (assoc best :at-edge? edge?))))

(defn pi-half-brackets
  "Brackets [E1 E2] where δ−π/2 changes sign (narrow resonance indicator for elastic δ)."
  [curve]
  (for [i (range (dec (count curve)))
        :let [a (nth curve i)
              b (nth curve (inc i))
              fa (- (:delta a) (* 0.5 Math/PI))
              fb (- (:delta b) (* 0.5 Math/PI))]
        :when (and (< (* fa fb) 0.0)
                   (every? #(Double/isFinite (:delta %)) [a b]))]
      [(:E a) (:E b)]))

(defn phase-shift-curve
  "Vector of {:E :delta} for E in [e-min, e-max], uniform grid."
  [l v0 r0 a0 e-min e-max n-steps h r-max r-match mass-f]
  (binding [mass-factor mass-f]
    (vec
     (for [i (range (inc n-steps))
           :let [e (+ e-min (* i (/ (- e-max e-min) (double n-steps))))
                 u (solve-numerov e l v0 r0 a0 h r-max)
                 d (phase-shift-from-numerov u h r-match e l)]
           :when (Double/isFinite d)]
       {:E e :delta d}))))

(defn phase-delta
  [E l v0 r0 a0 h r-max r-match mf]
  (binding [mass-factor mf]
    (let [u (solve-numerov E l v0 r0 a0 h r-max)]
      (phase-shift-from-numerov u h r-match E l))))

(defn first-pi-half-bracket-from-curve
  "First [E1 E2] with δ crossing π/2 (lowest energy), or nil."
  [curve]
  (first (sort-by first (pi-half-brackets curve))))

(defn refine-pi-half-crossing-E
  "Bisect δ(E) = π/2 on [e-lo e-hi] (MeV). Requires δ−π/2 to have opposite signs at endpoints."
  [l v0 r0 a0 e-lo e-hi mf h r-max r-match tol max-iter]
  (let [f (fn [e]
            (- (phase-delta e l v0 r0 a0 h r-max r-match mf) (* 0.5 Math/PI)))]
    (loop [lo (double e-lo) hi (double e-hi) i 0
           flo (f e-lo) fhi (f e-hi)]
      (cond
        (or (not (Double/isFinite flo)) (not (Double/isFinite fhi)) (>= (* flo fhi) 0.0))
        nil
        (>= i max-iter)
        (* 0.5 (+ lo hi))
        (< (- hi lo) tol)
        (* 0.5 (+ lo hi))
        :else
        (let [mid (* 0.5 (+ lo hi))
              fm (f mid)]
          (if (< (* flo fm) 0.0)
            (recur lo mid (inc i) flo fm)
            (recur mid hi (inc i) fm fhi)))))))

(defn wigner-width-from-slope
  "Γ ≈ 2/|dδ/dE| (MeV) with δ in rad, from symmetric finite difference at E_res."
  [l v0 r0 a0 E-res mf h r-max r-match dE]
  (let [d- (phase-delta (- E-res dE) l v0 r0 a0 h r-max r-match mf)
        d+ (phase-delta (+ E-res dE) l v0 r0 a0 h r-max r-match mf)
        slope (/ (- d+ d-) (* 2.0 dE))]
    (when (and (Double/isFinite slope) (pos? (Math/abs slope)))
      (/ 2.0 (Math/abs slope)))))

(defn closest-sin2-peak-to-target
  "Lowest-energy interior local max of (2ℓ+1)sin²δ on curve; nil if none."
  [l curve E-target]
  (let [peaks (resonance-candidates-from-sin2 l curve)]
    (when (seq peaks)
      (apply min-key (fn [p] (Math/abs (double (- (:E-MeV p) E-target)))) peaks))))

(defn scan-ws-for-p-resonance
  "Coarse grid on V0, R, a. Prefer lowest ℓ=1 **δ = π/2** crossing (refined); if none in [e-min,e-max],
   fall back to **closest interior (2ℓ+1)sin²δ peak** to E_target (Breit–Wigner peak proxy).
   Maps include :match :pi-half | :sin2peak and :cost = |E − E_target|."
  [{:keys [E-target Gamma-target l-scan e-min e-max n-scan h-coarse r-max v0-list r-list a-list mf
           refine-tol top-n width-dE]
    :or {l-scan 1 e-min 0.05 e-max 2.2 n-scan 56 refine-tol 2.0e-5 h-coarse 0.02 top-n 8
         width-dE 3.0e-4}}]
  (let [r-match (* 0.85 r-max)
        row-for
        (fn [v0 r0 a0 curve]
          (let [br (first-pi-half-bracket-from-curve curve)
                E-pi (when br
                       (refine-pi-half-crossing-E l-scan v0 r0 a0 (first br) (second br) mf h-coarse
                                                   r-max r-match refine-tol 80))
                sp (closest-sin2-peak-to-target l-scan curve E-target)
                pick (cond
                       E-pi {:match :pi-half :E-res E-pi
                             :gamma-est (wigner-width-from-slope l-scan v0 r0 a0 E-pi mf h-coarse r-max
                                                                  r-match width-dE)
                             :cost (Math/abs (- E-pi E-target))}
                       sp {:match :sin2peak :E-res (:E-MeV sp) :delta-at-peak (:delta-rad sp)
                           :peak-sin2wt (:sin2-weighted sp)
                           :gamma-est nil
                           :cost (Math/abs (- (:E-MeV sp) E-target))}
                       :else nil)]
            (when pick
              (merge {:v0 v0 :r0 r0 :a0 a0}
                     pick
                     {:gamma-cost (when (and (:gamma-est pick) Gamma-target)
                                     (Math/abs (- (:gamma-est pick) Gamma-target)))}))))
        candidates
        (for [v0 v0-list
              r0 r-list
              a0 a-list
              :let [curve (phase-shift-curve l-scan v0 r0 a0 e-min e-max n-scan h-coarse r-max r-match mf)
                    row (row-for v0 r0 a0 curve)]
              :when row]
          row)]
    (->> candidates
         (sort-by :cost)
         (take top-n)
         vec)))

(defn -main [& _args]
  (let [mu (reduced-mass-9Li-n)
        mf (mass-factor-from-mu mu)
        ;; Schematic volume WS for n–⁹Li relative motion (tune for your interaction).
        v0 50.0
        r0 3.0
        a0 0.65
        v-params [v0 r0 a0]
        r-max 30.0
        h 0.02
        r-match (* 0.85 r-max)
        l-min 1
        l-max 4
        e-pos-max 4.0
        pos-steps 120
        ;; Finer low-E grid for ℓ = 1 (¹⁰Li physics is near threshold).
        fine-e-min 0.005
        fine-e-max 2.5
        fine-steps 220
        ;; Shallow explorations of V0 for p-wave (shows resonance strength is interaction-dependent).
        v0-sensitivity [50.0 54.0 58.0 62.0]]
    (println "=== ⁹Li + n — bound states (E<0) and resonance candidates (phase-shift scan) ===")
    (println)
    (println (format "Reduced mass μ = %.4f MeV/c²  →  mass-factor 2μ/ℏ² = %.6f MeV⁻¹·fm⁻²"
                     mu mf))
    (println (format "Woods–Saxon: V0 = %.1f MeV, R = %.3f fm, a = %.3f fm"
                     v0 r0 a0))
    (println (format "Numerov: h = %.3f fm, r_max = %.1f fm, δ extracted at r = %.2f fm"
                     h r-max r-match))
    (println)
    (println "Bound-state search (ℓ ≥ 1 — Pauli/shell context in file header):")
    (doseq [l (range l-min (inc l-max))
            n [1 2]]
      (binding [mass-factor mf]
        (let [res (tr/find-bound-state-energy v-params l n r-max h)]
          (if (and res (every? some? [(:energy res) (:boundary-value res)])
                  (< (Math/abs (double (:boundary-value res))) 0.05))
            (println (format "  ℓ=%d  n=%d  E = %.4f MeV  nodes=%d  |u_end|=%.2e"
                             l n (:energy res) (:nodes res) (:boundary-value res)))
            (println (format "  ℓ=%d  n=%d  — no converged bound state in default search"
                             l n))))))
    (println)
    (println "Hypothetical — **no Pauli block**, include **ℓ = 0** (same V0, R, a, μ):")
    (doseq [n [1 2]]
      (binding [mass-factor mf]
        (let [res (tr/find-bound-state-energy v-params 0 n r-max h)]
          (if (and res (every? some? [(:energy res) (:boundary-value res)])
                  (< (Math/abs (double (:boundary-value res))) 0.05))
            (println (format "  ℓ=0  n=%d  E = %.4f MeV  nodes=%d  |u_end|=%.2e"
                             n (:energy res) (:nodes res) (:boundary-value res)))
            (do
              (println (format "  ℓ=0  n=%d  — no physical bound state (|u_end| not small)." n))
              (when res
                (println (format "    (search returned E ≈ %.3f MeV, |u_end| ≈ %.3g — not normalizable on this grid.)"
                                 (double (:energy res))
                                 (Math/abs (double (:boundary-value res)))))))))))
    (println)
    (println "Positive-energy scan (0.02–4 MeV, 120 steps): global sin²δ max, π/2 crossings, |dδ/dE| (unwrapped).")
    (println "  ‘No resonance’ here means: no δ≈π/2 crossing and no narrow sin²δ peak — expected if the well is")
    (println "  too shallow so |δ| stays ≪ π/2. Deeper V0 or a fitted optical potential can move δ through π/2.")
    (println)
    (doseq [l (range l-min (inc l-max))]
      (let [curve (phase-shift-curve l v0 r0 a0 0.02 e-pos-max pos-steps h r-max r-match mf)
            peaks (resonance-candidates-from-sin2 l curve)
            wigs (wigner-slope-peaks curve 4)
            gmax (global-sin2-peak l curve)
            x-half (pi-half-brackets curve)]
        (println (format "  ℓ = %d  (%d points):" l (count curve)))
        (when gmax
          (println (format "     global max (2ℓ+1)sin²δ = %.4f at E = %.4f MeV  (δ = %.4f rad)%s%s"
                           (:weight gmax) (:E gmax) (:delta gmax)
                           (if (:at-edge? gmax) "  [max at grid endpoint — not a localized resonance]" "")
                           (if (and (not (:at-edge? gmax)) (< (:weight gmax) 0.25))
                             "  → weak phase shift"
                             ""))))
        (if (seq x-half)
          (doseq [[e1 e2] x-half]
            (println (format "     δ crosses π/2 between E = %.4f and %.4f MeV (narrow-resonance indicator)"
                             e1 e2)))
          (println "     no bracket where δ crosses π/2 on this grid"))
        (if (seq peaks)
          (doseq [p (sort-by :E-MeV peaks)]
            (println (format "     interior sin²δ peak: E ≈ %.3f MeV  δ = %.3f rad  (2ℓ+1)sin²δ ≈ %.4f"
                             (:E-MeV p) (:delta-rad p) (:sin2-weighted p))))
          (println "     no interior local max of sin²δ (often: sin²δ rises/falls monotonically toward a grid edge)"))
        (when (seq wigs)
          (println "     largest |dδ_unwrap/dE| (rad/MeV):")
          (doseq [w wigs]
            (println (format "       E_mid ≈ %.3f  |dδ/dE| ≈ %.4f  (δ_wrap: %.3f → %.3f)"
                             (:E-mid w) (Math/abs (double (:slope w))) (:delta-a w) (:delta-b w)))))
        (when (= l l-min)
          (println (format "     δ(E) sample (same grid): %s"
                           (pr-str (for [e [0.1 0.5 1.0 2.0]]
                                     (let [u (binding [mass-factor mf]
                                               (solve-numerov e l v0 r0 a0 h r-max))
                                           d (binding [mass-factor mf]
                                               (phase-shift-from-numerov u h r-match e l))]
                                       [e d]))))))))
    (println)
    (println (format "ℓ = 1 finer grid [%.3f, %.3f] MeV (%d steps), same V0 = %.1f MeV:"
                     fine-e-min fine-e-max fine-steps v0))
    (let [curve-f (phase-shift-curve 1 v0 r0 a0 fine-e-min fine-e-max fine-steps h r-max r-match mf)
          gmax (global-sin2-peak 1 curve-f)
          x-half (pi-half-brackets curve-f)]
      (when gmax
        (println (format "  global max (2ℓ+1)sin²δ = %.4f at E = %.4f MeV  (δ = %.5f rad)%s"
                         (:weight gmax) (:E gmax) (:delta gmax)
                         (if (:at-edge? gmax) "  [at E_max — extend fine-e-max if probing higher E]" ""))))
      (if (seq x-half)
        (doseq [[e1 e2] x-half]
          (println (format "  δ crosses π/2 between E = %.4f and %.4f MeV" e1 e2)))
        (println "  still no π/2 crossing at this V0 — potential is too weak for a Wigner-style p resonance")))
    (println)
    (println "ℓ = 1: V0 sensitivity (fine grid, same R, a) — max (2ℓ+1)sin²δ and E at maximum:")
    (doseq [v0-try v0-sensitivity]
      (let [cur (phase-shift-curve 1 v0-try r0 a0 fine-e-min fine-e-max fine-steps h r-max r-match mf)
            g (global-sin2-peak 1 cur)
            x (pi-half-brackets cur)]
        (when g
          (println (format (str "  V0 = %5.1f MeV → max (2ℓ+1)sin²δ = %.4f at E = %.4f MeV"
                              "  |  π/2 brackets: %s%s")
                           v0-try (:weight g) (:E g)
                           (if (seq x) (str (count x)) "0")
                           (if (:at-edge? g) "  [E at grid edge]" ""))))))
    (println)
    (println "=== Schematic fit: scan Woods–Saxon (ℓ=1) near ¹⁰Li E_r = 0.62 MeV ===")
    (println "  Prefer δ=π/2 crossing (refined); else closest interior (2ℓ+1)sin²δ peak. Γ_est = 2/|dδ/dE| when π/2.")
    (let [v0-grid (mapv #(+ 28.0 (* 2.0 (double %))) (range 0 25))
          r-grid (mapv #(+ 2.35 (* 0.18 (double %))) (range 0 11))
          hits (scan-ws-for-p-resonance
                {:E-target 0.62 :Gamma-target 0.33 :mf mf :r-max r-max :v0-list v0-grid
                 :r-list r-grid :a-list [a0] :n-scan 64 :e-min 0.05 :e-max 5.5 :h-coarse 0.02
                 :refine-tol 1.5e-5 :top-n 10 :width-dE 2.5e-4})]
      (println (format "  Grid: V0 ∈ [%.1f, %.1f] MeV (%d), R ∈ [%.2f, %.2f] fm (%d), a = %.2f fm; E scan [%.2f, %.2f] MeV"
                       (first v0-grid) (peek v0-grid) (count v0-grid)
                       (first r-grid) (peek r-grid) (count r-grid)
                       a0 0.05 5.5))
      (if (empty? hits)
        (println "  No π/2 crossing and no interior sin²δ peak — widen :v0-list / :r-list / :e-max in -main.")
        (doseq [hrow hits]
          (println (format (str "  [%s] V0=%5.1f  R=%.2f  a=%.2f  →  E=%.4f MeV  ΔE=%.4f  "
                                "Γ_est=%s  ΔΓ=%s")
                           (name (:match hrow))
                           (:v0 hrow) (:r0 hrow) (:a0 hrow) (:E-res hrow) (:cost hrow)
                           (if-let [g (:gamma-est hrow)] (format "%.3f" g) "—")
                           (if (:gamma-cost hrow) (format "%.3f" (:gamma-cost hrow)) "—"))))))
    (println)
    (println "Done.")
    nil))

(-main)
