;; Example: 16O(p,d) Transfer Reaction
;; 
;; This example calculates the form factor, transfer amplitude (post formulation),
;; and differential cross section for the 16O(p,d) reaction.
;;
;; Reaction: 16O(p,d)15O
;; - Initial: neutron bound in 16O (l=1, E=-15.67 MeV)
;; - Final: neutron bound in deuteron (l=0, E=-2.214 MeV)
;;
;; Reference parameters:
;; - 16O: R0 = 2.7 fm, V0 = 62 MeV, a0 = 0.6 fm, l=1
;; - Deuteron: R0 = 1.5 fm, V0 = 50 MeV, a0 = 0.6 fm, l=0
;;
;; Sister example (stripping): `examples/example_16Odp.clj` — 16O(d,p)17O (**handbook** namespace `o16-dp-handbook`).
;;
;; **Web dashboard 16O(p,d):** `GET /api/transfer-default?target=16O&reaction_type=p-d` runs the same pipeline as
;; this script (`transfer-default-data` in `dwba_web/simple_core.clj`: post DWBA, L = 0…7, h = 0.01 fm, 20 MeV lab).
;; **`POST /api/transfer`** with `reaction_type: p-d` is a **different** schematic (zero-range overlap from the
;; shared WS form, global `functions/mass-factor`, `E_f ≈ 0.8 E_lab`, one partial wave) — numbers will not match below.
;;
;; Use (load-file "examples/example_16Opd.clj") or run from project root.
;; Own namespace avoids alias conflicts with dwba.core when loaded from REPL.
;; Plots: **output/16Opd_dcs.png** (linear y), **output/16Opd_dcs_log.png** (log₁₀ y).

(ns examples.example-16Opd
  (:require [dwba.transfer :as t]
            [dwba.form-factors :as ff]
            [dwba.inelastic :as inel]
            [functions :refer [solve-numerov mass-factor]]
            [fastmath.core :as m]
            [complex :as cx :refer [mag re im complex-cartesian add mul]]
            [incanter.core :as i]
            [incanter.charts :as c]
            [clojure.java.io :as io])
  (:import [org.jfree.chart.axis LogAxis]))

(defn- clamp-pos-for-log ^double [^double y]
  (if (or (Double/isNaN y) (Double/isInfinite y) (<= y 0.0))
    1e-30
    (max y 1e-30)))

(defn- series-for-log [ys]
  (mapv clamp-pos-for-log ys))

(defn- ->complex-16Opd
  "Coerce a real or complex T-amplitude to complex (mirrors dwba.transfer's private transfer-T->complex)."
  [T-L]
  (if (number? T-L)
    (complex-cartesian (double T-L) 0.0)
    T-L))

(defn- chart-set-log-range-y!
  [chart ^String y-label]
  (let [plot (.getPlot chart)
        axis (LogAxis. y-label)]
    (.setSmallestValue axis 1e-35)
    (.setRangeAxis plot 0 axis)
    chart))

(println "=== 16O(p,d) Transfer Reaction Calculation ===")
(println "")

;; ============================================================================
;; Bound State Wavefunctions
;; ============================================================================

(println "=== Step 1: Bound State Wavefunctions ===")

(let [;; Bound state geometry (well depths are *searched*, not fixed, so that each
      ;; wavefunction is a genuine eigenstate at its target separation energy below).
      R0-i 2.7 diff-i 0.6 r-max 20.0 h 0.01 l-i 1
      R0-f 1.5 diff-f 0.6 l-f 0
      Es-i -15.67  ; Neutron bound in 16O
      Es-f -2.214  ; Neutron bound in deuteron

      ;; Correct per-channel reduced-mass factor 2μ/ħ² (MeV⁻¹·fm⁻²) — NOT a shared
      ;; "typical" placeholder: phi-i is n+16O, phi-f is n+p (the deuteron).
      m-n 939.565
      m-p-bs 938.272
      m-16O-bs (* 16.0 931.494)
      m-f-i (/ (* 2.0 (/ (* m-n m-16O-bs) (+ m-n m-16O-bs))) (* 197.327 197.327))
      m-f-f (/ (* 2.0 (/ (* m-n m-p-bs) (+ m-n m-p-bs))) (* 197.327 197.327))

      ;; Bisect V0 so the Numerov eigenvalue at (l-i, m-f-i) matches Es-i exactly
      ;; (mirrors the search pattern in dwba.benchmark.o16-dp-handbook).
      find-v0 (fn [l mf Es v0-lo v0-hi]
                (let [eigen (fn [v0]
                              (binding [mass-factor mf]
                                (:energy (t/find-bound-state-energy [v0 R0-i diff-i] l r-max h
                                                                     (- v0 5.0) -0.05 0.005 nil))))]
                  (loop [a v0-lo b v0-hi n 0]
                    (let [mid (/ (+ a b) 2.0)
                          e (or (eigen mid) -0.001)]
                      (if (or (> n 30) (< (Math/abs (- e Es)) 0.005))
                        mid
                        (if (> e Es) (recur mid b (inc n)) (recur a mid (inc n))))))))
      v0-i (find-v0 l-i m-f-i Es-i 40.0 90.0)
      v0-f (find-v0 l-f m-f-f Es-f 30.0 70.0)

      ;; Calculate bound state wavefunctions at the calibrated depths
      phi-i-raw (t/solve-bound-state-numerov Es-i l-i v0-i R0-i diff-i m-f-i h r-max)
      phi-f-raw (t/solve-bound-state-numerov Es-f l-f v0-f R0-f diff-f m-f-f h r-max)
      ;; Normalize bound state wavefunctions
      phi-i (t/normalize-bound-state phi-i-raw h)
      phi-f (t/normalize-bound-state phi-f-raw h)

      ;; Calculate normalized overlap (form factor)
      overlap-norm (ff/normalized-overlap phi-i phi-f r-max h)
      
      ;; ============================================================================
      ;; Reaction Parameters
      ;; ============================================================================
      
      ;; Incident energy (typical for (p,d) reactions: 10-20 MeV)
      E-lab 20.0  ; Lab frame energy (MeV)
      
      ;; Masses (in MeV/c²)
      m-p  938.27     ; Proton mass
      m-16O 14899.0   ; 16O mass (approximate)
      m-d 1876.136      ; Deuteron mass
      m-15O 13975.0   ; 15O mass (approximate)
      
      ;; Reduced masses
      mu-i (/ (* m-p m-16O) (+ m-p m-16O))      ; Entrance channel: p+16O
      mu-f (/ (* m-d m-15O) (+ m-d m-15O))      ; Exit channel: d+15O
      
      ;; Mass factors
      mass-factor-i (/ (* 2.0 mu-i) (* 197.7 197.7))  ; 2μ/ħ² for entrance
      mass-factor-f (/ (* 2.0 mu-f) (* 197.7 197.7))  ; 2μ/ħ² for exit
      
      ;; CM frame energies
      E-CM-i (* E-lab (/ m-16O (+ m-16O m-p)))  ; Entrance CM energy
      ;; Q-value for 16O(p,d)15O reaction
      ;; Q = (m_p + m_16O - m_d - m_15O) * c²
      ;; Or using binding energies: Q = |Es-i| - |Es-f|
      ;; where Es-i is the neutron separation energy from 16O (negative, bound)
      ;; and Es-f is the deuteron binding energy (negative, bound)
      ;; For (p,d) pickup: Q = B_n(16O) - B_d = |Es-i| - |Es-f|
      ;; Note: Q is typically NEGATIVE for pickup reactions (endothermic)
      ;; The reaction is kinematically allowed if E_CM_i + Q > 0
      Q-value  (+ m-p m-16O (- m-d) (- m-15O))
      E-CM-f (+ E-CM-i Q-value)  ; Exit CM energy = E_i + Q (Q can be negative for endothermic reactions)
      
      ;; ============================================================================
      ;; Distorted Waves and Transfer Amplitudes (multiple L for angular dependence)
      ;; ============================================================================
      
      ;; Lab energies for optical potential calculation
      E-lab-i E-lab  ; Entrance channel lab energy
      E-lab-f (* E-CM-f (/ (+ m-d m-15O) m-15O))  ; Exit channel lab energy (approximate)
      
      ;; Include partial waves L = 0 .. L-max so dσ/dΩ depends on angle (e.g. drop at 70°)
      L-max 7
      D0 (t/zero-range-constant :p-d)  ; D₀ for (p,d) reaction
      ;; Build T-amplitudes map {L → T_L} for L = 0, 1, ..., L-max
      T-amplitudes (into {}
                        (for [L (range (inc L-max))]
                          (let [chi-i (inel/distorted-wave-entrance E-CM-i L nil h r-max
                                                                    :projectile-type :p
                                                                    :target-A 16
                                                                    :target-Z 8
                                                                    :E-lab E-lab-i
                                                                    :s 0.5
                                                                    :j (+ L 0.5)
                                                                    :mass-factor mass-factor-i)
                                chi-f (inel/distorted-wave-exit E-CM-i Q-value L nil h r-max
                                                                :outgoing-type :d
                                                                :residual-A 15
                                                                :residual-Z 8
                                                                :E-lab E-lab-f
                                                                :s 1
                                                                :j (inc L)
                                                                :mass-factor mass-factor-f)
                                T-L (t/transfer-amplitude-post chi-i chi-f phi-i phi-f r-max h
                                                               :zero-range D0)]
                            [L T-L])))
      
      ;; Keep L=0 waves for display (Step 3)
      L-i 0
      L-f 0
      chi-i (inel/distorted-wave-entrance E-CM-i L-i nil h r-max
                                          :projectile-type :p
                                          :target-A 16
                                          :target-Z 8
                                          :E-lab E-lab-i
                                          :s 0.5
                                          :j 0.5
                                          :mass-factor mass-factor-i)
      chi-f (inel/distorted-wave-exit E-CM-i Q-value L-f nil h r-max
                                      :outgoing-type :d
                                      :residual-A 15
                                      :residual-Z 8
                                      :E-lab E-lab-f
                                      :s 1
                                      :j 1
                                      :mass-factor mass-factor-f)
      T-post (get T-amplitudes 0)  ; Same as T for L=0
      
      ;; ============================================================================
      ;; Differential Cross Section
      ;; ============================================================================
      
      ;; Wavenumbers
      k-i (Math/sqrt (* mass-factor-i E-CM-i))
      k-f (Math/sqrt (* mass-factor-f E-CM-f))

      ;; ============================================================================
      ;; Austern (5.5) amplitude normalization (was MISSING — root cause of the
      ;; ~6000x undershoot vs experiment).
      ;;
      ;; `transfer-amplitude-post` returns the *bare* ZR radial overlap
      ;; T_raw = D0 * ∫ χ_f* F χ_i r² dr — this is the same "I"-type integral that
      ;; `dwba.benchmark.o16-dp-handbook` / `ca40-pd-handbook` compute via
      ;; `handbook-radial-integral-I-zr-from-neutron-bound`, BEFORE they multiply by
      ;; N. Austern's Eq. (5.5) prefactor (M_B/M_A)(4π/(k_α k_β)) and the
      ;; √(2ℓ+1) angular-momentum-coupling factor from T_m = D0·√(2ℓ+1)·β_m
      ;; (see `dwba.transfer` around `austern-radial-integral-prefactor-eq-5-5`,
      ;; whose docstring explicitly says: "I ... must use the same (5.5) prefactor
      ;; ... (LNPS / Austern angular reduction is derived for that I)").
      ;;
      ;; `transfer-differential-cross-section(-angular)` is that same downstream
      ;; formula, so it must be fed a T that already includes this normalization —
      ;; example_16Opd.clj was passing it the bare T_raw instead.
      ;;
      ;; M_A = M(16O) (target, no transferred nucleon in this :handbook-F-from :phi-i
      ;; convention), M_B = M(15O) (residual) — same (M-target, M-residual) labeling
      ;; used by o16-dp-handbook / ca40-pd-handbook.
      austern-P (t/austern-radial-integral-prefactor-eq-5-5 m-16O m-15O k-i k-f)
      ell-transfer (long l-i)  ; ℓ of the transferred neutron in 16O (:handbook-F-from :phi-i default)
      amp-scale (* austern-P (Math/sqrt (inc (* 2.0 (double ell-transfer)))))
      T-amplitudes-norm (into {}
                              (map (fn [[L T-L]]
                                     [L (mul (complex-cartesian amp-scale 0.0) (->complex-16Opd T-L))])
                                   T-amplitudes))
      T-post-norm (get T-amplitudes-norm 0)

      ;; Nuclear-spin statistical factor (same pattern as o16-dp-handbook, applied to
      ;; dσ/dΩ, not to T): (2J_f+1)/(2J_i+1) for J(16O)=0, J(15O g.s.)=1/2, times the
      ;; unpolarized-proton entrance-spin average 1/(2·1/2+1). For this transition it
      ;; evaluates to 2 * 1/2 = 1.0 (a no-op numerically) but is kept explicit/correct
      ;; rather than silently assumed.
      J-16O 0.0
      J-15O 0.5
      s-proton 0.5
      spin-factor (* (t/transfer-nuclear-spin-statistical-factor J-16O J-15O)
                     (/ 1.0 (inc (* 2.0 s-proton))))

      ;; Spectroscopic factor (typically 0 < S < 1, using 1.0 for this example)
      S-factor 1.0
      theta-deg 0.0  ; Forward angle (degrees)
      theta-rad (* theta-deg (/ Math/PI 180.0))
      ;; Pass l-i, l-f so only allowed L contribute (L selection: triangle + parity).
      ;; **Note:** This uses **|Y_L0|²** per L → CM shape is **σ(θ)=σ(180°−θ)** for parity-selected L.
      ;; Asymmetry needs full **m**-sum: **`transfer-differential-cross-section-angular-m-sum`** (meaningful
      ;; when several L contribute and **l_i≥1**; for **l_i=1,l_f=0** only **L=1** is allowed → still symmetric).
      dsigma (* spin-factor
                (t/transfer-differential-cross-section-angular T-amplitudes-norm S-factor k-i k-f
                                                             theta-rad mass-factor-i mass-factor-f 0.0 l-i l-f))
      ;; At 70° (compare to experiment: ~0.5 mb/sr at 20 MeV)
      theta-70-deg 70.0
      theta-70-rad (* theta-70-deg (/ Math/PI 180.0))
      dsigma-70 (* spin-factor
                   (t/transfer-differential-cross-section-angular T-amplitudes-norm S-factor k-i k-f
                                                                 theta-70-rad mass-factor-i mass-factor-f 0.0 l-i l-f))
      ;; Angles and DCS for plotting (0° to 180°, step 5°), dσ/dΩ in mb/sr
      angles-deg (range 0.0 181.0 5.0)
      angles-rad (mapv #(* % (/ Math/PI 180.0)) angles-deg)
      dsigma-mb-sr (mapv (fn [theta-rad]
                           (* spin-factor
                              (t/transfer-differential-cross-section-angular T-amplitudes-norm S-factor k-i k-f
                                                                          theta-rad mass-factor-i mass-factor-f 0.0 l-i l-f)))
                         angles-rad)
      ;; Same angle set as web dashboard: 20° to 160° step 20° (CM)
      angles-output-deg (range 20.0 181.0 20.0)
      dsigma-at-angles (mapv (fn [theta-deg]
                               (let [theta-rad (* theta-deg (/ Math/PI 180.0))]
                                 (* spin-factor
                                    (t/transfer-differential-cross-section-angular T-amplitudes-norm S-factor k-i k-f
                                                                                theta-rad mass-factor-i mass-factor-f 0.0 l-i l-f))))
                             angles-output-deg)]
  
  ;; ============================================================================
  ;; Output
  ;; ============================================================================
  
  (println (format "Initial state (neutron in 16O):"))
  (println (format "  Energy: E_i = %.2f MeV" Es-i))
  (println (format "  Angular momentum: l_i = %d" l-i))
  (println (format "  Wavefunction length: %d points" (count phi-i)))
  (println "")
  
  (println (format "Final state (neutron in deuteron):"))
  (println (format "  Energy: E_f = %.2f MeV" Es-f))
  (println (format "  Angular momentum: l_f = %d" l-f))
  (println (format "  Wavefunction length: %d points" (count phi-f)))
  (println "")
  
  (println (format "Normalized overlap (form factor): %.6f" overlap-norm))
  (println "")
  
  (println "=== Step 2: Reaction Parameters ===")
  (println (format "Incident energy: E_lab = %.2f MeV" E-lab))
  (println (format "Entrance CM energy: E_CM_i = %.2f MeV" E-CM-i))
  (println (format "Exit CM energy: E_CM_f = %.2f MeV" E-CM-f))
  (println (format "Q-value: Q = %.2f MeV" Q-value))
  (println "")
  
  (println "=== Step 3: Distorted Waves ===")
  (println (format "Entrance channel distorted wave:"))
  (println (format "  Energy: E_i = %.2f MeV (CM)" E-CM-i))
  (println (format "  Angular momentum: L_i = %d" L-i))
  (println (format "  Wavefunction length: %d points" (count chi-i)))
  (println "")
  
  (println (format "Exit channel distorted wave:"))
  (println (format "  Energy: E_f = %.2f MeV (CM)" E-CM-f))
  (println (format "  Angular momentum: L_f = %d" L-f))
  (println (format "  Wavefunction length: %d points" (count chi-f)))
  ;; Check magnitude of distorted waves
  (let [chi-i-max (apply max (map #(if (number? %) (Math/abs %) (mag %)) chi-i))
        chi-f-max (apply max (map #(if (number? %) (Math/abs %) (mag %)) chi-f))
        chi-i-avg (let [sum (reduce + (map #(if (number? %) (Math/abs %) (mag %)) chi-i))]
                    (/ sum (count chi-i)))
        chi-f-avg (let [sum (reduce + (map #(if (number? %) (Math/abs %) (mag %)) chi-f))]
                    (/ sum (count chi-f)))
        ;; Calculate integral of chi-i* · chi-f
        n-chi (min (count chi-i) (count chi-f))
        integrand-chi (mapv (fn [i]
                              (let [r (* i h)
                                    chi-i-val (get chi-i i)
                                    chi-f-val (get chi-f i)
                                    chi-i-conj (if (number? chi-i-val)
                                                chi-i-val
                                                (complex-cartesian (re chi-i-val) (- (im chi-i-val))))
                                    product (if (and (number? chi-i-conj) (number? chi-f-val))
                                             (* chi-i-conj chi-f-val)
                                             (mul chi-i-conj chi-f-val))]
                                product))
                            (range n-chi))
        simpson-sum-chi (loop [i 1 sum (complex-cartesian 0.0 0.0)]
                         (if (>= i (dec n-chi))
                           sum
                           (let [coeff (if (odd? i) 4.0 2.0)
                                 term-val (get integrand-chi i)
                                 coeff-complex (complex-cartesian coeff 0.0)
                                 term (mul coeff-complex term-val)]
                             (recur (inc i) (add sum term)))))
        first-term-chi (get integrand-chi 0)
        last-term-chi (get integrand-chi (dec n-chi))
        h-over-3-complex (complex-cartesian (/ h 3.0) 0.0)
        integral-chi (mul h-over-3-complex
                         (add first-term-chi last-term-chi simpson-sum-chi))]
    (println (format "  chi-i max magnitude: %.6e" chi-i-max))
    (println (format "  chi-f max magnitude: %.6e" chi-f-max))
    (println (format "  chi-i avg magnitude: %.6e" chi-i-avg))
    (println (format "  chi-f avg magnitude: %.6e" chi-f-avg))
    (println (format "  ∫ χ*_i · χ_f dr = %.6e + i%.6e" (re integral-chi) (im integral-chi)))
    (println (format "  |∫ χ*_i · χ_f dr| = %.6e" (mag integral-chi))))
  (println "")
  
  (println "=== Step 4: Transfer Amplitudes (Post Formulation) ===")
  (println (format "Zero-range constant: D₀ = %.2f MeV·fm^(3/2)" D0))
  (println (format "Partial waves included: L = 0 .. %d" L-max))
  (println "")
  (println (format "Austern (5.5) amplitude normalization: P = (M_B/M_A)(4π/(k_α k_β)) = %.4f, √(2ℓ+1) (ℓ=%d) = %.4f"
                  austern-P ell-transfer (Math/sqrt (inc (* 2.0 (double ell-transfer))))))
  (println (format "  → T-amplitude scale factor = %.4f (applied to every T_L below)" amp-scale))
  (println (format "Nuclear-spin statistical factor (applied to dσ/dΩ): %.4f" spin-factor))
  (println "")
  (println "Transfer amplitude L=0 (post formulation, bare ZR integral before Austern normalization):")
  (if (number? T-post)
    (println (format "  T_0 = %.6e" T-post))
    (println (format "  T_0 = %.6e + i%.6e" (re T-post) (im T-post))))
  (println (format "  |T_0| (bare) = %.6e" (if (number? T-post) (Math/abs T-post) (mag T-post))))
  (println (format "  |T_0| (Austern-normalized) = %.6e" (mag T-post-norm)))
  (println "")
  (println "|T_L| for all L (bare vs Austern-normalized; normalized values used in angular distribution):")
  (doseq [L (range (inc L-max))]
    (let [T-L (get T-amplitudes L)
          T-L-norm (get T-amplitudes-norm L)
          mag-T (if (number? T-L) (Math/abs T-L) (mag T-L))
          mag-T-norm (mag T-L-norm)]
      (println (format "  L = %2d: |T_L| bare = %.6e, normalized = %.6e" L mag-T mag-T-norm))))
  (println "")
  
  (println "=== Step 5: Differential Cross Section ===")
  (println (format "Wavenumbers:"))
  (println (format "  k_i (entrance): %.4f fm⁻¹" k-i))
  (println (format "  k_f (exit): %.4f fm⁻¹" k-f))
  (println (format "  Ratio k_f/k_i: %.4f" (/ k-f k-i)))
  (println (format "Spectroscopic factor: S = %.2f" S-factor))
  (println (format "Scattering angle: θ = %.1f° (%.4f rad)" theta-deg theta-rad))
  (println "")
  
  (println "Differential cross section (mb/sr):")
  (println (format "  dσ/dΩ(θ=%.1f°) = %.6e mb/sr" theta-deg dsigma))
  (println (format "  dσ/dΩ(θ=%.1f°) = %.6e mb/sr (experiment at 20 MeV, 70°: ~0.5 mb/sr)" theta-70-deg dsigma-70))
  (println "")
  (println "DCS at different angles (CM, same as web dashboard; mb/sr):")
  (println "  Angle (deg) | dσ/dΩ (mb/sr)")
  (doseq [[theta-deg dsigma-mb] (map vector angles-output-deg dsigma-at-angles)]
    (println (format "  θ = %5.1f°   | %.6e" theta-deg dsigma-mb)))
  (println "")

  (println "=== Summary ===")
  (println (format "Reaction: 16O(p,d)15O"))
  (println (format "Incident energy: E_lab = %.2f MeV" E-lab))
  (println (format "Normalized overlap: %.6f" overlap-norm))
  (println (format "Transfer amplitude: |T_post| (bare) = %.6e, (Austern-normalized, used below) = %.6e"
                  (if (number? T-post) (Math/abs T-post) (mag T-post))
                  (mag T-post-norm)))
  (println (format "Differential cross section: dσ/dΩ(0°) = %.6e mb/sr" dsigma))
  (println (format "Differential cross section: dσ/dΩ(70°) = %.6e mb/sr" dsigma-70))
  (println "")
  (println "=== Calculation Complete ===")

  ;; Plot DCS vs angle and save to file (linear + log y), plus a TSV of the
  ;; plotted curve so external plotting (e.g. the paper figures) can reuse it.
  (try
    (let [_ (io/make-parents (io/file "output/16Opd_dcs.png"))
          th (vec angles-deg)
          ys (vec dsigma-mb-sr)
          _ (do (spit "output/16Opd_dcs.tsv"
                      (apply str "# theta_deg\tdsigma_mb_sr\n"
                             (map (fn [t y] (format "%.1f\t%.10e\n" (double t) (double y))) th ys)))
                (println "Data saved: output/16Opd_dcs.tsv"))
          chart (c/xy-plot th ys
                           :title "16O(p,d)15O — Differential cross section"
                           :x-label "θ (deg)"
                           :y-label "dσ/dΩ (mb/sr)"
                           :series-label (format "E_lab = %.1f MeV" E-lab)
                           :legend true)]
      (i/save chart "output/16Opd_dcs.png" :width 800 :height 500)
      (println "Plot saved: output/16Opd_dcs.png")
      (let [chart-log (-> (c/xy-plot th (series-for-log ys)
                                     :title "16O(p,d)15O — Differential cross section (log y)"
                                     :x-label "θ (deg)"
                                     :y-label "dσ/dΩ (mb/sr), logarithmic"
                                     :series-label (format "E_lab = %.1f MeV" E-lab)
                                     :legend true)
                         (chart-set-log-range-y! "dσ/dΩ (mb/sr)"))]
        (i/save chart-log "output/16Opd_dcs_log.png" :width 800 :height 500)
        (println "Plot saved: output/16Opd_dcs_log.png (log y-axis)")))
    (catch Exception e
      (println (format "Note: Could not save plot (%s). Create output/ directory or check Incanter." (.getMessage e)))))

  ;; Return the differential cross section
  dsigma
  )

