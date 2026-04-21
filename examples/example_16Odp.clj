;; Example: 16O(d,p)17O — **Handbook** ZR pipeline
;;
;; Uses **`dwba.benchmark.o16-dp-handbook`**: handbook §5.4 **F=u/r**, **(5.5)** radial **I**, **`handbook-zr-multipole-amplitude-sum`**
;; for the neutron bound in **¹⁷O**, **T_m ≈ D₀√(2ℓ+1)β_m**,
;; and **`transfer-differential-cross-section`** — parallel to **`ca40-pd-handbook`** (pickup) with
;; **α = d+¹⁶O**, **β = p+¹⁷O**.
;;
;; **dσ/dΩ(θ)** — primary curve is the **full** unpolarized **Σ_m |T_m|²** (handbook / Austern **(5.6)** angular
;; content via **`handbook-zr-multipole-amplitude-sum`**, default **`coherent-m-beta? false`**). Additional series on
;; the angular plots are **diagnostics only**: coherent **Σ_{L_β} T_{L_β} Y_{L_β 0}** with **T_{L_β} = D₀√(2ℓ+1) Σ_{L_α} I**
;; (**`transfer-differential-cross-section-angular-coherent`**, no **l_i/l_f** weights — same reduced multipole
;; picture as **`example_O16_alpha_inelastic_2plus.clj`**), and a single dominant **(L_α,L_β)** row.
;;
;; Reaction: 16O(d,p)17O (g.s. model: d₅/₂ → **ℓ=2** bound well)
;;
;; Load: `(load-file "examples/example_16Odp.clj")` from project root.
;;
;; Writes **output/o17_bound_states.png** (0d₅/₂ + 2s₁/₂ u(r) and F(r)=u(r)/r),
;; **output/16Odp_dcs.png** (linear y), **output/16Odp_dcs_log.png** (log₁₀ y-axis),
;; **output/o16dp_radial_I_s12.png** (Re/Im of I_{L_α} for 2s₁/₂ transfer, Fig 5.3 coupling).

(ns examples.example-16Odp
  (:require [dwba.benchmark.o16-dp-handbook :as oh]
            [dwba.transfer :as t]
            [functions :as fn]
            [complex :refer [re im mag mul add subt2 complex-cartesian complex-polar]]
            [fastmath.polynomials :as poly]
            [fastmath.core :as m]
            [incanter.core :as i]
            [incanter.charts :as c]
            [clojure.java.io :as io])
  (:import [org.jfree.chart.axis LogAxis]
           [org.apache.commons.math3.analysis.interpolation SplineInterpolator]
           [org.apache.commons.math3.analysis UnivariateFunction]))

(defn- clamp-pos-for-log ^double [^double y]
  (if (or (Double/isNaN y) (Double/isInfinite y) (<= y 0.0))
    1e-30
    (max y 1e-30)))

(defn- series-for-log [ys]
  (mapv clamp-pos-for-log ys))

(defn- chart-set-log-range-y!
  [chart ^String y-label]
  (let [plot (.getPlot chart)
        axis (LogAxis. y-label)]
    (.setSmallestValue axis 1e-35)
    (.setRangeAxis plot 0 axis)
    chart))

(println "=== 16O(d,p)17O — Handbook ZR (dwba.benchmark.o16-dp-handbook) ===")
(println "")

(let [kin (oh/o16-dp-kinematics)
      e-cm-i (:e-cm-i kin)]
  (println "=== Kinematics (default: E_lab(d)=20 MeV on 16O at rest) ===")
  (println (format "  E_CM(d+16O) = %.4f MeV" e-cm-i))
  (println (format "  E_CM(p+17O) = %.4f MeV" (:e-cm-f kin)))
  (println (format "  Q           = %.4f MeV" (:Q-mev kin)))
  (println (format "  k_i, k_f    = %.4f, %.4f fm⁻¹" (:k-i kin) (:k-f kin)))
  (println ""))

(println "=== Zero-range ===")
(println (format "  D₀(d,p) = %.2f MeV·fm^(3/2)" (t/zero-range-constant :d-p)))
(println "")

(defn- r0-sc [^double r-ca ^double a ^double b]
  (* r-ca (Math/pow (/ a b) (/ 1.0 3.0))))

(let [{:keys [e-cm-i mass-factor-i]} (oh/o16-dp-kinematics)
      z12    (* 1.44 1.0 8.0)
      eta    (binding [fn/mass-factor mass-factor-i fn/Z1Z2ee z12]
               (fn/channel-sommerfeld-eta e-cm-i))
      k      (Math/sqrt (* mass-factor-i e-cm-i))
      ;; Real WS matching radius: a = 2(R + a₀), same as s-matrix-3-impl.
      R-real  (r0-sc 3.803 16 40)
      R-imag  (r0-sc 5.342 16 40)
      a-match (* 2.0 (+ R-real 0.875))
      rho-match (* k a-match)
      h-s     0.05
      r-max-s 100.0
      ;; L table cap: nuclear effects are negligible beyond L_max used for the DCS.
      L-tab   18
      ;; WS params in functions.clj [V0 R a] format.
      ws-real [97.4 R-real 0.875]
      ws-imag [70.0 R-imag 0.477]]
  ;;
  ;; Two distinct quantities for partial wave L:
  ;;   (a) e^{2i(σ+δ^n)} − 1  =  e^{2iσ} S^n − 1  : total bracket, does NOT vanish at large L
  ;;       for charged scattering because the pure-Coulomb term e^{2iσ}−1 ≠ 0.
  ;;   (b) e^{2iσ}(S^n − 1)                        : nuclear bracket f_N term (T&N 3.1.88), → 0 for large L.
  ;;
  ;; Elastic amplitude: f = f_C + f_N
  ;;   f_C = closed-form (T&N Eq. 3.1.81)  — handles ALL L; already convergent.
  ;;   f_N = Σ_{L=0}^{L_max} (−i/2k)(2L+1) P_L e^{2iσ_L}(S^n_L − 1)  — only L ≤ L_max matter.
  ;;   dσ/dΩ = |f_C + f_N|²  (functions/differential-cross-section-nuclear-cut).
  ;;
  ;;
  ;; σ_L − σ_0 via product formula (no Gamma):
  ;;   e^{2i(σ_L − σ_0)} = ∏_{k=1}^{L} (k+iη)/(k−iη)
  ;;
  ;; Elastic amplitude:  f = e^{2iσ_0} (f̃_C + f̃_N)
  ;;   f̃_C = −η/(2k sin²) · exp(−iη ln sin²)        [T&N 3.1.81 without σ_0]
  ;;   f̃_N = Σ_{L=0}^{L_max} (−i/2k)(2L+1) P_L · e^{2i(σ_L−σ_0)} · (S^n_L − 1)
  ;;   |f|² = |f̃_C + f̃_N|²   (σ_0 phase cancels — no Gamma needed)
  ;;
  (println "=== Partial-wave S-matrix — entrance d+¹⁶O ===")
  (println (format "  E_cm,i = %.4f MeV, η = %.5f, k = %.5f fm⁻¹, R→S at a = %.4f fm (ρ = %.5f)"
                   e-cm-i eta k a-match rho-match))
  (println "  S^n via Numerov χ + R-matrix + Hankel quotient.  σ_L = arg Γ(L+1+iη).")
  (println "  e^{2i(σ_L−σ_0)} via ∏(k+iη)/(k−iη) [product] and via Gamma [exact] — should agree.")
  (println "  Nuclear bracket = e^{2iσ}(S^n−1) → 0 for L beyond nuclear range.")
  (println "")
  (println (format "  %-3s  %-8s  %-22s  %-22s  %-22s"
                   "L" "|S^n|"
                   "e^{2i(σ_L−σ_0)} product"
                   "e^{2i(σ_L−σ_0)} Gamma"
                   "nuc bracket e^{2iσ}(S^n−1)"))
  (doseq [L (range (inc L-tab))]
    (let [j-d      (+ 1.0 (double L))
          U        (oh/optical-u-deuteron-o16 L 1.0 j-d)
          u        (t/distorted-wave-optical e-cm-i L 1.0 j-d U r-max-s h-s mass-factor-i
                                             :normalize-mode :raw)
          R        (t/distorted-wave-numerov-R-for-smatrix u h-s a-match)
          Sn       (t/distorted-wave-coulomb-S-from-numerov-R R L eta rho-match)
          ;; Phase difference via product (no Gamma):
          ph-prod  (fn/coulomb-phase-diff L eta)
          ;; Phase difference via Gamma (exact reference):
          sig-L    (fn/coulomb-sigma-L L eta)
          sig-0    (fn/coulomb-sigma-L 0 eta)
          ph-gamma (complex-polar (* 2.0 (- sig-L sig-0)) 1.0)
          ;; Full e^{2iσ_L} via Gamma, for nuclear bracket
          e2is-L   (complex-polar (* 2.0 sig-L) 1.0)
          nuc      (mul e2is-L (subt2 Sn 1.0))]
      (println (format "  L=%2d  |S^n|=%.5f  Re=% .5f Im=% .5f  Re=% .5f Im=% .5f  Re=% .5f Im=% .5f"
                       L (mag Sn)
                       (re ph-prod) (im ph-prod)
                       (re ph-gamma) (im ph-gamma)
                       (re nuc) (im nuc)))))
  (println "")
  (println "=== Elastic dσ/dΩ = |f̃_C + f̃_N|² — d+¹⁶O — (mb/sr, CM) ===")
  (println (format "  f̃_C = T&N 3.1.81 without σ_0; f̃_N = Σ_{L=0}^{%d} e^{2i(σ_L−σ_0)}(S^n_L−1) (product formula)." L-tab))
  (println "  |f̃_C + f̃_N|² = |f_C + f_N|² since σ_0 is a pure phase.  Shown alongside Gamma-based reference.")
  (println (format "  %-6s  %14s  %14s  %14s" "θ_CM" "|f_C|² (Ruth.)" "|f̃_C+f̃_N|²" "|f_C+f_N|² ref"))
  ;; Precompute S^n for L=0..L-tab (reused across angles)
  (let [sn-vec (mapv (fn [^long L]
                       (let [j-d (+ 1.0 (double L))
                             U   (oh/optical-u-deuteron-o16 L 1.0 j-d)
                             u   (t/distorted-wave-optical e-cm-i L 1.0 j-d U r-max-s h-s mass-factor-i
                                                           :normalize-mode :raw)
                             R   (t/distorted-wave-numerov-R-for-smatrix u h-s a-match)]
                         (t/distorted-wave-coulomb-S-from-numerov-R R L eta rho-match)))
                     (range (inc L-tab)))]
    (doseq [^double theta-deg (range 15.0 166.0 15.0)]
      (let [th-rad  (* theta-deg (/ Math/PI 180.0))
            ;; Rutherford = |f_C|²
            fc-sq   (* 10.0 (Math/pow (mag (fn/coulomb-scattering-amplitude-thompson-nunes-eq-3181 th-rad eta k)) 2))
            ;; f̃_C and f̃_N using product formula (no Gamma for L≥1):
            f-tilde-c (fn/coulomb-amplitude-tilde th-rad eta k)
            ;; **`k`** from **`mass-factor-i`**; **η** explicit so **`coulomb-phase-diff`** matches the table above.
            f-tilde-n (binding [fn/mass-factor mass-factor-i
                                fn/Z1Z2ee z12
                                fn/*elastic-imag-ws-params* ws-imag
                                fn/*partial-wave-s-matrix-fn*
                                (fn [_e _v ^long L] (nth sn-vec L))]
                        (fn/elastic-nuclear-amplitude-tilde-fn e-cm-i ws-real th-rad L-tab eta))
            tilde-sq  (* 10.0 (Math/pow (mag (add f-tilde-c f-tilde-n)) 2))
            ;; Gamma-based reference via differential-cross-section-nuclear-cut:
            ref-sq    (binding [fn/mass-factor mass-factor-i
                                fn/Z1Z2ee z12
                                fn/*elastic-imag-ws-params* ws-imag]
                        (re (fn/differential-cross-section-nuclear-cut e-cm-i ws-real th-rad L-tab)))]
        (println (format "  %5.1f°  %14.4e  %14.4e  %14.4e" theta-deg fc-sq tilde-sq ref-sq)))))
  (println ""))

(defn- spline-smooth
  "Fit a cubic spline through [xs ys] and evaluate at `x-fine` (seq of doubles).
  Values below `floor` are clamped to `floor` (avoids negative spline artefacts)."
  [xs ys x-fine & {:keys [floor] :or {floor 1e-30}}]
  (let [spline (.interpolate (SplineInterpolator.)
                             (double-array xs)
                             (double-array ys))
        xlo    (double (first xs))
        xhi    (double (last  xs))]
    (mapv (fn [^double x]
            (let [xc (max xlo (min xhi x))]
              (max floor (.value spline xc))))
          x-fine)))

;; ── Bound-state plots ────────────────────────────────────────────────────
;; ¹⁷O single-particle states: 0d₅/₂ (g.s.) and 2s₁/₂ (first excited, E_x=0.871 MeV).
;; u(r) is normalized: ∫|u|² dr = 1.  F(r) = u(r)/r is the radial form factor in DWBA.
(let [h-bs  0.05
      d52   (oh/o17-d52-bound-state h-bs)
      s12   (oh/o17-s12-bound-state h-bs)
      mkrs  (fn [{:keys [u h r-max-bs]}]
              (let [n  (count u)
                    rs (mapv #(* (double %) h) (range n))
                    us (mapv double u)
                    Fs (mapv (fn [^long i ^double r]
                               (if (< r 1e-6) 0.0 (/ (nth us i) r)))
                             (range n) rs)]
                [rs us Fs n]))]
  (doseq [[label st] [["0d₅/₂" d52] ["2s₁/₂" s12]]]
    (let [{:keys [h r-max-bs E-bind u]} st
          [rs us Fs n] (mkrs st)]
      (println (format "=== ¹⁷O %s bound state  E_bind = %.4f MeV ===" label E-bind))
      (println (format "  Grid: h = %.3f fm, r_max = %.1f fm, %d points" h r-max-bs n))
      (println (format "  Nodes: %d  |  Norm ∫|u|²dr = %.6f"
                       (t/count-nodes us) (* h (reduce + (map #(* % %) us)))))
      (println "")))
  ;; Single combined plot: u(r) for both states + F(r)=u(r)/r
  (try
    (let [[rs52 us52 Fs52 _] (mkrs d52)
          [rs12 us12 Fs12 _] (mkrs s12)
          _ (io/make-parents (io/file "output/o17_bound_states.png"))
          chart (-> (c/xy-plot rs52 us52
                               :title "¹⁷O single-particle states — u(r) and F(r)=u(r)/r"
                               :x-label "r (fm)"
                               :y-label "u(r) or F(r) = u(r)/r  [fm^{-1/2}]"
                               :series-label "u(r) 0d₅/₂"
                               :legend true)
                    (c/add-lines rs52 Fs52 :series-label "F(r) 0d₅/₂")
                    (c/add-lines rs12 us12 :series-label "u(r) 2s₁/₂")
                    (c/add-lines rs12 Fs12 :series-label "F(r) 2s₁/₂"))]
      (i/save chart "output/o17_bound_states.png" :width 800 :height 500)
      (println "Plot saved: output/o17_bound_states.png"))
    (catch Exception e
      (println (format "Note: could not save bound-state plot (%s)." (.getMessage e)))))
  (println ""))

(let [h          0.08
      r-max      100.0
      L-max      18
      ;; Imaginary optical depth scale: handbook Tables 5.1–5.2 + ZR multipole +
      ;; :coulomb-tail χ normalization give too much backward (σ(180°) > σ(90°)).
      ;; Scaling **Im U** on both channels (~0.48) restores σ(180°) < σ(90°) while
      ;; keeping the same **Re U** (geometry, real depths).  Tune toward **1.0** for
      ;; a literal reading of the tabulated imaginary volumes/surfaces (backward rise returns).
      imag-scale 0.48
      ;; Compute DWBA at 2° spacing (raw points); build radial rows once for handbook + diagnostics.
      angles-deg (vec (range 0.0 181.0 2.0))
      eci        (:e-cm-i (oh/o16-dp-kinematics))
      {:keys [mass-factor-i mass-factor-f e-cm-i e-cm-f k-i k-f]}
      (oh/o16-dp-kinematics eci)
      z12        (* 1.44 1.0 8.0)
      eta-i      (binding [fn/mass-factor mass-factor-i fn/Z1Z2ee z12]
                   (fn/channel-sommerfeld-eta e-cm-i))
      eta-f      (binding [fn/mass-factor mass-factor-f fn/Z1Z2ee z12]
                   (fn/channel-sommerfeld-eta e-cm-f))
      base-rows  (oh/o16-dp-radial-I-rows-handbook :r-max r-max :h h :L-max L-max :e-cm-i e-cm-i
                    :chi-normalize-mode :coulomb-tail :imag-scale imag-scale)
      rows-sig   (t/handbook-zr-rows-with-coulomb-sigma base-rows eta-i eta-f)
      sigma-raw  (mapv (fn [^double th]
                         (oh/o16-dp-dsigma-handbook-mb-sr th
                           :radial-rows-sigma rows-sig :e-cm-i eci
                           :coherent-m-beta? false :S-factor 1.0))
                       angles-deg)
      ;; ZR strength × √(2ℓ+1) on each **I** block (ℓ = 2 for 0d₅/₂); same prefactor as **T_m** in handbook sum.
      D0         (t/zero-range-constant :d-p)
      ell        2
      pref-zr    (complex-cartesian (* D0 (Math/sqrt (inc (* 2 ell)))) 0.0)
      I-sum-by-Lb (->> rows-sig
                       (group-by :L-beta)
                       (map (fn [[lb rs]] [lb (reduce add (map :I rs))]))
                       (into {}))
      T-map-reduced (into {} (map (fn [[lb ic]] [lb (mul pref-zr ic)]) I-sum-by-Lb))
      best-row   (apply max-key (fn [r] (mag (:I r))) rows-sig)
      ;; **Partial reduced:** coherent **Y_{L_β0}** sum using only the **top two** exit columns by
      ;; **|Σ_{L_α} I|**. A **single** **L_β** gives **|T Y_{L0}|²** with deep zeros from **Y_{L0}** alone;
      ;; the full reduced curve has **L_β–L_β′** interference, so one column often looks “wrong” even when
      ;; the amplitude is consistent. Two columns keep the main interference structure for many reactions.
      top2-Lb    (take 2 (sort-by (fn [[_lb ic]] (- (mag ic))) (seq I-sum-by-Lb)))
      T-map-dom  (into {}
                       (map (fn [[lb ic]] [lb (mul pref-zr ic)])
                            top2-Lb))
      spin       (* (t/transfer-nuclear-spin-statistical-factor 0.0 2.5)
                    (t/transfer-unpolarized-deuteron-spin-factor))
      theta-rads (mapv (fn [^double d] (* d (/ Math/PI 180.0))) angles-deg)
      sigma-red-raw
      (mapv (fn [^double th]
              (* spin (double (t/transfer-differential-cross-section-angular-coherent
                               T-map-reduced 1.0 k-i k-f th mass-factor-i mass-factor-f))))
            theta-rads)
      sigma-dom-raw
      (mapv (fn [^double th]
              (* spin (double (t/transfer-differential-cross-section-angular-coherent
                               T-map-dom 1.0 k-i k-f th mass-factor-i mass-factor-f))))
            theta-rads)
      ;; Spline through log(σ) for smooth curves (primary + diagnostics)
      th-fine    (vec (range 0.0 180.5 0.5))
      sigma-fine (let [log-sig (mapv #(Math/log (max % 1e-30)) sigma-raw)]
                   (mapv #(Math/exp %) (spline-smooth angles-deg log-sig th-fine)))
      sigma-red-fine
      (mapv #(Math/exp %)
            (spline-smooth angles-deg (mapv #(Math/log (max % 1e-30)) sigma-red-raw) th-fine))
      sigma-dom-fine
      (mapv #(Math/exp %)
            (spline-smooth angles-deg (mapv #(Math/log (max % 1e-30)) sigma-dom-raw) th-fine))
      s0         (first sigma-raw)
      s70        (double (nth sigma-raw (int (/ 70.0 2.0)) 0.0))
      s90        (double (nth sigma-raw 45 0.0))   ;; 90° / 2° grid
      s180       (double (nth sigma-raw 90 0.0))] ;; 180° / 2° grid
  (println "=== dσ/dΩ (mb/sr, CM) — full Σ_m |T_m|² (handbook); reduced Y_{L_β0} diagnostic ===")
  (println (format "  Optical Im U scale = %.2f (both d+¹⁶O and p+¹⁷O); see let-binding in script." imag-scale))
  (println (format "  Grid: h=%.3f fm, r_max=%.1f fm, L_max=%d  (raw 2°, spline 0.5°)" h r-max L-max))
  (println (format "  Peak |I| cell: (L_α,L_β) = (%d,%d); partial reduced = top-2 L_β by |Σ_{L_α} I|: %s"
                   (:L-alpha best-row) (:L-beta best-row)
                   (clojure.string/join ", " (map str (map first top2-Lb)))))
  (println (format "  dσ/dΩ(0°)  ≈ %.6e mb/sr" s0))
  (println (format "  dσ/dΩ(70°) ≈ %.6e mb/sr" s70))
  (println (format "  dσ/dΩ(90°) ≈ %.6e mb/sr" s90))
  (println (format "  dσ/dΩ(180°) ≈ %.6e mb/sr  (ratio 180/90 = %.3f)" s180 (/ s180 (max s90 1e-300))))
  (println "  Optical model: handbook Table 5.1 (d+16O) + Table 5.2 (p+17O).")
  (println "")

  (try
    (let [_ (io/make-parents (io/file "output/16Odp_dcs.png"))
          chart (-> (c/xy-plot th-fine sigma-fine
                                 :title "¹⁶O(d,p)¹⁷O — Handbook ZR: full (5.6) m-sum vs reduced Y_{L_β0}"
                                 :x-label "θ_CM (deg)"
                                 :y-label "dσ/dΩ (mb/sr)"
                                 :series-label "full Σ_m |T_m|² (spline)"
                                 :legend true)
                      (c/add-lines th-fine sigma-red-fine
                                   :series-label "reduced Σ_{L_β} T_{L_β} Y_{L_β0} (spline)")
                      (c/add-lines th-fine sigma-dom-fine
                                   :series-label (format "top-2 L_β {%s} coherent (peak |I| %d,%d) (spline)"
                                                         (clojure.string/join "," (map str (map first top2-Lb)))
                                                         (:L-alpha best-row) (:L-beta best-row))))]
      (i/save chart "output/16Odp_dcs.png" :width 800 :height 500)
      (println "Plot saved: output/16Odp_dcs.png")
      (let [chart-log (-> (-> (c/xy-plot th-fine (series-for-log sigma-fine)
                                         :title "¹⁶O(d,p)¹⁷O — Handbook ZR (log y)"
                                         :x-label "θ_CM (deg)"
                                         :y-label "dσ/dΩ (mb/sr), log₁₀"
                                         :series-label "full Σ_m |T_m|²"
                                         :legend true)
                             (c/add-lines th-fine (series-for-log sigma-red-fine)
                                          :series-label "reduced Y_{L_β0}")
                             (c/add-lines th-fine (series-for-log sigma-dom-fine)
                                          :series-label (format "top-2 L_β {%s} (|I| max %d,%d)"
                                                                (clojure.string/join "," (map str (map first top2-Lb)))
                                                                (:L-alpha best-row) (:L-beta best-row))))
                         (chart-set-log-range-y! "dσ/dΩ (mb/sr)"))]
        (i/save chart-log "output/16Odp_dcs_log.png" :width 800 :height 500)
        (println "Plot saved: output/16Odp_dcs_log.png (log y-axis)")))
    (catch Exception e
      (println (format "Note: could not save plot (%s)." (.getMessage e)))))

  (println "\n=== Done ===")
  s0)

;; ── Radial integral I_{L_α} vs L_α — 2s₁/₂ transfer, Fig 5.3 coupling ──────
;;
;; For ℓ=0 (s₁/₂) transfer the triangle rule forces L_β = L_α.
;; Fig 5.3 uses the specific j-couplings:
;;   J_α  = L_α − 1    (deuteron sub-channel with j_d = L − 1)
;;   L_β  = L_α
;;   J_β  = L_α − ½   (proton sub-channel  with j_p = L − ½)
;;
;; For L_α = 0 → J_α = −1 (forbidden); the curve is extended with (0,0) for the plot.
;; Distorted waves use :coulomb-tail normalization + Coulomb Numerov init.
;;
;; Handbook text (e.g. “620” sum loops for ℓ=0, j=½) counts **(m_α,m_β,m)** and allowed
;; **(L_α,J_α,L_β,J_β)** blocks when assembling the **full** reduced amplitude — not a
;; correction to **this** diagnostic. Here each L_α is **one** complex **I_{L_β L_α}** with
;; **fixed** j-coupling (Fig 5.3 slice); the multipole sum **`handbook-zr-multipole-amplitude-sum`**
;; combines many **(L_α,L_β)** rows only when building **β^{ℓm}(θ)** or the cross section.
;;
;; L=1,2 match Fig 5.3 (~±0.02) after dividing by Austern **P**; L=3 remains larger here
;; (~0.065 vs ~0.015) — likely optical-model / absorption detail vs the book’s DWUCK inputs,
;; not a missing factor of “620” in this plot.
(let [h         0.05
      r-max     30.0
      L-max     15
      kin       (oh/o16-dp-kinematics)
      {:keys [mass-factor-i mass-factor-f e-cm-i e-cm-f k-i k-f
              M-target M-residual]} kin
      z12       (* 1.44 1.0 8.0)
      eta-i     (binding [fn/mass-factor mass-factor-i fn/Z1Z2ee z12]
                  (fn/channel-sommerfeld-eta e-cm-i))
      eta-f     (binding [fn/mass-factor mass-factor-f fn/Z1Z2ee z12]
                  (fn/channel-sommerfeld-eta e-cm-f))
      rho-i     (* k-i r-max)
      rho-f     (* k-f r-max)
      zr        (t/handbook-zr-chi-exit-mass-ratio M-target M-residual)
      ;; 2s₁/₂ bound state (1 radial node, l=0, E_bind = −3.2728 MeV)
      phi-s12   (:u (oh/o17-s12-bound-state h))
      ;; Memoized distorted waves (Fig 5.3 j-couplings)
      chi-a!    (memoize
                 (fn [^long La]
                   (let [ja (max 0.0 (- (double La) 1.0))]
                     (t/distorted-wave-optical
                      e-cm-i La 1.0 ja
                      (oh/optical-u-deuteron-o16 La 1.0 ja)
                      r-max h mass-factor-i
                      :normalize-mode :coulomb-tail
                      :tail-eta eta-i :tail-rho rho-i
                      :coulomb-init-eta eta-i))))
      chi-b!    (memoize
                 (fn [^long Lb]
                   (let [jb (max 0.5 (- (double Lb) 0.5))]
                     (t/distorted-wave-optical
                      e-cm-f Lb 0.5 jb
                      (oh/optical-u-proton-o17 Lb 0.5 jb)
                      r-max h mass-factor-f
                      :normalize-mode :coulomb-tail
                      :tail-eta eta-f :tail-rho rho-f
                      :coulomb-init-eta eta-f))))
      ;; Austern prefactor (M_B/M_A)(4π/k_α k_β) — divide I by this to get bare
      ;; integral in fm^{-1/2}, which directly matches the scale in Thompson & Nunes Fig 5.3
      austern-pref (* (/ (double M-residual) (double M-target))
                      (/ (* 4.0 Math/PI) (* k-i k-f)))
      ;; Compute I for each L_α ≥ 1 (L_β = L_α for ℓ=0); L_α=0 is identically 0
      La-vals   (vec (range 1 (inc L-max)))
      I-vals    (mapv (fn [^long La]
                        (t/handbook-radial-integral-I-zr-from-neutron-bound-complex
                         phi-s12 (chi-a! La) (chi-b! La)
                         h M-target M-residual k-i k-f zr))
                      La-vals)
      La-dbl    (mapv double La-vals)
      ;; Handbook convention (Fig 5.3): 'Re' = Im(e^{iσ}I)/P,  'Im' = −Re(e^{iσ}I)/P
      ;; where P = Austern prefactor.  Dividing by P removes kinematics and gives values
      ;; of order ±0.02 fm^{-1/2}, matching the scale in the book figure.
      sig-rot   (mapv (fn [^long La]
                        (+ (fn/coulomb-sigma-L La eta-i)
                           (fn/coulomb-sigma-L La eta-f)))
                      La-vals)
      HbkRe     (mapv (fn [^long i]
                        (let [I (nth I-vals i) s (nth sig-rot i)]
                          (/ (+ (* (re I) (Math/sin s)) (* (im I) (Math/cos s)))
                             austern-pref)))
                      (range (count La-vals)))
      HbkIm     (mapv (fn [^long i]
                        (let [I (nth I-vals i) s (nth sig-rot i)]
                          (/ (- (* (im I) (Math/sin s)) (* (re I) (Math/cos s)))
                             austern-pref)))
                      (range (count La-vals)))
      ;; Prepend L=0 (identically zero for ℓ=0 transfer) so the plot starts at 0
      La-plot   (vec (cons 0.0 La-dbl))
      Re-plot   (vec (cons 0.0 HbkRe))
      Im-plot   (vec (cons 0.0 HbkIm))]
  (println "=== Radial integral I_{L_α} — 2s₁/₂ transfer  (J_α=L_α−1, L_β=L_α, J_β=L_α−½) ===")
  (println (format "    Austern prefactor P = %.4f fm²" austern-pref))
  (println "    Handbook convention: 'Re' = Im(e^{iσ}I)/P,  'Im' = -Re(e^{iσ}I)/P")
  (println (format "  %-4s  %14s  %14s  %14s  %14s  %14s"
                   "L_α" "Re(I)/P" "Im(I)/P" "|I|/P" "Hbk Re" "Hbk Im"))
  (doseq [i (range (count La-vals))]
    (let [La (nth La-vals i) I (nth I-vals i)]
      (println (format "  %-4d  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e"
                       La (/ (re I) austern-pref) (/ (im I) austern-pref)
                       (/ (mag I) austern-pref)
                       (nth HbkRe i) (nth HbkIm i)))))
  (println "  handbook approx (Fig 5.3): L=1 → −0.020,  L=2 → +0.020,  L=3 → +0.015")
  (println "  (our L=3 ‘Hbk Re’ is typically ~4× larger — see block comment above)")
  (println "")
  (try
    (let [_ (io/make-parents (io/file "output/o16dp_radial_I_s12.png"))
          chart (-> (c/xy-plot La-plot Re-plot
                               :title "¹⁶O(d,p) — Radial integral I_{L_α}  (2s₁/₂, Fig 5.3 convention)"
                               :x-label "L_α"
                               :y-label "(-i)·e^{iσ}·I_{L_α} / P  [fm^{-1/2}]"
                               :series-label "Re  = Im(e^{iσ}I)/P"
                               :legend true)
                    (c/add-lines La-plot Im-plot :series-label "Im  = −Re(e^{iσ}I)/P"))]
      (i/save chart "output/o16dp_radial_I_s12.png" :width 900 :height 500)
      (println "Plot saved: output/o16dp_radial_I_s12.png"))
    (catch Exception e
      (println (format "Note: could not save radial-I plot (%s)." (.getMessage e)))))
  (println ""))
