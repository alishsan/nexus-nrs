;; Example: 16O(d,p)17O — **Handbook** ZR pipeline
;;
;; Uses **`dwba.benchmark.o16-dp-handbook`**: handbook §5.4 **F=u/r**, **(5.5)** radial **I**, **`handbook-zr-multipole-amplitude-sum`**
;; for the neutron bound in **¹⁷O**, **T_m ≈ D₀√(2ℓ+1)β_m**,
;; and **`transfer-differential-cross-section`** — parallel to **`ca40-pd-handbook`** (pickup) with
;; **α = d+¹⁶O**, **β = p+¹⁷O**.
;;
;; Reaction: 16O(d,p)17O (g.s. model: d₅/₂ → **ℓ=2** bound well)
;;
;; Load: `(load-file "examples/example_16Odp.clj")` from project root.
;;
;; Writes **output/16Odp_dcs.png** (linear y) and **output/16Odp_dcs_log.png** (log₁₀ y-axis).

(ns examples.example-16Odp
  (:require [dwba.benchmark.o16-dp-handbook :as oh]
            [dwba.transfer :as t]
            [functions :as fn]
            [complex :refer [re im mag mul subt2 add2 complex-cartesian complex-polar]]
            [fastmath.polynomials :as poly]
            [fastmath.core :as m]
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
  ;;   f_N = Σ_{L=0}^{L_max} (−i/k)(2L+1) P_L e^{2iσ_L}(S^n_L − 1)  — only L ≤ L_max matter.
  ;;   dσ/dΩ = |f_C + f_N|²  (functions/differential-cross-section-nuclear-cut).
  ;;
  ;;
  ;; σ_L − σ_0 via product formula (no Gamma):
  ;;   e^{2i(σ_L − σ_0)} = ∏_{k=1}^{L} (k+iη)/(k−iη)
  ;;
  ;; Elastic amplitude:  f = e^{2iσ_0} (f̃_C + f̃_N)
  ;;   f̃_C = −η/(2k sin²) · exp(−iη ln sin²)        [T&N 3.1.81 without σ_0]
  ;;   f̃_N = Σ_{L=0}^{L_max} (−i/k)(2L+1) P_L · e^{2i(σ_L−σ_0)} · (S^n_L − 1)
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
            f-tilde-n (reduce (fn [acc ^long L]
                                (let [Sn      (nth sn-vec L)
                                      ph-prod (fn/coulomb-phase-diff L eta)
                                      pl      (double (poly/eval-legendre-P L (m/cos th-rad)))
                                      bracket (mul ph-prod (subt2 Sn 1.0))
                                      ;; (-i/k)(2L+1) P_L × bracket; -i = complex-polar(-π/2, 1)
                                      contrib (mul (complex-polar (* -0.5 Math/PI) (/ (inc (* 2 L)) k))
                                                   pl bracket)]
                                  (add2 acc contrib)))
                              (complex-cartesian 0.0 0.0)
                              (range (inc L-tab)))
            tilde-sq  (* 10.0 (Math/pow (mag (add2 f-tilde-c f-tilde-n)) 2))
            ;; Gamma-based reference via differential-cross-section-nuclear-cut:
            ref-sq    (binding [fn/mass-factor mass-factor-i
                                fn/Z1Z2ee z12
                                fn/*elastic-imag-ws-params* ws-imag]
                        (re (fn/differential-cross-section-nuclear-cut e-cm-i ws-real th-rad L-tab)))]
        (println (format "  %5.1f°  %14.4e  %14.4e  %14.4e" theta-deg fc-sq tilde-sq ref-sq)))))
  (println ""))

(let [h 0.08
      r-max 100.0
      L-max 18
      angles-deg (range 0.0 181.0 5.0)
      curve (oh/o16-dp-angular-curve-handbook-mb-sr
             angles-deg :h h :r-max r-max :L-max L-max)
      sigma-vec (mapv :differential_cross_section_mb_sr curve)
      s0 (:differential_cross_section_mb_sr
          (first (oh/o16-dp-angular-curve-handbook-mb-sr [0.0] :h h :r-max r-max :L-max L-max)))
      s70 (:differential_cross_section_mb_sr
           (first (oh/o16-dp-angular-curve-handbook-mb-sr [70.0] :h h :r-max r-max :L-max L-max)))]
  (println "=== dσ/dΩ (mb/sr, CM) — incoherent Σ_m |T_m|²; χ :coulomb-tail norm (handbook default) ===")
  (println (format "  Grid: h=%.3f fm, r_max=%.1f fm, L_max=%d" h r-max L-max))
  (println (format "  dσ/dΩ(0°)  ≈ %.6e mb/sr" s0))
  (println (format "  dσ/dΩ(70°) ≈ %.6e mb/sr" s70))
  (println "  Absolute mb/sr vs a handbook figure requires matching optics, CM energy, and conventions;")
  (println "  optional :chi-normalize-mode :coulomb-tail | :raw | :max on o16-dp-handbook (see transfer/distorted-wave-optical).")
  (println "")
  (println "Optical model: Ca40-listing-like depths, radii scaled to A≈16–17, Z=8 (`o16-dp-handbook` docstring).")
  (println "Bound neutron in 17O: l=2, E≈−4.14 MeV, illustrative Woods–Saxon [56, 2.85, 0.6] MeV/fm.")
  (println "")

  (try
    (let [_ (io/make-parents (io/file "output/16Odp_dcs.png"))
          th (vec angles-deg)
          chart (c/xy-plot th sigma-vec
                           :title "16O(d,p)17O — Handbook ZR multipole"
                           :x-label "θ_CM (deg)"
                           :y-label "dσ/dΩ (mb/sr)"
                           :series-label "o16-dp-handbook"
                           :legend true)]
      (i/save chart "output/16Odp_dcs.png" :width 800 :height 500)
      (println "Plot saved: output/16Odp_dcs.png")
      (let [chart-log (-> (c/xy-plot th (series-for-log sigma-vec)
                                     :title "16O(d,p)17O — Handbook ZR (log y)"
                                     :x-label "θ_CM (deg)"
                                     :y-label "dσ/dΩ (mb/sr), logarithmic"
                                     :series-label "o16-dp-handbook"
                                     :legend true)
                         (chart-set-log-range-y! "dσ/dΩ (mb/sr)"))]
        (i/save chart-log "output/16Odp_dcs_log.png" :width 800 :height 500)
        (println "Plot saved: output/16Odp_dcs_log.png (log y-axis)")))
    (catch Exception e
      (println (format "Note: could not save plot (%s)." (.getMessage e)))))

  (println "\n=== Done ===")
  s0)
