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

(let [h 0.08
      r-max 20.0
      L-max 6
      angles-deg (range 0.0 181.0 5.0)
      curve (oh/o16-dp-angular-curve-handbook-mb-sr
             angles-deg :h h :r-max r-max :L-max L-max)
      sigma-vec (mapv :differential_cross_section_mb_sr curve)
      s0 (:differential_cross_section_mb_sr
          (first (oh/o16-dp-angular-curve-handbook-mb-sr [0.0] :h h :r-max r-max :L-max L-max)))
      s70 (:differential_cross_section_mb_sr
           (first (oh/o16-dp-angular-curve-handbook-mb-sr [70.0] :h h :r-max r-max :L-max L-max)))]
  (println "=== dσ/dΩ (mb/sr, CM) — incoherent Σ_m |T_m|²; χ :coulomb-tail (default) ===")
  (println (format "  Grid: h=%.3f fm, r_max=%.1f fm, L_max=%d" h r-max L-max))
  (println (format "  dσ/dΩ(0°)  ≈ %.6e mb/sr" s0))
  (println (format "  dσ/dΩ(70°) ≈ %.6e mb/sr" s70))
  (println "  Absolute mb/sr vs a handbook figure requires matching optics, CM energy, and conventions;")
  (println "  optional :chi-normalize-mode :max | :coulomb-tail on o16-dp-handbook fns alters χ scaling.")
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
