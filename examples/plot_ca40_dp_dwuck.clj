;; Plot Ca40(d,p)Ca41 g.s. transfer — Nexus-NRS (`dwba.benchmark.ca40-dwuck`) vs embedded DWUCK4
;; reference σ(θ).  Listing column **Inelsig** is taken as **fm²/sr** and multiplied by
;; `dwba.benchmark.ca40-dwuck/dwuck-inelsig-fm2-sr->mb-sr` (10) for **mb/sr** curves. For (d,p)
;; stripping it is still **transfer** dσ/dΩ, not nuclear inelastic scattering (e.g. d,d′).
;;
;; From project root:
;;   lein repl
;;   (load-file "examples/plot_ca40_dp_dwuck.clj")
;;
;; Writes (gitignored):
;;   output/ca40_dp_dwuck_mb_sr.png      — both in mb/sr (linear y)
;;   output/ca40_dp_dwuck_mb_sr_log.png  — same data, **log₁₀ y** (easier when scales differ a lot)
;;   output/ca40_dp_dwuck_normalized.png — each curve divided by its own max(σ)
;;
;; **Scale vs shape:** Flux calibration matches σ at **30°** only (overall factor).
;; **Nexus curve** is pure DWBA (`ca40-dp-dsigma-mb-sr`, default **`:angular-mode :coherent`**): Coulomb in χ,
;; ZR POST, coherent multipole angular factor (**θ-dependent**, CM-symmetric for one **L**). **DWUCK4** listing
;; still has full **(1.9)–(1.11)** / asymmetry — not the same reduced model. Use **`:angular-mode :m-sum`** for
;; **l_i=0** isotropic orbital average (horizontal σ vs θ).

(ns examples.plot-ca40-dp-dwuck
  (:require [dwba.benchmark.ca40-dwuck :as ca40]
            [incanter.core :as i]
            [incanter.charts :as c]
            [clojure.java.io :as io])
  (:import [org.jfree.chart.axis LogAxis]))

(def ^:private dwuck-transfer-reference-inelsig-fm2-sr
  "Embedded myrun.lis angular table: **Inelsig** column as printed (**fm²/sr**), not mb/sr.
  Same digits as `dwba.ca40-dwuck-benchmark-test` sparse table / listing."
  [[0.0 5.4232e-2] [5.0 7.0429e-2] [10.0 1.1613e-1] [15.0 1.9254e-1] [20.0 3.0309e-1]
   [25.0 4.2502e-1] [30.0 5.0653e-1] [35.0 5.0825e-1] [40.0 4.3680e-1] [45.0 3.3440e-1]
   [50.0 2.4526e-1] [55.0 1.9219e-1] [60.0 1.7466e-1] [65.0 1.7977e-1] [70.0 1.9313e-1]
   [75.0 2.0362e-1] [80.0 2.0431e-1] [85.0 1.9295e-1] [90.0 1.7214e-1] [95.0 1.4758e-1]
   [100.0 1.2531e-1] [105.0 1.0947e-1] [110.0 1.0104e-1] [115.0 9.8157e-2] [120.0 9.7973e-2]
   [125.0 9.8494e-2] [130.0 9.8878e-2] [135.0 9.8858e-2] [140.0 9.8369e-2] [145.0 9.7310e-2]
   [150.0 9.5269e-2] [155.0 9.1769e-2] [160.0 8.6961e-2] [165.0 8.1744e-2] [170.0 7.7167e-2]
   [175.0 7.4035e-2] [180.0 7.2915e-2]])

(def ^:private dwuck-transfer-reference-mb-sr
  "Inelsig (**fm²/sr**) × `ca40/dwuck-inelsig-fm2-sr->mb-sr` → **mb/sr** for plots vs Nexus."
  (mapv (fn [[^double th ^double s]]
          [th (* (double ca40/dwuck-inelsig-fm2-sr->mb-sr) s)])
        dwuck-transfer-reference-inelsig-fm2-sr))

(defn- vec-max ^double [xs]
  (double (reduce max xs)))

(defn- normalize-by-max [xs]
  (let [m (vec-max xs)]
    (if (or (zero? m) (Double/isNaN m) (Double/isInfinite m))
      xs
      (mapv #(/ (double %) m) xs))))

(defn- clamp-pos-for-log ^double [^double y]
  "Floor for log-scale plotting only (σ must be > 0 for JFreeChart LogAxis)."
  (if (or (Double/isNaN y) (Double/isInfinite y) (<= y 0.0))
    1e-30
    (max y 1e-30)))

(defn- series-for-log [ys]
  (mapv clamp-pos-for-log ys))

(defn- chart-set-log-range-y!
  "Replace linear range axis with `LogAxis` (mutates chart). Base 10 ticks (JFreeChart default)."
  [chart ^String y-label]
  (let [plot (.getPlot chart)
        axis (LogAxis. y-label)]
    ;; Allow auto-range down to tiny σ (Nexus can be ~1e-30 after clamp at some angles)
    (.setSmallestValue axis 1e-35)
    (.setRangeAxis plot 0 axis)
    chart))

(defn- compute-nexus-curve-mb-sr
  "Nexus σ(θ) in **mb/sr** via `ca40-dp-dsigma-mb-sr` (DWBA, default **κ=0**)."
  [theta-deg-vec {:keys [h r-max] :or {h 0.08 r-max 20.0}}]
  (mapv (fn [^double th] (ca40/ca40-dp-dsigma-mb-sr th :h h :r-max r-max))
        theta-deg-vec))

(defn- mismatch-explanation!
  []
  (println "
=== Ca40(d,p) plot: Nexus vs DWUCK ===
  • **Nexus** (default **`:angular-mode :coherent`**) has **θ-dependent** dσ from reduced multipole interference
    (**|Σ T_L Y_{L0}|²**-style), CM-symmetric for one dominant **L** — not identical to the listing curve.
  • **`:angular-mode :m-sum`** + **l_i=0** → **θ-flat** dσ (orbital addition theorem).
  • **DWUCK4** uses full **T** / **S** / **β** / **P_L^m** (manual **1.9–1.11**) — can differ in shape and asymmetry.
  • **Similar y-scale** = one multiplicative factor at 30° (χ flux calibration).
  • **Units:** mb/sr; flux-matched Nexus on linear + log; **normalized** PNG = raw σ / each max.
================================================================
"))

(defn make-plots!
  "Generate PNGs under `output/`. Returns `{:dwuck-thetas ... :nexus-sigma ...}`.

  **Options:** **`:h`**, **`:r-max`** — passed through to `ca40-dp-dsigma-mb-sr`.

  **Note:** Nexus and DWUCK curves are **not** meant to agree; see `mismatch-explanation!` output."
  [& {:keys [h r-max] :or {h 0.08 r-max 20.0}}]
  (mismatch-explanation!)
  (let [curve-opts {:h h :r-max r-max}
        dwuck-thetas (mapv first dwuck-transfer-reference-mb-sr)
        dwuck-sigma-mb (mapv second dwuck-transfer-reference-mb-sr)
        ;; Finer grid for our curve (smooth line)
        fine-thetas (vec (range 0.0 181.0 3.0))
        nexus-fine-mb (compute-nexus-curve-mb-sr fine-thetas curve-opts)
        ;; Same angles as DWUCK for a pointwise table print
        nexus-at-dwuck-mb (compute-nexus-curve-mb-sr dwuck-thetas curve-opts)
        ;; ~10¹–10²: distorted-wave-optical max-normalization vs DWUCK incoming flux
        flux-scale (double (ca40/ca40-dp-flux-scale-to-embedded-dwuck :h h :r-max r-max))
        nexus-fine-matched (mapv #(* flux-scale ^double %) nexus-fine-mb)
        nexus-at-dwuck-matched (mapv #(* flux-scale ^double %) nexus-at-dwuck-mb)
        out-dir (io/file "output")
        _ (.mkdirs out-dir)]

    (println (format "Flux calibration (match @ 30° CM to embedded DWUCK): ×%.4g%n" flux-scale))
    (println "Ca40(d,p)Ca41 g.s. — mb/sr | DWUCK listing | Nexus matched | Nexus raw")
    (println "θ_cm (deg) | DWUCK4    | Nx matched | Nx raw")
    (doseq [[th sd sm sr] (map vector dwuck-thetas dwuck-sigma-mb nexus-at-dwuck-matched nexus-at-dwuck-mb)]
      (println (format " %6.1f    | %.4e | %.4e | %.4e" th sd sm sr)))
    (when-let [row90 (first (filter (fn [[^double th _ _ _]] (< (Math/abs (- th 90.0)) 0.1))
                                    (map vector dwuck-thetas dwuck-sigma-mb nexus-at-dwuck-matched nexus-at-dwuck-mb)))]
      (let [[_ dw sm] row90]
        (println (format "\n@ 90° CM: DWUCK ≈ %.4e  |  Nexus matched ≈ %.4e  (Nexus default: coherent angular)%n"
                         dw sm))))

    (try
      ;; Raw comparison: DWUCK (tabulated) + Nexus-NRS (smooth)
      (let [chart-raw (c/xy-plot dwuck-thetas dwuck-sigma-mb
                                 :title (format "Ca40(d,p) mb/sr — scale @30° (×%.3g); Nexus coherent angular vs DWUCK"
                                                flux-scale)
                                 :x-label "θ_cm (deg)"
                                 :y-label "dσ/dΩ (mb/sr)"
                                 :series-label "DWUCK4 (listing mb/sr)"
                                 :legend true)
            chart-raw (c/add-lines chart-raw fine-thetas nexus-fine-matched
                                   :series-label "Nexus-NRS (flux-matched)")]
        (i/save chart-raw "output/ca40_dp_dwuck_mb_sr.png" :width 900 :height 520)
        (println "\nSaved: output/ca40_dp_dwuck_mb_sr.png"))

      ;; Log y: clamp non-positive σ to a tiny floor for display only
      (let [dw-log (series-for-log dwuck-sigma-mb)
            nx-log (series-for-log nexus-fine-matched)
            chart-log (c/xy-plot dwuck-thetas dw-log
                                  :title (format "Ca40(d,p) log y — flux ×%.3g @30°; angular models differ"
                                                 flux-scale)
                                  :x-label "θ_cm (deg)"
                                  :y-label "dσ/dΩ (mb/sr), logarithmic"
                                  :series-label "DWUCK4 (mb/sr)"
                                  :legend true)
            chart-log (c/add-lines chart-log fine-thetas nx-log
                                   :series-label "Nexus-NRS (flux-matched)")
            chart-log (chart-set-log-range-y! chart-log "dσ/dΩ (mb/sr)")]
        (i/save chart-log "output/ca40_dp_dwuck_mb_sr_log.png" :width 900 :height 520)
        (println "Saved: output/ca40_dp_dwuck_mb_sr_log.png (log y-axis)"))

      (let [nd (normalize-by-max dwuck-sigma-mb)
            nn (normalize-by-max nexus-fine-mb)
            chart-n (c/xy-plot dwuck-thetas nd
                               :title "Ca40(d,p) σ/σ_max — intrinsic shapes (no flux factor); still incompatible models"
                               :x-label "θ_cm (deg)"
                               :y-label "σ(θ) / max σ"
                               :series-label "DWUCK4"
                               :legend true)
            chart-n (c/add-lines chart-n fine-thetas nn :series-label "Nexus-NRS")]
        (i/save chart-n "output/ca40_dp_dwuck_normalized.png" :width 900 :height 520)
        (println "Saved: output/ca40_dp_dwuck_normalized.png"))

      (catch Exception e
        (println "Plot failed:" (.getMessage e))
        (println "Ensure Incanter/JFreeChart can run headless (Java 2D).")))

    {:dwuck-thetas dwuck-thetas
     :dwuck-sigma-mb-sr dwuck-sigma-mb
     :flux-scale-dwuck-embedded-30deg flux-scale
     :nexus-at-dwuck-angles-mb-sr-raw nexus-at-dwuck-mb
     :nexus-at-dwuck-angles-mb-sr-matched nexus-at-dwuck-matched
     :fine-thetas fine-thetas
     :nexus-fine-mb-sr-raw nexus-fine-mb
     :nexus-fine-mb-sr-matched nexus-fine-matched}))

;; Run when loaded via load-file
(println "\n--- examples.plot-ca40-dp-dwuck ---")
(println "Invoking (make-plots!) — next block explains why Nexus ≠ DWUCK.\n")
(make-plots!)
