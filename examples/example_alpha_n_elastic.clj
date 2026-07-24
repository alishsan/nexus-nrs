;; Neutral **α + n** elastic (no Coulomb): Woods–Saxon **V(r)** only, **Z₁Z₂ = 0**.
;;
;; **dσ/dΩ** uses **`binding`** on **α–n** **`mass-factor`**, **`Z1Z2ee`**, and
;; **`functions/*partial-wave-s-matrix-fn*`** so each **S_L** is **`alpha-n-S-from-numerov-for-L`**
;; (**`distorted-wave-optical`**; no intermediate radial rescale — preserves **u(r)** scaling across **L**). Compare **dσ(90°)** to default **`s-matrix`**.
;;
;; Run from repository root:
;;   lein run -m clojure.main examples/example_alpha_n_elastic.clj
;;
;; Writes under **output/** (gitignored):
;;   **alpha_n_elastic_dsigma.png** — dσ/dΩ vs θ_CM (linear y)
;;   **alpha_n_elastic_dsigma_log.png** — same, log₁₀ y
;;
;; **Partial waves:** **η = 0** ⇒ **σ_L = 0**; sum **L = 0 … L_max** (**L_max** below). **Δθ** for smooth PNGs.

(import '[java.awt BasicStroke]
        '[org.jfree.chart.axis LogAxis])

(require '[dwba.benchmark.alpha-n-elastic-neutral :as an]
         '[functions :as fn :refer [differential-cross-section mass-factor Z1Z2ee
                                    *partial-wave-s-matrix-fn*]]
         '[complex :as cpx]
         '[incanter.core :as i]
         '[incanter.charts :as c]
         '[clojure.java.io :as io])

(defn- clamp-pos-for-log ^double [^double y]
  (if (or (Double/isNaN y) (Double/isInfinite y) (<= y 0.0))
    1e-30
    (max y 1e-30)))

(defn- chart-set-log-range-y!
  [chart ^String y-label]
  (let [plot (.getPlot chart)
        axis (LogAxis. y-label)]
    (.setSmallestValue axis 1e-12)
    (.setRangeAxis plot 0 axis)
    chart))

(defn- chart-smooth-series!
  "Anti-alias and a slightly thick rounded stroke so dense **(θ,σ)** points read as a smooth line."
  [chart]
  (try
    (let [plot (.getPlot chart)
          ^org.jfree.chart.renderer.xy.XYItemRenderer r (.getRenderer plot 0)]
      (when r
        (.setSeriesStroke r 0 (BasicStroke. 1.75
                                           BasicStroke/CAP_ROUND
                                           BasicStroke/JOIN_ROUND))))
    (catch Throwable _))
  (try (.setAntiAlias chart true) (catch Throwable _))
  chart)

(defn- S-dist [a b]
  (Math/sqrt (+ (Math/pow (- (double (cpx/re a)) (double (cpx/re b))) 2.0)
                (Math/pow (- (double (cpx/im a)) (double (cpx/im b))) 2.0))))

(def E-cm 10.0)
;; Same illustrative WS as **`alpha-n-elastic-neutral-test`** **[V₀ R a]** MeV, fm.
(def ws-params [22.0 3.4 0.62])
(def L-max-partial 16)
(def theta-step-deg 0.5)

(def r-max-numerov 50.0)
(def h-numerov 0.01)

(println "=== α + n elastic — neutral (no Coulomb), central Woods–Saxon ===")
(println "")
(println (format "  E_cm        = %.2f MeV" E-cm))
(println (format "  μ (reduced) = %.4f MeV/c²" (an/alpha-n-reduced-mass)))
(println (format "  mass-factor = %.6f MeV⁻¹·fm⁻² (2μ/(ℏc)²)" (an/alpha-n-mass-factor)))
(println (format "  Z₁Z₂ e²     = %.1f (off)" (an/alpha-n-z1z2ee)))
(println (format "  WS [V₀ R a] = [%.1f %.2f %.2f] MeV, fm, fm"
                 (first ws-params) (second ws-params) (nth ws-params 2)))
(println (format "  Matching **a** = 2(R+a) = %.3f fm" (an/alpha-n-matching-radius-fm ws-params)))
(println "")

(println "=== |S_L (Numerov R-match) − S_L (s-matrix)| — L = 0 … 10 ===")
(doseq [L (range 0 11)]
  (let [S-ref (an/alpha-n-sigma-L-ref E-cm ws-params L)
        S-num (an/alpha-n-S-from-numerov-for-L E-cm ws-params L r-max-numerov h-numerov)
        d (S-dist S-ref S-num)]
    (println (format "  L = %2d  ΔS = %.4e" L d))))
(println "")

(println "=== Partial-wave sum **L = 0 … L_max** ===")
(println (format "  L_max = %d" L-max-partial))
(println "")

(def theta-degrees (vec (range 0.0 (+ 180.01 (/ theta-step-deg 2)) theta-step-deg)))

;; One Numerov integration per L (**`distorted-wave-optical`**).
(def S-numerov-by-L
  (vec (for [L (range 0 (inc L-max-partial))]
         (an/alpha-n-S-from-numerov-for-L E-cm ws-params L r-max-numerov h-numerov))))

(let [s-fn-numerov (fn [E V L]
                     (if (and (= (double E) (double E-cm))
                              (= V ws-params)
                              (<= 0 (long L) (long L-max-partial)))
                       (nth S-numerov-by-L (long L))
                       (an/alpha-n-sigma-L-ref E V L)))
      th90 (* 90.0 (/ Math/PI 180.0))
      ;; Default **s-matrix** at 90° (no override).
      dcs90-s-mat-r (double (cpx/re (binding [mass-factor (an/alpha-n-mass-factor)
                                              Z1Z2ee (an/alpha-n-z1z2ee)
                                              *partial-wave-s-matrix-fn* nil]
                                     (differential-cross-section E-cm ws-params th90 L-max-partial))))
      sigma-mb-sr-vec
      (mapv (fn [^double deg]
              (let [th (* deg (/ Math/PI 180.0))
                    z (binding [mass-factor (an/alpha-n-mass-factor)
                                Z1Z2ee (an/alpha-n-z1z2ee)
                                *partial-wave-s-matrix-fn* s-fn-numerov]
                        (differential-cross-section E-cm ws-params th L-max-partial))]
                (double (cpx/re z))))
            theta-degrees)
      i90 (min (dec (count sigma-mb-sr-vec))
               (max 0 (int (Math/round (/ 90.0 theta-step-deg)))))
      s90-num (double (get sigma-mb-sr-vec i90))]
  (println "=== dσ/dΩ: **S_L** from Numerov (**distorted-wave-optical**) ===")
  (println (format "  dσ/dΩ(90° CM) ≈ %.6e mb/sr (L_max = %d)" s90-num L-max-partial))
  (println (format "  same at 90° with default **s-matrix** → %.6e mb/sr (ratio Numerov/def %.4g)"
                   dcs90-s-mat-r (if (zero? dcs90-s-mat-r) Double/NaN (/ s90-num dcs90-s-mat-r))))
  (println "")
  (try
  (io/make-parents (io/file "output/alpha_n_elastic_dsigma.png"))
  (-> (c/xy-plot theta-degrees sigma-mb-sr-vec
        :title "α + n elastic (neutral) — dσ/dΩ vs θ_CM"
        :x-label "θ_CM (deg)"
        :y-label "dσ/dΩ (mb/sr)"
        :series-label "partial-wave sum"
        :legend true)
      (chart-smooth-series!)
      (i/save "output/alpha_n_elastic_dsigma.png" :width 900 :height 520))
  (println "Saved: output/alpha_n_elastic_dsigma.png")
  (let [y-log (mapv clamp-pos-for-log sigma-mb-sr-vec)]
    (-> (c/xy-plot theta-degrees y-log
          :title "α + n elastic — dσ/dΩ vs θ_CM (log₁₀ y)"
          :x-label "θ_CM (deg)"
          :y-label "dσ/dΩ (mb/sr)"
          :series-label "partial-wave sum"
          :legend true)
        (chart-smooth-series!)
        (chart-set-log-range-y! "dσ/dΩ (mb/sr)")
        (i/save "output/alpha_n_elastic_dsigma_log.png" :width 900 :height 520)))
  (println "Saved: output/alpha_n_elastic_dsigma_log.png")
  (catch Exception e
    (println "Could not write plots:" (.getMessage e))))
  (println "")
  (println "Done."))
