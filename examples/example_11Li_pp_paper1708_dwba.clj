;; ¹¹Li**(p,p′)** — DWBA with **arXiv:1708.07719** Table **1** optics (**Set V** / **Set S**),
;; **Austern (5.5)+(5.6)** angular coupling (**ℓ = 1** dipole), **Coulomb-tail** **χ**,
;; and macroscopic transition form factor (**`:beta-scale`** tunes absolute **mb/sr** vs Fig. 5).
;;
;; **Set V** (volume **W** at small **r_w A^{1/3}**) + this **schematic F** typically needs **much** smaller **β_sc**
;; than **Set S** (surface **W**) for the same **mb/sr** — see ns doc **`dwba.benchmark.li11-pp-paper1708`**.
;;
;; Paper used **CHUCK3** + **H–D** / **Orlandini** E1(isoscalar) **FF** — see ns **`dwba.benchmark.li11-pp-paper1708`**.
;;
;; Run: `lein run -m clojure.main examples/example_11Li_pp_paper1708_dwba.clj`

(require '[dwba.benchmark.li11-pp-paper1708 :as li11]
         '[incanter.core :as i]
         '[incanter.charts :as c]
         '[clojure.java.io :as io])

(import '[org.jfree.chart.axis LogAxis])

(defn- clamp-pos-for-log ^double [^double y]
  (if (or (Double/isNaN y) (Double/isInfinite y) (<= y 0.0))
    1e-30
    (max y 1e-30)))

(defn- chart-set-log-range-y!
  [chart ^String y-label]
  (let [plot (.getPlot chart)
        axis (LogAxis. y-label)]
    (.setSmallestValue axis 1e-6)
    (.setRangeAxis plot 0 axis)
    chart))

(def E-p-lab 6.0)
(def E-ex 0.8)
;; Per-set **β_sc** targeting **~1 mb/sr** at **90° CM** (rough): **dσ ∝ β_sc²** for this pipeline.
;; Reference: single **β_sc = 0.08** gave **Set V ~ 3.6×10³** and **Set S ~ 7.3 mb/sr** here.
(def beta-sc-by-set
  {:V (* 0.08 (/ 1.0 (Math/sqrt 3635.0)))
   :S (* 0.08 (/ 1.0 (Math/sqrt 7.35)))})
(def theta-step 5.0)
(def L-max 20)
(def r-max 40.0)
(def h-step 0.02)

(println "=== ¹¹Li(p,p′) — paper optics + Austern (5.6), ℓ=1 ===")
(println (format "  E_p,lab = %.2f MeV   E_x = %.2f MeV   (per-set β_sc for ~1 mb/sr @ 90° with proxy F)"
                 E-p-lab E-ex))
(println "")

(doseq [oset [:V :S]]
  (println (str "  Optical Set " (name oset) ":"))
  (let [beta-scale (double (beta-sc-by-set oset))
        curve (li11/li11-pp-angular-curve-paper-mb-sr
                (vec (range 5.0 176.0 theta-step))
                :optical-set oset
                :e-p-lab-mev E-p-lab
                :e-ex-mev E-ex
                :beta-scale beta-scale
                :L-max L-max
                :r-max r-max
                :h h-step)
        mid (some #(when (< (Math/abs (- (:theta-deg %) 90.0)) 3.0) %) curve)
        s90 (or (:differential_cross_section_mb_sr mid)
                (li11/li11-pp-dsigma-paper-mb-sr 90.0 :optical-set oset :e-p-lab-mev E-p-lab :e-ex-mev E-ex
                        :beta-scale beta-scale :L-max L-max :r-max r-max :h h-step))
        out (format "output/11Li_pp_paper1708_Set_%s.png" (name oset))
        th (mapv :theta-deg curve)
        ys (mapv #(double (:differential_cross_section_mb_sr %)) curve)
        ylog (mapv clamp-pos-for-log ys)]
    (println (format "    β_sc = %.4e   dσ/dΩ ≈ %.3f mb/sr near 90° CM" beta-scale s90))
    (io/make-parents (io/file out))
    (-> (c/xy-plot th ys
                   :title (format "¹¹Li(p,p′) DWBA — paper Set %s (β_sc=%.2e)" (name oset) beta-scale)
                   :x-label "θ_CM (deg)"
                   :y-label "dσ/dΩ (mb/sr)"
                   :legend true)
        (i/save out :width 880 :height 520))
    (-> (c/xy-plot th ylog
                   :title (format "¹¹Li(p,p′) Set %s — log₁₀ y" (name oset))
                   :x-label "θ_CM (deg)"
                   :y-label "dσ/dΩ (mb/sr)"
                   :legend true)
        (chart-set-log-range-y! "dσ/dΩ (mb/sr)")
        (i/save (str out ".log.png") :width 880 :height 520))
    (println (format "    Wrote %s (+ .log.png)" out)))
  (println ""))

(println "Done.")
