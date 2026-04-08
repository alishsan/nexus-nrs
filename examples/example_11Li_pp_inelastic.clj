;; For **paper Table 1 optics**, **Austern (5.6)** multipole **θ** dependence, and **E_x = 0.8 MeV**, see
;; **`examples/example_11Li_pp_paper1708_dwba.clj`** / **`dwba.benchmark.li11-pp-paper1708`**.
;;
;; **¹¹Li(p,p′)** — schematic DWBA inelastic using **proton** distorted waves with **Coulomb + WS**
;; (and optional imaginary WS) via **`distorted-wave-entrance`** / **`distorted-wave-exit`**.
;;
;; **¹¹Li:** Z = 3, A = 11 (halo nucleus). Lab frame: proton beam on fixed ¹¹Li target.
;; Kinematics use non-relativistic reduced mass and CM energy for the Schrödinger equation.
;;
;; **Schematic physics:** quadrupole coupling (λ = 2), L_i = 0 → L_f = 2, excitation energy set to
;; **E_x ≈ 2.09 MeV** (same scale as the soft mode discussed for ¹¹Li in **`example_11Li_dd_inelastic.clj`**);
;; **β₂** and **V_params** are **not** fitted—adjust for serious comparisons.
;;
;; **dσ/dΩ scale:** **`distorted-wave-entrance`** / **`distorted-wave-exit`** default to **Coulomb-tail |u|** matching
;; on charged **WS + Coulomb** branches (**`:distorted-norm :coulomb-tail`**). Absolute **mb/sr** still need a full flux convention for fits.
;;
;; Writes (typically gitignored):
;;   **output/11Li_pp_inelastic_dsigma.png** — linear **y** (tiny **mb/sr** may look like a flat “zero” line).
;;   **output/11Li_pp_inelastic_dsigma_log.png** — log₁₀ **y** so the scale is visible.
;;
;; Run from repository root:
;;   lein run -m clojure.main examples/example_11Li_pp_inelastic.clj

(require '[dwba.inelastic :as inel]
         '[complex :as cpx]
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
    (.setSmallestValue axis 1e-40)
    (.setRangeAxis plot 0 axis)
    chart))

;; ---------------------------------------------------------------------------
;; Masses and kinematics (MeV/c² convention as in `example_11Li_dd_inelastic.clj`)
;; ---------------------------------------------------------------------------

(def ^:private m-p 938.272)
(def ^:private m-11Li (* 11.0 931.5))
(def ^:private A-11 11)
(def ^:private Z-11 3)

;; Proton lab kinetic energy (MeV)
(def E-lab 6.0)

;; CM kinetic energy: target at rest
(def E-CM (* E-lab (/ m-11Li (+ m-p m-11Li))))

(def mu-reduced (/ (* m-p m-11Li) (+ m-p m-11Li)))
(def mass-factor (/ (* 2.0 mu-reduced) (* 197.327 197.327)))

;; Excited state (schematic; multipolarity λ = 2)
(def E-ex 2.09)
(def lambda 2)
(def mu 0)
(def beta-2 0.12)

;; Woods–Saxon + optional absorption (same geometry for W as for V)
(def V-params [48.0 2.5 0.6])
(def W-params [12.0 2.5 0.6])

(def h 0.01)
(def r-max 30.0)

(def L-i 0)
(def L-f 2)

;; Crude exit lab energy for keyword (only matters for global CH89 paths; unused for WS+Coulomb U(r))
(def E-lab-exit (max 0.1 (- E-lab E-ex)))

;; ---------------------------------------------------------------------------
;; Distorted waves: **p + ¹¹Li** with **optical-potential-woods-saxon** (Coulomb + nuclear)
;; ---------------------------------------------------------------------------

(def chi-i
  (inel/distorted-wave-entrance E-CM L-i V-params h r-max
                                :projectile-type :p
                                :target-A A-11
                                :target-Z Z-11
                                :E-lab E-lab
                                :W-params W-params
                                :s 0.5
                                :j (+ L-i 0.5)
                                :mass-factor mass-factor))

(def chi-f
  (inel/distorted-wave-exit E-CM E-ex L-f V-params h r-max
                            :outgoing-type :p
                            :residual-A A-11
                            :residual-Z Z-11
                            :E-lab E-lab-exit
                            :W-params W-params
                            :s 0.5
                            :j (+ L-f 0.5)
                            :mass-factor mass-factor))

;; Transition potential on radial grid
(def n-steps (int (/ r-max h)))
(def V-trans-vec
  (mapv (fn [i] (inel/transition-potential-radial (* (double i) h) lambda mu beta-2 V-params))
        (range (inc n-steps))))

(def T-rad
  (inel/inelastic-amplitude-radial chi-i chi-f V-trans-vec r-max h))

(def angular-factor (* 4.0 Math/PI))
(def T-inel (if (number? T-rad)
              (* angular-factor T-rad)
              (cpx/mul angular-factor T-rad)))

(def k-i (Math/sqrt (* mass-factor E-CM)))
(def E-f (- E-CM E-ex))
(def k-f (Math/sqrt (* mass-factor E-f)))

(def dsigma-mb-sr
  (inel/inelastic-differential-cross-section T-inel k-i k-f E-CM E-ex mass-factor))

;; ---------------------------------------------------------------------------
;; Output
;; ---------------------------------------------------------------------------

(println "=== ¹¹Li(p,p′) — schematic DWBA (WS + Coulomb + imag WS) ===")
(println "")
(println (format "  E_lab(p) = %.2f MeV   E_CM = %.4f MeV   E_ex = %.2f MeV"
                 E-lab E-CM E-ex))
(println (format "  μ = %.4f MeV/c²   2μ/ℏ² = %.6e MeV⁻¹·fm⁻²" mu-reduced mass-factor))
(println (format "  λ = %d   β_λ = %.3f   L_i = %d → L_f = %d" lambda beta-2 L-i L-f))
(println (format "  V = [%.1f, %.2f, %.2f] MeV, fm   W = [%.1f, %.2f, %.2f] MeV, fm"
                 (first V-params) (second V-params) (nth V-params 2)
                 (first W-params) (second W-params) (nth W-params 2)))
(println "")
(println (format "  χ_i, χ_f grid: %d points, r_max = %.1f fm, h = %.3f fm"
                 (count chi-i) r-max h))
(println (format "  k_i = %.5f fm⁻¹   k_f = %.5f fm⁻¹" k-i k-f))
(println "")
(if (number? T-inel)
  (println (format "  T_inel (incl. 4π) = %.6e" T-inel))
  (println (format "  T_inel = %.6e + i %.6e"
                   (double (cpx/re T-inel)) (double (cpx/im T-inel)))))
(println (format "  dσ/dΩ (|T|² + kinematic prefactor; Coulomb-tail |u| norm) ≈ %.6e mb/sr"
                 dsigma-mb-sr))
(println "  (Values ~1e-20 mb/sr are not literal zero; open the _log.png for log10 y on this flat angular dist.)")
(println "")

;; Angular distribution **dσ/dΩ(θ)** (constant **T** from radial DWBA ⇒ flat curve; still defines scale)
(let [theta-step 2.0
      theta-rads (vec (map (fn [^double d] (* d (/ Math/PI 180.0)))
                           (range 0.0 (+ 181.0 (/ theta-step 2)) theta-step)))
      pairs (inel/inelastic-angular-distribution-function
             T-inel theta-rads k-i k-f E-CM E-ex mass-factor lambda mu)
      th-deg (mapv (fn [[th _]] (* (double th) (/ 180.0 Math/PI))) pairs)
      y-mb (mapv (fn [[_ sig]] (double sig)) pairs)
      out-path "output/11Li_pp_inelastic_dsigma.png"
      out-path-log "output/11Li_pp_inelastic_dsigma_log.png"]
  (try
    (io/make-parents (io/file out-path))
    (-> (c/xy-plot th-deg y-mb
                   :title "¹¹Li(p,p′) schematic — dσ/dΩ vs θ_CM (DWBA radial T; Coulomb-tail norm)"
                   :x-label "θ_CM (deg)"
                   :y-label "dσ/dΩ (mb/sr)"
                   :series-label "model"
                   :legend true)
        (i/save out-path :width 900 :height 520))
    (println (format "Plot saved: %s (%d points, Δθ = %.1f°, linear y)." out-path (count th-deg) theta-step))
    (let [y-log (mapv clamp-pos-for-log y-mb)]
      (-> (c/xy-plot th-deg y-log
                     :title "¹¹Li(p,p′) — dσ/dΩ vs θ_CM (log₁₀ y; same flat model)"
                     :x-label "θ_CM (deg)"
                     :y-label "dσ/dΩ (mb/sr), logarithmic"
                     :series-label "model"
                     :legend true)
          (chart-set-log-range-y! "dσ/dΩ (mb/sr)")
          (i/save out-path-log :width 900 :height 520)))
    (println (format "Plot saved: %s (log y)." out-path-log))
    (catch Exception e
      (println "Note: could not save plot:" (.getMessage e)))))

(println "")
(println "Done.")
