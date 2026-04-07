;; Elastic check: **p + ¹⁶O** (¹⁶O(p,p) elastic) at lab proton energy **E_lab(p)**.
;;
;; **Energy 35.2 MeV** matches **EXFOR O1198006** (Fabrici *et al.*, *Phys. Rev. C* **21**, **844** (1980));
;; saved page **`resources/nexus_nrs/X4sGetSubent.html`**. The vector **`exfor-o1198006-16opp-dsigma`** holds
;; **[θ_CM°, dσ/dΩ mb/sr]** for comparison.
;;
;; Optional **EXFOR fit:** with **`fit-to-exfor?` true**, minimize weighted **relative** MSE vs EXFOR, varying only
;; **W_r** = real WS depth **V₀** (MeV) and **W_i** = imag volume depth (MeV); **R = 2**, **a = 0.6** fm fixed.
;;
;; Mirrors **`examples/example_alpha_Sm148_elastic.clj`** (ratio plots, **s-matrix** table, Rutherford table).
;; Physics notes: **`differential-cross-section-nuclear-cut`** = **10 |f_C + f_N|²** mb/sr (T&N);
;; **`s-matrix`** = **S_L^n**; **σ** only in **e^{2iσ}(S^n−1)** for **f_N**. Optional **`*elastic-imag-ws-params*`**
;; matches the web dashboard elastic tab.
;;
;; **Partial-wave cutoff:** for p + light targets, **f_N** can grow erratically if **L_nuclear_cut** is large
;; (Coulomb–Hankel matching). Default **L_nuclear_cut = 22** aligns with the dashboard; raise only with care.
;;
;; Real Woods–Saxon **V = [V₀ R a]** (MeV, fm, fm). Point Coulomb only: **V₀ = R = 0**, **a > 0** for matching radius.
;;
;; Run from repository root:
;;   lein run -m clojure.main examples/example_16Opp_elastic.clj
;;
;; Writes (typically gitignored):
;;   **output/16Opp_elastic_sigma_ratio_Ruth.png** — dσ/σ_Ruth vs θ_CM: **model** + **EXFOR** (linear y)
;;   **output/16Opp_elastic_sigma_ratio_Ruth_log.png** — same, log₁₀ y-axis
;;   **output/16Opp_elastic_vs_EXFOR_dsigma.png** — dσ/dΩ (mb/sr): model curve + EXFOR points

(import '[org.jfree.chart.axis LogAxis])

(require '[functions :refer [differential-cross-section-nuclear-cut mass-factor mass-factor-from-mu
                             s-matrix phase-shift Z1Z2ee *elastic-imag-ws-params*]]
         :reload)
(require '[complex :as cpx]
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
    (.setSmallestValue axis 1e-35)
    (.setRangeAxis plot 0 axis)
    chart))

(def ^:private m-p 938.272)       ; MeV/c² (dashboard elastic handler)
(def ^:private A-16 16)
(def ^:private Z-16 8)
(def ^:private m-16O (* 931.5 A-16)) ; MeV/c² (same convention as dashboard)

;; **35.2 MeV** — EXFOR **O1198006** (¹⁶O(p,el), lab EN in subentry COMMON).
(def E-lab 35.2)

;; EXFOR **O1198006** — pairs **[θ_CM (deg), dσ/dΩ (mb/sr)]** from `DATA` block (ANG-CM, DATA-CM).
(def exfor-o1198006-16opp-dsigma
  [[21.2 771.0] [27.6 398.0] [30.1 272.0] [33.3 167.0] [36.0 108.0] [39.1 58.0]
   [42.3 36.7] [45.4 24.9] [51.7 21.5] [57.9 22.6] [64.1 18.7] [70.3 13.2]
   [73.4 9.61] [79.5 5.38] [85.0 3.25] [90.6 2.09] [96.6 1.57] [102.5 1.26]
   [108.4 0.941] [114.3 0.584] [120.2 0.341] [126.0 0.219] [131.8 0.212]
   [137.5 0.355] [143.3 0.496] [148.9 0.566] [154.6 0.540] [159.8 0.447]
   [164.9 0.303] [170.6 0.180]])

;; CM kinetic energy: target at rest, non-relativistic
(def E-CM (* E-lab (/ m-16O (+ m-p m-16O))))

(def mu-reduced (/ (* m-p m-16O) (+ m-p m-16O)))

;; Z₁ Z₂ e² with e² = 1.44 MeV·fm
(def z1z2ee (* 1 Z-16 1.44))

(def ^:private R-ws 2.0)
(def ^:private a-ws 0.6)

;; Set **true** to scan **W_r**, **W_i** vs EXFOR (slow); default **false** uses manual params below.
(def fit-to-exfor? false)

;; --- EXFOR fit: vary **W_r** (‑V real depth **V₀**) and **W_i** (imag volume **W₀**); geometry fixed. ---
(defn- v-params-from-wr [^double wr]
  [wr R-ws a-ws])

(defn- dsigma-model-mb
  "dσ/dΩ (mb/sr) at **theta-cm** rad for given **W_r**, **W_i**."
  [^double wr ^double wi ^double theta-rad ^long l-cut]
  (binding [mass-factor (mass-factor-from-mu mu-reduced)
            Z1Z2ee z1z2ee
            *elastic-imag-ws-params* (when (pos? wi) [wi R-ws a-ws])]
    (let [ds-c (differential-cross-section-nuclear-cut E-CM (v-params-from-wr wr) theta-rad l-cut)
          r (double (cpx/re ds-c))
          im (double (cpx/im ds-c))]
      (Math/sqrt (+ (* r r) (* im im))))))

(defn- exfor-relative-mse
  "Mean squared **relative** error vs EXFOR; **σ_i = max(5%·y_i, 0.05 mb/sr)** for weights."
  [^double wr ^double wi ^long l-cut]
  (let [ms (for [[^double th-deg ^double y-exp] exfor-o1198006-16opp-dsigma
                 :let [th (* th-deg (/ Math/PI 180.0))
                       y-mod (dsigma-model-mb wr wi th l-cut)
                       σ (max (* 0.05 y-exp) 0.05)]]
             (if (and (Double/isFinite y-mod) (pos? y-exp))
               (let [d (/ (- y-mod y-exp) σ)]
                 (* d d))
               1e12))]
    (/ (double (reduce + ms)) (count exfor-o1198006-16opp-dsigma))))

(defn- grid-min-2d
  "Evaluate **f [wr wi]** on a rectangular **wr** × **wi** grid; return `{:wr .. :wi .. :val ..}`."
  [f wrs wis]
  (reduce (fn [best [wr wi]]
            (let [v (f wr wi)]
              (if (< v (:val best)) {:wr wr :wi wi :val v} best)))
          {:wr (first wrs) :wi (first wis) :val Double/MAX_VALUE}
          (for [wr wrs wi wis] [wr wi])))

(defn- fit-wr-wi-to-exfor
  "Global medium-resolution grid on (**W_r**, **W_i**) — a tight local refinement can miss the physical basin
  if the coarse minimum lies in the wrong valley (e.g. low **W_r**). One **fine** pass around the best cell."
  [^long l-cut]
  (let [f (fn [wr wi] (exfor-relative-mse wr wi l-cut))
        ;; ~16×15 ≈ 240 global + ~100 local; **s-matrix** memoizes in **L** at fixed (**W_r**, **W_i**).
        phase1 (grid-min-2d f
                            (range 10.0 88.1 5.0)
                            (concat [0.0] (range 2.0 44.1 3.0)))
        {:keys [wr wi]} phase1
        wr-lo (max 8.0 (- wr 4.0)) wr-hi (min 95.0 (+ wr 4.0))
        wi-lo (max 0.0 (- wi 3.0)) wi-hi (min 48.0 (+ wi 3.0))
        phase2 (grid-min-2d f
                            (range wr-lo (+ wr-hi 0.001) 0.75)
                            (range wi-lo (+ wi-hi 0.001) 0.6))]
    phase2))

;; Manual fallbacks when **`fit-to-exfor?`** is false (**W_r**, **W_i** in MeV; same **R**, **a** as **`R-ws` / `a-ws`**):
;(def V-params [0.0 0.0 0.67])   ; point Coulomb only → dσ/σ_Ruth ≡ 1
(def V-params-manual [8.0 R-ws a-ws])
(def imag-ws-params-manual [30.0 R-ws a-ws])

(defn- imag-binding [params]
  (when (and (sequential? params) (pos? (double (first params))))
    params))

(def L-max 25)            ; print **s-matrix** for L = 0 … L-max
(def L-nuclear-cut 22)  ; **f_N** sum L = 0 … L_nuclear_cut (try 22–39; see header)

(def fit-wr-wi-result
  (when fit-to-exfor?
    (println "")
    (println "=== Fit to EXFOR O1198006: **W_r** (real depth), **W_i** (imag volume); R,a fixed ===")
    (flush)
    (let [r (fit-wr-wi-to-exfor L-nuclear-cut)]
      (println (format "  Best: W_r = %.3f MeV  W_i = %.3f MeV  (weighted relative MSE = %.5f)"
                       (:wr r) (:wi r) (:val r)))
      (println "")
      r)))

(def V-params
  (if fit-wr-wi-result
    [(double (:wr fit-wr-wi-result)) R-ws a-ws]
    V-params-manual))

(def imag-ws-params
  (if fit-wr-wi-result
    [(double (:wi fit-wr-wi-result)) R-ws a-ws]
    imag-ws-params-manual))

(defn rutherford-dsigma-mb-sr
  "Point-charge Rutherford dσ/dΩ (**mb/sr**), non-relativistic CM."
  [^double z1z2e2 ^double e-cm ^double theta-rad]
  (* 10.0
     (let [s (Math/sin (* 0.5 theta-rad))]
       (if (< (Math/abs s) 1e-15)
         Double/NaN
         (/ (Math/pow (/ z1z2e2 (* 4.0 e-cm)) 2.0)
            (Math/pow s 4.0))))))

(println "=== p + ¹⁶O elastic (¹⁶O(p,p)) ===")
(println "")
(println (format "  Projectile: p  (m = %.3f MeV/c², Z = 1)" m-p))
(println (format "  Target:     ¹⁶O (m ≈ %.1f MeV/c², Z = %d)" m-16O Z-16))
(println (format "  E_lab(p)    = %.2f MeV" E-lab))
(println (format "  E_CM        = %.4f MeV" E-CM))
(println (format "  μ (reduced) = %.4f MeV/c²" mu-reduced))
(println (format "  Z₁Z₂ e²     = %.4f MeV·fm" z1z2ee))
(println (format "  Real V:     [V₀ R a]     = [%.1f %.2f %.2f] (MeV, fm, fm)"
                 (first V-params) (second V-params) (nth V-params 2)))
(if-let [iw (imag-binding imag-ws-params)]
  (println (format "  Imag WS:    [Wᵢ Rᵢ aᵢ] = [%.1f %.2f %.2f] MeV, fm, fm"
                   (first iw) (second iw) (nth iw 2)))
  (println "  Imag WS:    (none — real potential only)"))
(println (format "  L_nuclear_cut = %d (partial waves in **f_N**)" L-nuclear-cut))
(println "")
(println "=== Point-charge Rutherford — spot checks (mb/sr) ===")
(println "")
(let [pref (/ z1z2ee (* 4.0 E-CM))
      pref2 (* pref pref)]
  (println (format "  (Z₁Z₂e² / (4 E_cm))² scale (×10 → mb/sr): %.6f (fm²)" pref2))
  (println "")
  (doseq [deg [10 90]]
    (let [th (* deg (/ Math/PI 180.0))
          s (Math/sin (* 0.5 th))
          s4 (Math/pow s 4)
          sig (/ pref2 s4)]
      (println (format "  θ_cm = %3d° → σ_Ruth = %.5f mb/sr" deg (* 10.0 sig))))))
(println "")

(binding [mass-factor (mass-factor-from-mu mu-reduced)
          Z1Z2ee z1z2ee
          *elastic-imag-ws-params* (imag-binding imag-ws-params)]
  (println "=== **`s-matrix`** = **S_L^n** (sample L = 0 … " L-max ") ===")
  (println "")
  (println "   L  |   Re(S_L^n)  |   Im(S_L^n)  | |S^n| | arg(S^n) rad |  δ (deg)")
  (println " -----+------------+------------+-------+-------------+----------")
  (doseq [L (range (inc L-max))]
    (try
      (let [S-val (s-matrix E-CM V-params L)
            reS (double (cpx/re S-val))
            imS (double (cpx/im S-val))]
        (if (or (Double/isNaN reS) (Double/isNaN imS)
                (Double/isInfinite reS) (Double/isInfinite imS))
          (println (format " %3d | (non-finite S^n)" L))
          (let [Sm (cpx/mag S-val)
                d-tot-rad (phase-shift E-CM V-params L)
                d-deg (* d-tot-rad (/ 180.0 Math/PI))]
            (println (format " %3d | %10.5f | %10.5f | %.5f | %11.5f | %9.3f"
                             L reS imS Sm (cpx/arg S-val) d-deg)))))
      (catch Exception e
        (println (format " %3d | ERROR: %s" L (.getMessage e))))))
  (println "")
  (println "=== Elastic dσ/dΩ vs Rutherford (**L ≤ " L-nuclear-cut "** in **f_N**) ===")
  (println "")
  (println "  θ_cm (deg) | dσ/dΩ / σ_Ruth | dσ/dΩ (mb/sr)")
  (println " ------------+------------------+---------------")
  (doseq [theta-deg [10 20 30 45 60 90 120 150 170]]
    (try
      (let [theta (* theta-deg (/ Math/PI 180.0))
            ds-c (differential-cross-section-nuclear-cut E-CM V-params theta L-nuclear-cut)
            r (double (cpx/re ds-c))
            i (double (cpx/im ds-c))
            ds-mb-sr (Math/sqrt (+ (* r r) (* i i)))
            ruth (rutherford-dsigma-mb-sr z1z2ee E-CM theta)
            ratio (if (and (Double/isFinite ruth) (> ruth 0.0))
                    (/ ds-mb-sr ruth)
                    Double/NaN)]
        (if (Double/isFinite ratio)
          (println (format "   %6.1f    |     %12.5f | %12.5e"
                           (double theta-deg) ratio ds-mb-sr))
          (println (format "   %6.1f    |     (n/a)      | %12.5e"
                           (double theta-deg) ds-mb-sr))))
      (catch Exception e
        (println (format "   %6.1f    | ERROR: %s" theta-deg (.getMessage e))))))
  (println "")
  (let [theta-step 0.25
        pairs (for [theta-deg (range 10.0 171.0 theta-step)]
                (try
                  (let [^double td theta-deg
                        theta (* td (/ Math/PI 180.0))
                        ds-c (differential-cross-section-nuclear-cut E-CM V-params theta L-nuclear-cut)
                        r (double (cpx/re ds-c))
                        im (double (cpx/im ds-c))
                        ds-mb-sr (Math/sqrt (+ (* r r) (* im im)))
                        ruth (rutherford-dsigma-mb-sr z1z2ee E-CM theta)
                        ratio (when (and (Double/isFinite ruth) (> ruth 0.0))
                                (/ ds-mb-sr ruth))]
                    (when (and ratio (Double/isFinite ratio))
                      [td (double ratio)]))
                  (catch Exception _ nil)))
        th-vec (vec (keep (fn [p] (when p (first p))) pairs))
        y-vec (vec (keep (fn [p] (when p (second p))) pairs))
        exfor-ratio-pts
        (vec
         (keep
          (fn [[^double th-deg ^double exp-mb]]
            (let [theta (* th-deg (/ Math/PI 180.0))
                  ruth (rutherford-dsigma-mb-sr z1z2ee E-CM theta)]
              (when (and (pos? exp-mb) (Double/isFinite ruth) (> ruth 0.0))
                [th-deg (/ exp-mb ruth)])))
          exfor-o1198006-16opp-dsigma))
        th-exp-r (mapv first exfor-ratio-pts)
        y-exp-r (mapv second exfor-ratio-pts)]
    (when (and (seq th-vec) (= (count th-vec) (count y-vec)))
      (try
        (io/make-parents (io/file "output/16Opp_elastic_sigma_ratio_Ruth.png"))
        (let [y-log (mapv clamp-pos-for-log y-vec)
              y-exp-r-log (mapv clamp-pos-for-log y-exp-r)]
          (-> (c/xy-plot th-vec y-vec
                         :title "p + ¹⁶O elastic — σ/σ_Ruth (CM); model + EXFOR O1198006"
                         :x-label "θ_CM (deg)"
                         :y-label "σ / σ_Ruth"
                         :series-label "model"
                         :legend true)
              (c/add-lines th-exp-r y-exp-r :series-label "EXFOR O1198006")
              (i/save "output/16Opp_elastic_sigma_ratio_Ruth.png" :width 800 :height 500))
          (-> (c/xy-plot th-vec y-log
                         :title "p + ¹⁶O elastic — σ/σ_Ruth (CM), log y; model + EXFOR"
                         :x-label "θ_CM (deg)"
                         :y-label "σ / σ_Ruth (log scale)"
                         :series-label "model"
                         :legend true)
              (c/add-lines th-exp-r y-exp-r-log :series-label "EXFOR O1198006")
              (chart-set-log-range-y! "σ / σ_Ruth (log scale)")
              (i/save "output/16Opp_elastic_sigma_ratio_Ruth_log.png" :width 800 :height 500)))
        (println (format (str "Plots saved: output/16Opp_elastic_sigma_ratio_Ruth.png, …_log.png "
                              "(%d model pts Δθ=%.2f° + %d EXFOR σ/σ_Ruth).")
                         (count th-vec) (double theta-step) (count th-exp-r)))
        (catch Exception e
          (println "Note: could not save ratio plot:" (.getMessage e))))))

  ;; --- EXFOR O1198006 (35.2 MeV lab): printed table + dσ overlay plot ---
  (println "")
  (println "=== Theory vs EXFOR **O1198006** (Fabrici et al., PRC 21 (1980) 844); E_lab = 35.2 MeV ===")
  (println (if fit-wr-wi-result
             (format "Model: **fitted** W_r=%.2f, W_i=%.2f MeV (R,a)=(%.1f,%.1f) fm for both wells."
                     (double (:wr fit-wr-wi-result)) (double (:wi fit-wr-wi-result)) R-ws a-ws)
             (format "Model: **manual** W_r=%.1f, W_i=%.1f MeV (R,a)=(%.1f,%.1f) fm for both wells."
                     (double (first V-params)) (double (first imag-ws-params)) R-ws a-ws)))
  (println "")
  (println "  θ_CM (deg) | EXFOR mb/sr | model mb/sr | model/EXFOR")
  (println " ------------+-------------+-------------+------------")
  (doseq [[^double th-deg ^double exp-mb] exfor-o1198006-16opp-dsigma]
    (try
      (let [theta (* th-deg (/ Math/PI 180.0))
            ds-c (differential-cross-section-nuclear-cut E-CM V-params theta L-nuclear-cut)
            r (double (cpx/re ds-c))
            i (double (cpx/im ds-c))
            model (Math/sqrt (+ (* r r) (* i i)))
            q (if (pos? exp-mb) (/ model exp-mb) Double/NaN)]
        (println (format "    %6.2f    | %11.3f | %11.3f | %10.3f"
                         th-deg exp-mb model q)))
      (catch Exception e
        (println (format "    %6.2f    | ERROR: %s" th-deg (.getMessage e))))))
  (println "")
  (try
    (let [step 0.5
          th-dense (vec (range 10.0 171.0 step))
          model-dense
          (mapv (fn [^double td]
                  (let [theta (* td (/ Math/PI 180.0))
                        ds-c (differential-cross-section-nuclear-cut E-CM V-params theta L-nuclear-cut)
                        r (double (cpx/re ds-c))
                        im (double (cpx/im ds-c))]
                    (Math/sqrt (+ (* r r) (* im im)))))
                th-dense)
          th-exp (mapv (fn [[a _]] (double a)) exfor-o1198006-16opp-dsigma)
          y-exp (mapv (fn [[_ y]] (double y)) exfor-o1198006-16opp-dsigma)
          chart (-> (c/xy-plot th-dense model-dense
                             :title "¹⁶O(p,p) elastic @ 35.2 MeV lab — dσ/dΩ (CM); model vs EXFOR O1198006"
                             :x-label "θ_CM (deg)"
                             :y-label "dσ/dΩ (mb/sr)"
                             :series-label (format "model W_r=%.1f W_i=%.1f (R,a)=(%.1f,%.1f), L≤%d"
                                                   (first V-params) (first imag-ws-params)
                                                   R-ws a-ws L-nuclear-cut)
                             :legend true)
                    (c/add-lines th-exp y-exp :series-label "EXFOR O1198006"))]
      (io/make-parents (io/file "output/16Opp_elastic_vs_EXFOR_dsigma.png"))
      (i/save chart "output/16Opp_elastic_vs_EXFOR_dsigma.png" :width 900 :height 520)
      (println (format "Plot saved: output/16Opp_elastic_vs_EXFOR_dsigma.png (model Δθ = %.1f°)." step)))
    (catch Exception e
      (println "Note: could not save EXFOR comparison plot:" (.getMessage e))))

  (println "")
  (println "Done."))
