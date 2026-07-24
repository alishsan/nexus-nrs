;; ¹⁶O(α,α′)¹⁶O* — **6.917 MeV 2⁺** (handbook ~p. 74, **Table 7.1**). Lab energy **40 MeV** (= **10 MeV/u** for α).
;; Collective **λ = 2** coupling: outgoing orbital **L_β** satisfies **|L_α − 2| ≤ L_β ≤ L_α + 2** with
;; **L_α + L_β** even (parity). **χ_α(L_α)** uses **incoming** optical **U**; **χ_β(L_β)** uses **outgoing** **U**
;; (different kinematic energy in Table 7.1). Transition form factor **F_λ ∝ β_λ R dU/dr** uses the **incoming**
;; channel **real** Woods–Saxon only (**dU/dr** from entrance well), not the exit potential.
;;
;; Intended benchmark: **Ptolemy** (collective inelastic). This script is Nexus radial DWBA + a reduced angular
;; factor (below); tune **β₂** to data / **B(E2)**.
;;
;; Distorted waves: **`optical-potential-woods-saxon`** (volume + derivative surface imag), **α**, **s = 0**,
;; **r_c A^(1/3)** Coulomb. Partial waves use **`:distorted-norm :coulomb-tail`** + **`:coulomb-Z-pair [2 8]`** so
;; **L**-dependent norms are not scrambled (**`:max`** breaks coherent **(5.6)** sums and lifts large-angle **dσ/dΩ**).
;;
;; **dσ/dΩ(θ)** — two levels:
;; 1. **Full (handbook / Austern Eq. (5.6)):** **`handbook-zr-multipole-amplitude-sum`** on rows
;; **`{:L-alpha :L-beta :I}`** with **:I** = collective radial **T_{L_α L_β}** from **`inelastic-amplitude-radial`**
;; (**∫ u_f^* F_λ u_i dr**, **u = r χ**). We **do not** call **`handbook-zr-rows-with-coulomb-sigma`**: explicit
;; **e^{i(σ_α+σ_β)}** in the multipole sum is for **bare ZR **I** from (5.5)**; our **T** already carries Coulomb
;; phases through the distorted **u** (same issue as double-phasing if **σ** were folded twice). Multipole **ℓ = λ**,
;; **Σ_m |β^{ℓm}|²**, then **`inelastic-differential-cross-section`** on **√(Σ|β|²)**. **Absolute scale / phases** vs
;; **Ptolemy** may still need **χ** normalization (**`:distorted-norm`**) tuning.
;; 2. **Reduced diagnostic:** coherent **Y_{L_β 0}** sum (**`transfer-differential-cross-section-angular-coherent`**)
;; for shape comparison only.
;;
;; Writes **output/o16_alpha_inelastic_2plus_TL.tsv** (`L_alpha`, `L_beta`, …) and PNGs:
;; **…_dU_dr.png** — **Re dU/dr** and **Im dU/dr** for the **full** incoming optical **U(r)** (Table 7.1:
;; real WS + volume/surface imaginary + Coulomb). Collective **F_λ** in this script still uses **`transition-form-factor`**
;; (**real** WS **dU/dr** only), not this full complex derivative.
;; **…_TL_vs_L.png**, **…_dsigma_vs_L.png**, **…_dsigma_theta.png**, **…_dsigma_theta_log.png**.
;;
;; Run: `lein run -m clojure.main examples/example_O16_alpha_inelastic_2plus.clj`

(require '[clojure.java.io :as io]
         '[dwba.inelastic :as inel]
         '[dwba.transfer :as t]
         '[functions :as fn :refer [mass-factor-from-mu]]
         '[complex :refer [mag re im add]]
         '[incanter.core :as i]
         '[incanter.charts :as c])

(import '[org.jfree.chart.axis LogAxis])

(def ^:private m-alpha 3727.379) ; MeV/c²
(def ^:private m-16O (* 16.0 931.494)) ; MeV/c² (approx.)

(def E-lab-alpha 40.0) ; MeV — Table 7.1 (≈ 10 MeV/u for α)
(def lambda-ex 2)      ; 2⁺ excitation, collective quadrupole
(def E-ex 6.917)       ; MeV — 2⁺ state

;; Table 7.1 — entrance (α + ¹⁶O g.s.)
(def V0-ent 150.355)
(def Wi-ent 1.619)   ; volume imag depth (MeV), matches −Vᵢ i in handbook notation
(def Ws-ent 24.394)  ; surface imag strength (MeV), derivative WS

;; Table 7.1 — exit (α + ¹⁶O* at E_i − E*)
(def V0-ex 154.365)
(def Wi-ex 0.644)
(def Ws-ex 24.945)

;; Geometry: r₀, r₀ᵢ, rₛ, r_c in fm as **coefficients × A^(1/3)** (handbook style)
(def r0-real 1.342)
(def r0-imag 1.426)
(def r-surf 1.293)
(def r-coul 1.350)

(def a0 0.658)
(def a0-imag 0.658)
(def a-surf 0.636)

;; Collective deformation **β₂** — not in your paste; adjust to match handbook Fig. / B(E2).
(def beta-2 0.295)

(def A 16)
(def Z-target 8)

(def h-step 0.05)
(def r-max 40.0)
;; Partial-wave cutoff: **dσ/dΩ(θ)** from **(5.6)** needs many **(L_α, L_β)** pairs; **L ≲ 24** truncates
;; the sum too early and shifts/peaks the distribution (e.g. forward angle); **30** matches converged shape.
(def L-max 30)

(defn- A13 ^double [^double a] (Math/pow a (/ 1.0 3.0)))

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

(defn- allowed-L-betas
  "For multipole **λ**, exit partial waves **L_β** with **|L_α−λ| ≤ L_β ≤ L_α+λ** and **L_α+L_β** even."
  [^long L-alpha ^long lambda]
  (let [lo (Math/abs (- L-alpha lambda))
        hi (+ L-alpha lambda)]
    (filter (fn [^long Lb]
              (and (<= lo Lb hi)
                   (zero? (mod (+ L-alpha Lb) 2))))
            (range lo (inc hi)))))

(defn- T-map-coherent-by-L-beta
  "Complex sum of **T_{L_α L_β}** sharing the same exit **L_β** (for reduced **Y_{L_β 0}** angular sum)."
  [rows]
  (->> rows
       (group-by :L-beta)
       (map (fn [[Lb pair-rows]]
              [Lb (reduce add (map :T pair-rows))]))
       (into {})))

(defn- dsigma-austern-56-m-sum
  "dσ/dΩ (**mb/sr**) at CM **theta-rad**: **Σ_{m=-ℓ}^ℓ |β^{ℓm}(θ)|²** with **β** from **`handbook-zr-multipole-amplitude-sum`**,
  then **`inelastic-differential-cross-section`** on **√(Σ|β|²)** (unpolarized **m** sum)."
  [theta-rad ell radial-rows-sigma k-in k-out E-cm E-ex mf]
  (let [th (double theta-rad)
        l (long ell)
        T-sq (reduce (fn [^double acc ^long mm]
                       (+ acc (let [b (t/handbook-zr-multipole-amplitude-sum l mm th radial-rows-sigma)
                                     m (mag b)]
                                (* m m))))
                     0.0
                     (range (- l) (inc l)))]
    (double (inel/inelastic-differential-cross-section (Math/sqrt (max T-sq 1e-300))
                                                       k-in k-out E-cm E-ex mf))))

(def R-V (* r0-real (A13 A)))
(def R-Wi (* r0-imag (A13 A)))
(def R-S (* r-surf (A13 A)))
(def R-C (* r-coul (A13 A)))

(def v-params-real-ent [V0-ent R-V a0])
(def w-vol-ent [Wi-ent R-Wi a0-imag])
(def w-surf-ent [Ws-ent R-S a-surf])

(def v-params-real-ex [V0-ex R-V a0])
(def w-vol-ex [Wi-ex R-Wi a0-imag])
(def w-surf-ex [Ws-ex R-S a-surf])

(defn- U-ent ^Object [^double r]
  (t/optical-potential-woods-saxon r v-params-real-ent w-vol-ent w-surf-ent
                                   nil nil nil 0 0.0 0.0 2 Z-target R-C))

(defn- U-ex ^Object [^double r]
  (t/optical-potential-woods-saxon r v-params-real-ex w-vol-ex w-surf-ex
                                   nil nil nil 0 0.0 0.0 2 Z-target R-C))

(defn- ws-df-dr
  "df/dr for **f = 1/(1 + exp((r−R)/a))** = **−f(1−f)/a** (surface peaked toward **r ≈ R**)."
  ^double [^double r ^double R ^double a]
  (let [x (/ (- r R) a)
        ex (Math/exp x)
        den (+ 1.0 ex)
        f (/ 1.0 den)]
    (- (/ (* f (- 1.0 f)) a))))

(defn- d-re-U-entrance-dr
  "d/dr Re **U_ent(r)** — real Woods–Saxon + Coulomb (same branch as **`optical-potential-woods-saxon`**; α, no l·s)."
  ^double [^double r]
  (let [[V0 Rv av] v-params-real-ent
        ;; V_real = −V0 f_V
        dVreal (- (* V0 (ws-df-dr r Rv av)))
        Zeff (* 2.0 (double Z-target) 1.44)
        dVC (if (> r R-C)
              (/ (- Zeff) (* r r))
              (/ (* (- Zeff) r) (* R-C R-C R-C)))]
    (+ dVreal dVC)))

(defn- d-im-U-entrance-dr
  "d/dr Im **U_ent(r)** — volume + derivative surface imaginary (same as **`optical-potential-woods-saxon`**)."
  ^double [^double r]
  (let [[W0 Rw aw] w-vol-ent
        [Wd Rs as-surf] w-surf-ent
        ;; W_vol = −W0 f_W  ⇒  d/dr = −W0 · df_W/dr
        dWvol (- (* W0 (ws-df-dr r Rw aw)))
        fS (/ 1.0 (+ 1.0 (Math/exp (/ (- r Rs) as-surf))))
        ;; W_surf = −4 W_D f(1−f);  d/dr = −4 W_D · d[f(1−f)]/dr = −4 W_D (1−2f) df/dr
        dWsurf (* -4.0 Wd (ws-df-dr r Rs as-surf) (- 1.0 (* 2.0 fS)))]
    (+ dWvol dWsurf)))

(let [mu (/ (* m-alpha m-16O) (+ m-alpha m-16O))
      mf (mass-factor-from-mu mu)
      E-cm (* E-lab-alpha (/ m-16O (+ m-alpha m-16O)))
      E-out (- E-cm E-ex)
      k-in (Math/sqrt (* mf E-cm))
      k-out (Math/sqrt (* mf E-out))]
  (println "=== ¹⁶O(α,α′)¹⁶O* 2⁺ (6.917 MeV) — 10 MeV/u / 40 MeV lab, Table 7.1 ===")
  (println "")
  (println (format "  E_lab(α)   = %.3f MeV" E-lab-alpha))
  (println (format "  E_cm       = %.4f MeV" E-cm))
  (println (format "  E_f = E_cm − E_x = %.4f MeV" E-out))
  (println (format "  μ          = %.4f MeV/c²" mu))
  (println (format "  mass-factor (2μ/ℏ²) = %.6f MeV⁻¹·fm⁻²" mf))
  (println (format "  k_in, k_out = %.4f, %.4f fm⁻¹" k-in k-out))
  (println (format "  R_V, R_C   = %.4f, %.4f fm" R-V R-C))
  (println (format "  β₂ (tune!) = %.4f" beta-2))
  (println "")
  (binding [fn/mass-factor mf]
    (let [rows
          (for [L-a (range 0 (inc L-max))
                L-b (allowed-L-betas L-a lambda-ex)
                :let [                      chi-i (inel/distorted-wave-entrance E-cm L-a nil h-step r-max
                                                          :optical-potential-fn U-ent
                                                          :mass-factor mf
                                                          :s 0
                                                          ;; **:coulomb-tail** keeps correct relative **L** weights for coherent
                                                          ;; multipole sums; **:max** rescales each **L** separately and inflates tails.
                                                          :distorted-norm :coulomb-tail
                                                          :coulomb-Z-pair [2 Z-target])
                      chi-f (inel/distorted-wave-exit E-cm E-ex L-b nil h-step r-max
                                                      :optical-potential-fn U-ex
                                                      :mass-factor mf
                                                      :s 0
                                                      :distorted-norm :coulomb-tail
                                                      :coulomb-Z-pair [2 Z-target])
                      n (min (count chi-i) (count chi-f))
                      ;; **dU/dr** from **incoming** real WS only (Table 7.1 entrance **V**).
                      Vtrans (mapv (fn [i]
                                     (inel/transition-form-factor (* (double i) h-step)
                                                                  lambda-ex beta-2 v-params-real-ent))
                                   (range n))
                      T (inel/inelastic-amplitude-radial chi-i chi-f Vtrans r-max h-step)
                      Tmag (mag T)]]
            {:L-alpha L-a :L-beta L-b :T T :mag Tmag})
          L-nums (vec (range 0 (inc L-max)))
          ;; largest |T| among (L_α, L_β) pairs
          best (apply max-key :mag rows)
          mag-max-by-La
          (mapv (fn [^long la]
                  (let [ms (->> rows (filter #(= la (:L-alpha %))) (map :mag))]
                    (if (seq ms) (apply max ms) 0.0)))
                L-nums)
          ds-best (inel/inelastic-differential-cross-section (:T best) k-in k-out E-cm E-ex mf)
          ;; naive incoherent sum of |T|² over pairs (diagnostic only)
          sum-T2 (reduce (fn [^double acc r] (+ acc (* (:mag r) (:mag r)))) 0.0 rows)
          ds-naive (inel/inelastic-differential-cross-section (Math/sqrt sum-T2) k-in k-out E-cm E-ex mf)]
      (println (format "  λ = %d: for each L_α, L_β ∈ [|L_α−λ|, L_α+λ] with L_α+L_β even." lambda-ex))
      (println "  Partial waves (L_α → χ_in, L_β → χ_out); F_λ from incoming dU/dr:")
      (doseq [r (take 16 rows)]
        (let [arg-deg (* (/ 180.0 Math/PI)
                         (Math/atan2 (im (:T r)) (re (:T r))))]
          (println (format "  L_α=%2d  L_β=%2d  |T|=% .4e  arg/deg=% .2f"
                           (:L-alpha r) (:L-beta r) (:mag r)
                           arg-deg))))
      (when (> (count rows) 16)
        (println "  ..."))
      (println "")
      (println (format "  Largest |T| at (L_α,L_β) = (%d,%d)  (|T| = %.4e)"
                       (:L-alpha best) (:L-beta best) (:mag best)))
      (println (format "  dσ/dΩ from that L only (illustrative): %.4e mb/sr" (double ds-best)))
      (println (format "  dσ/dΩ from √(Σ |T_{L_α L_β}|²) over pairs (naive, not physical angular dist): %.4e mb/sr"
                       (double ds-naive)))
      (println "")
      (println "  dσ/dΩ(θ): primary curve = Austern (5.6) Σ_m |β^{ℓm}|² + inelastic prefactor; reduced Y_{L_β 0} = diagnostic.")
      (println "  dU/dr plot: Re and Im for full incoming optical U (WS + Im + Coulomb); F_λ still uses real WS dU/dr only.")
      (println "  Export TSV for DWUCK/ECIS partial-wave listings.")
      (try
        (let [T-map (T-map-coherent-by-L-beta rows)
              ;; **Partial reduced:** top **two** **L_β** columns by **|Σ_{L_α} T|** (coherent **Y_{L_β0}**).
              ;; One column alone often looks pathological (zeros from a single **Y_{L0}**; no **L_β–L_β′** terms).
              top2-Lb (take 2 (sort-by (fn [[_lb tc]] (- (mag tc))) (seq T-map)))
              T-map-dom (into {} top2-Lb)
              ds-per-La (mapv (fn [^long la]
                                (let [ds (->> rows
                                              (filter #(= la (:L-alpha %)))
                                              (map (fn [r]
                                                     (double (inel/inelastic-differential-cross-section
                                                              (:T r) k-in k-out E-cm E-ex mf)))))]
                                  (if (seq ds) (apply max ds) 0.0)))
                              L-nums)
              theta-step 2.0
              theta-rads (vec (map (fn [^double d] (* d (/ Math/PI 180.0)))
                                   (range 0.0 (+ 181.0 (/ theta-step 2)) theta-step)))
              th-deg (mapv (fn [^double th] (* th (/ 180.0 Math/PI))) theta-rads)
              y-theta (mapv (fn [^double th]
                              (double (t/transfer-differential-cross-section-angular-coherent
                                       T-map 1.0 k-in k-out th mf mf)))
                            theta-rads)
              y-theta-dominant-L
              (mapv (fn [^double th]
                      (double (t/transfer-differential-cross-section-angular-coherent
                               T-map-dom 1.0 k-in k-out th mf mf)))
                    theta-rads)
              ;; No **`handbook-zr-rows-with-coulomb-sigma`** — see file header (**T** vs bare ZR **I**).
              radial-rows-austern
              (mapv (fn [row] {:L-alpha (:L-alpha row) :L-beta (:L-beta row) :I (:T row)}) rows)
              ell (long lambda-ex)
              y-theta-austern
              (mapv (fn [^double th]
                      (dsigma-austern-56-m-sum th ell radial-rows-austern k-in k-out E-cm E-ex mf))
                    theta-rads)
              ;; **dU/dr** for full incoming optical **U** = Re U + i Im U (Table 7.1 entrance).
              dr-plot 0.06
              r-plot (vec (take-while #(<= % (+ r-max 1e-9))
                                      (iterate (fn [^double x] (+ x dr-plot)) 0.0)))
              dU-dr-re (mapv (fn [^double r] (double (d-re-U-entrance-dr r))) r-plot)
              dU-dr-im (mapv (fn [^double r] (double (d-im-U-entrance-dr r))) r-plot)]
          (io/make-parents (io/file "output/o16_alpha_inelastic_2plus_dU_dr.png"))
          (-> (c/xy-plot r-plot dU-dr-re
                         :title "¹⁶O(α,α′) entrance — Re dU/dr and Im dU/dr (Table 7.1 full U; α, s=0)"
                         :x-label "r (fm)"
                         :y-label "dU/dr (MeV/fm)"
                         :series-label "Re dU/dr"
                         :legend true)
              (c/add-lines r-plot dU-dr-im :series-label "Im dU/dr")
              (i/save "output/o16_alpha_inelastic_2plus_dU_dr.png" :width 900 :height 520))
          (println "Plot saved: output/o16_alpha_inelastic_2plus_dU_dr.png (Re and Im dU/dr)")
          (io/make-parents (io/file "output/o16_alpha_inelastic_2plus_TL_vs_L.png"))
          (-> (c/xy-plot L-nums mag-max-by-La
                         :title "¹⁶O(α,α′) 2⁺ — max_{L_β} |T_{L_α L_β}| vs L_α (λ=2 coupling, Table 7.1)"
                         :x-label "L_α"
                         :y-label "max |T|"
                         :legend false)
              (i/save "output/o16_alpha_inelastic_2plus_TL_vs_L.png" :width 900 :height 520))
          (println "Plot saved: output/o16_alpha_inelastic_2plus_TL_vs_L.png")
          (-> (c/xy-plot L-nums ds-per-La
                         :title "¹⁶O(α,α′) 2⁺ — max_{L_β} dσ/dΩ from single pair (illustrative) vs L_α"
                         :x-label "L_α"
                         :y-label "dσ/dΩ (mb/sr)"
                         :legend false)
              (i/save "output/o16_alpha_inelastic_2plus_dsigma_vs_L.png" :width 900 :height 520))
          (println "Plot saved: output/o16_alpha_inelastic_2plus_dsigma_vs_L.png")
          (-> (c/xy-plot th-deg y-theta-austern
                         :title "¹⁶O(α,α′) 2⁺ — dσ/dΩ(θ): Austern (5.6) Σ_m|β^m|² (ℓ=2) vs reduced Y_{L_β0}"
                         :x-label "θ_CM (deg)"
                         :y-label "dσ/dΩ (mb/sr)"
                         :series-label "Austern (5.6) full m-sum"
                         :legend true)
              (c/add-lines th-deg y-theta :series-label "reduced Σ L_β Y_{L_β0}")
              (c/add-lines th-deg y-theta-dominant-L
                           :series-label (format "top-2 L_β {%s} coherent (peak |T| %d,%d)"
                                                 (clojure.string/join "," (map str (map first top2-Lb)))
                                                 (:L-alpha best) (:L-beta best)))
              (i/save "output/o16_alpha_inelastic_2plus_dsigma_theta.png" :width 900 :height 520))
          (println "Plot saved: output/o16_alpha_inelastic_2plus_dsigma_theta.png")
          (let [chart-log (-> (c/xy-plot th-deg (series-for-log y-theta-austern)
                                         :title "¹⁶O(α,α′) 2⁺ — dσ/dΩ(θ) Austern (5.6) + reduced (log y)"
                                         :x-label "θ_CM (deg)"
                                         :y-label "dσ/dΩ (mb/sr), logarithmic"
                                         :series-label "Austern (5.6)"
                                         :legend true)
                             (c/add-lines th-deg (series-for-log y-theta) :series-label "reduced Y_{L_β0}")
                             (c/add-lines th-deg (series-for-log y-theta-dominant-L)
                                          :series-label (format "top-2 L_β {%s} (|T| max %d,%d)"
                                                                (clojure.string/join "," (map str (map first top2-Lb)))
                                                                (:L-alpha best) (:L-beta best)))
                             (chart-set-log-range-y! "dσ/dΩ (mb/sr)"))]
            (i/save chart-log "output/o16_alpha_inelastic_2plus_dsigma_theta_log.png"
                    :width 900 :height 520))
          (println "Plot saved: output/o16_alpha_inelastic_2plus_dsigma_theta_log.png"))
        (catch Exception e
          (println "Note: could not save plots:" (.getMessage e))))
      (try
        (io/make-parents (io/file "output/o16_alpha_inelastic_2plus_TL.tsv"))
        (spit "output/o16_alpha_inelastic_2plus_TL.tsv"
              (str "L_alpha\tL_beta\treT\timT\tmagT\n"
                   (apply str (map (fn [r]
                                     (format "%d\t%d\t%.8e\t%.8e\t%.8e\n"
                                             (:L-alpha r) (:L-beta r)
                                             (re (:T r)) (im (:T r)) (:mag r)))
                                   rows))))
        (println "Wrote output/o16_alpha_inelastic_2plus_TL.tsv")
        (catch Exception e
          (println "Note: could not write TSV:" (.getMessage e)))))))

(println "\n=== Done ===")
