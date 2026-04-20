;; Elastic check: α + ¹⁴⁸Sm at E_lab(α) = 50 MeV
;;
;; **Point Coulomb** in the radial equation: **V = [0  R_C  a_WS]** with **V₀ = 0** (no nuclear
;; Woods–Saxon) and **R_C = 0** so **`Coulomb-pot`** uses **Z₁Z₂e² / r** for **r > 0**. The third
;; entry **a_WS** is only used for the **matching radius** **2(R_C + a_WS)** in **`s-matrix`** when
;; **V₀ = 0**; it should stay **> 0**. Uniform/volume charge is **R_C > 0** — enable that later.
;;
;; **R-matrix match** uses **`s-matrix`** default **a = 2(R_C + a₀)** (MeV·fm potential **V = [V₀ R_C a]**).
;; **Do not** bind **`*elastic-match-radius-fm*`** / **`*r-matrix-numerov-dr-fm*`** here for Fig. 2.2-style
;; handbook plots — large **a** (e.g. 300 fm) or coarse **dr** changes **S_L^n** and ruins agreement with the book.
;; To experiment, uncomment inside **`binding`** below: **`*elastic-match-radius-fm*`**, **`*r-matrix-numerov-dr-fm*`**
;; (and **`require`** those vars from **`functions`**).
;;
;; **`s-matrix`** = **S_L^n** (same R-matrix / Hankel quotient as neutral **`s-matrix0`**, **Hankel±** for **η ≠ 0**; **no σ** in **`s-matrix`**). **σ** only in **e^{2iσ}(S^n−1)** for **f_N**. Elastic **dσ/dΩ** via
;; **`differential-cross-section-nuclear-cut`**: **dσ = 10 |f_C + f_N|² mb/sr** with
;; **f_C** = Thompson & Nunes **Eq. (3.181)** (`coulomb-scattering-amplitude-thompson-nunes-eq-3181`),
;; **f_N** = partial waves **L ≤ L_nuclear_cut** (vanishes when **`s-matrix` = 1**, pure point Coulomb in **(3.1.84)**).
;; Same convention as the web dashboard elastic API (*elastic-imag-ws-params*).
;;
;; A printed table compares that reference to **|f̃_C + f̃_N|²** (**`coulomb-amplitude-tilde`** + **`coulomb-phase-diff`** partial sum), same as **`examples/example_16Odp.clj`** — **dσ** is unchanged (**σ_0** is a global phase).
;;
;; Optional full optical + finite charge radius (swap **V-params** / **imag-ws-params**):
;;   Real:    V₀ = 65 MeV,  R = 7.5 fm, a = 0.67 fm
;;   Imag:   Wᵢ = 30 MeV,  Rᵢ = 7.5 fm, aᵢ = 0.67 fm
;;
;; Run from repository root:
;;   lein run -m clojure.main examples/example_alpha_Sm148_elastic.clj
;;
;; Writes (gitignored):
;;   **output/alpha_Sm148_elastic_sigma_ratio_Ruth.png** — σ/σ_Ruth vs θ_CM (linear y), **L** cut = **L_nuclear_cut**
;;   **output/alpha_Sm148_elastic_sigma_ratio_Ruth_log.png** — same, log₁₀ y-axis
;;   **output/alpha_Sm148_elastic_fig22b_Lconvergence.png** — Fig. 2.2(b)-style: **L ≤ 20, 35, 45** (convergence)
;;   **output/alpha_Sm148_elastic_fig22b_Lconvergence_log.png** — same, log y
;;   **output/alpha_Sm148_elastic_SL_re_im.png** — **Re(S_L^n)**, **Im(S_L^n)** vs **L** — **nuclear** **S_L^n** from **`s-matrix`**
;;     (R-matrix / Coulomb **Hankel±** match; **Coulomb σ_L** is **not** folded into **S**, T&N **§3.1**).
;;   **output/alpha_Sm148_elastic_SL_exp2sigma_re_im.png** — **Re/Im** of **e^{2iσ_L} S_L^n** (same **S^n**, multiplied by Coulomb phase)
;;     for comparison to tables that absorb **σ** into the partial-wave **S**.
;;
;; Handbook-style plots (e.g. Fig. 2.2(b), p.11) compare **f_N** sums truncated at different **L_max**; the curves
;; merge where the partial-wave series has converged. Sample **θ** densely — the physics is not piecewise-linear.

(import '[org.jfree.chart.axis LogAxis])

;; :reload so REPL picks up src/functions.clj changes (avoids stale dσ logic + s-matrix-3-memo confusion).
(require '[functions :refer [differential-cross-section-nuclear-cut mass-factor mass-factor-from-mu
                             s-matrix phase-shift Z1Z2ee *elastic-imag-ws-params*
                             coulomb-sigma-L channel-sommerfeld-eta
                             coulomb-amplitude-tilde coulomb-phase-diff
                             coulomb-scattering-amplitude-thompson-nunes-eq-3181]]
         :reload)
(require '[complex :as cpx]
         '[fastmath.polynomials :as poly]
         '[fastmath.core :as m]
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

(def ^:private m-alpha 3727.379)   ; MeV/c² (same as web-dashboard elastic handler)
(def ^:private A-sm 148)
(def ^:private Z-sm 62)
(def ^:private m-sm (* 931.5 A-sm)) ; MeV/c² (same mass convention as dashboard)

(def E-lab 50.0)  ; Lab kinetic energy of alpha (MeV)

;; CM kinetic energy: target at rest, non-relativistic
(def E-CM (* E-lab (/ m-sm (+ m-alpha m-sm))))

(def mu-reduced (/ (* m-alpha m-sm) (+ m-alpha m-sm)))

;; Z₁ Z₂ e² with e² = 1.44 MeV·fm (nuclear convention, same as simple_core)
(def z1z2ee (* 2 Z-sm 1.44))

;; Real channel: [V₀ (MeV), R_C or R_WS (fm), a (fm)]. Point Coulomb: V₀=0, R_C=0, a>0 for matching.
;(def V-params [0.0 0.0 0.67])
(def V-params [65. 7.5 0.67])
;; Imaginary WS: W₀ (MeV), R_W (fm), a_W (fm) — absorption (bound to *elastic-imag-ws-params*)
;(def imag-ws-params [0.0 0. 0.0])
(def imag-ws-params [30.0 7.5 0.67])

;; **`f_N`** uses **(S_L^n − 1)**; only **S_L^n → +1** (complex **1 + 0i**) kills that wave.  **|S_L^n| ≈ 1**
;; (unitarity) is **not** enough — e.g. **S_L^n ≈ −1** still has **|S_L^n − 1| ≈ 2**.  For this 65 MeV real WS,
;; **`s-matrix`** at high **L** stays on the unit circle but **not** at **+1**, so **L_cut = 41** can still change
;; **dσ**.  Point Coulomb (**V₀ = 0**, **R_C = 0**) gives **S_L^n = 1** and then **L_cut** beyond that is harmless.
;; If plots looked **jagged** with imag **WS**, check **non-finite `s-matrix`** (fixed path: Numerov in
;; **`r-matrix-complex-imag-ws`**).
(def L-max 45)
(def L-nuclear-cut 41)

;; Partial-wave cutoffs for **Fig. 2.2(b)**-style convergence plot (**f_N** sum **L = 0 … L_cut**).
(def ^:private fig22b-L-cuts [20 35 45])

(defn rutherford-dsigma-mb-sr
  "Point-charge Rutherford dσ/dΩ (**mb/sr**), non-relativistic CM. Z1Z2e² in MeV·fm, E_cm in MeV, θ rad.
  Leading fraction is **fm²/sr**; **×10** is **1 fm² = 10 mb** (same as **`functions/differential-cross-section`**)."
  [^double z1z2e2 ^double e-cm ^double theta-rad]
  (* 10.0 ;; fm²/sr → mb/sr
     (let [s (Math/sin (* 0.5 theta-rad))]
       (if (< (Math/abs s) 1e-15)
         Double/NaN
         (/ (Math/pow (/ z1z2e2 (* 4.0 e-cm)) 2.0)
            (Math/pow s 4.0))))))

(defn- elastic-dsigma-tilde-mb-sr
  "**10 |f̃_C + f̃_N|²** (**mb/sr**): **f̃_C** = **`coulomb-amplitude-tilde`**, **f̃_N** = partial waves **L ≤ L-cut** with
  **`coulomb-phase-diff`** (same construction as **`example_16Odp`**). Requires **`mass-factor`**, **`Z1Z2ee`**, **`*elastic-imag-ws-params*`** bound like **`differential-cross-section-nuclear-cut`**."
  [^double e-cm V-params ^double theta-rad ^long L-cut]
  (let [eta (channel-sommerfeld-eta e-cm)
        k (Math/sqrt (* (double mass-factor) e-cm))
        f-tilde-c (coulomb-amplitude-tilde theta-rad eta k)
        f-tilde-n (reduce (fn [acc ^long L]
                            (try
                              (let [Sn (s-matrix e-cm V-params L)
                                    reS (double (cpx/re Sn))
                                    imS (double (cpx/im Sn))]
                                (if (or (Double/isNaN reS) (Double/isNaN imS)
                                        (Double/isInfinite reS) (Double/isInfinite imS))
                                  acc
                                  (let [ph-prod (coulomb-phase-diff L eta)
                                        pl (double (poly/eval-legendre-P L (m/cos theta-rad)))
                                        bracket (cpx/mul ph-prod (cpx/subt2 Sn 1.0))
                                        contrib (cpx/mul (cpx/complex-polar (* -0.5 Math/PI) (/ (inc (* 2 L)) (* 2.0 k)))
                                                         pl bracket)]
                                    (cpx/add2 acc contrib))))
                              (catch Exception _ acc)))
                          (cpx/complex-cartesian 0.0 0.0)
                          (range (inc L-cut)))]
    (* 10.0 (Math/pow (cpx/mag (cpx/add2 f-tilde-c f-tilde-n)) 2))))

(defn- sigma-ratio-vectors-for-L-cut
  "Vectors **[θ_deg]**, **[σ/σ_Ruth]** for **differential-cross-section-nuclear-cut** with nuclear sum **L ≤ L-cut**."
  [e-cm V-params z12ee theta-step-deg L-cut]
  (let [e-cm (double e-cm)
        z12ee (double z12ee)
        theta-step-deg (double theta-step-deg)
        Lc (long L-cut)
        pairs (for [theta-deg (range 10.0 (+ 171.0 (/ theta-step-deg 2.0)) theta-step-deg)]
                (try
                  (let [^double td theta-deg
                        theta (* td (/ Math/PI 180.0))
                        ds-c (differential-cross-section-nuclear-cut e-cm V-params theta Lc)
                        r (double (cpx/re ds-c))
                        im (double (cpx/im ds-c))
                        ds-mb-sr (Math/sqrt (+ (* r r) (* im im)))
                        ruth (rutherford-dsigma-mb-sr z12ee e-cm theta)
                        ratio (when (and (Double/isFinite ruth) (> ruth 0.0))
                                (/ ds-mb-sr ruth))]
                    (when (and ratio (Double/isFinite ratio))
                      [td (double ratio)]))
                  (catch Exception _ nil)))
        th-vec (vec (keep (fn [p] (when p (first p))) pairs))
        y-vec (vec (keep (fn [p] (when p (second p))) pairs))]
    [th-vec y-vec]))

(println "=== α + ¹⁴⁸Sm elastic — point Coulomb in V(r), nuclear WS off (sanity check) ===")
(println "")
(println (format "  Projectile: α  (m = %.3f MeV/c², Z = 2)" m-alpha))
(println (format "  Target:     ¹⁴⁸Sm (m ≈ %.1f MeV/c², Z = %d)" m-sm Z-sm))
(println (format "  E_lab(α)    = %.2f MeV" E-lab))
(println (format "  E_CM        = %.4f MeV" E-CM))
(println (format "  μ (reduced) = %.4f MeV/c²" mu-reduced))
(println (format "  Z₁Z₂ e²     = %.4f MeV·fm" z1z2ee))
(println (format "  Real V:     [V₀ R a]     = [%.1f %.2f %.2f] (MeV, fm — R=charge/nuclear radius, fm)"
                 (first V-params) (second V-params) (nth V-params 2)))
(println (format "  Imag WS:    [Wᵢ Rᵢ aᵢ] = [%.1f %.2f %.2f] MeV, fm, fm"
                 (first imag-ws-params) (second imag-ws-params) (nth imag-ws-params 2)))
(println "")
(println "=== Point-charge Rutherford (Coulomb only, non-relativistic CM) — verify the numbers ===")
(println "")
(let [pref (/ z1z2ee (* 4.0 E-CM))
      pref2 (* pref pref)]
  (println "  dσ/dΩ|_Ruth = (Z₁Z₂e² / (4 E_cm))² / sin⁴(θ/2)  →  **mb/sr** here (×10 from fm²/sr).")
  (println "")
  (println (format "  Z₁Z₂e² / (4 E_cm) = %.4f / (4 × %.4f) = %.6f fm" z1z2ee E-CM pref))
  (println (format "  (Z₁Z₂e² / (4 E_cm))² = %.6f (fm²); ×10 ⇒ mb/sr prefactor scale" pref2))
  (println "")
  (println "  Forward peaking comes entirely from 1/sin⁴(θ/2): sin(5°) ≈ 0.0872, sin⁴(5°) ≈ 5.77×10⁻⁵,")
  (println "  so σ_Ruth(10°) is ~1/sin⁴ larger than at θ where sin⁴ ~ 𝒪(1) (e.g. near 90°).")
  (println "")
  (doseq [deg [10 90]]
    (let [th (* deg (/ Math/PI 180.0))
          s (Math/sin (* 0.5 th))
          s4 (Math/pow s 4)
          sig (/ pref2 s4)]
      (println (format "  θ_cm = %3d°: sin(θ/2) = %.6f, sin⁴(θ/2) = %.5e  →  σ_Ruth = %12.1f mb/sr"
                       deg s s4 (* 10.0 sig))))))
(println "")
(println "  Elastic **dσ** here: **`differential-cross-section-nuclear-cut`** ⇒ **10 |f_C + f_N|²** with")
(println "    **f_C** — Thompson & Nunes **(3.181)**; **f_N** — **(3.1.88)**, **L ≤ L_nuclear_cut**, each term")
(println "    **−i/k · (2L+1) P_L · e^{2iσ_L} (S_L^n − 1)** with **S_L^n = `s-matrix`** (R-matrix / **Hankel±**). **(3.1.84):** incoming/outgoing Coulomb match; **S_L^n = 1** ⇒ no **f_N**.")
(println "    **`phase-shift`** = **½ arg(S_L^n)** — **not** **σ + …**; **σ** enters only the **e^{2iσ}** in **f_N**.")
(println "  Tiny **S_L^n − 1** from numerics can shift **dσ** when **L_cut** changes.")
(println "  Ratio **dσ/σ_Ruth** uses point Rutherford; at small θ **f_C** dominates.")
(println "")

(binding [mass-factor (mass-factor-from-mu mu-reduced)
          Z1Z2ee z1z2ee
          *elastic-imag-ws-params* imag-ws-params
          ;; *elastic-match-radius-fm* 300.0
          ;; *r-matrix-numerov-dr-fm* 0.1
          ]
  (println "=== **`s-matrix`** = **S_L^n** (T&N **§3.1** / **(3.1.84)**) ===")
  (println (format "   (**f_N** uses **L ≤ %d**; **phase-shift** = ½ arg(**S_L^n**).)" L-nuclear-cut))
  (println "")
  (println "   L  |   Re(S_L^n)  |   Im(S_L^n)  | |S^n| | arg(S^n) rad |  δ^n (deg)  (= phase-shift)")
  (println " -----+------------+------------+-------+-------------+---------------------------")
  (let [sl-rows (vec (for [L (range (inc L-max))]
                       (try
                         (let [S-val (s-matrix E-CM V-params L)
                               reS (double (cpx/re S-val))
                               imS (double (cpx/im S-val))]
                           {:L L :S-val S-val :re reS :im imS})
                         (catch Exception e
                           {:L L :error e}))))]
    (doseq [row sl-rows]
      (if-let [e (:error row)]
        (println (format " %3d | ERROR: %s" (:L row) (.getMessage ^Exception e)))
        (let [{:keys [L S-val re im]} row]
          (if (or (Double/isNaN re) (Double/isNaN im)
                  (Double/isInfinite re) (Double/isInfinite im))
            (println (format " %3d | (non-finite S^n — Coulomb/Hankel U numerics at this L, η, ρ)" L))
            (let [Sm (cpx/mag S-val)
                  d-tot-rad (phase-shift E-CM V-params L)
                  d-deg (* d-tot-rad (/ 180.0 Math/PI))]
              (println (format " %3d | %10.5f | %10.5f | %.5f | %11.5f | %11.5f"
                               L re im Sm (cpx/arg S-val) d-deg)))))))
    (let [finite-rows (filterv (fn [row]
                                  (and (not (:error row))
                                       (let [re (:re row) im (:im row)]
                                         (and (Double/isFinite re) (Double/isFinite im)))))
                                sl-rows)
          L-plot (mapv :L finite-rows)
          re-plot (mapv :re finite-rows)
          im-plot (mapv :im finite-rows)
          eta-cm (channel-sommerfeld-eta E-CM)
          dressed (mapv (fn [{:keys [L S-val]}]
                          (let [sig (coulomb-sigma-L L eta-cm)
                                e2is (cpx/complex-polar (* 2.0 sig) 1.0)
                                ds (cpx/mul e2is S-val)]
                            {:L L :re (double (cpx/re ds)) :im (double (cpx/im ds))}))
                        finite-rows)
          re-d (mapv :re dressed)
          im-d (mapv :im dressed)]
      (println "  **S_L^n** plots: values are **`s-matrix`** = nuclear **S** (Coulomb-matched); **σ_L** enters **f_N** as **e^{2iσ_L}(S^n−1)**, not inside **S^n**.")
      (when (seq L-plot)
        (try
          (io/make-parents (io/file "output/alpha_Sm148_elastic_SL_re_im.png"))
          (-> (c/xy-plot L-plot re-plot
                         :title "α + ¹⁴⁸Sm elastic — Re/Im S_L^n (nuclear; σ not in S, T&N)"
                         :x-label "L"
                         :y-label "S_L^n"
                         :series-label "Re S_L^n"
                         :legend true)
              (c/add-lines L-plot im-plot :series-label "Im S_L^n")
              (i/save "output/alpha_Sm148_elastic_SL_re_im.png" :width 900 :height 520))
          (-> (c/xy-plot L-plot re-d
                         :title "α + ¹⁴⁸Sm elastic — Re/Im e^{2iσ_L} S_L^n (Coulomb-phased)"
                         :x-label "L"
                         :y-label "e^{2iσ} S^n"
                         :series-label "Re e^{2iσ} S^n"
                         :legend true)
              (c/add-lines L-plot im-d :series-label "Im e^{2iσ} S^n")
              (i/save "output/alpha_Sm148_elastic_SL_exp2sigma_re_im.png" :width 900 :height 520))
          (println (format "S-matrix plots: …SL_re_im.png and …SL_exp2sigma_re_im.png (%d finite L, L_max = %d)."
                           (count L-plot) L-max))
          (catch Exception e
            (println "Note: could not save S_L Re/Im plot:" (.getMessage e)))))))
  (println "")
  (println "=== Elastic dσ/dΩ (|f_C + f_N|²); **f_N** partial waves **L = 0 … " L-nuclear-cut " ===")
  (println "(Partial waves with non-finite **S_L^n** / match are omitted — see differential-cross-section doc.)")
  (println "σ_Ruth: same point-Coulomb formula as the web dashboard (e² = 1.44 MeV·fm).")
  (println "")
  (println "  θ_cm (deg) | dσ/dΩ / σ_Ruth | dσ/dΩ (mb/sr)")
  (println " ------------+------------------+---------------")
  (doseq [theta-deg [10 20 30 45 60 90 120 150 170]]
    (try
      (let [theta (* theta-deg (/ Math/PI 180.0))
            ds-c (differential-cross-section-nuclear-cut E-CM V-params theta L-nuclear-cut)
            ;; |f|² returned as near-real complex; Cartesian magnitude avoids refer/shadow bugs.
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
          (println (format "   %6.1f    |     (n/a θ→0°) | %12.5e"
                           (double theta-deg) ds-mb-sr))))
      (catch Exception e
        (println (format "   %6.1f    | ERROR: %s" theta-deg (.getMessage e))))))
  (println "")
  (println "=== **|f̃_C + f̃_N|²** vs **`differential-cross-section-nuclear-cut`** (same **dσ**; **σ_0** phase) ===")
  (println (format "  Partial waves **L = 0 … %d**; **|f_C|²** column = **10 |f_C(3.181)|²** (not Rutherford at finite **η**)." L-nuclear-cut))
  (println (format "  %-6s  %14s  %14s  %14s" "θ_CM" "|f_C|²" "|f̃_C+f̃_N|²" "|f_C+f_N|² ref"))
  (doseq [theta-deg [10 20 30 45 60 90 120 150 170]]
    (try
      (let [td (double theta-deg)
            theta (* td (/ Math/PI 180.0))
            eta (channel-sommerfeld-eta E-CM)
            k (Math/sqrt (* (double mass-factor) E-CM))
            fc-sq (* 10.0 (Math/pow (cpx/mag (coulomb-scattering-amplitude-thompson-nunes-eq-3181 theta eta k)) 2))
            tilde-sq (elastic-dsigma-tilde-mb-sr E-CM V-params theta L-nuclear-cut)
            ref-c (differential-cross-section-nuclear-cut E-CM V-params theta L-nuclear-cut)
            ref-sq (double (cpx/re ref-c))]
        (println (format "  %5.1f°  %14.4e  %14.4e  %14.4e"
                         td fc-sq tilde-sq ref-sq)))
      (catch Exception e
        (println (format "  %5.1f°  ERROR: %s" (double theta-deg) (.getMessage e))))))
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
        y-vec (vec (keep (fn [p] (when p (second p))) pairs))]
    (when (and (seq th-vec) (= (count th-vec) (count y-vec)))
      (try
        (io/make-parents (io/file "output/alpha_Sm148_elastic_sigma_ratio_Ruth.png"))
        (let [y-log (mapv clamp-pos-for-log y-vec)]
          (-> (c/xy-plot th-vec y-vec
                         :title "α + ¹⁴⁸Sm elastic (point V_C) — |f_C+f_N|² / σ_Ruth (CM)"
                         :x-label "θ_CM (deg)"
                         :y-label "σ / σ_Ruth"
                         :series-label "σ/σ_Ruth"
                         :legend true)
              (i/save "output/alpha_Sm148_elastic_sigma_ratio_Ruth.png" :width 800 :height 500))
          (-> (c/xy-plot th-vec y-log
                         :title "α + ¹⁴⁸Sm elastic (point V_C) — |f_C+f_N|² / σ_Ruth (CM), log y"
                         :x-label "θ_CM (deg)"
                         :y-label "σ / σ_Ruth (log scale)"
                         :series-label "σ/σ_Ruth"
                         :legend true)
              (chart-set-log-range-y! "σ / σ_Ruth (log scale)")
              (i/save "output/alpha_Sm148_elastic_sigma_ratio_Ruth_log.png" :width 800 :height 500)))
        (println (format "Plots saved: output/alpha_Sm148_elastic_sigma_ratio_Ruth.png, …_log.png (%d points, Δθ = %.1f°)."
                         (count th-vec) theta-step))
        (catch Exception e
          (println "Note: could not save ratio plot:" (.getMessage e)))))
    ;; Fig. 2.2(b): same ratio, three **L** truncations (convergence of **f_N**).
    (try
      (let [theta-step-22b 0.25
            series (for [^long Lc fig22b-L-cuts]
                     (let [[thv yv] (sigma-ratio-vectors-for-L-cut E-CM V-params z1z2ee theta-step-22b Lc)]
                       {:L-cut Lc :th thv :y yv :label (format "L ≤ %d" Lc)}))
            first-s (first series)
            rest-s (rest series)]
        (when (and first-s (seq (:th first-s)))
          (io/make-parents (io/file "output/alpha_Sm148_elastic_fig22b_Lconvergence.png"))
          (let [chart-base (c/xy-plot (:th first-s) (:y first-s)
                                      :title "α + ¹⁴⁸Sm elastic — σ/σ_Ruth vs L_cut (Fig. 2.2(b) style)"
                                      :x-label "θ_CM (deg)"
                                      :y-label "σ / σ_Ruth"
                                      :series-label (:label first-s)
                                      :legend true)
                chart (reduce (fn [ch s]
                                (c/add-lines ch (:th s) (:y s) :series-label (:label s)))
                              chart-base
                              rest-s)]
            (i/save chart "output/alpha_Sm148_elastic_fig22b_Lconvergence.png" :width 900 :height 520))
          (let [chart-log-base (c/xy-plot (:th first-s) (mapv clamp-pos-for-log (:y first-s))
                                           :title "α + ¹⁴⁸Sm elastic — σ/σ_Ruth vs L_cut (log y)"
                                           :x-label "θ_CM (deg)"
                                           :y-label "σ / σ_Ruth (log scale)"
                                           :series-label (:label first-s)
                                           :legend true)
                chart-log (reduce (fn [ch s]
                                    (c/add-lines ch (:th s) (mapv clamp-pos-for-log (:y s))
                                                 :series-label (:label s)))
                                  chart-log-base
                                  rest-s)]
            (-> chart-log
                (chart-set-log-range-y! "σ / σ_Ruth (log scale)")
                (i/save "output/alpha_Sm148_elastic_fig22b_Lconvergence_log.png" :width 900 :height 520)))
          (println (format "Fig. 2.2(b) plots: output/alpha_Sm148_elastic_fig22b_Lconvergence.png, …_log.png (L_cut ∈ %s, Δθ = %.2f°)."
                           (pr-str (vec fig22b-L-cuts)) theta-step-22b))))
      (catch Exception e
        (println "Note: could not save Fig. 2.2(b) convergence plot:" (.getMessage e)))))
  (println "")
  (println "Note: **f_C** = T&N (3.181); **f_N** sums L ≤ L_nuclear_cut (non-finite **S_L^n** / match omitted); σ_Ruth = point Rutherford.")
  (println "Done."))
