(ns dwba.benchmark.ca40-pd-austern
  "⁴⁰Ca(p,d)³⁹Ca g.s. pickup — **dσ/dΩ** from **N. Austern** radial **(5.5)**, angular **(5.6)**, amplitude **T_m ≈ D₀√(2ℓ+1) β_m**
  (**(4.60)** link), and **`transfer-differential-cross-section`**.

  Contrasts with **`dwba.benchmark.ca40-dwuck`**, which is **(d,p)** stripping and by default uses the
  reduced multipole angular map instead of the explicit **β(θ)** sum.

  **Conventions (same MeV/c² mass table as the DWUCK Ca40 block)**
  - Entrance **α**: **p + ⁴⁰Ca**; exit **β**: **d + ³⁹Ca**.
  - **Q** (nonrelativistic): **Q = m(⁴⁰Ca)+m(p) − m(³⁹Ca)−m(d)** ⇒ **e-cm-f = e-cm-i + Q**.
  - **ZR** sampling: **`austern-zr-chi-exit-mass-ratio` = M(⁴⁰Ca)/M(³⁹Ca)** in **(5.5)** / **`f-betaL`**.
  - Bound overlap: pickup **φ_i** = neutron-like **l=3** well in the target; **φ_f** = **l=0** deuteron
    (**`solve-bound-state-numerov`** parameters aligned with the **(d,p)** benchmark for the shell/deuteron).
  - Partial waves in **(5.5)/(5.6):** only **(L_α,L_β)** allowed by coupling **L⃗_β+ℓ⃗→L⃗_α** are built (**`austern-eq-5-6-admissible-L-beta-values`**), not a dense **L_α×L_β** square.

  **Angular / spin**
  - **Amplitude after (5.6):** we use **T_m ≈ D₀ √(2ℓ+1) β_m** from the **(4.60)** side note **(2ℓ+1)^{1/2} β = I/i^ℓ**
    (same **β_m** as **(5.6)**). **Full Austern (4.59)** needs a **sum over (m_a, m_b)**; a single **m_a = m_b = ½**
    forces **⟨ℓ 0; ½ 0 | j 0⟩ ≡ 0** for half-integer **j** (3-j parity), so **`T_m ≡ 0`** if you plug only that pair into **(4.59)**.
  - **Default angular weight:** **incoherent** **Σ_m |T_m|²** — often **≈ symmetric** about **90°** CM for Coulomb **σ_L(η)** in **(5.6)**.
    Nuclear OM sits in **χ** and in **(5.5)**. **`:coherent-m-beta? true`** uses **|Σ_m T_m|²**. **`:cm-asymmetry-kappa`** — **1 + κ cos θ** like **`ca40-dp-dsigma-mb-sr`**.
  - **Not elastic Rutherford:** **(5.6)** builds **transfer** **β_{sj}^{ℓm}** from radial **I_{L_β L_α}** and spherical harmonics — not the **elastic** amplitude **f_el(θ)**. See **`dwba.austern-eq56-and-rutherford-test`**.
  - Spin weight like **`ca40-dp-dsigma-mb-sr`**: **(2J_f+1)/(2J_i+1)**, unpolarized **d** **× 1/3**, entrance **p** **× ½**.

  **REPL:** `(require '[dwba.benchmark.ca40-pd-austern :as pd])` then `(pd/ca40-pd-dsigma-austern-eq-56-mb-sr 30.0)`.
  One entrance wave: `(pd/ca40-pd-dsigma-austern-eq-56-mb-sr 30.0 :L-alpha-only 5)` or **`ca40-pd-angular-curve-austern-eq-56-mb-sr`** with the same key (sums all **L_β** coupled to that **L_α**).

  Optics are **literature-style placeholders** (no DWUCK listing for this channel) — interpret shapes and
  relative **θ** dependence before absolute comparison to experiment."
  (:require [dwba.transfer :as t]
            [functions :as fn :refer [mass-factor-from-mu channel-sommerfeld-eta]]
            [complex :refer [mag mul add complex-cartesian]]))

(def ^:private ca40-pd-bound-ell 3)
(def ^:private ca40-pd-bound-j 3.5)
(def ^:private ca40-pd-J-i 0.0)
(def ^:private ca40-pd-J-f 3.5)

(defn ca40-pd-kinematics
  "Map **:mass-factor-i :mass-factor-f :e-cm-i :e-cm-f :k-i :k-f :Q-mev** for **p + ⁴⁰Ca → d + ³⁹Ca**.

  **e-cm-i** — CM kinetic energy (MeV) in the **entrance** channel (default **18.0**, open pickup).
  Uses **Q ≈ −5.85 MeV** with integer **A×931.494** masses + tabulated **m_p**, **m_d**."
  ([] (ca40-pd-kinematics 18.0))
  ([^double e-cm-i]
   (let [m-p 938.272
         m-d 1875.613
         m40 (* 40.0 931.494)
         m39 (* 39.0 931.494)
         Q (- (+ m40 m-p) (+ m39 m-d))
         e-cm-f (+ e-cm-i Q)
         _ (when (< e-cm-f 0.5)
             (throw (ex-info "ca40-pd-kinematics: e-cm-f too low — raise e-cm-i"
                             {:e-cm-i e-cm-i :Q Q :e-cm-f e-cm-f})))
         mu-i (/ (* m-p m40) (+ m-p m40))
         mu-f (/ (* m-d m39) (+ m-d m39))
         mfi (mass-factor-from-mu mu-i)
         mff (mass-factor-from-mu mu-f)]
     {:mass-factor-i mfi :mass-factor-f mff
      :e-cm-i e-cm-i :e-cm-f e-cm-f
      :k-i (Math/sqrt (* mfi e-cm-i))
      :k-f (Math/sqrt (* mff e-cm-f))
      :Q-mev Q
      :M-target m40 :M-residual m39})))

(defn optical-u-proton-ca40
  "Woods–Saxon **p + ⁴⁰Ca** (same functional form as **`optical-u-proton-ca41`** in `ca40-dwuck`, **RC** from **A=40**)."
  [L s j]
  (let [z1 1 z2 20
        rc (* 4.3103 (Math/pow (/ 40.0 41.0) (/ 1.0 3.0)))]
    (fn [^double r]
      (t/optical-potential-woods-saxon r [49.47 4.0689 0.70] [19.8 4.3172 0.75]
                                       24.2 4.0689 0.70 L s j z1 z2 rc))))

(defn optical-u-deuteron-ca39
  "Woods–Saxon **d + ³⁹Ca** — **`optical-u-deuteron-ca40`** from the listing, **R_C** scaled by **(39/40)^{1/3}**."
  [L s j]
  (let [z1 1 z2 20
        rc (* 4.7879 (Math/pow (/ 39.0 40.0) (/ 1.0 3.0)))]
    (fn [^double r]
      (t/optical-potential-woods-saxon r [97.4 3.803 0.875] [70.0 5.342 0.477]
                                       nil nil nil L s j z1 z2 rc))))

(defn- proton-j-for-partial-wave
  ^double [^long L]
  (+ 0.5 (double L)))

(defn- deuteron-j-for-partial-wave
  ^double [^long L]
  (if (= L 3)
    3.0
    (+ 1.0 (double L))))

(defn- sommerfeld-eta-channel
  ^double [^double e-cm ^double mfactor ^double z1z2ee]
  (binding [fn/mass-factor mfactor
            fn/Z1Z2ee z1z2ee]
    (channel-sommerfeld-eta e-cm)))

(defn- ca40-pd-cm-asymmetry-factor
  "Same exploratory factor as **`ca40-dp-cm-asymmetry-factor`**: **1 + κ cos θ** (θ rad)."
  ^double [^double theta-rad ^double kappa]
  (if (< (Math/abs kappa) 1e-15)
    1.0
    (max 1e-300 (+ 1.0 (* kappa (Math/cos theta-rad))))))

(defn- ca40-pd-rows-coulomb-sigma
  "Attach Coulomb **σ_L(η)** for **(5.6)** via **`austern-radial-rows-with-coulomb-sigma`**."
  [base-rows eta-i eta-f]
  (t/austern-radial-rows-with-coulomb-sigma base-rows eta-i eta-f))

(defn- ca40-pd-filter-rows-by-L-alpha
  "**(5.6)** sum restricted to **L_α = Lα** (all admissible **L_β** for that **L_α** retained). **nil** ⇒ no filter."
  [rows L-alpha-only]
  (if (nil? L-alpha-only)
    rows
    (let [La (long L-alpha-only)]
      (filterv #(= (long (:L-alpha %)) La) rows))))

(defn ca40-pd-radial-I-rows
  "Build **{:L-alpha :L-beta :I}** from **(5.5)** for **L_α ≤ L-max** and **L_β** coupled to **L_α** by transferred **ℓ**
  (**`austern-eq-5-6-admissible-L-beta-values`**), i.e. the same pairs retained in **(5.6)**.

  **`:transfer-ell`** — bound orbital **ℓ** in **(5.6)** (default **`ca40-pd-bound-ell`**, **3** here).

  Distorted **χ** use **`:coulomb-tail`** normalization (**`distorted-wave-optical`**) — **|u(r_max)|** matched to
  **|H_L^+(η, k r_max)|** per **L** so partial waves are not reweighted only by interior maxima (**mitigates
  artificial **θ → π** growth when **L_max** is large and (5.6) multiplies by **√(2L_β+1) Y_{L_β}^{m}(π)**)."
  [& {:keys [r-max h L-max e-cm-i transfer-ell]
      :or {r-max 100.0 h 0.05 L-max 20 e-cm-i 18.0 transfer-ell ca40-pd-bound-ell}}]
  (let [{:keys [mass-factor-i mass-factor-f e-cm-f k-i k-f M-target M-residual]}
        (ca40-pd-kinematics e-cm-i)
        zr (t/austern-zr-chi-exit-mass-ratio M-target M-residual)
        ell (long transfer-ell)
        phi-n (t/normalize-bound-state
               (t/solve-bound-state-numerov -8.364 3 58.4538 4.0355 0.7 0.048 h r-max) h)
        phi-d (t/normalize-bound-state
               (t/solve-bound-state-numerov -2.224 0 42.0 3.9 0.65 0.048 h r-max) h)
        z12 (* 1.44 1.0 20.0)
        eta-i (sommerfeld-eta-channel e-cm-i mass-factor-i z12)
        eta-f (sommerfeld-eta-channel e-cm-f mass-factor-f z12)
        rho-i (* k-i r-max)
        rho-f (* k-f r-max)
        chi-alpha!
        (memoize
         (fn [^long La]
           (t/distorted-wave-optical e-cm-i La 0.5 (proton-j-for-partial-wave La)
                                     (optical-u-proton-ca40 La 0.5 (proton-j-for-partial-wave La))
                                     r-max h mass-factor-i
                                     :normalize-mode :coulomb-tail
                                     :tail-eta eta-i :tail-rho rho-i)))
        chi-beta!
        (memoize
         (fn [^long Lb]
           (t/distorted-wave-optical e-cm-f Lb 1.0 (deuteron-j-for-partial-wave Lb)
                                     (optical-u-deuteron-ca39 Lb 1.0 (deuteron-j-for-partial-wave Lb))
                                     r-max h mass-factor-f
                                     :normalize-mode :coulomb-tail
                                     :tail-eta eta-f :tail-rho rho-f)))]
    (vec
     (for [La (range 0 (inc (long L-max)))
           Lb (t/austern-eq-5-6-admissible-L-beta-values La ell (long L-max))
           :let [Ireal (t/austern-radial-integral-I-eq-5-5-from-F-lsj
                        phi-n phi-d (chi-alpha! La) (chi-beta! Lb) h
                        M-target M-residual k-i k-f zr)]
           :when (> (Math/abs (double Ireal)) 1e-30)]
       {:L-alpha La :L-beta Lb :I Ireal}))))

(defn ca40-pd-T-squared-sum-from-austern-eq-56
  "With **`:coherent-m-beta? false`** (default): **Σ_m |T_m|²**. With **`:coherent-m-beta? true`**: **|Σ_m T_m|²**
  (**m-channel interference** — use sparingly).

  **T_m = D₀ √(2ℓ+1) β_m**; **β_m** from **`austern-reduced-amplitude-beta-sum-eq-5-6`**."
  [theta-rad radial-rows-sigma D0 & {:keys [coherent-m-beta?] :or {coherent-m-beta? false}}]
  (let [ell (long ca40-pd-bound-ell)
        sqrt2l1 (Math/sqrt (inc (* 2.0 (double ell))))
        pref (complex-cartesian (* (double D0) sqrt2l1) 0.0)
        ms (range (- ell) (inc ell))]
    (if coherent-m-beta?
      (let [sum-b (reduce (fn [acc ^long m-ell]
                            (add acc (mul pref (t/austern-reduced-amplitude-beta-sum-eq-5-6
                                                  ell m-ell theta-rad radial-rows-sigma))))
                          (complex-cartesian 0.0 0.0)
                          ms)
            Tmag (mag sum-b)]
        (* Tmag Tmag))
      (double
       (reduce
        (fn [^double acc ^long m-ell]
          (let [beta (t/austern-reduced-amplitude-beta-sum-eq-5-6 ell m-ell theta-rad radial-rows-sigma)
                Tm (mul pref beta)
                Tmag (mag Tm)]
            (+ acc (* Tmag Tmag))))
        0.0
        ms)))))

(defn ca40-pd-dsigma-austern-eq-56-mb-sr
  "**dσ/dΩ (mb/sr)** at CM angle **theta-deg** via **(5.5)**, **(5.6)**, **T_m ≈ D₀√(2ℓ+1) β_m** (**(4.60)** link), and **`transfer-differential-cross-section`** (**S=1**).
  Expensive: builds partial-wave **I** table each call unless you pass **`:radial-rows-sigma`**.

  **Options:** **`:e-cm-i`**, **`:r-max`**, **`:h`**, **`:L-max`**, **`:S-factor`** (default **1**),
  **`:radial-rows-sigma`**, **`:coherent-m-beta?`**, **`:cm-asymmetry-kappa`**,
  **`:L-alpha-only`** — keep only rows with this entrance **L_α** (slice **(5.6)** / **(5.5)** to one partial wave)."
  [theta-deg & {:keys [e-cm-i r-max h L-max S-factor radial-rows-sigma
                       coherent-m-beta? cm-asymmetry-kappa L-alpha-only]
                :or {e-cm-i 18.0 r-max 100.0 h 0.05 L-max 20 S-factor 1.0
                     coherent-m-beta? false
                     cm-asymmetry-kappa 0.0}}]
  (let [{:keys [mass-factor-i mass-factor-f k-i k-f e-cm-i e-cm-f]} (ca40-pd-kinematics e-cm-i)
        z12 (* 1.44 1.0 20.0)
        eta-i (sommerfeld-eta-channel e-cm-i mass-factor-i z12)
        eta-f (sommerfeld-eta-channel e-cm-f mass-factor-f z12)
        rows-sig (ca40-pd-filter-rows-by-L-alpha
                  (or radial-rows-sigma
                      (ca40-pd-rows-coulomb-sigma
                        (ca40-pd-radial-I-rows :r-max r-max :h h :L-max L-max :e-cm-i e-cm-i)
                        eta-i eta-f))
                  L-alpha-only)
        D0 (t/zero-range-constant :p-d)
        theta-rad (* (double theta-deg) (/ Math/PI 180.0))
        T-sq (ca40-pd-T-squared-sum-from-austern-eq-56 theta-rad rows-sig D0
               :coherent-m-beta? (boolean coherent-m-beta?))
        ;; Same weight recipe as `ca40-dp-dsigma-mb-sr`: stat × (1/3) d; × (1/2) for p m_a (pickup).
        spin (* (t/transfer-nuclear-spin-statistical-factor ca40-pd-J-i ca40-pd-J-f)
                (t/transfer-unpolarized-deuteron-spin-factor)
                0.5)
        ds (t/transfer-differential-cross-section (Math/sqrt (max T-sq 0.0)) S-factor k-i k-f
                                                   mass-factor-i mass-factor-f)
        asym (ca40-pd-cm-asymmetry-factor theta-rad (double cm-asymmetry-kappa))]
    (* (double ds) spin asym)))

(defn ca40-pd-angular-curve-austern-eq-56-mb-sr
  "Seq **`{:theta-deg θ :differential_cross_section_mb_sr σ}`**; optional **`:e-cm-i`**, **`:r-max`**, **`:h`**, **`:L-max`**, **`:S-factor`**,
  **`:coherent-m-beta?`**, **`:cm-asymmetry-kappa`**, **`:L-alpha-only`**. Precomputes Coulomb **σ** rows once per curve."
  [theta-degrees & {:keys [e-cm-i r-max h L-max S-factor
                           coherent-m-beta? cm-asymmetry-kappa L-alpha-only]
                    :or {e-cm-i 18.0 r-max 100.0 h 0.05 L-max 20 S-factor 1.0
                         coherent-m-beta? false
                         cm-asymmetry-kappa 0.0}}]
  (let [{:keys [mass-factor-i mass-factor-f e-cm-i e-cm-f]} (ca40-pd-kinematics e-cm-i)
        z12 (* 1.44 1.0 20.0)
        eta-i (sommerfeld-eta-channel e-cm-i mass-factor-i z12)
        eta-f (sommerfeld-eta-channel e-cm-f mass-factor-f z12)
        base-rows (ca40-pd-radial-I-rows :r-max r-max :h h :L-max L-max :e-cm-i e-cm-i)
        rows-sig (ca40-pd-filter-rows-by-L-alpha
                  (ca40-pd-rows-coulomb-sigma base-rows eta-i eta-f)
                  L-alpha-only)]
    (mapv (fn [^double th]
            {:theta-deg th
             :differential_cross_section_mb_sr
             (ca40-pd-dsigma-austern-eq-56-mb-sr th
               :radial-rows-sigma rows-sig
               :e-cm-i e-cm-i :r-max r-max :h h :L-max L-max :S-factor S-factor
               :coherent-m-beta? coherent-m-beta?
               :cm-asymmetry-kappa cm-asymmetry-kappa)})
          theta-degrees)))
