(ns dwba.benchmark.o16-dp-handbook
  "¹⁶O(d,p)¹⁷O g.s. stripping — **handbook** (*Handbook of direct nuclear reaction for retarded theorist*)
  ZR conventions, parallel to **`dwba.benchmark.ca40-pd-handbook`** with **α = d+¹⁶O**, **β = p+¹⁷O**.

  **Radial:** **F_{ℓsj} = R_n = u/r** (bound neutron in **¹⁷O**); **`handbook-radial-integral-I-zr-from-neutron-bound`**
  uses **(5.5)** **(M_B/M_A)(4π/(k_α k_β))** on **∫ F R_α R_β r² dr** (same **I** as Austern ZR for this **β** / **dσ** chain). **M_A = M(¹⁶O)**, **M_B = M(¹⁷O)**.

  **Angular:** **`handbook-zr-multipole-amplitude-sum`** + **`handbook-zr-rows-with-coulomb-sigma`**.

  **Amplitude:** **T_m ≈ D₀ √(2ℓ+1) β_m** and **`transfer-differential-cross-section`** → **mb/sr** (same dimensional unit as typical handbook **mb/sr** plots).

  **Distorted waves:** default **`:chi-normalize-mode` `:coulomb-tail`** (**|**u**| ∝ |H_L^+|** at **ρ = k r_max** per channel). Optional **`:raw`** or **`:max`**. With strong absorption, tail match can mis-scale **χ**; use **`:raw`** for relative **L** weights. Elastic **R → S** uses **`distorted-wave-numerov-R-for-smatrix`** + **`distorted-wave-coulomb-S-from-numerov-R`**. Match **potentials, energies, and conventions** from any reference figure. For **Ca40(d,p)** vs **DWUCK4** listing scale, see **`ca40-dp-flux-scale-to-embedded-dwuck`**.

  **Spin:** **(2J_f+1)/(2J_i+1)** for **J(¹⁶O)=0**, **J(¹⁷O)=5/2** and unpolarized deuteron **× 1/3** (same pattern as **`ca40-dp-dsigma-mb-sr`**, no extra **½** — that factor is for **(p,d)** entrance proton).

  Optics are **Ca40 listing depths/radii scaled** to **A≈16–17**, **Z=8** — illustrative; tune before comparing to experiment or a specific handbook figure."
  (:require [dwba.transfer :as t]
            [functions :as fn :refer [mass-factor-from-mu channel-sommerfeld-eta lab-to-cm-energy]]
            [complex :refer [mag mul add complex-cartesian]]))

(def ^:private o16-dp-bound-ell 2)
(def ^:private o16-dp-J-i 0.0)
(def ^:private o16-dp-J-f 2.5)

(defn- m16 [] (* 16.0 931.494))
(defn- m17 [] (* 17.0 931.494))

(defn o16-dp-kinematics
  "Map **:mass-factor-i :mass-factor-f :e-cm-i :e-cm-f :k-i :k-f :Q-mev :M-target :M-residual**
  for **d + ¹⁶O → p + ¹⁷O**.

  **e-cm-i** — CM kinetic (MeV) in **entrance** (**d+¹⁶O**). Nullary **`()`** uses
  **E_lab(d)=20 MeV** on **¹⁶O** at rest: **`lab-to-cm-energy 20 m_d m_16O`**."
  ([]
   (o16-dp-kinematics (lab-to-cm-energy 20.0 1875.613 (m16))))
  ([^double e-cm-i]
   (let [m-p 938.272
         m-d 1875.613
         m16 (m16)
         m17 (m17)
         Q (- (+ m16 m-d) (+ m17 m-p))
         e-cm-f (+ e-cm-i Q)
         _ (when (< e-cm-f 0.5)
             (throw (ex-info "o16-dp-kinematics: e-cm-f too low — raise e-cm-i or check masses"
                             {:e-cm-i e-cm-i :Q Q :e-cm-f e-cm-f})))
         mu-i (/ (* m-d m16) (+ m-d m16))
         mu-f (/ (* m-p m17) (+ m-p m17))
         mfi (mass-factor-from-mu mu-i)
         mff (mass-factor-from-mu mu-f)]
     {:mass-factor-i mfi :mass-factor-f mff
      :e-cm-i e-cm-i :e-cm-f e-cm-f
      :k-i (Math/sqrt (* mfi e-cm-i))
      :k-f (Math/sqrt (* mff e-cm-f))
      :Q-mev Q
      :M-target m16 :M-residual m17})))

(defn- r0-sc ^double [^double r-ca ^double a ^double b]
  (* r-ca (Math/pow (/ a b) (/ 1.0 3.0))))

(defn- a16 ^double [] (Math/pow 16.0 (/ 1.0 3.0)))

(defn optical-u-deuteron-o16
  "Woods–Saxon **d + ¹⁶O**: Table 5.1 from the handbook, **10 MeV/u** (E_lab = 20 MeV deuteron).
  Radii **R = r₀ × 16^{1/3}**: real V = 88.955, r₀ = 1.149, a = 0.751;
  volume imaginary Wi = 2.348, r₀i = 1.345, ai = 0.603;
  surface imaginary Ws = 10.218, r₀s = 1.397, as = 0.687 (derivative WS: −4Wsf(1−f));
  spin–orbit Vso = 3.557, r₀so = 0.972, aso = 1.011; Coulomb Rc = 1.303 × 16^{1/3}."
  [L s j]
  (let [z1 1 z2 8
        a16 (a16)
        rc (* 1.303 a16)]
    (fn [^double r]
      (t/optical-potential-woods-saxon r
                                       [88.955  (* 1.149 a16)  0.751]   ; real WS
                                       [2.348   (* 1.345 a16)  0.603]   ; volume imaginary
                                       [10.218  (* 1.397 a16)  0.687]   ; surface imaginary
                                       3.557 (* 0.972 a16) 1.011       ; spin-orbit
                                       L s j z1 z2 rc))))

(defn- a17 ^double [] (Math/pow 17.0 (/ 1.0 3.0)))

(defn optical-u-proton-o17
  "Woods–Saxon **p + ¹⁷O**: Table 5.2 from the handbook, **10 MeV/u** context.
  Radii **R = r₀ × 17^{1/3}**: real V = 49.544, r₀ = 1.146, a = 0.675;
  volume imaginary Vi = 2.061, same r₀ and a as real;
  surface imaginary Vs = 7.670, r₀s = 1.302, as = 0.528 (derivative WS: −4Vsf(1−f));
  spin–orbit Vso = 5.296 (real part; Im = 0.106 neglected), r₀so = 0.934, aso = 0.590;
  Coulomb Rc = 1.419 × 17^{1/3}."
  [L s j]
  (let [z1 1 z2 8
        a17 (a17)
        rc  (* 1.419 a17)]
    (fn [^double r]
      (t/optical-potential-woods-saxon r
                                       [49.544 (* 1.146 a17) 0.675]  ; real WS
                                       [2.061  (* 1.146 a17) 0.675]  ; volume imaginary
                                       [7.670  (* 1.302 a17) 0.528]  ; surface imaginary
                                       5.296 (* 0.934 a17) 0.590    ; spin-orbit (real)
                                       L s j z1 z2 rc))))

(defn- deuteron-j-for-partial-wave
  ^double [^long L]
  (+ 1.0 (double L)))

(defn- proton-j-for-partial-wave
  ^double [^long L]
  (+ 0.5 (double L)))

(defn- sommerfeld-eta-channel
  ^double [^double e-cm ^double mfactor ^double z1z2ee]
  (binding [fn/mass-factor mfactor
            fn/Z1Z2ee z1z2ee]
    (channel-sommerfeld-eta e-cm)))

(defn- o16-dp-cm-asymmetry-factor
  ^double [^double theta-rad ^double kappa]
  (if (< (Math/abs kappa) 1e-15)
    1.0
    (max 1e-300 (+ 1.0 (* kappa (Math/cos theta-rad))))))

(defn- o16-dp-rows-coulomb-sigma
  [base-rows eta-i eta-f]
  (t/handbook-zr-rows-with-coulomb-sigma base-rows eta-i eta-f))

(defn- o16-dp-filter-rows-by-L-alpha
  [rows L-alpha-only]
  (if (nil? L-alpha-only)
    rows
    (let [La (long L-alpha-only)]
      (filterv #(= (long (:L-alpha %)) La) rows))))

(defn o16-dp-radial-I-rows-handbook
  "Build **{:L-alpha :L-beta :I}** with **`handbook-radial-integral-I-zr-from-neutron-bound`**.
  **φ_n** — neutron in **¹⁷O** (**l=2**, illustrative well); **χ_α** — **d** on **¹⁶O**; **χ_β** — **p** on **¹⁷O**.

  **`:chi-normalize-mode`** — **`:coulomb-tail`** (default), **`:raw`**, or **`:max`** (**`distorted-wave-optical`**)."
  [& {:keys [r-max h L-max e-cm-i transfer-ell chi-normalize-mode]
      :or {r-max 100.0 h 0.05 L-max 20 transfer-ell o16-dp-bound-ell
           chi-normalize-mode :coulomb-tail}}]
  (when-not (#{:raw :max :coulomb-tail} chi-normalize-mode)
    (throw (ex-info "o16-dp-radial-I-rows-handbook: :chi-normalize-mode must be :raw, :max, or :coulomb-tail"
                    {:chi-normalize-mode chi-normalize-mode})))
  (let [e-cm-i (double (or e-cm-i (:e-cm-i (o16-dp-kinematics))))
        {:keys [mass-factor-i mass-factor-f e-cm-f k-i k-f M-target M-residual]}
        (o16-dp-kinematics e-cm-i)
        zr (t/handbook-zr-chi-exit-mass-ratio M-target M-residual)
        ell (long transfer-ell)
        phi-n (t/normalize-bound-state
               (t/solve-bound-state-numerov -4.1438 2 56.0 2.85 0.6 0.048 h r-max {:no-spin-orbit true}) h)
        z12 (* 1.44 1.0 8.0)
        eta-i (sommerfeld-eta-channel e-cm-i mass-factor-i z12)
        eta-f (sommerfeld-eta-channel e-cm-f mass-factor-f z12)
        rho-i (* k-i r-max)
        rho-f (* k-f r-max)
        chi-alpha!
        (memoize
         (fn [^long La]
           (if (= :raw chi-normalize-mode)
             (t/distorted-wave-optical e-cm-i La 1.0 (deuteron-j-for-partial-wave La)
                                       (optical-u-deuteron-o16 La 1.0 (deuteron-j-for-partial-wave La))
                                       r-max h mass-factor-i)
             (t/distorted-wave-optical e-cm-i La 1.0 (deuteron-j-for-partial-wave La)
                                       (optical-u-deuteron-o16 La 1.0 (deuteron-j-for-partial-wave La))
                                       r-max h mass-factor-i
                                       :normalize-mode chi-normalize-mode
                                       :tail-eta eta-i :tail-rho rho-i))))
        chi-beta!
        (memoize
         (fn [^long Lb]
           (if (= :raw chi-normalize-mode)
             (t/distorted-wave-optical e-cm-f Lb 0.5 (proton-j-for-partial-wave Lb)
                                       (optical-u-proton-o17 Lb 0.5 (proton-j-for-partial-wave Lb))
                                       r-max h mass-factor-f)
             (t/distorted-wave-optical e-cm-f Lb 0.5 (proton-j-for-partial-wave Lb)
                                       (optical-u-proton-o17 Lb 0.5 (proton-j-for-partial-wave Lb))
                                       r-max h mass-factor-f
                                       :normalize-mode chi-normalize-mode
                                       :tail-eta eta-f :tail-rho rho-f))))]
    (vec
     (for [La (range 0 (inc (long L-max)))
           Lb (t/handbook-zr-partial-wave-L-beta-values La ell (long L-max))
           :let [Ireal (t/handbook-radial-integral-I-zr-from-neutron-bound
                        phi-n (chi-alpha! La) (chi-beta! Lb) h
                        M-target M-residual k-i k-f zr)]
           :when (> (Math/abs (double Ireal)) 1e-30)]
       {:L-alpha La :L-beta Lb :I Ireal}))))

(defn o16-dp-T-squared-sum-handbook
  [theta-rad radial-rows-sigma D0 & {:keys [coherent-m-beta?] :or {coherent-m-beta? false}}]
  (let [ell (long o16-dp-bound-ell)
        sqrt2l1 (Math/sqrt (inc (* 2.0 (double ell))))
        pref (complex-cartesian (* (double D0) sqrt2l1) 0.0)
        ms (range (- ell) (inc ell))]
    (if coherent-m-beta?
      (let [sum-b (reduce (fn [acc ^long m-ell]
                            (add acc (mul pref (t/handbook-zr-multipole-amplitude-sum
                                                  ell m-ell theta-rad radial-rows-sigma))))
                          (complex-cartesian 0.0 0.0)
                          ms)
            Tmag (mag sum-b)]
        (* Tmag Tmag))
      (double
       (reduce
        (fn [^double acc ^long m-ell]
          (let [beta (t/handbook-zr-multipole-amplitude-sum ell m-ell theta-rad radial-rows-sigma)
                Tm (mul pref beta)
                Tmag (mag Tm)]
            (+ acc (* Tmag Tmag))))
        0.0
        ms)))))

(defn o16-dp-dsigma-handbook-mb-sr
  "**dσ/dΩ (mb/sr)** — handbook **F_n**, **(5.5)** radial **I**, **`handbook-zr-multipole-amplitude-sum`**, **`transfer-differential-cross-section`**.

  **`:chi-normalize-mode`** — passed to **`o16-dp-radial-I-rows-handbook`** when **`:radial-rows-sigma`** is omitted (**`:coulomb-tail`** default)."
  [theta-deg & {:keys [e-cm-i r-max h L-max S-factor radial-rows-sigma
                       coherent-m-beta? cm-asymmetry-kappa L-alpha-only chi-normalize-mode]
                :or {r-max 100.0 h 0.05 L-max 20 S-factor 1.0
                     coherent-m-beta? false
                     cm-asymmetry-kappa 0.0
                     chi-normalize-mode :coulomb-tail}}]
  (let [eci (if (some? e-cm-i) (double e-cm-i) (:e-cm-i (o16-dp-kinematics)))
        {:keys [mass-factor-i mass-factor-f e-cm-f k-i k-f]} (o16-dp-kinematics eci)
        z12 (* 1.44 1.0 8.0)
        eta-i (sommerfeld-eta-channel eci mass-factor-i z12)
        eta-f (sommerfeld-eta-channel e-cm-f mass-factor-f z12)
        rows-sig (o16-dp-filter-rows-by-L-alpha
                  (or radial-rows-sigma
                      (o16-dp-rows-coulomb-sigma
                       (o16-dp-radial-I-rows-handbook :r-max r-max :h h :L-max L-max :e-cm-i eci
                         :chi-normalize-mode chi-normalize-mode)
                       eta-i eta-f))
                  L-alpha-only)
        D0 (t/zero-range-constant :d-p)
        theta-rad (* (double theta-deg) (/ Math/PI 180.0))
        T-sq (o16-dp-T-squared-sum-handbook theta-rad rows-sig D0
               :coherent-m-beta? (boolean coherent-m-beta?))
        spin (* (t/transfer-nuclear-spin-statistical-factor o16-dp-J-i o16-dp-J-f)
                (t/transfer-unpolarized-deuteron-spin-factor))
        ds (t/transfer-differential-cross-section (Math/sqrt (max T-sq 0.0)) S-factor k-i k-f
                                                  mass-factor-i mass-factor-f)
        asym (o16-dp-cm-asymmetry-factor theta-rad (double cm-asymmetry-kappa))]
    (* (double ds) spin asym)))

(defn o16-dp-angular-curve-handbook-mb-sr
  [theta-degrees & {:keys [e-cm-i r-max h L-max S-factor
                           coherent-m-beta? cm-asymmetry-kappa L-alpha-only chi-normalize-mode]
                    :or {r-max 100.0 h 0.05 L-max 20 S-factor 1.0
                         coherent-m-beta? false
                         cm-asymmetry-kappa 0.0
                         chi-normalize-mode :coulomb-tail}}]
  (let [eci (if (some? e-cm-i) (double e-cm-i) (:e-cm-i (o16-dp-kinematics)))
        {:keys [mass-factor-i mass-factor-f e-cm-i e-cm-f]} (o16-dp-kinematics eci)
        z12 (* 1.44 1.0 8.0)
        eta-i (sommerfeld-eta-channel e-cm-i mass-factor-i z12)
        eta-f (sommerfeld-eta-channel e-cm-f mass-factor-f z12)
        base-rows (o16-dp-radial-I-rows-handbook :r-max r-max :h h :L-max L-max :e-cm-i e-cm-i
                   :chi-normalize-mode chi-normalize-mode)
        rows-sig (o16-dp-filter-rows-by-L-alpha
                  (o16-dp-rows-coulomb-sigma base-rows eta-i eta-f)
                  L-alpha-only)]
    (mapv (fn [^double th]
            {:theta-deg th
             :differential_cross_section_mb_sr
             (o16-dp-dsigma-handbook-mb-sr th
               :radial-rows-sigma rows-sig
               :e-cm-i eci :r-max r-max :h h :L-max L-max :S-factor S-factor
               :coherent-m-beta? coherent-m-beta?
               :cm-asymmetry-kappa cm-asymmetry-kappa)})
          theta-degrees)))
