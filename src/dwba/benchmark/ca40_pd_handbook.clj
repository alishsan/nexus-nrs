(ns dwba.benchmark.ca40-pd-handbook
  "⁴⁰Ca(p,d)³⁹Ca g.s. pickup — **handbook** (*Handbook of direct nuclear reaction for retarded theorist*)
  ZR pipeline only (no separate Austern-named entry points).

  **Radial:** **F_{ℓsj} = R_n = u/r**; **I** from **`handbook-radial-integral-I-zr-from-neutron-bound`** (**(5.5)** **(M_B/M_A)(4π/(k_α k_β))** on **∫ F R_α R_β r² dr**).

  **Angular:** **`handbook-zr-multipole-amplitude-sum`** + **`handbook-zr-rows-with-coulomb-sigma`** (**σ_L** on rows); **`handbook-zr-partial-wave-L-beta-values`**
  for **L_β + ℓ ↔ L_α** (same triangle as the Austern benchmark, different public API).

  **Amplitude:** **T_m ≈ D₀ √(2ℓ+1) β_m** and **`transfer-differential-cross-section`**. **D₀**: cluster / ZR strength. Spin weights match **`ca40-pd-austern`**."
  (:require [dwba.transfer :as t]
            [functions :as fn :refer [mass-factor-from-mu channel-sommerfeld-eta]]
            [complex :refer [mag mul add complex-cartesian]]))

(def ^:private ca40-pd-bound-ell 3)
(def ^:private ca40-pd-J-i 0.0)
(def ^:private ca40-pd-J-f 3.5)

(defn ca40-pd-kinematics
  "Same as **`dwba.benchmark.ca40-pd-austern/ca40-pd-kinematics`**."
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
  [L s j]
  (let [z1 1 z2 20
        rc (* 4.3103 (Math/pow (/ 40.0 41.0) (/ 1.0 3.0)))]
    (fn [^double r]
      (t/optical-potential-woods-saxon r [49.47 4.0689 0.70] [19.8 4.3172 0.75]
                                       24.2 4.0689 0.70 L s j z1 z2 rc))))

(defn optical-u-deuteron-ca39
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
  ^double [^double theta-rad ^double kappa]
  (if (< (Math/abs kappa) 1e-15)
    1.0
    (max 1e-300 (+ 1.0 (* kappa (Math/cos theta-rad))))))

(defn- ca40-pd-rows-coulomb-sigma
  [base-rows eta-i eta-f]
  (t/handbook-zr-rows-with-coulomb-sigma base-rows eta-i eta-f))

(defn- ca40-pd-filter-rows-by-L-alpha
  [rows L-alpha-only]
  (if (nil? L-alpha-only)
    rows
    (let [La (long L-alpha-only)]
      (filterv #(= (long (:L-alpha %)) La) rows))))

(defn ca40-pd-radial-I-rows-handbook
  "Build **{:L-alpha :L-beta :I}** with **`t/handbook-radial-integral-I-zr-from-neutron-bound`** (nucleon-only **F**)."
  [& {:keys [r-max h L-max e-cm-i transfer-ell]
      :or {r-max 100.0 h 0.05 L-max 20 e-cm-i 18.0 transfer-ell ca40-pd-bound-ell}}]
  (let [{:keys [mass-factor-i mass-factor-f e-cm-f k-i k-f M-target M-residual]}
        (ca40-pd-kinematics e-cm-i)
        zr (t/handbook-zr-chi-exit-mass-ratio M-target M-residual)
        ell (long transfer-ell)
        phi-n (t/normalize-bound-state
               (t/solve-bound-state-numerov -8.364 3 58.4538 4.0355 0.7 0.048 h r-max) h)
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
           Lb (t/handbook-zr-partial-wave-L-beta-values La ell (long L-max))
           :let [Ireal (t/handbook-radial-integral-I-zr-from-neutron-bound
                        phi-n (chi-alpha! La) (chi-beta! Lb) h
                        M-target M-residual k-i k-f zr)]
           :when (> (Math/abs (double Ireal)) 1e-30)]
       {:L-alpha La :L-beta Lb :I Ireal}))))

(defn ca40-pd-T-squared-sum-handbook
  [theta-rad radial-rows-sigma D0 & {:keys [coherent-m-beta?] :or {coherent-m-beta? false}}]
  (let [ell (long ca40-pd-bound-ell)
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

(defn ca40-pd-dsigma-handbook-mb-sr
  "**dσ/dΩ (mb/sr)** — handbook **F_n**, **(5.5)** radial **I**, **`handbook-zr-multipole-amplitude-sum`**, **`transfer-differential-cross-section`**."
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
                        (ca40-pd-radial-I-rows-handbook :r-max r-max :h h :L-max L-max :e-cm-i e-cm-i)
                        eta-i eta-f))
                  L-alpha-only)
        D0 (t/zero-range-constant :p-d)
        theta-rad (* (double theta-deg) (/ Math/PI 180.0))
        T-sq (ca40-pd-T-squared-sum-handbook theta-rad rows-sig D0
               :coherent-m-beta? (boolean coherent-m-beta?))
        spin (* (t/transfer-nuclear-spin-statistical-factor ca40-pd-J-i ca40-pd-J-f)
                (t/transfer-unpolarized-deuteron-spin-factor)
                0.5)
        ds (t/transfer-differential-cross-section (Math/sqrt (max T-sq 0.0)) S-factor k-i k-f
                                                   mass-factor-i mass-factor-f)
        asym (ca40-pd-cm-asymmetry-factor theta-rad (double cm-asymmetry-kappa))]
    (* (double ds) spin asym)))

(defn ca40-pd-angular-curve-handbook-mb-sr
  [theta-degrees & {:keys [e-cm-i r-max h L-max S-factor
                           coherent-m-beta? cm-asymmetry-kappa L-alpha-only]
                    :or {e-cm-i 18.0 r-max 100.0 h 0.05 L-max 20 S-factor 1.0
                         coherent-m-beta? false
                         cm-asymmetry-kappa 0.0}}]
  (let [{:keys [mass-factor-i mass-factor-f e-cm-i e-cm-f]} (ca40-pd-kinematics e-cm-i)
        z12 (* 1.44 1.0 20.0)
        eta-i (sommerfeld-eta-channel e-cm-i mass-factor-i z12)
        eta-f (sommerfeld-eta-channel e-cm-f mass-factor-f z12)
        base-rows (ca40-pd-radial-I-rows-handbook :r-max r-max :h h :L-max L-max :e-cm-i e-cm-i)
        rows-sig (ca40-pd-filter-rows-by-L-alpha
                  (ca40-pd-rows-coulomb-sigma base-rows eta-i eta-f)
                  L-alpha-only)]
    (mapv (fn [^double th]
            {:theta-deg th
             :differential_cross_section_mb_sr
             (ca40-pd-dsigma-handbook-mb-sr th
               :radial-rows-sigma rows-sig
               :e-cm-i e-cm-i :r-max r-max :h h :L-max L-max :S-factor S-factor
               :coherent-m-beta? coherent-m-beta?
               :cm-asymmetry-kappa cm-asymmetry-kappa)})
          theta-degrees)))
