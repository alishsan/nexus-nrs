(ns dwba.benchmark.li11-pp-paper1708
  "¹¹Li**(p,p′)** DWBA aligned with **Tanaka et al., arXiv:1708.07719** (Table 1 optical **Set S** / **Set V**).

  **Scope**
  - **Optics:** Woods–Saxon + Thomas SO + Coulomb + imaginary (**volume** for **V**, **derivative surface** for **S**)
    as in the paper’s Eq. (2). **CHUCK3** micro-differences may remain (e.g. SO **4(ħ/m_π c)²** folded into **V_so**).
  - **Kinematics:** normal kinematics **p** lab on ¹¹Li at rest (**E_p,lab ≈ 6 MeV** ↔ **66 MeV** ¹¹Li inverse kinematics
    on **H₂**). CM kinetic **E_cm**, exit **E_cm − E_x** with **E_x = 0.80 MeV** (paper peak).
  - **Multipole:** **Austern (5.5)** radial **I_{L_β L_α}** with **complex** **u** (**`austern-radial-integral-I-zr-eq-5-5-from-u-complex`**) and
    **(5.6)** angular sum (**`austern-reduced-amplitude-beta-sum-eq-5-6`**, Coulomb **σ_L** only on rows — same pattern as
    **`dwba.benchmark.ca40-pd-austern`**). **ℓ = 1** (dipole / **ΔL = 1** coupling).
  - **Form factor:** macroscopic **F(r) = β_{sc} R_V (−d/dr)(V_v f_v)** (surface vibrational proxy). Paper **Harakeh–Dieperink**
    / **Orlandini** E1(isoscalar) shapes are **not** encoded here — tune **`:beta-scale`** for absolute comparison to Fig. 5.

  **Set V vs Set S (same β_sc):** **Set V** uses **volume** **W** with **r_w A^{1/3} ≈ 0.61 A^{1/3} fm** (inward); the schematic **FF** peaks in the **surface**. **χ** can stay large where **F** is large, so **dσ** may exceed **Set S** (surface derivative **W**) by **orders** — mostly **geometry + proxy FF**. **dσ ∝ β_sc²** (roughly); rescale **β_sc** by **√(dσ_ref/dσ_target)** when matching a figure.

  **REPL:** `(require '[dwba.benchmark.li11-pp-paper1708 :as li11])`
  then `(li11/li11-pp-dsigma-paper-mb-sr 90.0 :optical-set :V)`."
  (:require [dwba.transfer :as t]
            [dwba.inelastic :as inel]
            [functions :as fn :refer [mass-factor-from-mu channel-sommerfeld-eta]]
            [complex :refer [mag mul add complex-cartesian]]))

(def ^:private m-p 938.272)
(def ^:private m-11Li (* 11.0 931.494))
(def ^:private Z-target 3)
(def ^:private A-target 11)

(defn- A13 []
  (Math/pow A-target (/ 1.0 3.0)))

;; ---------------------------------------------------------------------------
;; Table 1 — arXiv:1708.07719 (parameter *r_i* multiplies **A^{1/3}** in **f(r,r_i,a_i)**).
;; ---------------------------------------------------------------------------

(def ^:private paper-set-V
  "Volume imaginary **W_w**; no surface **W_s**."
  {:Vv 54.2 :rv 1.16 :av 0.75
   :Vso 6.23 :rso 1.16 :aso 0.75
   :rC 1.16
   :Ww 14.3 :rw 0.61 :aw 1.98
   :Ws 0.0 :rs 0.0 :as 0.0})

(def ^:private paper-set-S
  "Surface imaginary **4 a_s d/dr(W_s f)**; **W_w = 0**."
  {:Vv 35.2 :rv 1.71 :av 0.92
   :Vso 7.95 :rso 1.71 :aso 0.92
   :rC 1.71
   :Ww 0.0 :rw 0.0 :aw 0.0
   :Ws 14.4 :rs 1.49 :as 0.54})

(defn- ws-f ^double [^double r ^double R ^double a]
  (/ 1.0 (+ 1.0 (Math/exp (/ (- r R) a)))))

(defn- ws-df-dr ^double [^double r ^double R ^double a]
  (let [f (ws-f r R a)]
    (- (/ (* f (- 1.0 f)) a))))

(defn- l-dot-s ^double [^double l ^double s ^double j]
  (/ (- (* j (+ j 1.0)) (* l (+ l 1.0)) (* s (+ s 1.0))) 2.0))

(defn- coulomb-v ^double [^double r ^double Rc ^double z1z2ee]
  (if (> r Rc)
    (/ z1z2ee r)
    (* z1z2ee (/ (- 3.0 (/ (* r r) (* Rc Rc))) (* 2.0 Rc)))))

(defn li11-paper-optical-fn
  "Complex **U(r)** for **p + ¹¹Li** — **`:optical-set` **`:V`** | **`:S`** (Table 1).
  **l**, **s**, **j** — same partial-wave labels as **`distorted-wave-optical`**."
  [optical-set ^double l ^double s ^double j]
  (let [pm (case optical-set :V paper-set-V :S paper-set-S
                            (throw (ex-info "li11-paper-optical-fn: :optical-set must be :V or :S"
                                            {:optical-set optical-set})))
        a13 (A13)
        Rv (* (:rv pm) a13)
        av (:av pm)
        Rso (* (:rso pm) a13)
        aso (:aso pm)
        Rc (* (:rC pm) a13)
        z1z2ee (* 1.44 1.0 Z-target)
        Ww (:Ww pm)
        Rw (* (:rw pm) a13)
        aw (:aw pm)
        Ws (:Ws pm)
        Rs (* (:rs pm) a13)
        asurf (:as pm)
        Vv (:Vv pm)
        Vso (:Vso pm)]
    (fn [^double r]
      (let [fV (ws-f r Rv av)
            Vcentr (* -1.0 Vv fV)
            fso (ws-f r Rso aso)
            dfso (* -1.0 (/ (* fso (- 1.0 fso)) aso))
            Vls (if (< r 1e-12)
                  0.0
                  (* Vso (l-dot-s l s j) dfso (/ 1.0 r)))
            Vc (coulomb-v r Rc z1z2ee)
            Ureal (+ Vcentr Vls Vc)
            ;; Imaginary: volume − W_w f_w (code convention: Cartesian imag part = W_im)
            W-vol (- (* Ww (ws-f r Rw aw)))
            ;; Surface: + 4 a_s W_s (d/dr)f_s  — paper **i 4 a_s d/dr(W_s f)**
            W-surf (* 4.0 asurf Ws (ws-df-dr r Rs asurf))
            W-tot (+ W-vol W-surf)]
        (complex-cartesian Ureal W-tot)))))

(defn li11-pp-kinematics
  "**E_p,lab** (MeV) proton on ¹¹Li at rest; **E_ex** (MeV) excitation (**0.80** in paper).
  Returns **:e-cm-i :e-cm-f :mass-factor-i :mass-factor-f :k-i :k-f :M-entrance :M-exit**."
  [^double e-p-lab-mev ^double e-ex-mev]
  (let [e-cm (* e-p-lab-mev (/ m-11Li (+ m-p m-11Li)))
        e-cm-f (- e-cm e-ex-mev)
        _ (when (< e-cm-f 0.05)
            (throw (ex-info "li11-pp-kinematics: exit CM energy too low"
                            {:e-p-lab e-p-lab-mev :e-ex e-ex-mev :e-cm-f e-cm-f})))
        mu (/ (* m-p m-11Li) (+ m-p m-11Li))
        mf (mass-factor-from-mu mu)
        Mi (+ m-p m-11Li)
        ;; Neglect **MeV**–level mass shift of ¹¹Li* on **M_B/M_A** ratio.
        Mf Mi]
    {:e-cm-i e-cm :e-cm-f e-cm-f
     :mass-factor-i mf :mass-factor-f mf
     :k-i (Math/sqrt (* mf e-cm))
     :k-f (Math/sqrt (* mf e-cm-f))
     :M-entrance Mi :M-exit Mf}))

(defn- proton-j ^double [^long L]
  (+ 0.5 (double L)))

(defn- sommerfeld-eta ^double [^double e-cm ^double mfactor]
  (binding [fn/mass-factor mfactor
            fn/Z1Z2ee (* 1.44 1.0 Z-target)]
    (channel-sommerfeld-eta e-cm)))

(defn li11-paper-collective-F-vec
  "Macroscopic transition factor **F(r_i) ≈ β_sc R_V × (−d/dr)(V_v f_V)** on **r = i h**
  (scales with **β_sc**). Uses the **real volume** geometry of the chosen optical set."
  [optical-set ^double h ^long n-points beta-scale]
  (let [pm (case optical-set :V paper-set-V :S paper-set-S)
        Rv (* (:rv pm) (A13))
        Vv (:Vv pm)
        av (:av pm)]
    (mapv (fn [^long i]
            (let [r (* (double i) h)
                  ;; **U = −V_v f**  ⇒  **dU/dr = − d(V_v f)/dr**
                  dU-dr (- (inel/woods-saxon-derivative r [Vv Rv av]))]
              (* (double beta-scale) Rv dU-dr)))
          (range n-points))))

(defn li11-paper-radial-I-rows
  "Build **{:L-alpha :L-beta :I}** with **complex** **I** (absorbing optics). **ℓ** = **1** (**ΔL = 1** dipole coupling)."
  [& {:keys [optical-set e-p-lab-mev e-ex-mev r-max h L-max beta-scale transfer-ell]
      :or {optical-set :V e-p-lab-mev 6.0 e-ex-mev 0.8 r-max 45.0 h 0.02 L-max 22
           beta-scale 1.0 transfer-ell 1}}]
  (let [{:keys [e-cm-i e-cm-f mass-factor-i mass-factor-f k-i k-f M-entrance M-exit]}
        (li11-pp-kinematics e-p-lab-mev e-ex-mev)
        zr (t/austern-zr-chi-exit-mass-ratio M-entrance M-exit)
        ell (long transfer-ell)
        n-pts (inc (int (/ (double r-max) (double h))))
        F-vec (li11-paper-collective-F-vec optical-set (double h) n-pts (double beta-scale))
        eta-i (sommerfeld-eta e-cm-i mass-factor-i)
        eta-f (sommerfeld-eta e-cm-f mass-factor-f)
        rho-i (* k-i r-max)
        rho-f (* k-f r-max)
        chi-a! (memoize
                (fn [^long La]
                  (t/distorted-wave-optical
                    e-cm-i La 0.5 (proton-j La)
                    (li11-paper-optical-fn optical-set (double La) 0.5 (proton-j La))
                    r-max h mass-factor-i
                    :normalize-mode :coulomb-tail
                    :tail-eta eta-i
                    :tail-rho rho-i)))
        chi-b! (memoize
                (fn [^long Lb]
                  (t/distorted-wave-optical
                    e-cm-f Lb 0.5 (proton-j Lb)
                    (li11-paper-optical-fn optical-set (double Lb) 0.5 (proton-j Lb))
                    r-max h mass-factor-f
                    :normalize-mode :coulomb-tail
                    :tail-eta eta-f
                    :tail-rho rho-f)))]
    (vec
      (for [La (range 0 (inc (long L-max)))
            Lb (t/austern-eq-5-6-admissible-L-beta-values La ell (long L-max))
            :let [I (t/austern-radial-integral-I-zr-eq-5-5-from-u-complex
                      F-vec (chi-a! La) (chi-b! Lb) h
                      M-entrance M-exit k-i k-f zr)]
            :when (> (mag I) 1e-30)]
        {:L-alpha La :L-beta Lb :I I}))))

(defn- li11-pp-T-squared-sum-from-beta
  "Incoherent **Σ_m |√(2ℓ+1) β_m|²** (no **D₀** — collective strength lives in **F** and **β_sc**)."
  [theta-rad radial-rows-sigma ^long transfer-ell]
  (let [ell (long transfer-ell)
        sqrt2l1 (Math/sqrt (inc (* 2.0 (double ell))))
        pref (complex-cartesian sqrt2l1 0.0)
        ms (range (- ell) (inc ell))]
    (double
     (reduce
      (fn [^double acc ^long m-ell]
        (let [beta (t/austern-reduced-amplitude-beta-sum-eq-5-6 ell m-ell theta-rad radial-rows-sigma)
              Tm (mul pref beta)
              Tmag (mag Tm)]
          (+ acc (* Tmag Tmag))))
      0.0
      ms))))

(defn li11-pp-dsigma-paper-mb-sr
  "**dσ/dΩ** (**mb/sr**, CM) at lab angle **θ_cm** (degrees). Uses **`transfer-differential-cross-section`**
  with **μ_i ≈ μ_f** (same reduced mass). Tune **`:beta-scale`** for absolute agreement with Fig. 5.

  **Options:** **`:optical-set`** **`:V`** or **`:S`**, **`:e-p-lab-mev`** (default **6**), **`:e-ex-mev`** (**0.80**),
  **`:r-max`**, **`:h`**, **`:L-max`**, **`:beta-scale`**, **`:transfer-ell`** (default **1**), **`:radial-rows-sigma`**, **`:S-factor`**."
  [theta-deg & {:keys [optical-set e-p-lab-mev e-ex-mev r-max h L-max beta-scale transfer-ell
                       radial-rows-sigma S-factor]
                :or {optical-set :V e-p-lab-mev 6.0 e-ex-mev 0.8 r-max 45.0 h 0.02 L-max 22
                     beta-scale 1.0 transfer-ell 1 S-factor 1.0}}]
  (let [{:keys [mass-factor-i mass-factor-f k-i k-f e-cm-i e-cm-f]}
        (li11-pp-kinematics e-p-lab-mev e-ex-mev)
        eta-i (sommerfeld-eta e-cm-i mass-factor-i)
        eta-f (sommerfeld-eta e-cm-f mass-factor-f)
        base-rows (or radial-rows-sigma
                      (li11-paper-radial-I-rows :optical-set optical-set :e-p-lab-mev e-p-lab-mev
                                                 :e-ex-mev e-ex-mev :r-max r-max :h h :L-max L-max
                                                 :beta-scale beta-scale :transfer-ell transfer-ell))
        rows-sig (t/austern-radial-rows-with-coulomb-sigma base-rows eta-i eta-f)
        theta-rad (* (double theta-deg) (/ Math/PI 180.0))
        T-sq (li11-pp-T-squared-sum-from-beta theta-rad rows-sig (long transfer-ell))
        ds (t/transfer-differential-cross-section (Math/sqrt (max T-sq 0.0)) (double S-factor) k-i k-f
                                                  mass-factor-i mass-factor-f)]
    (double ds)))

(defn li11-pp-angular-curve-paper-mb-sr
  "Vector **{:theta-deg … :differential_cross_section_mb_sr …}**; precomputes **σ** rows once."
  [theta-degrees & {:keys [optical-set e-p-lab-mev e-ex-mev r-max h L-max beta-scale transfer-ell S-factor]
                    :or {optical-set :V e-p-lab-mev 6.0 e-ex-mev 0.8 r-max 45.0 h 0.02 L-max 22
                         beta-scale 1.0 transfer-ell 1 S-factor 1.0}}]
  (let [{:keys [e-cm-i e-cm-f mass-factor-i mass-factor-f]} (li11-pp-kinematics e-p-lab-mev e-ex-mev)
        eta-i (sommerfeld-eta e-cm-i mass-factor-i)
        eta-f (sommerfeld-eta e-cm-f mass-factor-f)
        base-rows (li11-paper-radial-I-rows :optical-set optical-set :e-p-lab-mev e-p-lab-mev
                                             :e-ex-mev e-ex-mev :r-max r-max :h h :L-max L-max
                                             :beta-scale beta-scale :transfer-ell transfer-ell)
        rows-sig (t/austern-radial-rows-with-coulomb-sigma base-rows eta-i eta-f)]
    (mapv (fn [^double th]
            {:theta-deg th
             :differential_cross_section_mb_sr
             (li11-pp-dsigma-paper-mb-sr th :optical-set optical-set :e-p-lab-mev e-p-lab-mev :e-ex-mev e-ex-mev
                                         :r-max r-max :h h :L-max L-max :beta-scale beta-scale
                                         :transfer-ell transfer-ell :radial-rows-sigma rows-sig
                                         :S-factor S-factor)})
          theta-degrees)))