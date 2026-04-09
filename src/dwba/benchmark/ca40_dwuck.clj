(ns dwba.benchmark.ca40-dwuck
  "Ca40(d,p)Ca41 g.s. (L=3, f7/2 stripping) using the same kinematics and **fitted** optical
  parameters printed in a DWUCK4 listing (first case in DW4TST.DAT).

  **Angular model — `:angular-mode` on `ca40-dp-dsigma-mb-sr`**

  - **Default `:coherent`** — **`transfer-differential-cross-section-angular-coherent`**
    (**`|Σ_L w_L T_L Y_{L0}|²`** / multipole interference). Gives **θ-dependent** dσ/dΩ for this benchmark
    (still **CM symmetric** for a single dominant **L**). This is the usual **reduced** DWBA angular factor
    when you are not summing over bound **m_i**, **m_f** explicitly.
  - **`:m-sum`** — **`transfer-differential-cross-section-angular-m-sum`**
    (**`transfer-angular-distribution-m-sum-unpolarized`**). For **l_i=0**, **l_f=3** the orbital piece
    collapses to **|wT|² (2l_f+1)/(4π)** — **no θ** (addition theorem). dσ/dΩ is **constant in θ**; use
    this only when you want that isotropic orbital average.

  Then **(2J_f+1)/(2J_i+1)** and unpolarized deuteron **1/3**. Full **6j** / DWUCK **(1.9)–(1.11)** is still
  **not** folded here.

  Default **`:Ls`** is **`[3]`** — triangle+parity for **l_i=0 → l_f=3** allows **only L=3**.

  **DWBA content (no phenomenological Coulomb-on-σ hacks)**

  - **Coulomb** enters through **`optical-potential-woods-saxon`** in **`distorted-wave-optical`** (χ_d, χ_p).
  - **Radial** zero-range POST: **`transfer-amplitude-post`** → **T_L** per partial wave **L**, with
    **handbook** single-nucleon **F_{ℓsj}=R** from **`:handbook-F-from :phi-f`** (captured nucleon in **Ca41**)
    and **Austern Eq. (5.3)** exit-radius sampling **`:zr-chi-exit-mass-ratio` = M_Ca40/M_Ca41**.
  - **Angular:** coherent or m-sum per **`:angular-mode`**, **l_i=0**, **l_f=3**, **`:Ls`**.

  **vs DWUCK4 listing:** still a **reduced** model (ZR POST, **`:bind-flux`** χ outer fix from **|**H⁻|** at **k r_max** — **BIND**-style **k**/flux). Listing **asymmetry** needs full
  **T** / **S** / **β** / **P_L^m** — not this benchmark alone.

  Optional **`:cm-asymmetry-kappa`** (default **0**) is **non-DWBA** (**1 + κ cos θ_cm**) for exploratory
  plots only (e.g. rough visual match to listing shape).

  **Absolute scale (mb/sr)** — `transfer-amplitude-post` converts **u(r)=r·R(r)** to **R(r)** in the
  radial overlap. Distorted waves use **:coulomb-tail** matching at **r_max**, not unit asymptotic flux — use
  **`ca40-dp-flux-scale-to-embedded-dwuck`** when comparing to listing **Inelsig** (fm²/sr ×
  **`dwuck-inelsig-fm2-sr->mb-sr`**); **`:bind-flux`** narrows the gap vs **`:coulomb-tail`** / **`:raw`**. Spectroscopic **C²S** not folded (**S=1**).

  Reference listing values are embedded in `dwba.ca40-dwuck-benchmark-test`.

  **REPL:** `(require '[dwba.benchmark.ca40-dwuck :as c])` then `(c/ca40-dp-kinematics)` —
  not `(dwba/benchmark/...)` (slashes are division, not namespaces)."
  (:require [dwba.transfer :as t]
            [functions :as fn :refer [mass-factor-from-mu channel-sommerfeld-eta]]
            [complex :refer [mag]]))

(def dwuck-inelsig-fm2-sr->mb-sr
  "DWUCK listing column **Inelsig** is differential cross-section in **fm²/sr** (not mb/sr).
  Multiply by this to compare to Nexus / plots in **mb/sr** (1 fm² = 10 mb)."
  10.0)

(def ^:private dwuck-embedded-inelsig-fm2-sr-at-30deg
  "Raw **Inelsig** at θ_cm = 30° from myrun.lis (**fm²/sr**). Same digit string as the listing;
  do not treat as mb/sr without multiplying by `dwuck-inelsig-fm2-sr->mb-sr`."
  0.50653)

(def ^:private dwuck-embedded-sigma-mb-sr-at-30deg
  "Inelsig at 30° converted to mb/sr for `ca40-dp-flux-scale-to-embedded-dwuck`."
  (* dwuck-inelsig-fm2-sr->mb-sr dwuck-embedded-inelsig-fm2-sr-at-30deg))

;; Bound orbitals for this benchmark: deuteron s-state (l=0) → Ca41 f7/2-like (l=3).
(def ^:private ca40-dp-bound-l-i 0)
(def ^:private ca40-dp-bound-l-f 3)

(defn- ca40-dp-cm-asymmetry-factor
  "Exploratory CM factor **1 + κ cos θ** (θ radians). κ=0 → 1. Not from first-principles DWBA here."
  ^double [^double theta-rad ^double kappa]
  (if (< (Math/abs kappa) 1e-15)
    1.0
    (max 1e-300 (+ 1.0 (* kappa (Math/cos theta-rad))))))
(defn ca40-dp-kinematics
  "Returns map with :mass-factor-i :mass-factor-f :e-cm-i :e-cm-f :k-i :k-f from DWUCK listing.
  Uses `mass-factor-from-mu`: k² = (2μ/ħ²)·E with μ in MeV/c² (same convention as `functions`)."
  []
  (let [m-d 1875.613
        m-p 938.272
        m-t (* 40.0 931.494)
        m-r (* 41.0 931.494)
        mu-i (/ (* m-d m-t) (+ m-d m-t))
        mu-f (/ (* m-p m-r) (+ m-p m-r))
        mfi (mass-factor-from-mu mu-i)
        mff (mass-factor-from-mu mu-f)
        e-cm-i 13.0476
        e-cm-f 19.1886]
    {:mass-factor-i mfi :mass-factor-f mff
     :e-cm-i e-cm-i :e-cm-f e-cm-f
     :k-i (Math/sqrt (* mfi e-cm-i))
     :k-f (Math/sqrt (* mff e-cm-f))}))

(defn- ca40-dp-tail-eta ^double [^double e-cm ^double mf ^double z12ee]
  (binding [fn/mass-factor mf fn/Z1Z2ee z12ee]
    (channel-sommerfeld-eta e-cm)))

(defn optical-u-deuteron-ca40
  "Volume + imaginary Woods–Saxon from DWUCK particle-1 block (d + Ca40), R_C from listing RC."
  [L s j]
  (let [z1 1 z2 20 rc 4.7879]
    (fn [^double r]
      (t/optical-potential-woods-saxon r [97.4 3.803 0.875] [70.0 5.342 0.477]
                                       nil nil nil L s j z1 z2 rc))))

(defn optical-u-proton-ca41
  "Volume + imag + spin–orbit from DWUCK particle-2 block (p + Ca41)."
  [L s j]
  (let [z1 1 z2 20 rc 4.3103]
    (fn [^double r]
      (t/optical-potential-woods-saxon r [49.47 4.0689 0.70] [19.8 4.3172 0.75]
                                       24.2 4.0689 0.70 L s j z1 z2 rc))))

(defn- ca40-dp-deuteron-j-for-partial-wave
  "Total j for d-channel partial wave with orbital L and s=1 (L=3→j=3; else j=L+1)."
  ^double [^long L]
  (if (= L 3)
    3.0
    (+ 1.0 (double L))))

(defn ca40-dp-transfer-amplitude-for-multipole-L
  "Radial POST zero-range **T_L** for transfer multipole **L** (orbital partial wave index in χ_d, χ_p).
  Same φ_i (s), φ_f (f7/2-like l=3 well) as the L=3 benchmark; only χ_i, χ_f vary with L.

  Uses **Austern (5.3)**: exit distorted wave **χ_p** sampled at **(M_target/M_residual) r** along the
  entrance grid (Ca40/Ca41 mass ratio, same nucleon-mass convention as `ca40-dp-kinematics`)."
  [L {:keys [r-max h]
      :or {r-max 22.0 h 0.05}}]
  (let [L (long L)
        {:keys [mass-factor-i mass-factor-f e-cm-i e-cm-f k-i k-f]} (ca40-dp-kinematics)
        j-d (ca40-dp-deuteron-j-for-partial-wave L)
        j-p (+ 0.5 (double L))
        m-t (* 40.0 931.494)
        m-r (* 41.0 931.494)
        zr-ratio (t/austern-zr-chi-exit-mass-ratio m-t m-r)
        z12 (* 1.44 1.0 20.0)
        eta-i (ca40-dp-tail-eta e-cm-i mass-factor-i z12)
        eta-f (ca40-dp-tail-eta e-cm-f mass-factor-f z12)
        phi-f (t/normalize-bound-state
               (t/solve-bound-state-numerov -8.364 3 58.4538 4.0355 0.7 0.048 h r-max {:no-spin-orbit true}) h)
        phi-i (t/normalize-bound-state
               (t/solve-bound-state-numerov -2.224 0 42.0 3.9 0.65 0.048 h r-max {:no-spin-orbit true}) h)
        chi-i (t/distorted-wave-optical e-cm-i L 1.0 j-d
                                        (optical-u-deuteron-ca40 L 1.0 j-d)
                                        r-max h mass-factor-i
                                        :normalize-mode :bind-flux
                                        :bind-eta eta-i)
        chi-f (t/distorted-wave-optical e-cm-f L 0.5 j-p
                                        (optical-u-proton-ca41 L 0.5 j-p)
                                        r-max h mass-factor-f
                                        :normalize-mode :bind-flux
                                        :bind-eta eta-f)
        D0 (t/zero-range-constant :d-p)]
    (t/transfer-amplitude-post chi-i chi-f phi-i phi-f r-max h :zero-range D0
                               {:zr-chi-exit-mass-ratio zr-ratio
                                :handbook-F-from :phi-f})))


(defn ca40-dp-transfer-T-map
  "Map **{L → T_L}** for each integer L in `Ls` (default **`[3]`** — sole allowed multipole for l_i=0, l_f=3).
  Expensive: one Numerov+overlap per L."
  [& {:keys [r-max h Ls]
      :or {r-max 22.0 h 0.05 Ls [3]}}]
  (into {}
        (for [L Ls
              :let [L (long L)]]
          [L (ca40-dp-transfer-amplitude-for-multipole-L L {:r-max r-max :h h})])))

(defn ca40-dp-transfer-amplitude-L3
  "Same as `(ca40-dp-transfer-amplitude-for-multipole-L 3 opts)`."
  [opts]
  (ca40-dp-transfer-amplitude-for-multipole-L 3 opts))

(defn ca40-dp-dsigma-mb-sr
  "dσ/dΩ (**mb/sr**) at CM angle `theta-deg` for Ca40(d,p), S=1 — **DWBA only** (Coulomb in χ, ZR POST).

  **`:angular-mode`** (default **`:coherent`**): **`:coherent`** →
  **`transfer-differential-cross-section-angular-coherent`** (θ-dependent reduced multipole angular factor);
  **`:m-sum`** → **`transfer-differential-cross-section-angular-m-sum`**. For **l_i=0**, **`:m-sum`** makes the
  orbital angular factor **θ-independent** (constant dσ vs θ — not a numerical bug).

  Then **(2J_f+1)/(2J_i+1)** for **J_i=0**, **J_f=7/2** and unpolarized deuteron **1/3**.

  **Options**
  - `:angular-mode` — **`:coherent`** | **`:m-sum`** (default **`:coherent`**).
  - `:Ls` — partial waves in the T-map (default **`[3]`**). Forbidden L for (0→3) get zero weight.
  - `:cm-asymmetry-kappa` — optional exploratory **1 + κ cos θ_cm** (default **0**); not part of DWBA."
  [theta-deg & {:keys [r-max h Ls cm-asymmetry-kappa angular-mode]
                :or {Ls [3] cm-asymmetry-kappa 0.0 angular-mode :coherent}}]
  (let [r (or r-max 22.0)
        hh (or h 0.05)
        kappa (double (or cm-asymmetry-kappa 0.0))
        mode (or angular-mode :coherent)
        T-map (ca40-dp-transfer-T-map :r-max r :h hh :Ls Ls)
        {:keys [mass-factor-i mass-factor-f k-i k-f]} (ca40-dp-kinematics)
        theta-rad (* (double theta-deg) (/ Math/PI 180.0))
        base (case mode
               :m-sum (t/transfer-differential-cross-section-angular-m-sum
                       T-map 1.0 k-i k-f theta-rad
                       mass-factor-i mass-factor-f 0.0
                       ca40-dp-bound-l-i ca40-dp-bound-l-f)
               :coherent (t/transfer-differential-cross-section-angular-coherent
                          T-map 1.0 k-i k-f theta-rad
                          mass-factor-i mass-factor-f 0.0
                          ca40-dp-bound-l-i ca40-dp-bound-l-f)
               (throw (ex-info "ca40-dp-dsigma-mb-sr: :angular-mode must be :coherent or :m-sum"
                               {:angular-mode mode})))
        base-d (double (if (number? base) base (mag base)))
        spin-scale (* (t/transfer-nuclear-spin-statistical-factor 0 3.5)
                      (t/transfer-unpolarized-deuteron-spin-factor))
        asym (ca40-dp-cm-asymmetry-factor theta-rad kappa)]
    (* base-d spin-scale asym)))

(defn ca40-dp-flux-scale-to-embedded-dwuck
  "Return **σ_DWUCK(30°) / σ_Nexus(30°)**. Numerator: embedded **Inelsig** at 30° (**fm²/sr**) ×
  `dwuck-inelsig-fm2-sr->mb-sr`; denominator: `ca40-dp-dsigma-mb-sr` at **30°** (**mb/sr**), with same
  **`:angular-mode`** / **`:cm-asymmetry-kappa`** / **`:h`** / **`:r-max`** as passed. Typical **10¹–10²**
  from max-norm χ vs DWUCK flux. Does **not** fix angular shape."
  [& {:keys [h r-max cm-asymmetry-kappa angular-mode]
      :or {h 0.08 r-max 20.0 cm-asymmetry-kappa 0.0 angular-mode :coherent}}]
  (/ (double dwuck-embedded-sigma-mb-sr-at-30deg)
     (max 1e-300
          (double (ca40-dp-dsigma-mb-sr 30.0 :h h :r-max r-max :cm-asymmetry-kappa cm-asymmetry-kappa
                  :angular-mode angular-mode)))))

(defn ca40-dp-dsigma-mb-sr-dwuck-matched
  "`ca40-dp-dsigma-mb-sr` × flux scale to embedded DWUCK at 30°. For many angles, compute the scale
  once with `ca40-dp-flux-scale-to-embedded-dwuck` and multiply in a `let` (avoids recomputing T).
  Pass **`:angular-mode`** and **`:cm-asymmetry-kappa`** through to match the raw curve."
  [theta-deg & {:keys [h r-max cm-asymmetry-kappa angular-mode]
                :or {h 0.08 r-max 20.0 cm-asymmetry-kappa 0.0 angular-mode :coherent}}]
  (* (double (ca40-dp-flux-scale-to-embedded-dwuck :h h :r-max r-max :cm-asymmetry-kappa cm-asymmetry-kappa
                                                   :angular-mode angular-mode))
     (double (ca40-dp-dsigma-mb-sr theta-deg :h h :r-max r-max :cm-asymmetry-kappa cm-asymmetry-kappa
                                   :angular-mode angular-mode))))

(defn ca40-dp-angular-curve-mb-sr
  "Seq of `{:theta-deg θ :differential_cross_section_mb_sr σ}` for each angle in `theta-degrees`.
  Optional map keys are passed through to `ca40-dp-dsigma-mb-sr` (e.g. `:angular-mode`, `:Ls`, `:h`, `:r-max`, `:cm-asymmetry-kappa`)."
  [theta-degrees & {:as opts}]
  (mapv (fn [^double th]
          {:theta-deg th
           :differential_cross_section_mb_sr
           (if (seq opts)
             (apply ca40-dp-dsigma-mb-sr th (mapcat identity (seq opts)))
             (ca40-dp-dsigma-mb-sr th))})
        theta-degrees))
