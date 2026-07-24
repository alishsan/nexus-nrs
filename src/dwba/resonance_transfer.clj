(ns dwba.resonance-transfer
  "Single-neutron transfer DWBA to an UNBOUND (resonant) final state.

   Worked case: **¹¹Li(p,d)¹⁰Li*** — the ¹⁰Li **1p₁/₂** resonance at **E_r ≈ 0.62 MeV**
   above the n+⁹Li threshold (Sanetullaev *et al.*, *Phys. Lett. B* **755**, 481 (2016),
   EXFOR C2217001). ¹⁰Li is particle-unbound, so the ordinary bound-state machinery in
   **`dwba.transfer`** (`solve-bound-state-numerov` / `find-bound-state-energy` /
   `normalize-bound-state`, all of which require **E < 0**) cannot supply the neutron
   overlap wavefunction needed for the DWBA form factor.

   **Physics convention — 'quasi-bound single-particle resonance' (standard, not novel):**
   For an ℓ=1 (p-wave) resonance the centrifugal barrier plus the real nuclear
   Woods–Saxon well is a genuine potential barrier, so a positive-energy scattering
   solution at E = E_r is *metastable*: it has large interior amplitude near the nuclear
   surface and only slowly leaks outward (this is exactly what a finite width Γ means).
   The standard, well-established trick used throughout the DWBA transfer-to-unbound-states
   literature (e.g. C.M. Vincent & H.T. Fortune, Phys. Rev. C 2, 782 (1970), for treating a
   real-potential single-particle resonance as an effective bound state in a reaction form
   factor) is to:
     1. Tune the real Woods–Saxon well depth V0 so the ℓ=1 phase shift δ₁(E) passes through
        π/2 (the Breit–Wigner / unitarity-limit resonance condition) exactly at E = E_r.
     2. Truncate the E = +E_r radial solution at a **non-arbitrary** finite radius — here,
        the first zero crossing beyond the nuclear surface, i.e. the outer edge of the
        single interior 'quasi-bound lobe' — and normalize it there exactly as if it were
        a genuine bound state (`dwba.transfer/normalize-bound-state`).
     3. Use that normalized, truncated wavefunction as the ordinary bound-state overlap
        φ in the existing zero-range POST DWBA machinery
        (`dwba.transfer/transfer-amplitude-post`) — no change to the transfer-amplitude
        integral itself.

   **What this is NOT:** a genuine three-body (n+n+⁹Li) or R-matrix/Gamow continuum
   treatment of ¹⁰Li. It is a single-channel, single-particle approximation, standard for
   a first DWBA estimate of transfer to an isolated, moderately narrow resonance, and its
   limitations (see `resonance-width-wigner`, whose output should be compared to the
   measured Γ as a sanity check — not force-fit to it) are reported explicitly rather than
   hidden.

   **ℓ, j convention:** the n+⁹Li Woods–Saxon here is CENTRAL ONLY (no Thomas spin–orbit),
   same simplification already used in `examples/example_9Li_n_bound_resonances.clj`. We
   therefore resolve the ℓ=1 channel but not the j=1/2 vs j=3/2 splitting; the 'p₁/₂' label
   follows the literature assignment of the observed resonance's dominant configuration,
   not an explicit spin–orbit-split calculation here.

   **Overlap role (φ_i vs φ_f):** for (p,d) pickup, `transfer-amplitude-post`'s default
   `:handbook-F-from :phi-i` builds the zero-range form factor from the neutron's
   wavefunction *within the target nucleus* (¹¹Li), relative to the ⁹Li core — exactly the
   object this namespace builds (see `examples/example_16Opd.clj`, the existing (p,d)
   pickup precedent in this codebase, for the same `:phi-i` role and the same Austern
   (5.5)/(5.3) mass-ratio convention: M_A = target nucleus mass (the nucleus described by
   φ, here M(¹¹Li)), M_B = residual nucleus mass (M_A + m_n, here M(¹⁰Li*)))."
  (:require [functions :as fn :refer [mass-factor-from-mu solve-numerov phase-shift-from-numerov
                                       channel-sommerfeld-eta]]
            [dwba.transfer :as t]
            [complex :refer [mag mul complex-cartesian]]))

;; ============================================================================
;; Masses (MeV/c²) — M = A·u + Δ, u = 931.494 MeV, Δ from AME2016-style mass excess
;; tables (rounded), same convention/values as examples/example_11Li_pd_10Li.clj.
;; ============================================================================

(def m-p 938.2720813)
(def m-d 1875.612762)
(def m-n 939.5654205)
(def m-11Li (+ (* 11.0 931.494) 40.789))
(def m-9Li (+ (* 9.0 931.494) 12.415))

(def E-r-10Li
  "¹⁰Li 1p₁/₂ resonance centroid (MeV) above the n+⁹Li threshold (Sanetullaev et al. 2016)."
  0.62)

(def Gamma-lit-10Li
  "Literature total width (MeV) of the ¹⁰Li 1p₁/₂ resonance (Sanetullaev et al. 2016),
   for SANITY COMPARISON only — never used as a fit target in this namespace."
  0.33)

(def m-10Li-res
  "¹⁰Li* mass at the resonance centroid: m(⁹Li) + m(n) + E_r."
  (+ m-9Li m-n E-r-10Li))

;; ============================================================================
;; n+⁹Li resonance search: standard global Woods–Saxon geometry (r0=1.25 fm,
;; a0=0.65 fm), central only (no Thomas spin–orbit — see namespace docstring).
;; ============================================================================

(def n9Li-r0 1.25)
(def n9Li-a0 0.65)
(def n9Li-R0 (* n9Li-r0 (Math/pow 9.0 (/ 1.0 3.0))))
(def n9Li-l 1)

(defn n9Li-mu
  "Reduced mass (MeV/c²) of the n+⁹Li system."
  ^double []
  (/ (* m-n m-9Li) (+ m-n m-9Li)))

(defn n9Li-mass-factor
  "2μ/ħ² (MeV⁻¹·fm⁻²) for the n+⁹Li system."
  ^double []
  (mass-factor-from-mu (n9Li-mu)))

(defn resonance-phase-delta
  "ℓ=1 phase shift δ(E) (rad, principal value in (−π/2, π/2]) for the n+⁹Li central
   Woods–Saxon [V0, n9Li-R0, n9Li-a0] at energy E (MeV), via `functions/solve-numerov` +
   `functions/phase-shift-from-numerov` (same S-matrix/Hankel quotient used throughout
   `dwba.transfer` — see `examples/example_9Li_n_bound_resonances.clj` for the same
   pattern applied as an E-scan instead of a V0-scan)."
  [E l v0 h r-max r-match]
  (binding [fn/mass-factor (n9Li-mass-factor)]
    (let [u (solve-numerov E l v0 n9Li-R0 n9Li-a0 h r-max)]
      (phase-shift-from-numerov u h r-match E l))))

(defn find-resonance-v0
  "Bisect the well depth V0 (MeV) so that δ_ℓ(E-r) = π/2 exactly — the standard
   Breit–Wigner / unitarity-limit resonance condition for a real potential.

   Because `phase-shift-from-numerov` returns the PRINCIPAL value of δ in (−π/2, π/2],
   δ(V0) does not cross π/2 continuously as V0 increases through resonance: it rises
   smoothly toward π/2 and then jumps (branch cut) to near −π/2. We bracket that jump
   with a coarse V0 scan (`high-branch?` = δ > 1.0 rad, i.e. still on the pre-crossing
   rising branch) and bisect on the BOOLEAN branch indicator, which is well-defined even
   though the underlying real-valued (δ − π/2) is not continuous across the jump.

   Returns `{:converged? true :v0 ... :delta ... :l :E-r :h :r-max :r-match}` on success,
   or `{:converged? false :error \"...\"}` if no crossing is found in [v0-lo, v0-hi] — a
   genuine dead end that callers should handle rather than silently ignore."
  [& {:keys [l E-r h r-max r-match v0-lo v0-hi v0-coarse-step tol max-iter]
      :or {l n9Li-l E-r E-r-10Li h 0.02 r-max 40.0 v0-lo 5.0 v0-hi 100.0
           v0-coarse-step 1.0 tol 1.0e-5 max-iter 60}}]
  (let [r-match (double (or r-match (* 0.85 r-max)))
        delta-at (fn [v0] (resonance-phase-delta E-r l v0 h r-max r-match))
        high-branch? (fn [v0] (> (delta-at v0) 1.0))
        coarse (vec (range v0-lo v0-hi v0-coarse-step))
        bracket (first (for [i (range (dec (count coarse)))
                              :let [a (nth coarse i) b (nth coarse (inc i))]
                              :when (and (high-branch? a) (not (high-branch? b)))]
                          [a b]))]
    (if-not bracket
      {:converged? false
       :error (format "no delta=pi/2 crossing found in V0 in [%.1f,%.1f] MeV at E=%.3f MeV -- widen v0-lo/v0-hi"
                       (double v0-lo) (double v0-hi) (double E-r))}
      (let [[a0 b0] bracket
            v0-found (loop [a a0 b b0 i 0]
                       (if (or (>= i max-iter) (< (- b a) tol))
                         (/ (+ a b) 2.0)
                         (let [mid (/ (+ a b) 2.0)]
                           (if (high-branch? mid)
                             (recur mid b (inc i))
                             (recur a mid (inc i))))))]
        {:converged? true :v0 v0-found :delta (delta-at v0-found)
         :l l :E-r E-r :h h :r-max r-max :r-match r-match}))))

(defn- unwrap-toward
  "Add/subtract π from x until it lies within π/2 of ref (undoes the branch-cut jump
   in the PRINCIPAL δ so a finite-difference slope across the jump is meaningful)."
  ^double [^double x ^double ref]
  (loop [x x]
    (cond
      (> (- x ref) (/ Math/PI 2)) (recur (- x Math/PI))
      (< (- x ref) (- (/ Math/PI 2))) (recur (+ x Math/PI))
      :else x)))

(defn resonance-width-wigner
  "Wigner estimate of the total width: Γ = 2/|dδ/dE| (MeV), evaluated at fixed V0 by a
   central finite difference of δ(E) around E-r, unwrapping the principal-value branch
   cut (see `find-resonance-v0`) before differencing. Returns nil if the slope is not
   finite/positive."
  [& {:keys [l v0 E-r h r-max r-match dE]
      :or {l n9Li-l E-r E-r-10Li h 0.02 r-max 40.0 dE 0.002}}]
  (let [r-match (double (or r-match (* 0.85 r-max)))
        delta-at (fn [E] (resonance-phase-delta E l v0 h r-max r-match))
        d- (delta-at (- E-r dE))
        d+ (unwrap-toward (delta-at (+ E-r dE)) d-)
        slope (/ (- d+ d-) (* 2.0 dE))]
    (when (and (Double/isFinite slope) (pos? (Math/abs slope)))
      (/ 2.0 (Math/abs slope)))))

(defn- first-zero-crossing-index
  "First index i >= start-idx such that u changes sign between i and i+1; nil if none."
  [u start-idx]
  (loop [i (long start-idx)]
    (cond
      (>= i (dec (count u))) nil
      (let [a (double (nth u i)) b (double (nth u (inc i)))]
        (and (not (zero? a)) (not= (Math/signum a) (Math/signum b))))
      i
      :else (recur (inc i)))))

(defn quasi-bound-resonance-wavefunction
  "Build the 'quasi-bound' n+⁹Li p₁/₂ overlap wavefunction (see namespace docstring).

   Integrates the E = +E-r Numerov solution (`functions/solve-numerov`, same central WS
   as `resonance-phase-delta`) out to r-max, truncates at the FIRST zero crossing beyond
   ~1.5×n9Li-R0 (the outer edge of the single interior lobe — a non-arbitrary cutoff, not
   a free parameter), then normalizes over [0, r-cut] exactly like a genuine bound state
   (`dwba.transfer/normalize-bound-state`).

   Returns `{:u-raw :u :r-cut :node-idx :h :v0 :E-r :l}`; `:u` is the normalized,
   truncated wavefunction to pass as `phi-i` to `dwba.transfer/transfer-amplitude-post`."
  [& {:keys [l v0 E-r h r-max]
      :or {l n9Li-l E-r E-r-10Li h 0.02 r-max 40.0}}]
  (let [u-raw-full (binding [fn/mass-factor (n9Li-mass-factor)]
                      (solve-numerov E-r l v0 n9Li-R0 n9Li-a0 h r-max))
        start-idx (int (/ (* 1.5 n9Li-R0) h))
        node-idx (or (first-zero-crossing-index u-raw-full start-idx)
                      (dec (count u-raw-full)))
        r-cut (* node-idx h)
        u-raw (subvec u-raw-full 0 (inc node-idx))
        u-norm (t/normalize-bound-state u-raw h)]
    {:u-raw u-raw :u u-norm :r-cut r-cut :node-idx node-idx :h h :v0 v0 :E-r E-r :l l}))

(defn build-resonance
  "Convenience: find V0, estimate Γ, and build the quasi-bound wavefunction in one call.
   Returns `(assoc (find-resonance-v0 ...) :gamma-wigner ... :wavefunction ...)`, or the
   `{:converged? false ...}` map unchanged if no resonance was found."
  [& {:keys [l E-r h r-max] :or {l n9Li-l E-r E-r-10Li h 0.02 r-max 40.0}}]
  (let [search (find-resonance-v0 :l l :E-r E-r :h h :r-max r-max)]
    (if-not (:converged? search)
      search
      (let [v0 (:v0 search)
            gamma (resonance-width-wigner :l l :v0 v0 :E-r E-r :h h :r-max r-max)
            wf (quasi-bound-resonance-wavefunction :l l :v0 v0 :E-r E-r :h h :r-max r-max)]
        (assoc search :gamma-wigner gamma :wavefunction wf)))))

;; ============================================================================
;; p + ¹¹Li -> d + ¹⁰Li* kinematics
;; ============================================================================

(defn li11-pd-kinematics
  "Kinematics for p + ¹¹Li -> d + ¹⁰Li*(E_r), for the SAME inverse-kinematics beam
   documented in `examples/example_11Li_pd_10Li.clj` (¹¹Li at K/A = 5.7 MeV/u on a
   hydrogen target at rest), reduced to the frame-independent relative (CM) kinetic
   energy, then re-expressed as the EQUIVALENT 'light particle beam on the heavy
   partner at rest' lab energy needed only to look up CH89 / Daehnick80 global optical
   potential DEPTHS (the Numerov dynamics themselves always use the CM energy, which
   is frame-independent).

   Daehnick80 is documented (see `dwba.global-potentials`) as parameterized by deuteron
   lab energy PER NUCLEON, so `:E-lab-d-equiv-per-nucleon` divides the total equivalent
   deuteron lab energy by 2."
  [& {:keys [K-per-u] :or {K-per-u 5.7}}]
  (let [K-11Li (* 11.0 (double K-per-u))
        E-cm-i (* K-11Li (/ m-p (+ m-p m-11Li)))
        Q (- (+ m-p m-11Li) (+ m-d m-10Li-res))
        E-cm-f (+ E-cm-i Q)
        _ (when (< E-cm-f 0.1)
            (throw (ex-info "li11-pd-kinematics: E-cm-f too low -- check masses/Q"
                            {:E-cm-i E-cm-i :Q Q :E-cm-f E-cm-f})))
        mu-i (/ (* m-p m-11Li) (+ m-p m-11Li))
        mu-f (/ (* m-d m-10Li-res) (+ m-d m-10Li-res))
        mf-i (mass-factor-from-mu mu-i)
        mf-f (mass-factor-from-mu mu-f)
        E-lab-p-equiv (* E-cm-i (/ (+ m-p m-11Li) m-11Li))
        E-lab-d-equiv-total (* E-cm-f (/ (+ m-d m-10Li-res) m-10Li-res))]
    {:K-per-u K-per-u :K-11Li K-11Li
     :E-cm-i E-cm-i :E-cm-f E-cm-f :Q Q
     :mass-factor-i mf-i :mass-factor-f mf-f
     :k-i (Math/sqrt (* mf-i E-cm-i)) :k-f (Math/sqrt (* mf-f E-cm-f))
     :E-lab-p-equiv E-lab-p-equiv
     :E-lab-d-equiv-total E-lab-d-equiv-total
     :E-lab-d-equiv-per-nucleon (/ E-lab-d-equiv-total 2.0)}))

;; ============================================================================
;; Optical potentials (global CH89 / Daehnick80 — see caveat in namespace docstring:
;; both are extrapolated well below their fitted mass range, A=11/A=10 vs the fitted
;; A=40-209 / A=12-208; this is the same "use the nearest global parameterization,
;; flagged as illustrative" convention already used in
;; `dwba.benchmark.o16-dp-handbook` for A=16-17) and distorted waves.
;; ============================================================================

(def li11-Z 3) (def li11-A 11)
(def li10-Z 3) (def li10-A 10)

(defn- chi-i-potential-fn [L kin]
  (fn [r] (t/optical-potential-entrance-channel r :p li11-A li11-Z (:E-lab-p-equiv kin) L 0.5 (+ L 0.5))))

(defn- chi-f-potential-fn [L kin]
  (fn [r] (t/optical-potential-exit-channel r :d li10-A li10-Z (:E-lab-d-equiv-per-nucleon kin) L 1.0 (inc L))))

(defn distorted-wave-i
  "Entrance-channel distorted wave χ_i (reduced u = rχ) for p+¹¹Li, partial wave L,
   global CH89 optical potential, `:coulomb-tail` normalization."
  [L kin r-max h]
  (let [z12 (* 1.44 1.0 (double li11-Z))
        eta (binding [fn/mass-factor (:mass-factor-i kin) fn/Z1Z2ee z12]
              (channel-sommerfeld-eta (:E-cm-i kin)))
        rho (* (:k-i kin) r-max)]
    (t/distorted-wave-optical (:E-cm-i kin) L 0.5 (+ L 0.5) (chi-i-potential-fn L kin) r-max h (:mass-factor-i kin)
                               :normalize-mode :coulomb-tail :tail-eta eta :tail-rho rho :coulomb-init-eta eta)))

(defn distorted-wave-f
  "Exit-channel distorted wave χ_f (reduced u = rχ) for d+¹⁰Li*, partial wave L, global
   Daehnick80 optical potential, `:coulomb-tail` normalization."
  [L kin r-max h]
  (let [z12 (* 1.44 1.0 (double li10-Z))
        eta (binding [fn/mass-factor (:mass-factor-f kin) fn/Z1Z2ee z12]
              (channel-sommerfeld-eta (:E-cm-f kin)))
        rho (* (:k-f kin) r-max)]
    (t/distorted-wave-optical (:E-cm-f kin) L 1.0 (inc L) (chi-f-potential-fn L kin) r-max h (:mass-factor-f kin)
                               :normalize-mode :coulomb-tail :tail-eta eta :tail-rho rho :coulomb-init-eta eta)))

;; ============================================================================
;; Transfer amplitude / differential cross section
;; ============================================================================

(def li11-pd-l-i
  "Orbital angular momentum of the transferred neutron (p₁/₂ resonance): ℓ=1."
  n9Li-l)

(def li11-pd-l-f
  "Orbital angular momentum of the neutron within the deuteron (dominant S-state): ℓ=0.
   Used only for the triangle+parity L-selection in
   `dwba.transfer/transfer-differential-cross-section-angular` — with l_i=1, l_f=0 only
   the total transferred L=1 partial wave survives, same structure as
   `examples/example_16Opd.clj`."
  0)

(defn li11-pd-spin-factor
  "Spin-averaging factor applied to dσ/dΩ.

   We apply only the UNAMBIGUOUS unpolarized-proton entrance average, 1/(2·½+1) = 1/2
   (same factor `examples/example_16Opd.clj` applies for its (p,d) entrance channel).

   We deliberately OMIT the explicit nuclear (2J_f+1)/(2J_i+1) statistical factor used in
   `dwba.benchmark.o16-dp-handbook` (there J_i=0, J_f=1/2 are both firmly known). Here
   J(¹¹Li)=3/2⁻ is well established but J^π of the ¹⁰Li p₁/₂ resonance is NOT
   unambiguously resolved in the literature (candidate unresolved 1⁻/2⁻ doublet); rather
   than guess, we set that factor to 1 and flag it as an open normalization uncertainty
   (of order (2J_f+1)/(2J_i+1) ~ 0.5-1.3 depending on the true J_f) in the example
   report rather than silently absorbing it."
  ^double []
  (/ 1.0 (inc (* 2.0 0.5))))

(defn li11-pd-dsigma-mb-sr
  "dσ/dΩ (mb/sr, CM frame) for p+¹¹Li -> d+¹⁰Li*(E_r) at `theta-deg`, via
   `dwba.transfer/transfer-amplitude-post` (zero-range POST, D0 = `zero-range-constant
   :p-d`) summed over partial waves L=0..L-max, normalized by the Austern (5.5)/(5.3)
   prefactor `(M_B/M_A)(4π/(k_α k_β)) sqrt(2ℓ_i+1)` exactly as in
   `examples/example_16Opd.clj` (M_A = M(¹¹Li), M_B = M(¹⁰Li*) — see namespace docstring).

   `resonance` must be a map with `:u` (the quasi-bound wavefunction from
   `quasi-bound-resonance-wavefunction`/`build-resonance`) and `:h` matching the `h`
   passed here (all radial arrays in `transfer-amplitude-post` share one grid).

   `S-factor` is the spectroscopic factor (default 1.0 — NOT fit to data; see
   namespace docstring and the example script for the honest ratio-to-data report)."
  [theta-deg & {:keys [kin resonance r-max h L-max S-factor]
                :or {r-max 60.0 h 0.02 L-max 4 S-factor 1.0}}]
  (let [kin (or kin (li11-pd-kinematics))
        phi-i (:u resonance)
        _ (when (not= (double h) (double (:h resonance)))
            (throw (ex-info "li11-pd-dsigma-mb-sr: h must match resonance :h (shared radial grid)"
                            {:h h :resonance-h (:h resonance)})))
        D0 (t/zero-range-constant :p-d)
        austern-P (t/austern-radial-integral-prefactor-eq-5-5 m-11Li m-10Li-res (:k-i kin) (:k-f kin))
        amp-scale (* austern-P (Math/sqrt (inc (* 2.0 (double li11-pd-l-i)))))
        T-amplitudes (into {}
                            (for [L (range (inc (long L-max)))]
                              (let [chi-i (distorted-wave-i L kin r-max h)
                                    chi-f (distorted-wave-f L kin r-max h)
                                    T-raw (t/transfer-amplitude-post chi-i chi-f phi-i phi-i r-max h :zero-range D0)
                                    T-c (if (number? T-raw) (complex-cartesian T-raw 0.0) T-raw)]
                                [L (mul (complex-cartesian amp-scale 0.0) T-c)])))
        theta-rad (* (double theta-deg) (/ Math/PI 180.0))]
    (* (li11-pd-spin-factor)
       (t/transfer-differential-cross-section-angular T-amplitudes S-factor (:k-i kin) (:k-f kin)
                                                        theta-rad (:mass-factor-i kin) (:mass-factor-f kin)
                                                        0.0 li11-pd-l-i li11-pd-l-f))))

(def exfor-angles-deg
  "The 7 EXFOR C2217001 angles (deg) at E_r = 0.62 MeV from paper/experimental_data.json."
  [19.6 29.6 39.6 64.9 74.8 84.9 94.9])

(def exfor-dsigma-mb-sr
  "dσ/dΩ (mb/sr) at `exfor-angles-deg`, from paper/experimental_data.json
   (Sanetullaev et al. 2016, EXFOR C2217001)."
  [32.20 22.59 11.30 3.82 2.69 1.79 1.38])

(defn li11-pd-angular-comparison
  "dσ/dΩ(θ) at the 7 EXFOR angles plus the ratio to the measured values.
   Returns a vector of `{:theta-deg :dsigma-mb-sr :exp-mb-sr :ratio}`."
  [& {:keys [kin resonance r-max h L-max S-factor]
      :or {r-max 60.0 h 0.02 L-max 4 S-factor 1.0}}]
  (let [kin (or kin (li11-pd-kinematics))]
    (mapv (fn [th exp]
            (let [d (li11-pd-dsigma-mb-sr th :kin kin :resonance resonance :r-max r-max :h h
                                           :L-max L-max :S-factor S-factor)]
              {:theta-deg th :dsigma-mb-sr d :exp-mb-sr exp :ratio (/ d exp)}))
          exfor-angles-deg exfor-dsigma-mb-sr)))
