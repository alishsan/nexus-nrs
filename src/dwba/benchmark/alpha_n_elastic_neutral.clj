(ns dwba.benchmark.alpha-n-elastic-neutral
  "Neutral **α + n** elastic sandbox: **no Coulomb** (**Z₁Z₂ = 0**), central Woods–Saxon only.

  **Purpose:** Check that **`distorted-wave-optical`** (**Numerov** for **u = rχ**) plus
  **`distorted-wave-numerov-R-for-smatrix`** and **`distorted-wave-coulomb-S-from-numerov-R`** with **η = 0**
  reproduces **`functions/s-matrix`** (same **Hankel±** ratio as elastic partial waves).

  **Note:** **S** from **R = u/(a u′)** is **independent** of any global scale on **u**, so **`:normalize-mode`**
  does not affect **`alpha-n-S-from-numerov-for-L`**."

  (:require [dwba.transfer :as t]
            [functions :as fn :refer [mass-factor-from-mu]]))

(def ^:private m-alpha 3727.379)
(def ^:private m-neutron 939.565)

(defn alpha-n-reduced-mass
  ^double []
  (/ (* m-alpha m-neutron) (+ m-alpha m-neutron)))

(defn alpha-n-mass-factor
  ^double []
  (mass-factor-from-mu (alpha-n-reduced-mass)))

(defn alpha-n-z1z2ee
  "Coulomb strength **Z₁Z₂ e²** (MeV·fm); **0** for n + α."
  ^double []
  0.0)

(defn alpha-n-matching-radius-fm
  "**a = 2(R + a₀)** — same construction as **`functions/s-matrix-3-impl`** for **V = [V₀ R a₀]**."
  ^double [[_v0 ^double r0 ^double a0]]
  (* 2.0 (+ r0 a0)))

(defn alpha-n-optical-fn-central
  "Central WS + explicit **nil** Coulomb (no **l·s**). **j = l + ½** for chosen **l**."
  [V-params ^long L]
  (let [l (double L)
        s 0.5
        j (+ s l)]
    (fn [^double r]
      (t/optical-potential-woods-saxon r V-params nil nil nil nil l s j nil nil nil))))

(defn alpha-n-S-from-numerov-for-L
  "**S_L** from **`distorted-wave-optical`** grid → **R** → **`distorted-wave-coulomb-S-from-numerov-R`** (**η = 0**)."
  [e-cm V-params L r-max h]
  (let [e-cm (double e-cm)
        L (long L)
        r-max (double r-max)
        h (double h)
        mf (alpha-n-mass-factor)
        k (Math/sqrt (* mf e-cm))
        a (double (alpha-n-matching-radius-fm V-params))
        rho (* k a)
        U (alpha-n-optical-fn-central V-params L)
        u-vec (t/distorted-wave-optical e-cm L 0.5 (+ 0.5 (double L)) U r-max h mf
                                        :normalize-mode :raw)
        Rmat (t/distorted-wave-numerov-R-for-smatrix u-vec h a)]
    (t/distorted-wave-coulomb-S-from-numerov-R Rmat L 0.0 rho)))

(defn alpha-n-sigma-L-ref
  "Reference **S** from **`fn/s-matrix`** with **Z1Z2ee = 0** and α–n **mass-factor**."
  [e-cm V-params L]
  (binding [fn/mass-factor (alpha-n-mass-factor)
            fn/Z1Z2ee (alpha-n-z1z2ee)]
    (fn/s-matrix (double e-cm) V-params (long L))))
