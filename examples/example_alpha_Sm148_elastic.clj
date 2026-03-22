;; Elastic check: α + ¹⁴⁸Sm at E_lab(α) = 50 MeV
;;
;; Prints S-matrix for L = 0 … 50 and sample dσ/dΩ (partial-wave sum 0…50) using
;; functions.clj (Coulomb + real Woods–Saxon + imaginary Woods–Saxon absorption), same
;; convention as the web dashboard elastic API (*elastic-imag-ws-params*).
;;
;; Woods–Saxon (user-specified):
;;   Real:    V₀ = 65 MeV,  R = 7.5 fm, a = 0.67 fm
;;   Imag:   Wᵢ = 30 MeV,  Rᵢ = 7.5 fm, aᵢ = 0.67 fm
;;
;; Run from repository root:
;;   lein run -m clojure.main examples/example_alpha_Sm148_elastic.clj

;; :reload so REPL picks up src/functions.clj changes (avoids stale dσ logic + s-matrix-3-memo confusion).
(require '[functions :refer [differential-cross-section mass-factor mass-factor-from-mu
                             phase-shift s-matrix Z1Z2ee *elastic-imag-ws-params*]]
         :reload)
(require '[complex :as cpx])

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

;; Real Woods–Saxon: V₀ (MeV), R₀ (fm), a (fm)
(def V-params [65.0 7.5 0.67])
;; Imaginary WS: W₀ (MeV), R_W (fm), a_W (fm) — absorption (bound to *elastic-imag-ws-params*)
(def imag-ws-params [30.0 7.5 0.67])

(def L-max 50)

(defn rutherford-dsigma-fm2
  "Point-charge Rutherford dσ/dΩ (fm²/sr), non-relativistic CM. Z1Z2e² in MeV·fm, E_cm in MeV, θ rad."
  [^double z1z2e2 ^double e-cm ^double theta-rad]
  (let [s (Math/sin (* 0.5 theta-rad))]
    (if (< (Math/abs s) 1e-15)
      Double/NaN
      (/ (Math/pow (/ z1z2e2 (* 4.0 e-cm)) 2.0)
         (Math/pow s 4.0)))))

(println "=== α + ¹⁴⁸Sm elastic (check script) ===")
(println "")
(println (format "  Projectile: α  (m = %.3f MeV/c², Z = 2)" m-alpha))
(println (format "  Target:     ¹⁴⁸Sm (m ≈ %.1f MeV/c², Z = %d)" m-sm Z-sm))
(println (format "  E_lab(α)    = %.2f MeV" E-lab))
(println (format "  E_CM        = %.4f MeV" E-CM))
(println (format "  μ (reduced) = %.4f MeV/c²" mu-reduced))
(println (format "  Z₁Z₂ e²     = %.4f MeV·fm" z1z2ee))
(println (format "  Real WS:    [V₀ R₀ a]   = [%.1f %.2f %.2f] MeV, fm, fm"
                 (first V-params) (second V-params) (nth V-params 2)))
(println (format "  Imag WS:    [Wᵢ Rᵢ aᵢ] = [%.1f %.2f %.2f] MeV, fm, fm"
                 (first imag-ws-params) (second imag-ws-params) (nth imag-ws-params 2)))
(println "")
(println "=== Point-charge Rutherford (Coulomb only, non-relativistic CM) — verify the numbers ===")
(println "")
(let [pref (/ z1z2ee (* 4.0 E-CM))
      pref2 (* pref pref)]
  (println "  dσ/dΩ|_Ruth = (Z₁Z₂e² / (4 E_cm))² / sin⁴(θ/2)")
  (println "  Units: Z₁Z₂e² in MeV·fm, E_cm in MeV  →  (…)² has dimensions fm²  →  dσ/dΩ in fm²/sr.")
  (println "  Conversion printed elsewhere: 1 fm² = 10 mb  →  multiply fm²/sr by 10 for mb/sr.")
  (println "")
  (println (format "  Z₁Z₂e² / (4 E_cm) = %.4f / (4 × %.4f) = %.6f fm" z1z2ee E-CM pref))
  (println (format "  (Z₁Z₂e² / (4 E_cm))² = %.6f fm²  (angle-independent prefactor)" pref2))
  (println "")
  (println "  Forward peaking comes entirely from 1/sin⁴(θ/2): sin(5°) ≈ 0.0872, sin⁴(5°) ≈ 5.77×10⁻⁵,")
  (println "  so σ_Ruth(10°) is ~1/sin⁴ larger than at θ where sin⁴ ~ 𝒪(1) (e.g. near 90°).")
  (println "")
  (doseq [deg [10 90]]
    (let [th (* deg (/ Math/PI 180.0))
          s (Math/sin (* 0.5 th))
          s4 (Math/pow s 4)
          sig (/ pref2 s4)]
      (println (format "  θ_cm = %3d°: sin(θ/2) = %.6f, sin⁴(θ/2) = %.5e  →  σ_Ruth = %10.2f fm²/sr = %12.1f mb/sr"
                       deg s s4 sig (* 10 sig))))))
(println "")
(println "  What the dσ table below is NOT:")
(println "  • functions/differential-cross-section is |Σ_L f_L|² with f_L ∝ [(2L+1)/k] P_L(cos θ) (S_L − 1).")
(println "    No separate Coulomb (Rutherford) amplitude is added to that sum.")
(println "  • So dσ_code / σ_Ruth is usually ≪ 1 (especially forward): numerator ≠ full elastic |f_C + f_N|².")
(println "    The ratio is useful as a scale check on σ_Ruth and on how big the (S−1) partial-wave piece is,")
(println "    not as “data should go to 1 at forward angles” unless your S-matrix convention embeds all Coulomb physics.")
(println "")

(binding [mass-factor (mass-factor-from-mu mu-reduced)
          Z1Z2ee z1z2ee
          *elastic-imag-ws-params* imag-ws-params]
  (println "=== S-matrix and phase shifts (L = 0 … L_max) ===")
  (println "")
  (println "   L  |    Re(S)   |    Im(S)   |  |S|   | arg(S) rad |  δ_L (deg)  (δ = arg(S)/2)")
  (println " -----+------------+------------+-------+-------------+---------------------------")
  (doseq [L (range (inc L-max))]
    (try
      (let [S-val (s-matrix E-CM V-params L)
            reS (double (cpx/re S-val))
            imS (double (cpx/im S-val))]
        (if (or (Double/isNaN reS) (Double/isNaN imS)
                (Double/isInfinite reS) (Double/isInfinite imS))
          (println (format " %3d | (non-finite S — Coulomb/Hankel U numerics at this L, η, ρ)" L))
          (let [d (phase-shift E-CM V-params L)
                Sm (cpx/mag S-val)]
            (println (format " %3d | %10.5f | %10.5f | %.5f | %11.5f | %11.5f"
                               L reS imS Sm (cpx/arg S-val) (* d (/ 180.0 Math/PI)))))))
      (catch Exception e
        (println (format " %3d | ERROR: %s" L (.getMessage e))))))
  (println "")
  (println "=== Elastic dσ/dΩ (|f|²), sum over L = 0 … " L-max " ===")
  (println "(Partial waves with non-finite S are omitted in the sum — see differential-cross-section doc.)")
  (println "σ_Ruth: same point-Coulomb formula as the web dashboard (e² = 1.44 MeV·fm).")
  (println "")
  (println "  θ_cm (deg) | dσ/dΩ / σ_Ruth | dσ/dΩ (mb/sr)")
  (println " ------------+------------------+---------------")
  (doseq [theta-deg [10 20 30 45 60 90 120 150 170]]
    (try
      (let [theta (* theta-deg (/ Math/PI 180.0))
            ds-c (differential-cross-section E-CM V-params theta L-max)
            ;; |f|² returned as near-real complex; Cartesian magnitude avoids refer/shadow bugs.
            r (double (cpx/re ds-c))
            i (double (cpx/im ds-c))
            ds-fm2 (Math/sqrt (+ (* r r) (* i i)))
            ruth (rutherford-dsigma-fm2 z1z2ee E-CM theta)
            ratio (if (and (Double/isFinite ruth) (> ruth 0.0))
                    (/ ds-fm2 ruth)
                    Double/NaN)
            ds-mb (* ds-fm2 10.0)]
        (if (Double/isFinite ratio)
          (println (format "   %6.1f    |     %12.5f | %12.5e"
                           (double theta-deg) ratio ds-mb))
          (println (format "   %6.1f    |     (n/a θ→0°) | %12.5e"
                           (double theta-deg) ds-mb))))
      (catch Exception e
        (println (format "   %6.1f    | ERROR: %s" theta-deg (.getMessage e))))))
  (println "")
  (println "Note: numeric L_max ⇒ sum L = 0…L_max; a collection sums only those L.")
  (println "Done."))
