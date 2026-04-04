;; Example: **¹¹Li(p,d)¹⁰Li** one-neutron pickup — kinematics and literature context
;;
;; **¹⁰Li is unbound**; experiment sees a low-lying resonance (width Γ) rather than a true bound
;; final state. The built-in **`dwba.transfer`** zero-range post-form path uses **bound** neutron
;; overlaps in **`solve-bound-state-numerov`** (**φ_f** must be normalizable). A full comparison
;; to measured **dσ/dΩ** therefore needs a **continuum / Gamow** (or energy-averaged) treatment for
;; **¹⁰Li**; this script does **not** claim to reproduce the published DWBA curve—it collects
;; **two-body kinematics** and citations so production runs are easier to wire later.
;;
;; **Experiment (summary):** The first **¹¹Li(p,d)¹⁰Li** one-neutron transfer measurement with the
;; **IRIS** facility at **TRIUMF** used a **5.7 MeV/u** **¹¹Li** beam on a **solid H₂** target.
;; **¹⁰Li** appeared as a strong resonance at **E_r = 0.62 ± 0.04 MeV** with total width
;; **Γ = 0.33 ± 0.07 MeV**. The angular distribution is consistent with the valence neutron in the
;; **1p₁/₂** orbital. A **DWBA** analysis gives spectroscopic factor **S = 0.67 ± 0.12** for **p₁/₂**
;; removal from the **¹¹Li** ground state to the peak region.
;;
;; **Reference:** A. Sanetullaev *et al.*, *Phys. Lett. B* **755**, 481–485 (2016),
;; [doi:10.1016/j.physletb.2016.02.060](https://doi.org/10.1016/j.physletb.2016.02.060).
;;
;; **Run (repo root):** `lein run -m clojure.main examples/example_11Li_pd_10Li.clj`

(ns examples.example-11Li-pd-10Li)

(def ^:private m-p 938.2720813)
(def ^:private m-d 1875.612762)
(def ^:private m-n 939.5654205)
;; Masses M = A·u + Δ (MeV/c²), with u = 931.494 MeV and mass excess Δ from AME2016-style tables (rounded).
(def ^:private m-11Li (+ (* 11.0 931.494) 40.789))
(def ^:private m-9Li (+ (* 9.0 931.494) 12.415))
;; Residue mass at the **measured resonance centroid** above n + ⁹Li (paper E_r ≈ 0.62 MeV).
(def ^:private E-r-10Li 0.62)
(def ^:private m-10Li-res (+ m-9Li m-n E-r-10Li))

(defn Q-pd-10Li
  "Q-value (MeV) for p + ¹¹Li → d + ¹⁰Li* at resonance centroid mass."
  ^double []
  (+ m-p m-11Li (- m-d) (- m-10Li-res)))

(defn T-cm-from-K-projectile
  "Entrance-channel CM kinetic (MeV), target at rest, non-relativistic: T_cm = K·m_t/(m_p + m_t)."
  ^double [^double K-lab-projectile ^double m-projectile ^double m-target]
  (* K-lab-projectile (/ m-target (+ m-projectile m-target))))

(defn sqrt-s-nonrel
  "Invariant √s (MeV) in lab: projectile kinetic K, target at rest (non-relativistic)."
  ^double [^double m-proj ^double m-targ ^double K-proj]
  (Math/sqrt (+ (* m-proj m-proj) (* m-targ m-targ) (* 2.0 m-targ (+ m-proj K-proj)))))

(defn -main [& _args]
  (println "=== ¹¹Li(p,d)¹⁰Li — kinematics (¹⁰Li unbound; see file header) ===")
  (println)
  (let [K-per-u 5.7
        K-Li (* 11.0 K-per-u)
        Q (Q-pd-10Li)
        Tcm (T-cm-from-K-projectile K-Li m-11Li m-p)
        W (sqrt-s-nonrel m-11Li m-p K-Li)]
    (println "Paper-style inverse kinematics (heavy projectile on proton at rest):")
    (println (format "  Beam: ¹¹Li, K/A = %.2f MeV/u  →  K_lab(¹¹Li) = %.3f MeV" K-per-u K-Li))
    (println (format "  Masses used: m(¹¹Li) = %.3f, m(¹⁰Li*) ≈ m(⁹Li)+m_n+E_r = %.3f MeV/c² (E_r=%.2f MeV)"
                     m-11Li m-10Li-res E-r-10Li))
    (println (format "  Q = m_p + m(¹¹Li) - m_d - m(¹⁰Li*) = %.3f MeV" Q))
    (println (format "  T_cm (entrance, non-rel.) ≈ %.4f MeV" Tcm))
    (println (format "  √s (non-rel. lab invariant) ≈ %.4f MeV" W))
    (println)
    (println "For a bound-state **schematic** pickup tutorial (full angular distribution) see")
    (println "  `examples/example_16Opd.clj` (¹⁶O(p,d)¹⁵O).")
    (println)
    (println "Citation: Sanetullaev A. et al., Phys. Lett. B 755 (2016) 481–485.")
    nil))

(-main)
