#!/usr/bin/env clojure
;; Benchmark Calculations Script
;;
;; This script runs all benchmark reactions and generates plots.
;; Usage: clj -M:run-benchmarks or (load-file "paper/run_benchmarks.clj")
;;
;; Output:
;; - Individual plots for each reaction in paper/plots/
;; - Summary plot with all reactions in paper/plots/summary.png
;; - Results data in paper/benchmark_results.edn

(ns paper.run-benchmarks
  (:require [dwba.transfer :as t]
            [dwba.form-factors :as ff]
            [dwba.inelastic :as inel]
            [dwba.halo-nuclei :as halo]
            [functions :refer [differential-cross-section mass-factor mass-factor-from-mu Z1Z2ee]]
            [complex :as cx :refer [mag]]
            [incanter.core :as i]
            [incanter.charts :as c]
            [clojure.java.io :as io]
            [clojure.data.json :as json]))

;; ============================================================================
;; Configuration
;; ============================================================================

(def output-dir "paper/plots")
(def results-file "paper/benchmark_results.edn")
(def experimental-data-file "paper/experimental_data.json")
(def h 0.01)  ; Step size (fm)
(def r-max 20.0)  ; Maximum radius (fm)

;; ============================================================================
;; Load Experimental Data
;; ============================================================================

(defn load-experimental-data
  "Load experimental data from JSON file"
  []
  (try
    (if (.exists (io/file experimental-data-file))
      (let [data (json/read-str (slurp experimental-data-file) :key-fn keyword)]
        (println (format "  ✓ Loaded experimental data from %s" experimental-data-file))
        data)
      (do
        (println (format "  ⚠ No experimental data file found at %s" experimental-data-file))
        {}))
    (catch Exception e
      (println (format "  ✗ Error loading experimental data: %s" (.getMessage e)))
      {})))

(def experimental-data (atom nil))

;; ============================================================================
;; Helper Functions
;; ============================================================================

(defn ensure-dir [path]
  (io/make-parents (io/file path)))

(defn save-plot [chart filename _title]
  (try
    (ensure-dir filename)
    (i/save chart filename :width 1000 :height 600)
    (println (format "  ✓ Plot saved: %s" filename))
    true
    (catch Exception e
      (println (format "  ✗ Could not save plot %s: %s" filename (.getMessage e)))
      false)))

(defn deg->rad [deg]
  (* deg (/ Math/PI 180.0)))

(defn rad->deg [rad]
  (* rad (/ 180.0 Math/PI)))

;; ============================================================================
;; Benchmark 1: 16O(p,d)15O Transfer Reaction
;; ============================================================================

(defn benchmark-16Opd []
  (println "\n=== Benchmark 1: 16O(p,d)15O Transfer Reaction ===")
  
  (let [;; Bound state parameters
        v0-i 62.0 R0-i 2.7 diff-i 0.6 r-max-bs 20.0 l-i 1
        v0-f 50.0 R0-f 1.5 diff-f 0.6 l-f 0
        Es-i -15.67  ; Neutron bound in 16O
        Es-f -2.214  ; Neutron bound in deuteron
        m-f 0.048  ; Mass factor
        
        ;; Calculate bound states
        phi-i-raw (t/solve-bound-state-numerov Es-i l-i v0-i R0-i diff-i m-f h r-max-bs)
        phi-f-raw (t/solve-bound-state-numerov Es-f l-f v0-f R0-f diff-f m-f h r-max-bs)
        phi-i (t/normalize-bound-state phi-i-raw h)
        phi-f (t/normalize-bound-state phi-f-raw h)
        overlap-norm (ff/normalized-overlap phi-i phi-f r-max-bs h)
        
        ;; Reaction parameters
        E-lab 20.0
        m-p 938.27
        m-16O 14899.0
        m-d 1876.136
        m-15O 13975.0
        mu-i (/ (* m-p m-16O) (+ m-p m-16O))
        mu-f (/ (* m-d m-15O) (+ m-d m-15O))
        mass-factor-i (/ (* 2.0 mu-i) (* 197.7 197.7))
        mass-factor-f (/ (* 2.0 mu-f) (* 197.7 197.7))
        E-CM-i (* E-lab (/ m-16O (+ m-16O m-p)))
        Q-value (+ m-p m-16O (- m-d) (- m-15O))
        E-CM-f (+ E-CM-i Q-value)
        E-lab-i E-lab
        E-lab-f (* E-CM-f (/ (+ m-d m-15O) m-15O))
        
        ;; Transfer amplitudes
        L-max 7
        D0 (t/zero-range-constant :p-d)
        T-amplitudes (into {}
                          (for [L (range (inc L-max))]
                            (let [chi-i (inel/distorted-wave-entrance E-CM-i L nil h r-max
                                                                      :projectile-type :p
                                                                      :target-A 16
                                                                      :target-Z 8
                                                                      :E-lab E-lab-i
                                                                      :s 0.5
                                                                      :j (+ L 0.5)
                                                                      :mass-factor mass-factor-i)
                                  chi-f (inel/distorted-wave-exit E-CM-i Q-value L nil h r-max
                                                                  :outgoing-type :d
                                                                  :residual-A 15
                                                                  :residual-Z 8
                                                                  :E-lab E-lab-f
                                                                  :s 1
                                                                  :j (inc L)
                                                                  :mass-factor mass-factor-f)
                                  T-L (t/transfer-amplitude-post chi-i chi-f phi-i phi-f r-max h
                                                                 :zero-range D0)]
                              [L T-L])))
        
        ;; Angular distribution
        k-i (Math/sqrt (* mass-factor-i E-CM-i))
        k-f (Math/sqrt (* mass-factor-f E-CM-f))
        S-factor 1.0
        angles-deg (range 0.0 181.0 5.0)
        angles-rad (mapv deg->rad angles-deg)
        dsigma-mb-sr (mapv (fn [theta-rad]
                             (* 10.0 (t/transfer-differential-cross-section-angular T-amplitudes S-factor k-i k-f
                                                                                    theta-rad mass-factor-i mass-factor-f 0.0 l-i l-f)))
                           angles-rad)
        
        ;; Key values
        dsigma-0 (* 10.0 (t/transfer-differential-cross-section-angular T-amplitudes S-factor k-i k-f 0.0 mass-factor-i mass-factor-f 0.0 l-i l-f))
        dsigma-70 (* 10.0 (t/transfer-differential-cross-section-angular T-amplitudes S-factor k-i k-f (deg->rad 70.0) mass-factor-i mass-factor-f 0.0 l-i l-f))
        dsigma-90 (* 10.0 (t/transfer-differential-cross-section-angular T-amplitudes S-factor k-i k-f (deg->rad 90.0) mass-factor-i mass-factor-f 0.0 l-i l-f))
        
        ;; Load experimental data
        exp-data (get @experimental-data "16O(p,d)15O")
        exp-angles (if exp-data
                     (mapv :angle_deg (:data exp-data))
                     [])
        exp-dsigma (if exp-data
                     (mapv :dsigma_mb_sr (:data exp-data))
                     [])
        exp-errors (if exp-data
                     (mapv :error_mb_sr (:data exp-data))
                     [])
        
        ;; Plot - theoretical curve
        exp-ref (if exp-data (:reference exp-data) nil)
        title-str (if exp-ref
                    (format "16O(p,d)15O — Differential Cross Section (Exp: %s)" exp-ref)
                    "16O(p,d)15O — Differential Cross Section")
        chart-base (c/xy-plot (vec angles-deg) (vec dsigma-mb-sr)
                              :title title-str
                              :x-label "θ (deg)"
                              :y-label "dσ/dΩ (mb/sr)"
                              :series-label (format "Theoretical (E_lab = %.1f MeV)" E-lab)
                              :legend true)
        
        ;; Add experimental data points
        chart (if (seq exp-angles)
                (c/add-lines chart-base (vec exp-angles) (vec exp-dsigma)
                            :series-label (if exp-ref 
                                            (format "Experimental (%s)" exp-ref)
                                            "Experimental")
                            :points true
                            :point-size 10
                            :point-type :circle
                            :color :red)
                chart-base)]
    
    (println (format "  E_lab = %.1f MeV" E-lab))
    (println (format "  dσ/dΩ(0°) = %.4e mb/sr" dsigma-0))
    (println (format "  dσ/dΩ(70°) = %.4e mb/sr" dsigma-70))
    (println (format "  dσ/dΩ(90°) = %.4e mb/sr" dsigma-90))
    (when (seq exp-angles)
      (println (format "  Experimental data points: %d" (count exp-angles)))
      (doseq [[angle dsigma error] (map vector exp-angles exp-dsigma exp-errors)]
        (println (format "    θ = %.1f°: %.2f ± %.2f mb/sr" angle dsigma error))))
    (println (format "  Normalized overlap = %.6f" overlap-norm))
    
    (save-plot chart (str output-dir "/16Opd_dcs.png") "16O(p,d)15O")
    
    {:reaction "16O(p,d)15O"
     :type "transfer"
     :E_lab E-lab
     :angles angles-deg
     :dsigma_mb_sr dsigma-mb-sr
     :key_values {:dsigma_0 dsigma-0
                  :dsigma_70 dsigma-70
                  :dsigma_90 dsigma-90
                  :overlap_norm overlap-norm}}))

;; ============================================================================
;; Benchmark 2: 11Li(p,d)10Li Transfer Reaction
;; ============================================================================

(defn benchmark-11Li-pd []
  (println "\n=== Benchmark 2: 11Li(p,d)10Li Transfer Reaction ===")
  
  (let [;; Reaction parameters from EXFOR
        E-lab 5.7  ; MeV/A, so E_lab = 5.7 MeV for proton
        E-relative 0.62  ; Relative energy between neutron and Li-9 in Li-10 (MeV)
        
        ;; Bound state parameters
        ;; Initial: neutron bound in 11Li (halo nucleus)
        ;; 11Li has a very weakly bound neutron (halo structure)
        ;; For 11Li, the neutron separation energy is very small (~0.3 MeV)
        ;; We'll use a simplified approach with appropriate parameters
        Es-i -0.3  ; Neutron separation energy from 11Li (approximate, very small for halo)
        l-i 0  ; s-wave halo
        v0-i 50.0 R0-i 2.5 diff-i 0.6 r-max-bs 30.0  ; Larger radius for halo
        
        ;; Final: neutron bound in deuteron
        Es-f -2.225  ; Deuteron binding energy
        l-f 0
        v0-f 50.0 R0-f 1.5 diff-f 0.6
        m-f 0.048  ; Mass factor
        
        ;; Calculate bound states
        phi-i-raw (t/solve-bound-state-numerov Es-i l-i v0-i R0-i diff-i m-f h r-max-bs)
        phi-f-raw (t/solve-bound-state-numerov Es-f l-f v0-f R0-f diff-f m-f h r-max-bs)
        phi-i (t/normalize-bound-state phi-i-raw h)
        phi-f (t/normalize-bound-state phi-f-raw h)
        overlap-norm (ff/normalized-overlap phi-i phi-f r-max-bs h)
        
        ;; Masses (in MeV/c²)
        m-p 938.27
        m-11Li 10252.0  ; Approximate
        m-d 1876.136
        m-10Li 9395.0  ; Approximate (11Li - neutron)
        
        ;; Reduced masses
        mu-i (/ (* m-p m-11Li) (+ m-p m-11Li))  ; Entrance: p+11Li
        mu-f (/ (* m-d m-10Li) (+ m-d m-10Li))  ; Exit: d+10Li
        
        ;; Mass factors
        mass-factor-i (/ (* 2.0 mu-i) (* 197.7 197.7))
        mass-factor-f (/ (* 2.0 mu-f) (* 197.7 197.7))
        
        ;; CM frame energies
        E-CM-i (* E-lab (/ m-11Li (+ m-11Li m-p)))
        ;; Q-value: Q = m_p + m_11Li - m_d - m_10Li
        Q-value (+ m-p m-11Li (- m-d) (- m-10Li))
        E-CM-f (+ E-CM-i Q-value)
        E-lab-i E-lab
        E-lab-f (* E-CM-f (/ (+ m-d m-10Li) m-10Li))
        
        ;; Transfer amplitudes
        L-max 7
        D0 (t/zero-range-constant :p-d)  ; (p,d) reaction
        T-amplitudes (into {}
                          (for [L (range (inc L-max))]
                            (let [chi-i (inel/distorted-wave-entrance E-CM-i L nil h r-max
                                                                      :projectile-type :p
                                                                      :target-A 11
                                                                      :target-Z 3
                                                                      :E-lab E-lab-i
                                                                      :s 0.5
                                                                      :j (+ L 0.5)
                                                                      :mass-factor mass-factor-i)
                                  chi-f (inel/distorted-wave-exit E-CM-i Q-value L nil h r-max
                                                                  :outgoing-type :d
                                                                  :residual-A 10
                                                                  :residual-Z 3
                                                                  :E-lab E-lab-f
                                                                  :s 1
                                                                  :j (inc L)
                                                                  :mass-factor mass-factor-f)
                                  T-L (t/transfer-amplitude-post chi-i chi-f phi-i phi-f r-max h
                                                                 :zero-range D0)]
                              [L T-L])))
        
        ;; Angular distribution
        k-i (Math/sqrt (* mass-factor-i E-CM-i))
        k-f (Math/sqrt (* mass-factor-f E-CM-f))
        S-factor 1.0
        angles-deg (range 0.0 181.0 5.0)
        angles-rad (mapv deg->rad angles-deg)
        dsigma-mb-sr (mapv (fn [theta-rad]
                             (* 10.0 (t/transfer-differential-cross-section-angular T-amplitudes S-factor k-i k-f
                                                                                    theta-rad mass-factor-i mass-factor-f 0.0 l-i l-f)))
                           angles-rad)
        
        ;; Key values
        dsigma-0 (* 10.0 (t/transfer-differential-cross-section-angular T-amplitudes S-factor k-i k-f 0.0 mass-factor-i mass-factor-f 0.0 l-i l-f))
        dsigma-70 (* 10.0 (t/transfer-differential-cross-section-angular T-amplitudes S-factor k-i k-f (deg->rad 70.0) mass-factor-i mass-factor-f 0.0 l-i l-f))
        dsigma-90 (* 10.0 (t/transfer-differential-cross-section-angular T-amplitudes S-factor k-i k-f (deg->rad 90.0) mass-factor-i mass-factor-f 0.0 l-i l-f))
        
        ;; Load experimental data
        exp-data (get @experimental-data "11Li(p,d)10Li")
        exp-angles (if exp-data
                     (mapv :angle_deg (:data exp-data))
                     [])
        exp-dsigma (if exp-data
                     (mapv :dsigma_mb_sr (:data exp-data))
                     [])
        exp-errors (if exp-data
                     (mapv :error_mb_sr (:data exp-data))
                     [])
        
        ;; Plot - theoretical curve
        exp-ref (if exp-data (:reference exp-data) nil)
        title-str (if exp-ref
                    (format "11Li(p,d)10Li — Differential Cross Section (Exp: %s)" exp-ref)
                    "11Li(p,d)10Li — Differential Cross Section")
        chart-base (c/xy-plot (vec angles-deg) (vec dsigma-mb-sr)
                              :title title-str
                              :x-label "θ (deg)"
                              :y-label "dσ/dΩ (mb/sr)"
                              :series-label (format "Theoretical (E_lab = %.1f MeV)" E-lab)
                              :legend true)
        
        ;; Add experimental data points
        chart (if (seq exp-angles)
                (c/add-lines chart-base (vec exp-angles) (vec exp-dsigma)
                            :series-label (if exp-ref 
                                            (format "Experimental (%s)" exp-ref)
                                            "Experimental")
                            :points true
                            :point-size 10
                            :point-type :circle
                            :color :red)
                chart-base)]
    
    (println (format "  E_lab = %.1f MeV (from EXFOR: %.1f MeV/A)" E-lab E-lab))
    (println (format "  E_relative = %.2f MeV" E-relative))
    (println (format "  dσ/dΩ(0°) = %.4e mb/sr" dsigma-0))
    (println (format "  dσ/dΩ(70°) = %.4e mb/sr" dsigma-70))
    (println (format "  dσ/dΩ(90°) = %.4e mb/sr" dsigma-90))
    (when (seq exp-angles)
      (println (format "  Experimental data points: %d" (count exp-angles)))
      (doseq [[angle dsigma error] (map vector exp-angles exp-dsigma exp-errors)]
        (println (format "    θ = %.1f°: %.2f ± %.2f mb/sr" angle dsigma error))))
    (println (format "  Normalized overlap = %.6f" overlap-norm))
    
    (save-plot chart (str output-dir "/11Li_pd_dcs.png") "11Li(p,d)10Li")
    
    {:reaction "11Li(p,d)10Li"
     :type "transfer"
     :E_lab E-lab
     :E_relative E-relative
     :angles angles-deg
     :dsigma_mb_sr dsigma-mb-sr
     :key_values {:dsigma_0 dsigma-0
                  :dsigma_70 dsigma-70
                  :dsigma_90 dsigma-90
                  :overlap_norm overlap-norm}}))

;; ============================================================================
;; Benchmark 3: 11Li(d,d) Elastic Scattering
;; ============================================================================

(defn benchmark-11Li-dd-elastic []
  (println "\n=== Benchmark 2: 11Li(d,d) Elastic Scattering ===")
  
  (let [E-lab 14.2  ; E/A = 7.1 MeV → E_d = 14.2 MeV
        V-params [50.0 2.5 0.6]
        m-d 1875.6
        m-Li11 10252.0
        mu-reduced (/ (* m-d m-Li11) (+ m-d m-Li11))
        E-CM (* E-lab (/ m-Li11 (+ m-Li11 m-d)))
        mass-f (mass-factor-from-mu mu-reduced)
        L-max 10
        
        ;; Calculate angular distribution
        angles-deg [0 15 30 45 60 75 90 105 120 135 150 165 180]
        angles-rad (mapv deg->rad angles-deg)
        dsigma-mb-sr (mapv (fn [theta]
                             (binding [mass-factor mass-f
                                       Z1Z2ee (* 3 1.44)]
                               (let [dsigma-complex (differential-cross-section E-CM V-params theta L-max)
                                     dsigma (mag dsigma-complex)]
                                 (* dsigma 10.0))))
                           angles-rad)
        
        ;; Key values
        dsigma-0 (binding [mass-factor mass-f Z1Z2ee (* 3 1.44)]
                   (* 10.0 (mag (differential-cross-section E-CM V-params 0.0 L-max))))
        dsigma-90 (binding [mass-factor mass-f Z1Z2ee (* 3 1.44)]
                    (* 10.0 (mag (differential-cross-section E-CM V-params (deg->rad 90.0) L-max))))
        dsigma-180 (binding [mass-factor mass-f Z1Z2ee (* 3 1.44)]
                     (* 10.0 (mag (differential-cross-section E-CM V-params Math/PI L-max))))
        
        ;; Experimental data (none available yet, but structure ready)
        exp-angles []
        exp-dsigma []
        
        ;; Plot - theoretical curve
        chart-base (c/xy-plot (vec angles-deg) (vec dsigma-mb-sr)
                              :title "11Li(d,d) — Elastic Scattering"
                              :x-label "θ (deg)"
                              :y-label "dσ/dΩ (mb/sr)"
                              :series-label (format "Theoretical (E_lab = %.1f MeV)" E-lab)
                              :legend true)
        
        ;; Add experimental data if available
        chart (if (seq exp-angles)
                (c/add-lines chart-base (vec exp-angles) (vec exp-dsigma)
                            :series-label "Experimental"
                            :points true
                            :point-size 8)
                chart-base)]
    
    (println (format "  E_lab = %.1f MeV (E/A = 7.1 MeV)" E-lab))
    (println (format "  E_CM = %.2f MeV" E-CM))
    (println (format "  dσ/dΩ(0°) = %.4e mb/sr" dsigma-0))
    (println (format "  dσ/dΩ(90°) = %.4e mb/sr" dsigma-90))
    (println (format "  dσ/dΩ(180°) = %.4e mb/sr" dsigma-180))
    
    (save-plot chart (str output-dir "/11Li_dd_elastic.png") "11Li(d,d)")
    
    {:reaction "11Li(d,d)"
     :type "elastic"
     :E_lab E-lab
     :angles angles-deg
     :dsigma_mb_sr dsigma-mb-sr
     :key_values {:dsigma_0 dsigma-0
                  :dsigma_90 dsigma-90
                  :dsigma_180 dsigma-180}}))

;; ============================================================================
;; Benchmark 4: 11Li(d,d') Inelastic Scattering (Monopole L=0)
;; ============================================================================

(defn benchmark-11Li-dd-inelastic []
  (println "\n=== Benchmark 3: 11Li(d,d') Inelastic Scattering (Monopole) ===")
  
  (let [E-lab 14.2
        E-ex 2.09  ; Monopole state
        lambda 0
        mu 0
        beta-0 0.10
        V-params [50.0 2.5 0.6]
        m-d 1875.6
        m-Li11 10252.0
        mu-reduced (/ (* m-d m-Li11) (+ m-d m-Li11))
        mass-f (/ (* 2.0 mu-reduced) (* 197.7 197.7))
        E-CM (* E-lab (/ m-Li11 (+ m-Li11 m-d)))
        L-i 0
        L-f 0
        
        ;; Distorted waves
        chi-i (inel/distorted-wave-entrance E-CM L-i V-params h r-max)
        chi-f (inel/distorted-wave-exit E-CM E-ex L-f V-params h r-max)
        
        ;; Inelastic amplitude
        T-inel (inel/inelastic-amplitude chi-i chi-f lambda mu beta-0 V-params r-max h)
        
        ;; Angular distribution
        k-i (Math/sqrt (* mass-f E-CM))
        E-f (- E-CM E-ex)
        k-f (Math/sqrt (* mass-f E-f))
        angles-deg [0 15 30 45 60 75 90 105 120 135 150 165 180]
        angles-rad (mapv deg->rad angles-deg)
        dsigma-mb-sr (mapv (fn [theta]
                             (let [dsigma (inel/inelastic-angular-distribution T-inel theta k-i k-f
                                                                               E-CM E-ex mass-f lambda mu)]
                               (* dsigma 10.0)))
                           angles-rad)
        
        ;; Key values
        dsigma-0 (* 10.0 (inel/inelastic-angular-distribution T-inel 0.0 k-i k-f E-CM E-ex mass-f lambda mu))
        dsigma-90 (* 10.0 (inel/inelastic-angular-distribution T-inel (deg->rad 90.0) k-i k-f E-CM E-ex mass-f lambda mu))
        
        ;; Experimental data (none available yet, but structure ready)
        exp-angles []
        exp-dsigma []
        
        ;; Plot - theoretical curve
        chart-base (c/xy-plot (vec angles-deg) (vec dsigma-mb-sr)
                              :title "11Li(d,d') — Inelastic Scattering (L=0 Monopole)"
                              :x-label "θ (deg)"
                              :y-label "dσ/dΩ (mb/sr)"
                              :series-label (format "Theoretical (E_ex = %.2f MeV)" E-ex)
                              :legend true)
        
        ;; Add experimental data if available
        chart (if (seq exp-angles)
                (c/add-lines chart-base (vec exp-angles) (vec exp-dsigma)
                            :series-label "Experimental"
                            :points true
                            :point-size 8)
                chart-base)]
    
    (println (format "  E_lab = %.1f MeV" E-lab))
    (println (format "  E_ex = %.2f MeV (monopole)" E-ex))
    (println (format "  β₀ = %.2f" beta-0))
    (println (format "  dσ/dΩ(0°) = %.4e mb/sr" dsigma-0))
    (println (format "  dσ/dΩ(90°) = %.4e mb/sr" dsigma-90))
    (println (format "  |T_inel| = %.4e" (mag T-inel)))
    
    (save-plot chart (str output-dir "/11Li_dd_inelastic.png") "11Li(d,d')")
    
    {:reaction "11Li(d,d')"
     :type "inelastic"
     :E_lab E-lab
     :E_ex E-ex
     :angles angles-deg
     :dsigma_mb_sr dsigma-mb-sr
     :key_values {:dsigma_0 dsigma-0
                  :dsigma_90 dsigma-90
                  :T_inel_mag (mag T-inel)}}))

;; ============================================================================
;; Benchmark 5: ¹²C(α,α')¹²C* (2⁺ State)
;; ============================================================================

(defn benchmark-C12-alpha-alpha-prime []
  (println "\n=== Benchmark 4: ¹²C(α,α')¹²C* (2⁺ State) ===")
  
  (let [E-lab 10.0
        E-ex 4.44  ; First 2⁺ state
        lambda 2
        mu 0
        beta-2 0.25
        V-params [50.0 2.0 0.6]
        
        ;; Simplified: use same mass factor as alpha
        mass-f 0.048  ; Approximate for alpha
        E-CM E-lab  ; Simplified
        L-i 0
        L-f 0
        
        ;; Distorted waves
        chi-i (inel/distorted-wave-entrance E-CM L-i V-params h r-max)
        chi-f (inel/distorted-wave-exit E-CM E-ex L-f V-params h r-max)
        
        ;; Inelastic amplitude
        T-inel (inel/inelastic-amplitude chi-i chi-f lambda mu beta-2 V-params r-max h)
        
        ;; Angular distribution
        k-i (Math/sqrt (* mass-f E-CM))
        E-f (- E-CM E-ex)
        k-f (Math/sqrt (* mass-f E-f))
        angles-deg [0 15 30 45 60 75 90 105 120 135 150 165 180]
        angles-rad (mapv deg->rad angles-deg)
        dsigma-mb-sr (mapv (fn [theta]
                             (let [dsigma (inel/inelastic-angular-distribution T-inel theta k-i k-f
                                                                               E-CM E-ex mass-f lambda mu)]
                               (* dsigma 10.0)))
                           angles-rad)
        
        ;; Key values
        dsigma-0 (* 10.0 (inel/inelastic-angular-distribution T-inel 0.0 k-i k-f E-CM E-ex mass-f lambda mu))
        dsigma-90 (* 10.0 (inel/inelastic-angular-distribution T-inel (deg->rad 90.0) k-i k-f E-CM E-ex mass-f lambda mu))
        
        ;; Experimental data (none available yet, but structure ready)
        exp-angles []
        exp-dsigma []
        
        ;; Plot - theoretical curve
        chart-base (c/xy-plot (vec angles-deg) (vec dsigma-mb-sr)
                              :title "¹²C(α,α')¹²C* — Inelastic Scattering (2⁺)"
                              :x-label "θ (deg)"
                              :y-label "dσ/dΩ (mb/sr)"
                              :series-label (format "Theoretical (E_ex = %.2f MeV)" E-ex)
                              :legend true)
        
        ;; Add experimental data if available
        chart (if (seq exp-angles)
                (c/add-lines chart-base (vec exp-angles) (vec exp-dsigma)
                            :series-label "Experimental"
                            :points true
                            :point-size 8)
                chart-base)]
    
    (println (format "  E_lab = %.1f MeV" E-lab))
    (println (format "  E_ex = %.2f MeV (2⁺ state)" E-ex))
    (println (format "  β₂ = %.2f" beta-2))
    (println (format "  dσ/dΩ(0°) = %.4e mb/sr" dsigma-0))
    (println (format "  dσ/dΩ(90°) = %.4e mb/sr" dsigma-90))
    
    (save-plot chart (str output-dir "/C12_alpha_alpha_prime.png") "¹²C(α,α')¹²C*")
    
    {:reaction "¹²C(α,α')¹²C*"
     :type "inelastic"
     :E_lab E-lab
     :E_ex E-ex
     :angles angles-deg
     :dsigma_mb_sr dsigma-mb-sr
     :key_values {:dsigma_0 dsigma-0
                  :dsigma_90 dsigma-90}}))

;; ============================================================================
;; Benchmark 6: Halo Nuclei Properties
;; ============================================================================

(defn benchmark-halo-nuclei []
  (println "\n=== Benchmark 5: Halo Nuclei Properties ===")
  
  (let [;; ¹¹Be
        result-11be (halo/example-11be)
        E-b-11be (:binding-energy result-11be)
        r-match-11be (:matching-radius result-11be)
        
        ;; ⁸B
        result-8b (halo/example-8b)
        E-b-8b (:binding-energy result-8b)
        r-match-8b (:matching-radius result-8b)]
    
    (println "  ¹¹Be (Neutron Halo):")
    (println (format "    E_b = %.3f MeV" E-b-11be))
    (println (format "    r_match = %.2f fm" r-match-11be))
    (println "")
    (println "  ⁸B (Proton Halo):")
    (println (format "    E_b = %.3f MeV" E-b-8b))
    (println (format "    r_match = %.2f fm" r-match-8b))
    
    {:reaction "Halo Nuclei"
     :type "properties"
     :key_values {:Be11 {:E_b E-b-11be
                         :r_match r-match-11be}
                  :B8 {:E_b E-b-8b
                       :r_match r-match-8b}}}))

;; ============================================================================
;; Summary Plot
;; ============================================================================

(defn create-summary-plot [results]
  (println "\n=== Creating Summary Plot ===")
  
  (let [;; Filter reactions with angular distributions
        reactions-with-data (filter #(contains? % :dsigma_mb_sr) results)
        
        ;; Create multi-series plot
        chart (c/xy-plot [] []
                        :title "Benchmark Reactions — Differential Cross Sections"
                        :x-label "θ (deg)"
                        :y-label "dσ/dΩ (mb/sr)"
                        :legend true)]
    
    ;; Add each reaction as a series (theoretical)
    (doseq [result reactions-with-data]
      (let [angles (:angles result)
            dsigma (:dsigma_mb_sr result)
            reaction-name (:reaction result)]
        (c/add-lines chart (vec angles) (vec dsigma) 
                     :series-label (str reaction-name " (Theoretical)"))))
    
    (save-plot chart (str output-dir "/summary.png") "Summary")
    (println "  ✓ Summary plot saved")))

;; ============================================================================
;; Main Execution
;; ============================================================================

(defn run-all-benchmarks []
  (println "=" 80)
  (println "Nexus-NRS Benchmark Calculations")
  (println "=" 80)
  
  ;; Load experimental data
  (reset! experimental-data (load-experimental-data))
  
  (ensure-dir (str output-dir "/dummy"))
  (ensure-dir results-file)
  
  (let [results (atom [])]
    
    ;; Run all benchmarks
    (try
      (swap! results conj (benchmark-16Opd))
      (catch Exception e
        (println (format "  ✗ Error in 16O(p,d)15O: %s" (.getMessage e)))))
    
    (try
      (swap! results conj (benchmark-11Li-pd))
      (catch Exception e
        (println (format "  ✗ Error in 11Li(p,d)10Li: %s" (.getMessage e)))))
    
    (try
      (swap! results conj (benchmark-11Li-dd-elastic))
      (catch Exception e
        (println (format "  ✗ Error in 11Li(d,d): %s" (.getMessage e)))))
    
    (try
      (swap! results conj (benchmark-11Li-dd-inelastic))
      (catch Exception e
        (println (format "  ✗ Error in 11Li(d,d'): %s" (.getMessage e)))))
    
    (try
      (swap! results conj (benchmark-C12-alpha-alpha-prime))
      (catch Exception e
        (println (format "  ✗ Error in ¹²C(α,α')¹²C*: %s" (.getMessage e)))))
    
    (try
      (swap! results conj (benchmark-halo-nuclei))
      (catch Exception e
        (println (format "  ✗ Error in Halo Nuclei: %s" (.getMessage e)))))
    
    ;; Create summary plot
    (try
      (create-summary-plot @results)
      (catch Exception e
        (println (format "  ✗ Error creating summary plot: %s" (.getMessage e)))))
    
    ;; Save results
    (try
      (spit results-file (pr-str @results))
      (println (format "\n✓ Results saved to %s" results-file))
      (catch Exception e
        (println (format "  ✗ Error saving results: %s" (.getMessage e)))))
    
    (println "\n" "=" 80)
    (println "Benchmark calculations complete!")
    (println (format "  - %d reactions calculated" (count @results)))
    (println (format "  - Plots saved to %s/" output-dir))
    (println (format "  - Results saved to %s" results-file))
    (println "=" 80)
    
    @results))

;; Main entry point
(defn -main [& _args]
  (run-all-benchmarks))

;; Run if loaded directly
(comment
  (run-all-benchmarks))
