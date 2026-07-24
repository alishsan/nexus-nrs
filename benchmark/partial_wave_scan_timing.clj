;; Single-threaded wall-clock timing of a representative partial-wave scan:
;; complex Numerov integration of the reduced radial equation for a proton
;; on a Woods-Saxon + Coulomb + spin-orbit optical potential (p+16O, E_lab=20 MeV),
;; L = 0..30, r_max = 20 fm, h = 0.01 fm (2000 steps per partial wave).
;;
;; Companion to benchmark/partial_wave_scan_timing.py (equivalent NumPy/Python
;; implementation of the same algorithm) for a Clojure/JVM vs NumPy/CPython
;; single-thread wall-clock comparison. Run from project root:
;;   lein run -m clojure.main benchmark/partial_wave_scan_timing.clj

(require '[dwba.transfer :as t])

(def E-lab 20.0)
(def r-max 20.0)
(def h 0.01)
(def L-max 30)
(def mass-factor (/ (* 2.0 (/ (* 938.272 (* 16.0 931.494)) (+ 938.272 (* 16.0 931.494)))) (* 197.327 197.327)))

(defn optical-fn [L s j]
  (fn [^double r]
    (t/optical-potential-woods-saxon r
      [50.0 (* 1.25 (Math/pow 16.0 (/ 1.0 3.0))) 0.65]     ; real WS
      [8.0  (* 1.25 (Math/pow 16.0 (/ 1.0 3.0))) 0.47]     ; volume imaginary
      nil
      6.0 (* 1.01 (Math/pow 16.0 (/ 1.0 3.0))) 0.75         ; spin-orbit
      L 0.5 j 1 8 (* 1.25 (Math/pow 16.0 (/ 1.0 3.0))))))

(defn partial-wave-scan []
  (doseq [L (range (inc L-max))]
    (let [j (+ L 0.5)]
      (t/distorted-wave-optical E-lab L 0.5 j (optical-fn L 0.5 j) r-max h mass-factor))))

(defn time-scan []
  (let [t0 (System/nanoTime)]
    (partial-wave-scan)
    (/ (- (System/nanoTime) t0) 1.0e6)))

(println "=== Nexus-NRS (Clojure/JVM) partial-wave scan timing ===")
(println (format "L = 0..%d, r_max = %.1f fm, h = %.3f fm (%d steps/partial wave)"
                  L-max r-max h (int (/ r-max h))))
(println "")

(let [cold-ms (time-scan)]
  (println (format "First call (JIT cold):  %.1f ms" cold-ms)))

(let [warm-times (mapv (fn [_] (time-scan)) (range 5))
      best (apply min warm-times)
      avg (/ (reduce + warm-times) (count warm-times))]
  (println (format "5 warm runs (ms):        %s" (mapv #(format "%.1f" %) warm-times)))
  (println (format "Best warm run:           %.1f ms" best))
  (println (format "Mean warm run:           %.1f ms" avg))
  (println (format "Per-partial-wave (warm): %.3f ms" (/ best (inc L-max)))))

(println "")
(println "Done.")
