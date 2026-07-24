;; Export ¹⁶O(d,p) ZR radial **I_{L_β L_α}** (Nexus) to TSV and diff vs a DWUCK reference file if present.
;;
;;   lein run -m clojure.main examples/o16_dp_radial_I_diff.clj
;;
;; Writes **`output/o16_dp_radial_I_nexus.tsv`**. If **`output/o16_dp_radial_I_dwuck_ref.tsv`** exists
;; (copy from **`benchmark/o16_dp_radial_i_dwuck_ref.example.tsv`** and replace numbers from DWUCK),
;; prints relative |I| error and Δφ per (Lα,Lβ).

(require '[dwba.benchmark.o16-dp-radial-i-benchmark :as rib])

(def nexus-out "output/o16_dp_radial_I_nexus.tsv")
(def ref-path "output/o16_dp_radial_I_dwuck_ref.tsv")

;; Match grid to your DWUCK run / example_16Odp.clj as needed.
(rib/run-o16-dp-radial-i-diff-files! nexus-out ref-path
  :h 0.08 :r-max 100.0 :L-max 18
  :chi-normalize-mode :coulomb-tail
  :imag-scale 1.0
  :attach-coulomb-sigma? true
  :nexus-field :I-eff)

(println (format "Nexus table written: %s" nexus-out))
