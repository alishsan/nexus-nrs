(ns dwba.benchmark.o16-dp-radial-i-benchmark
  "Export and **diff** ¹⁶O(d,p) ZR radial integrals **I_{L_β L_α}** (Nexus handbook pipeline)
  against a tab-separated reference (e.g. values transcribed from **DWUCK4** listings).

  Reference file: lines **`L_alpha L_beta ReI ImI`** (whitespace or tab). Lines starting with **`#`** are ignored.

  Nexus exports **Re/I Im/I** from **`handbook-radial-integral-I-zr-from-neutron-bound`**, Coulomb **σ_L** on rows when
  requested, and **I_eff = I · e^{i(σ_α+σ_β)}** — the same phase attached in **`handbook-zr-multipole-amplitude-sum`**.
  Compare **I_eff** to a DWUCK line that already includes Coulomb phase in **I**, or **I** to a bare radial integral."
  (:require [clojure.java.io :as io]
            [clojure.set :as set]
            [clojure.string :as str]
            [dwba.benchmark.o16-dp-handbook :as oh]
            [dwba.transfer :as t]
            [functions :as fn]
            [complex :refer [mag mul re im complex-polar]]))

(def ^:private benchmark-opt-keys
  #{:r-max :h :L-max :e-cm-i :chi-normalize-mode :imag-scale :attach-coulomb-sigma? :no-spin-orbit})

(defn o16-dp-radial-i-rows-benchmark
  "Radial **`{:L-alpha :L-beta :I ...}`** rows for ¹⁶O(d,p); optionally **`handbook-zr-rows-with-coulomb-sigma`**.
  **`:no-spin-orbit true`** ⇒ **Vso = 0** on both distorted-wave optical potentials
  (see **`dwba.benchmark.o16-dp-handbook/o16-dp-radial-I-rows-handbook`**)."
  [& {:keys [r-max h L-max e-cm-i chi-normalize-mode imag-scale attach-coulomb-sigma? no-spin-orbit]
      :or {r-max 100.0 h 0.08 L-max 18 chi-normalize-mode :coulomb-tail
           imag-scale 1.0 attach-coulomb-sigma? true no-spin-orbit false}}]
  (let [eci (double (or e-cm-i (:e-cm-i (oh/o16-dp-kinematics))))
        {:keys [mass-factor-i mass-factor-f e-cm-f]} (oh/o16-dp-kinematics eci)
        z12 (* 1.44 1.0 8.0)
        eta-i (binding [fn/mass-factor mass-factor-i fn/Z1Z2ee z12]
                (fn/channel-sommerfeld-eta eci))
        eta-f (binding [fn/mass-factor mass-factor-f fn/Z1Z2ee z12]
                (fn/channel-sommerfeld-eta e-cm-f))
        base (oh/o16-dp-radial-I-rows-handbook :r-max r-max :h h :L-max L-max :e-cm-i eci
               :chi-normalize-mode chi-normalize-mode :imag-scale imag-scale :no-spin-orbit no-spin-orbit)]
    (vec (if attach-coulomb-sigma?
           (t/handbook-zr-rows-with-coulomb-sigma base eta-i eta-f)
           base))))

(defn o16-dp-row-I-eff
  "**I · e^{i(σ_α+σ_β)}** for a **`{:I :sigma-alpha :sigma-beta}`** row (same phase as the multipole sum)."
  [row]
  (let [Iv (:I row)
        sa (double (or (:sigma-alpha row) 0.0))
        sb (double (or (:sigma-beta row) 0.0))]
    (if (< (+ (Math/abs sa) (Math/abs sb)) 1e-20)
      Iv
      (mul Iv (complex-polar (+ sa sb) 1.0)))))

(defn o16-dp-radial-i-rows->tsv-string
  "Sorted TSV with header; includes **I_eff** for phase-matched comparison to DWUCK."
  [rows]
  (let [sb (StringBuilder.)
        sorted (sort-by (juxt :L-alpha :L-beta) rows)]
    (.append sb "# Nexus ¹⁶O(d,p) — ZR I from handbook (5.5); I_eff = I * exp(i(σ_α+σ_β)) when σ on rows\n")
    (.append sb "L_alpha\tL_beta\tRe_I\tIm_I\tsigma_alpha\tsigma_beta\tRe_I_eff\tIm_I_eff\n")
    (doseq [row sorted]
      (let [Iv (:I row)
            sig-a (double (or (:sigma-alpha row) 0.0))
            sig-b (double (or (:sigma-beta row) 0.0))
            Ie (o16-dp-row-I-eff row)]
        (.append sb (format "%d\t%d\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n"
                            (long (:L-alpha row)) (long (:L-beta row))
                            (re Iv) (im Iv) sig-a sig-b (re Ie) (im Ie)))))
    (str sb)))

(defn o16-dp-write-nexus-radial-i-tsv!
  [^String path & {:as opts}]
  (io/make-parents (io/file path))
  (let [bench-opts (select-keys opts benchmark-opt-keys)]
    (spit path (o16-dp-radial-i-rows->tsv-string (apply o16-dp-radial-i-rows-benchmark (mapcat identity bench-opts)))))
  path)

(defn parse-o16-dp-radial-i-ref-tsv
  "Parse reference text → **`{[Lα Lβ] {:re .. :im ..}}`**. Only first four columns used (**Lα Lβ Re Im**)."
  [^String text]
  (into {}
        (keep (fn [line]
                (let [t (str/trim line)]
                  (when-not (or (str/blank? t) (str/starts-with? t "#"))
                    (let [parts (str/split t #"\s+")]
                      (when-not (re-matches #"(?i)l_alpha.*" (first parts))
                        (when (< (count parts) 4)
                          (throw (ex-info "ref line needs L_alpha L_beta ReI ImI"
                                          {:line line :parts parts})))
                        (let [La (Long/parseLong (nth parts 0))
                              Lb (Long/parseLong (nth parts 1))
                              rre (Double/parseDouble (nth parts 2))
                              rim (Double/parseDouble (nth parts 3))]
                          [[La Lb] {:re rre :im rim}]))))))
              (str/split-lines text))))

(defn load-o16-dp-radial-i-ref-tsv-file
  [^String path]
  (parse-o16-dp-radial-i-ref-tsv (slurp path)))

(defn- complex-from-row-part
  [row field]
  (case field
    :I-raw (:I row)
    :I-eff (o16-dp-row-I-eff row)))

(defn report-o16-dp-radial-i-diff
  "Print comparison of **nexus** rows to **ref** map **`{[La Lb] {:re :im}}`**. Returns a summary map."
  [nexus-rows ref-map & {:keys [nexus-field print?]
                          :or {nexus-field :I-eff print? true}}]
  (when print?
    (println "=== ¹⁶O(d,p) radial I — Nexus vs reference ===")
    (println (format "  Nexus field: %s  (%s)"
                     (name nexus-field)
                     (if (= nexus-field :I-eff)
                       "compare to DWUCK I that already includes Coulomb phase"
                       "compare to bare integral (no σ on rows)")))
    (println ""))
  (let [sorted (sort-by (juxt :L-alpha :L-beta) nexus-rows)
        tol-mag 1e-9]
    (loop [rows sorted
           max-rel 0.0
           max-row nil
           ncmp 0]
      (if (empty? rows)
        (let [ref-keys (set (keys ref-map))
              nx-keys  (set (map (fn [r] [(long (:L-alpha r)) (long (:L-beta r))]) sorted))
              only-ref (set/difference ref-keys nx-keys)
              only-nx  (set/difference nx-keys ref-keys)]
          (when print?
            (when (seq only-ref)
              (println "  Only in reference (no Nexus row):")
              (doseq [k (sort only-ref)] (println "   " k)))
            (when (seq only-nx)
              (println "  Only in Nexus (no ref line):")
              (doseq [k (sort only-nx)] (println "   " k)))
            (println (format "%n  Compared pairs: %d  max rel |Δ|I||: %.4e  worst keys: %s%n"
                             ncmp max-rel (pr-str max-row))))
          {:compared-pairs ncmp
           :max-rel-mag-err max-rel
           :worst-keys max-row
           :only-in-reference (vec (sort only-ref))
           :only-in-nexus (vec (sort only-nx))})
        (let [row (first rows)
              k   [(long (:L-alpha row)) (long (:L-beta row))]
              ref (get ref-map k)]
          (if (nil? ref)
            (recur (rest rows) max-rel max-row ncmp)
            (let [Iv (complex-from-row-part row nexus-field)
                  nr (mag Iv)
                  rr (Math/hypot (:re ref) (:im ref))
                  rel (if (< rr tol-mag)
                        (if (< nr tol-mag) 0.0 Double/POSITIVE_INFINITY)
                        (/ (Math/abs (- nr rr)) rr))
                  dph (- (Math/atan2 (im Iv) (re Iv))
                         (Math/atan2 (:im ref) (:re ref)))
                  dph (if (> dph Math/PI) (- dph (* 2.0 Math/PI))
                          (if (< dph (- Math/PI)) (+ dph (* 2.0 Math/PI)) dph))]
              (when print?
                (println (format "  La=%d Lb=%d  |N|=%.6e |R|=%.6e  rel|Δmag|=%.4e  Δφ/deg=%.3f"
                                 (first k) (second k) nr rr rel (* (/ 180.0 Math/PI) dph))))
              (recur (rest rows)
                     (double (max max-rel rel))
                     (if (> rel max-rel) k max-row)
                     (inc ncmp)))))))))

(defn run-o16-dp-radial-i-diff-files!
  "Write Nexus TSV to **nexus-out**, then if **ref-path** exists load ref and **`report-o16-dp-radial-i-diff`**."
  [^String nexus-out ^String ref-path & {:as opts}]
  (let [bench-opts (select-keys opts benchmark-opt-keys)
        rep-opts (select-keys opts [:nexus-field :print?])]
    (apply o16-dp-write-nexus-radial-i-tsv! nexus-out (mapcat identity bench-opts))
    (if (.exists (io/file ref-path))
      (report-o16-dp-radial-i-diff (apply o16-dp-radial-i-rows-benchmark (mapcat identity bench-opts))
                                   (load-o16-dp-radial-i-ref-tsv-file ref-path)
                                   (merge {:nexus-field :I-eff :print? true} rep-opts))
      (do (println (format "No reference file %s — paste DWUCK I(Lα,Lβ) there and re-run." ref-path))
          nil))))
