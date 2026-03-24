(ns dwba.ca40-dwuck-benchmark-test
  "Benchmark Ca40(d,p)Ca41 g.s. against DWUCK4 listing `myrun.lis` (first reaction block).

  Embedded **Inelsig** column from **myrun.lis** is stored here as **fm²/sr** (listing digits).
  Plots and flux scaling use `ca40/dwuck-inelsig-fm2-sr->mb-sr` (×10) to match `ca40-dp-dsigma-mb-sr`.
  Default **`ca40-dp-dsigma-mb-sr`** uses **`:angular-mode :coherent`** (θ-dependent). **`:m-sum`** gives
  **l_i=0** θ-flat σ — see `ca40-dp-angular-m-sum-isotropic-for-l0`."
  (:require [clojure.test :refer :all]
            [dwba.benchmark.ca40-dwuck :as ca40]))

;; Sparse angles from the same listing as `examples/plot_ca40_dp_dwuck` — **Inelsig** in **fm²/sr**.
;; (Multiply by `ca40/dwuck-inelsig-fm2-sr->mb-sr` for mb/sr.)
(def ^:private dwuck-transfer-reference-inelsig-fm2-sr
  [[0.0 5.4232e-2] [10.0 1.1613e-1] [20.0 3.0309e-1] [30.0 5.0653e-1]
   [40.0 4.3680e-1] [50.0 2.4526e-1] [60.0 1.7466e-1] [70.0 1.9313e-1]
   [80.0 2.0431e-1] [90.0 1.7214e-1] [100.0 1.2531e-1] [120.0 9.7973e-2]
   [150.0 9.5269e-2] [180.0 7.2915e-2]])

(deftest ca40-dp-dwuck-angular-sanity-benchmark
  "Default **:coherent** angular mode: θ-dependent dσ/dΩ (L=3 **Y_{30}**-type factor), not constant vs θ."
  (let [angles (mapv first dwuck-transfer-reference-inelsig-fm2-sr)
        ours (mapv (fn [^double th] [th (ca40/ca40-dp-dsigma-mb-sr th :h 0.08 :r-max 20.0)])
                   angles)
        s30 (ca40/ca40-dp-dsigma-mb-sr 30.0 :h 0.08 :r-max 20.0)
        s90 (ca40/ca40-dp-dsigma-mb-sr 90.0 :h 0.08 :r-max 20.0)
        sigmas (mapv second ours)]
    (is (every? #(Double/isFinite (double %)) sigmas))
    (is (every? #(>= % 0.0) sigmas))
    (is (pos? s30))
    ;; Coherent **L=3**: **P_3(0)=0** at **90°** → σ(90°) is tiny; dσ is **not** θ-flat (mb/sr magnitudes ≪1).
    (is (< (double s90) (* 1e-10 (double s30)))
        "coherent L=3: σ(90°) ≪ σ(30°)")
    (is (> (double (reduce max sigmas)) (* 1e6 (double (reduce min sigmas))))
        "σ min (near 90°) ≪ σ max across sampled angles")))

(deftest ca40-dp-L-selection-ignores-forbidden-multipoles
  "l_i=0, l_f=3: only odd L in [0..4] contribute; padding :Ls with forbidden L must not change σ."
  (let [h 0.08 rmax 20.0 th 25.0
        s-default (ca40/ca40-dp-dsigma-mb-sr th :h h :r-max rmax)
        s-padded (ca40/ca40-dp-dsigma-mb-sr th :h h :r-max rmax :Ls [0 1 2 3 4])]
    (is (< (Math/abs (- (double s-default) (double s-padded)))
           (* 1e-9 (max 1.0 (Math/abs (double s-default)))))
        "Default [3] same as [0 1 2 3 4] with DWBA L weights")))

(deftest ca40-dp-explicit-Ls-matches-default
  (let [h 0.08 rmax 20.0 th 40.0
        s-def (ca40/ca40-dp-dsigma-mb-sr th :h h :r-max rmax)
        s3 (ca40/ca40-dp-dsigma-mb-sr th :h h :r-max rmax :Ls [3])]
    (is (pos? (double s-def)))
    (is (< (Math/abs (- (double s-def) (double s3))) (* 1e-9 (max 1.0 (Math/abs (double s-def))))))))

(deftest ca40-dp-angular-m-sum-isotropic-for-l0
  "**:angular-mode :m-sum** with **l_i=0**: orbital factor has no θ — σ flat (κ=0)."
  (let [h 0.08 rmax 20.0
        opts [:h h :r-max rmax :angular-mode :m-sum]
        s30 (apply ca40/ca40-dp-dsigma-mb-sr 30.0 opts)
        s90 (apply ca40/ca40-dp-dsigma-mb-sr 90.0 opts)
        s150 (apply ca40/ca40-dp-dsigma-mb-sr 150.0 opts)]
    (is (< (Math/abs (- (double s30) (double s90)))
           (* 1e-9 (max 1.0 (double s30))))
        "m-sum l_i=0: σ(30°) = σ(90°)")
    (is (< (Math/abs (- (double s30) (double s150)))
           (* 1e-9 (max 1.0 (double s30) (double s150))))
        "σ(30°) = σ(150°) up to round-off")))

(deftest ca40-dsigma-symmetric-cm-default-coherent-l3
  "Default coherent **L=3**: σ(θ)=σ(π−θ) at 30° and 150° (same |P_3(cos θ)|)."
  (let [h 0.08 rmax 20.0
        s30 (ca40/ca40-dp-dsigma-mb-sr 30.0 :h h :r-max rmax)
        s150 (ca40/ca40-dp-dsigma-mb-sr 150.0 :h h :r-max rmax)]
    (is (< (Math/abs (- (double s30) (double s150)))
           (* 1e-9 (max 1.0 (double s30) (double s150))))
        "σ(30°) = σ(150°) for symmetric |P_L|²")))

(deftest ca40-dp-kinematics-match-listing
  (let [{:keys [k-i k-f e-cm-i e-cm-f]} (ca40/ca40-dp-kinematics)]
    (is (< (Math/abs (- k-i 1.0904)) 0.02) "k_i ≈ 1.0904 fm⁻¹ from DWUCK listing")
    (is (< (Math/abs (- k-f 0.9467)) 0.02) "k_f ≈ 0.9467 fm⁻¹ from DWUCK listing")
    (is (< (Math/abs (- e-cm-i 13.0476)) 0.01))
    (is (< (Math/abs (- e-cm-f 19.1886)) 0.01))))
