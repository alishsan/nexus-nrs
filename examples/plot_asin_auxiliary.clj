;; Unrelated diagnostic: y = π + 2x_rad − 4 asin((sin x_rad)/n);
;; x and y on the chart are in degrees (x: 0–90°).
;; Run: lein run -m clojure.main examples/plot_asin_auxiliary.clj

(import '[java.awt Color])

(require '[clojure.java.io :as io]
         '[incanter.core :as i]
         '[incanter.charts :as c])

(def ^:private n-red 1.331)
(def ^:private n-blue 1.340)
(def ^:private rad->deg (/ 180.0 Math/PI))

(defn- y-rad [^double n ^double x-rad]
  (+ Math/PI (* 2.0 x-rad)
     (* -4.0 (Math/asin (/ (Math/sin x-rad) n)))))

(defn- set-line-colors!
  "Incanter `add-lines` often uses a second dataset (series 0 on renderer 1).
  If both series share one dataset, set series 0 and 1 on renderer 0."
  [chart colors]
  (let [plot (.getPlot chart)
        r0   (when (pos? (.getDatasetCount plot)) (.getRenderer plot 0))
        ds0  (when (pos? (.getDatasetCount plot)) (.getDataset plot 0))]
    (if (and ds0 r0 (>= (.getSeriesCount ds0) 2))
      (dotimes [i (min (count colors) (.getSeriesCount ds0))]
        (.setSeriesPaint r0 i (nth colors i)))
      (dotimes [i (min (count colors) (.getDatasetCount plot))]
        (when-let [r (.getRenderer plot i)]
          (.setSeriesPaint r 0 (nth colors i)))))
    chart))

(defn- chart-set-y-range!
  [chart ^double lo ^double hi]
  (let [plot   (.getPlot chart)
        y-axis (.getRangeAxis plot)]
    (.setRange y-axis lo hi))
  chart)

(let [x-deg (mapv (fn [^long i] (* (/ i 500.0) 90.0)) (range 501))
      x-rad (mapv #(* (double %) (/ Math/PI 180.0)) x-deg)
      ys-red-deg  (mapv #(* rad->deg (y-rad n-red %)) x-rad)
      ys-blue-deg (mapv #(* rad->deg (y-rad n-blue %)) x-rad)
      out         "output/asin_auxiliary_y_vs_x.png"]
  (io/make-parents (io/file out))
  (-> (c/xy-plot x-deg ys-red-deg
                 :title "y = π + 2x − 4 asin((sin x)/n) — x and y in degrees"
                 :x-label "x (deg)"
                 :y-label "y (deg)"
                 :series-label (format "n = %.3f" n-red)
                 :legend true)
      (c/add-lines x-deg ys-blue-deg :series-label (format "n = %.3f" n-blue))
      (set-line-colors! [Color/red Color/blue])
      (chart-set-y-range! 130.0 180.0)
      (i/save out :width 800 :height 500))
  (println "Saved:" out))
