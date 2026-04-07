(ns nexus-nrs.optical-potentials
  "Load **`resources/nexus_nrs/standard_optical_potentials.json`** — schematic / tutorial OM and bound-state wells.

  Parsed map keys are keywords (`:entries`, `:schema_version`, …). Each entry has `:ws_vector_real` ready for `functions/s-matrix` style **[V0 R0 a0]**.")

(require '[clojure.java.io :as io]
         '[clojure.data.json :as json])

(defn load-standard-optical-potentials
  "Return the JSON root as a Clojure map (`:key-fn keyword`)."
  []
  (with-open [r (io/reader (io/resource "nexus_nrs/standard_optical_potentials.json"))]
    (json/read r :key-fn keyword)))

(defn entry-by-id
  "First entry whose **:id** equals **id-kw** or string; else nil."
  [id]
  (let [id-str (name (if (keyword? id) id (symbol (str id))))
        xs (:entries (load-standard-optical-potentials))]
    (some (fn [e] (when (= id-str (name (:id e))) e)) xs)))
