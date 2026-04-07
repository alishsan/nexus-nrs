(ns nexus-nrs.dashboard-runner
  "Starts the web UI by running `lein run` inside web-dashboard/. Use `lein run` or `lein dashboard` from repo root."
  (:require [clojure.java.shell :as sh]))

(defn -main [& _args]
  (let [root (System/getProperty "user.dir")
        wd (java.io.File. root "web-dashboard")]
    (when-not (.isDirectory wd)
      (binding [*out* *err*]
        (println "Expected web-dashboard/ directory in" root))
      (System/exit 1))
    (println "Starting web dashboard (equivalent: cd web-dashboard && lein run)…")
    (let [{:keys [out err exit]} (sh/sh "lein" "run" :dir (.getPath wd))]
      (print out)
      (binding [*out* *err*] (print err))
      (System/exit (long (or exit 0))))))
