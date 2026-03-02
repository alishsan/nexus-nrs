;; Helper script to load and run benchmarks from REPL
;; Usage: (load-file "paper/load_benchmarks.clj")

;; Get the project root directory
(def project-root 
  (let [file-path *file*]
    (if file-path
      (.getParentFile (.getParentFile (java.io.File. file-path)))
      (java.io.File. "."))))

(def benchmark-file (java.io.File. project-root "paper" "run_benchmarks.clj"))

(if (.exists benchmark-file)
  (do
    (println (format "Loading benchmark script from: %s" (.getAbsolutePath benchmark-file)))
    (load-file (.getAbsolutePath benchmark-file))
    (println "Benchmark script loaded! Run with:")
    (println "  (paper.run-benchmarks/run-all-benchmarks)")
    (println "")
    (println "Or use the shell script:")
    (println "  ./paper/run_benchmarks.sh"))
  (do
    (println "ERROR: Could not find benchmark script!")
    (println (format "  Looking for: %s" (.getAbsolutePath benchmark-file)))
    (println (format "  Current directory: %s" (System/getProperty "user.dir")))))
