(defproject nexus-nrs "0.1.0-SNAPSHOT"
  :global-vars {*warn-on-reflection* false}
  :description "Nexus-NRS: Nuclear Reaction Suite - scattering, transfer (DWBA), inelastic, halo nuclei"
  :url "https://github.com/alishsan/nexus-nrs"
  :license {:name "EPL-2.0 OR GPL-2.0-or-later WITH Classpath-exception-2.0"
            :url "https://www.eclipse.org/legal/epl-2.0/"}
  ;; **`lein run`** delegates to **web-dashboard/** (Jetty + `/api/*`). Core library has no server here.
  :main nexus-nrs.dashboard-runner
  :aliases {"dashboard" ["run" "-m" "nexus-nrs.dashboard-runner"]}
  :dependencies [[org.clojure/clojure "1.12.0"]
                 ;; `fastmath.special` is shadowed by `src/fastmath/special.clj` (vendored alpha4)
                 ;; so `Si` / `si` do not compile to colliding `special$Si` vs `special$si` on APFS.
                 [generateme/fastmath "3.0.0-alpha4" :exclusions [com.github.haifengl/smile-mkl]]
                 [incanter/incanter-core "1.9.3"]
                 [incanter/incanter-charts "1.9.3"]
                 [com.hypirion/clj-xchart "0.2.0"]
                 [org.clojure/data.json "2.4.0"]]

  :repl-options {:init-ns dwba.core})
