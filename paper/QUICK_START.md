# Quick Start - Running Benchmarks

## Problem: FileNotFoundException

If you get `FileNotFoundException` when trying to load the benchmark script, it's because the REPL's working directory is not the project root.

## Solutions

### Solution 1: Use Absolute Path (Easiest)

```clojure
;; Replace with your actual project path
(load-file "/Users/a.sanetullaev/Development/clojure/nexus-nrs/paper/run_benchmarks.clj")
(paper.run-benchmarks/run-all-benchmarks)
```

### Solution 2: Use Helper Script

```clojure
(load-file "paper/load_benchmarks.clj")
;; Then follow the instructions it prints
```

### Solution 3: Check and Change Directory

```clojure
;; Check current directory
(System/getProperty "user.dir")

;; If not in project root, change it (replace with your path)
(System/setProperty "user.dir" "/Users/a.sanetullaev/Development/clojure/nexus-nrs")

;; Then load
(load-file "paper/run_benchmarks.clj")
(paper.run-benchmarks/run-all-benchmarks)
```

### Solution 4: Use Shell Script (Recommended)

Just run from terminal:
```bash
cd /Users/a.sanetullaev/Development/clojure/nexus-nrs
./paper/run_benchmarks.sh
```

### Solution 5: Start REPL from Project Root

```bash
cd /Users/a.sanetullaev/Development/clojure/nexus-nrs
lein repl
# or
clj -M
```

Then in REPL:
```clojure
(load-file "paper/run_benchmarks.clj")
(paper.run-benchmarks/run-all-benchmarks)
```

## Recommended Approach

**Use the shell script** - it's the simplest and handles all path issues automatically:

```bash
./paper/run_benchmarks.sh
```
