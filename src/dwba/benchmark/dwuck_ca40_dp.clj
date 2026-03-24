(ns dwba.benchmark.dwuck-ca40-dp
  "Deprecated alias — prefer `dwba.benchmark.ca40-dwuck` (clearer name, same API)."
  (:require [dwba.benchmark.ca40-dwuck :as impl]))

(def dwuck-inelsig-fm2-sr->mb-sr impl/dwuck-inelsig-fm2-sr->mb-sr)
(def ca40-dp-kinematics impl/ca40-dp-kinematics)
(def optical-u-deuteron-ca40 impl/optical-u-deuteron-ca40)
(def optical-u-proton-ca41 impl/optical-u-proton-ca41)
(def ca40-dp-transfer-amplitude-L3 impl/ca40-dp-transfer-amplitude-L3)
(def ca40-dp-transfer-amplitude-for-multipole-L impl/ca40-dp-transfer-amplitude-for-multipole-L)
(def ca40-dp-transfer-T-map impl/ca40-dp-transfer-T-map)
(def ca40-dp-dsigma-mb-sr impl/ca40-dp-dsigma-mb-sr)
(def ca40-dp-flux-scale-to-embedded-dwuck impl/ca40-dp-flux-scale-to-embedded-dwuck)
(def ca40-dp-dsigma-mb-sr-dwuck-matched impl/ca40-dp-dsigma-mb-sr-dwuck-matched)
(def ca40-dp-angular-curve-mb-sr impl/ca40-dp-angular-curve-mb-sr)
