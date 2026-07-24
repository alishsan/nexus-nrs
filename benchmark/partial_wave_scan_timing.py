"""Single-threaded wall-clock timing of a representative partial-wave scan:
complex Numerov integration of the reduced radial equation for a proton
on a Woods-Saxon + Coulomb + spin-orbit optical potential (p+16O, E_lab=20 MeV),
L = 0..30, r_max = 20 fm, h = 0.01 fm (2000 steps per partial wave).

Companion to benchmark/partial_wave_scan_timing.clj (Clojure/JVM Nexus-NRS
implementation of the same algorithm) for a single-thread wall-clock comparison.
Uses NumPy for scalar/complex arithmetic (the natural way most physicists
would prototype this in Python); the Numerov recursion itself is inherently
sequential and is not vectorized across r, in either implementation.

Run: python3 benchmark/partial_wave_scan_timing.py
"""
import time
import numpy as np

E_lab = 20.0
r_max = 20.0
h = 0.01
L_max = 30

m_p = 938.272
m_16O = 16.0 * 931.494
mu = (m_p * m_16O) / (m_p + m_16O)
hbarc = 197.327
mass_factor = 2.0 * mu / (hbarc * hbarc)

A16_THIRD = 16.0 ** (1.0 / 3.0)
R_real = 1.25 * A16_THIRD
R_imag = 1.25 * A16_THIRD
R_so = 1.01 * A16_THIRD
R_c = 1.25 * A16_THIRD
V0, a_real = 50.0, 0.65
W0, a_imag = 8.0, 0.47
Vso, a_so = 6.0, 0.75
Z1Z2e2 = 1.44 * 1 * 8


def coulomb_potential(r):
    if r <= R_c:
        return Z1Z2e2 / (2.0 * R_c) * (3.0 - (r * r) / (R_c * R_c))
    return Z1Z2e2 / r


def optical_potential(r, L, j):
    f_real = 1.0 / (1.0 + np.exp((r - R_real) / a_real))
    V_real = -V0 * f_real

    f_imag = 1.0 / (1.0 + np.exp((r - R_imag) / a_imag))
    W_imag = -W0 * f_imag

    s = 0.5
    ls = 0.5 * (j * (j + 1.0) - L * (L + 1.0) - s * (s + 1.0))
    if r == 0.0:
        V_so = 0.0
    else:
        x = (r - R_so) / a_so
        f_so = 1.0 / (1.0 + np.exp(x))
        dfso_dr = -f_so * (1.0 - f_so) / a_so
        V_so = Vso * ls * (1.0 / r) * dfso_dr

    Vc = coulomb_potential(r)
    return complex(V_real + V_so + Vc, W_imag)


def f_r_numerov_complex(r, E, L, U, mfactor):
    centrifugal = 0.0 if r == 0.0 else L * (L + 1) / (r * r)
    return complex(mfactor * (U.real - E) + centrifugal, mfactor * U.imag)


def distorted_wave_numerov(E, L, j, r_max, h, mfactor):
    steps = int(r_max / h)
    u = np.zeros(steps + 2, dtype=complex)
    u[1] = h ** (L + 1)
    rs = [i * h for i in range(steps + 2)]
    fs = [f_r_numerov_complex(r, E, L, optical_potential(r, L, j), mfactor) for r in rs]

    h2_12 = h * h / 12.0
    for n in range(1, steps - 1):
        fn_1, fn, fnp1 = fs[n - 1], fs[n], fs[n + 1]
        rhs = (2.0 + (5.0 / 6.0) * h * h * fn) * u[n] - (1.0 - h2_12 * fn_1) * u[n - 1]
        u[n + 1] = rhs / (1.0 - h2_12 * fnp1)
    return u


def partial_wave_scan():
    for L in range(L_max + 1):
        j = L + 0.5
        distorted_wave_numerov(E_lab, L, j, r_max, h, mass_factor)


def time_scan():
    t0 = time.perf_counter()
    partial_wave_scan()
    return (time.perf_counter() - t0) * 1000.0


if __name__ == "__main__":
    print("=== NumPy/CPython partial-wave scan timing ===")
    print(f"L = 0..{L_max}, r_max = {r_max:.1f} fm, h = {h:.3f} fm ({int(r_max / h)} steps/partial wave)")
    print()

    cold_ms = time_scan()
    print(f"First call:              {cold_ms:.1f} ms")

    warm_times = [time_scan() for _ in range(5)]
    best = min(warm_times)
    avg = sum(warm_times) / len(warm_times)
    print(f"5 runs (ms):              {[f'{t:.1f}' for t in warm_times]}")
    print(f"Best run:                 {best:.1f} ms")
    print(f"Mean run:                 {avg:.1f} ms")
    print(f"Per-partial-wave:         {best / (L_max + 1):.3f} ms")
    print()
    print("Done.")
