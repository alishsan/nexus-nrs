#!/usr/bin/env python3
"""Parse DWUCK4 O16(d,p)O17 RADIAL MATRIX ELEMENTS listing -> complex I(La,Lb) TSV.

Confirmed format (reverse-engineered + validated against output/o16_dp_radial_I_dwuck_ref.tsv):
  Block header "0 RADIAL MATRIX ELEMENTS,  J1=L1+ 2/2,  J2=L2+ 1/2" (aligned/stretched
  j_d=La+1, j_p=Lb+1/2, matching Nexus's deuteron-j-for-partial-wave / proton-j-for-partial-wave).
  Each data row: "L2  <6 floats>" = 3 complex (Re,Im) pairs for L1 = |L2-LTR|, |L2-LTR|+2, |L2-LTR|+4
  (LTR = transferred ell = 2 here). I = F_complex * 4*pi/(k_alpha * k_beta).
"""
import re, sys

def parse_derived_k(text, occurrence):
    ks = re.findall(r'K\s*=\s*([\-0-9.]+)\s+ETA', text)
    return float(ks[occurrence])

def parse_block(text, ltr, header="RADIAL MATRIX ELEMENTS,  J1=L1+ 2/2,  J2=L2+ 1/2"):
    idx = text.find(header)
    if idx < 0:
        raise RuntimeError("block not found: " + header)
    chunk = text[idx:]
    lines = chunk.split("\n")
    rows = {}
    # data lines start after the two header lines; format: "%4d" L2 then 6 numbers each field width 11
    expected_L2 = 0
    for ln in lines[2:2+30]:
        if not ln.strip():
            break
        m = re.match(r'^\s*(\d+)\s*(.*)$', ln)
        if not m:
            break
        L2 = int(m.group(1))
        if L2 != expected_L2:
            # table ended (page break / next block header parsed as a bogus row) -- stop
            break
        rest = m.group(2)
        # numbers are fixed-width 11-char fields with no guaranteed spaces between negatives
        # e.g. "-8.992E-02 3.524E-02 9.025E-02-8.936E-02 9.327E-02 2.228E-01"
        nums = re.findall(r'-?\d\.\d{3}E[+-]\d{2}', rest)
        if len(nums) != 6:
            # a genuine data row always has exactly 6 numeric tokens; anything else
            # (page-break / header text swallowed by the regex) means the table ended
            break
        vals = [float(x) for x in nums[:6]]
        base = abs(L2 - ltr)
        for k in range(3):
            L1 = base + 2*k
            re_v = vals[2*k]
            im_v = vals[2*k+1]
            rows[(L1, L2)] = complex(re_v, im_v)
        expected_L2 += 1
    return rows

def main(lis_path, ltr=2, header=None):
    text = open(lis_path).read()
    k_alpha = parse_derived_k(text, 0)  # particle 1 = deuteron entrance
    k_beta  = parse_derived_k(text, 1)  # particle 2 = proton exit
    conv = 4*3.14159265358979/(k_alpha*k_beta)
    if header is None:
        # spin-orbit ON deck prints 6 j-split blocks; use the stretched/aligned one
        # (j_d=La+1, j_p=Lb+1/2) matching Nexus's convention. Spin-orbit OFF deck
        # collapses to a single unsplit block since U no longer depends on j.
        if "J1=L1+ 2/2,  J2=L2+ 1/2" in text:
            header = "RADIAL MATRIX ELEMENTS,  J1=L1+ 2/2,  J2=L2+ 1/2"
        elif "J1=L1+ 0/2,  J2=L2+ 0/2" in text:
            header = "RADIAL MATRIX ELEMENTS,  J1=L1+ 0/2,  J2=L2+ 0/2"
        else:
            raise RuntimeError("no recognized RADIAL MATRIX ELEMENTS block header found")
    rows = parse_block(text, ltr, header)
    print(f"# header='{header}'")
    print(f"# k_alpha={k_alpha} k_beta={k_beta} conv=4pi/(ka kb)={conv}")
    print("L_alpha\tL_beta\tRe_I\tIm_I")
    for (L1, L2) in sorted(rows.keys()):
        F = rows[(L1, L2)]
        if abs(F) < 1e-12:
            continue
        I = F * conv
        print(f"{L1}\t{L2}\t{I.real:.10e}\t{I.imag:.10e}")

if __name__ == "__main__":
    main(sys.argv[1])
