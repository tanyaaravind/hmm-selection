# Benchmark summary

## FST results (summary)

**YRI_vs_CEU** — count=1831, mean=0.061, median=0.000, max=1.663

Top sites (pos, fst):

- 133245983: 1.6628

- 133245991: 1.6628

- 133246032: 1.6417

- 133246083: 1.6002

- 133281870: 1.4597

- 133281691: 1.4342

- 133258618: 1.3951

- 133258655: 1.3302

- 133280762: 1.2900

- 133247231: 1.2631



**YRI_vs_CHB** — count=1831, mean=0.074, median=0.000, max=2.255

Top sites (pos, fst):

- 133270337: 2.2553

- 133271928: 2.2494

- 133258618: 2.0989

- 133258655: 2.0989

- 133258835: 1.9928

- 133251128: 1.9201

- 133259404: 1.8770

- 133260437: 1.8104

- 133252969: 1.7173

- 133247823: 1.4784



**CEU_vs_CHB** — count=1831, mean=0.022, median=0.000, max=1.045

Top sites (pos, fst):

- 133251128: 1.0451

- 133283519: 0.7747

- 133279776: 0.7533

- 133262748: 0.7150

- 133281870: 0.6864

- 133280580: 0.6812

- 133280762: 0.6666

- 133279766: 0.6572

- 133281691: 0.6521

- 133260437: 0.6029



## Tajima's D (YRI) — sliding windows

Window size/step: see file; top deltaAF threshold (>= 0e+00) used to mark windows: 0.441561

count    179.000000
mean      -0.154881
std        0.580553
min       -1.308597
25%       -0.536169
50%       -0.208330
75%        0.168477
max        1.700531



## Notes for Sheki Okwayo (lso24) 

- These benchmarks used the allele-frequency table in `results/abo_expanded/allele_frequencies.csv` and the filtered delta AF table in `results/abo_expanded/delta_af_nonzero.csv`.

- I flagged windows overlapping positions with DeltaAF >= the 95th percentile (shown as red-shaded regions in the Tajima's D plot).

- If you want different population pairs, window sizes, or more annotation (e.g., overlapping genes), tell me and I’ll extend this.
