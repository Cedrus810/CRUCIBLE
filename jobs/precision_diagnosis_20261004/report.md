# Float32 arithmetic probe — diagnostic only

Input: the previous 48-frame independent GROMACS audit in
`../gromacs_engine_audit_20261004_3/`. The probe used its preserved translated
coordinates and the source AMBER parameters. It computed numpy-IR energies on
float32-rounded coordinates, plus a simple float32 harmonic-bond/acos-angle
arithmetic model. Per-frame results are in `results.json`.

For acetonitrile's original near-linear angle:

| Quantity | Delta from original numpy reference, kJ/mol |
|---|---:|
| Real GROMACS | approximately -0.001122927 |
| Float32 coordinates, float64 IR arithmetic | -0.000000606 |
| Simple float32 acos arithmetic model | -0.002692980 |

Coordinate casting alone does not explain the observed angle error; the simple
float32 arithmetic model is also not an exact model of GROMACS. Across all
frames its largest difference from GROMACS is 0.000121144 for bonds and
0.001570053 for angles. These observations motivate a double-precision control,
not a tolerance waiver or a claim that all residuals are already explained.
