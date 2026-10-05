# Gate G0 报告

## ethanol
### config 0  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 0.477286675 | 0.477286675 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | angle | 1.032522937 | 1.032522937 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | proper | 0.008657437374 | 0.008657437374 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | improper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | -0.01682627822 | -0.01682627822 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb | 24.22182188 | 24.22182188 | 0 | 1.24e-06 | True |
| numpy-ir | reference-itself | lj14 | 0.7145709439 | 0.7145709439 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -34.88660758 | -34.88660758 | 0 | 1.35e-06 | True |
| numpy-ir | reference-itself | total | -8.448573984 | -8.448573984 | 0 | 1.08e-06 | True |
| openmm-ir | numpy-ir | bond | 0.477286675 | 0.477286675 | 5.55e-16 | 1e-06 | True |
| openmm-ir | numpy-ir | angle | 1.032522937 | 1.032522937 | 2.89e-15 | 1.01e-06 | True |
| openmm-ir | numpy-ir | proper | 0.008657437374 | 0.008657437374 | 4.16e-17 | 1e-06 | True |
| openmm-ir | numpy-ir | improper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | -0.01682627822 | -0.01682627822 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb | 24.22182188 | 24.22182188 | 0 | 1.24e-06 | True |
| openmm-ir | numpy-ir | lj14 | 0.7145709439 | 0.7145709439 | 4.44e-16 | 1.01e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -34.88660758 | -34.88660758 | 0 | 1.35e-06 | True |
| openmm-ir | numpy-ir | total | -8.448573984 | -8.448573984 | 3.55e-15 | 1.35e-06 | True |
| openmm-xml | numpy-ir | bond | 0.477286675 | 0.477286675 | 5.55e-16 | 1e-06 | True |
| openmm-xml | numpy-ir | angle | 1.032522937 | 1.032522937 | 2.89e-15 | 1.01e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 0.008657437374 | 0.008657437374 | 4.16e-17 | 1e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | -9.967041159 | -9.967041033 | 1.26e-07 | 1.1e-06 | True |
| openmm-xml | numpy-ir | total | -8.44857411 | -8.448573984 | 1.26e-07 | 1.35e-06 | True |
| openmm-prmtop | numpy-ir | bond | 0.477286675 | 0.477286675 | 5.55e-16 | 0.0001 | True |
| openmm-prmtop | numpy-ir | angle | 1.032522937 | 1.032522937 | 8.88e-16 | 0.0001 | True |
| openmm-prmtop | numpy-ir | torsion_all | 0.008657437374 | 0.008657437374 | 4.16e-17 | 0.0001 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | -9.967041163 | -9.967041033 | 1.3e-07 | 0.000101 | True |
| openmm-prmtop | numpy-ir | total | -8.448574114 | -8.448573984 | 1.3e-07 | 0.000103 | True |
| openmm-top | numpy-ir | bond | 0.477286675 | 0.477286675 | 5.16e-15 | 0.0001 | True |
| openmm-top | numpy-ir | angle | 1.032522956 | 1.032522937 | 1.87e-08 | 0.0001 | True |
| openmm-top | numpy-ir | torsion_all | 0.008657437374 | 0.008657437374 | 2.27e-14 | 0.0001 | True |
| openmm-top | numpy-ir | nonbonded_all | -9.967041176 | -9.967041033 | 1.43e-07 | 0.000101 | True |
| openmm-top | numpy-ir | total | -8.448574108 | -8.448573984 | 1.25e-07 | 0.000103 | True |
| sander-ours | numpy-ir | bond | 0.1141 | 0.1140742531 | 2.57e-05 | 0.0002 | True |
| sander-ours | numpy-ir | angle | 0.2468 | 0.2467789046 | 2.11e-05 | 0.0002 | True |
| sander-ours | numpy-ir | dihed | 0.0021 | 0.002069177193 | 3.08e-05 | 0.0002 | True |
| sander-ours | numpy-ir | vdw | -0.004 | -0.004021577012 | 2.16e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 5.789 | 5.78915437 | 0.000154 | 0.000489 | True |
| sander-ours | numpy-ir | vdw14 | 0.1708 | 0.1707865545 | 1.34e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -8.3378 | -8.338099326 | 0.000299 | 0.000617 | True |
| sander-tleap | sander-ours | bond | 0.1141 | 0.1141 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | angle | 0.2468 | 0.2468 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | dihed | 0.0021 | 0.0021 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw | -0.004 | -0.004 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 5.789 | 5.789 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | vdw14 | 0.1708 | 0.1708 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -8.3378 | -8.3378 | 0 | 0.000201 | True |
| gmx | numpy-ir | bond | 0.4772803783 | 0.477286675 | 6.3e-06 | 0.0001 | True |
| gmx | numpy-ir | angle | 1.032525063 | 1.032522937 | 2.13e-06 | 0.0001 | True |
| gmx | numpy-ir | proper | 0.008657487109 | 0.008657437374 | 4.97e-08 | 0.0001 | True |
| gmx | numpy-ir | improper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj | -0.01682626456 | -0.01682627822 | 1.37e-08 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 24.22181892 | 24.22182188 | 2.96e-06 | 0.000102 | True |
| gmx | numpy-ir | lj14 | 0.7145693898 | 0.7145709439 | 1.55e-06 | 0.0001 | True |
| gmx | numpy-ir | coulomb14 | -34.88660049 | -34.88660758 | 7.09e-06 | 0.000103 | True |
| gmx | numpy-ir | total | -8.448575974 | -8.448573984 | 1.99e-06 | 0.000101 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 2.61e-08 | 0.001 | True |
| openmm-top | openmm-ir | 8.16e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 1  passed=False

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 100.2243053 | 100.2243053 | 0 | 2e-06 | True |
| numpy-ir | reference-itself | angle | 21.71297483 | 21.71297483 | 0 | 1.22e-06 | True |
| numpy-ir | reference-itself | proper | 0.2951387194 | 0.2951387194 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | improper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | -0.01425454276 | -0.01425454276 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb | 23.77185701 | 23.77185701 | 0 | 1.24e-06 | True |
| numpy-ir | reference-itself | lj14 | 0.169373514 | 0.169373514 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -34.54968119 | -34.54968119 | 0 | 1.35e-06 | True |
| numpy-ir | reference-itself | total | 111.6097137 | 111.6097137 | 0 | 2.12e-06 | True |
| openmm-ir | numpy-ir | bond | 100.2243053 | 100.2243053 | 2.84e-14 | 2e-06 | True |
| openmm-ir | numpy-ir | angle | 21.71297483 | 21.71297483 | 3.91e-14 | 1.22e-06 | True |
| openmm-ir | numpy-ir | proper | 0.2951387194 | 0.2951387194 | 3.33e-16 | 1e-06 | True |
| openmm-ir | numpy-ir | improper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | -0.01425454276 | -0.01425454276 | 1.73e-18 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb | 23.77185701 | 23.77185701 | 0 | 1.24e-06 | True |
| openmm-ir | numpy-ir | lj14 | 0.169373514 | 0.169373514 | 5.55e-16 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -34.54968119 | -34.54968119 | 7.11e-15 | 1.35e-06 | True |
| openmm-ir | numpy-ir | total | 111.6097137 | 111.6097137 | 7.11e-14 | 2e-06 | True |
| openmm-xml | numpy-ir | bond | 100.2243053 | 100.2243053 | 2.84e-14 | 2e-06 | True |
| openmm-xml | numpy-ir | angle | 21.71297483 | 21.71297483 | 3.91e-14 | 1.22e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 0.2951387194 | 0.2951387194 | 3.33e-16 | 1e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | -10.62270534 | -10.62270521 | 1.28e-07 | 1.11e-06 | True |
| openmm-xml | numpy-ir | total | 111.6097135 | 111.6097137 | 1.28e-07 | 2e-06 | True |
| openmm-prmtop | numpy-ir | bond | 100.2243053 | 100.2243053 | 2.84e-14 | 0.00011 | True |
| openmm-prmtop | numpy-ir | angle | 21.71297483 | 21.71297483 | 3.2e-14 | 0.000102 | True |
| openmm-prmtop | numpy-ir | torsion_all | 0.2951387194 | 0.2951387194 | 3.33e-16 | 0.0001 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | -10.62270534 | -10.62270521 | 1.31e-07 | 0.000101 | True |
| openmm-prmtop | numpy-ir | total | 111.6097135 | 111.6097137 | 1.31e-07 | 0.00011 | True |
| openmm-top | numpy-ir | bond | 100.2243053 | 100.2243053 | 9.95e-14 | 0.00011 | True |
| openmm-top | numpy-ir | angle | 21.71297487 | 21.71297483 | 4.24e-08 | 0.000102 | True |
| openmm-top | numpy-ir | torsion_all | 0.2951387194 | 0.2951387194 | 1.55e-12 | 0.0001 | True |
| openmm-top | numpy-ir | nonbonded_all | -10.62270536 | -10.62270521 | 1.45e-07 | 0.000101 | True |
| openmm-top | numpy-ir | total | 111.6097136 | 111.6097137 | 1.03e-07 | 0.00011 | True |
| sander-ours | numpy-ir | bond | 23.9542 | 23.95418387 | 1.61e-05 | 0.000202 | True |
| sander-ours | numpy-ir | angle | 5.1895 | 5.189525533 | 2.55e-05 | 0.000201 | True |
| sander-ours | numpy-ir | dihed | 0.0705 | 0.0705398469 | 3.98e-05 | 0.0002 | True |
| sander-ours | numpy-ir | vdw | -0.0034 | -0.003406917487 | 6.92e-06 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 5.6814 | 5.681610184 | 0.00021 | 0.000484 | True |
| sander-ours | numpy-ir | vdw14 | 0.0405 | 0.0404812414 | 1.88e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -8.2573 | -8.257571986 | 0.000272 | 0.000613 | True |
| sander-tleap | sander-ours | bond | 23.9542 | 23.9542 | 0 | 0.000202 | True |
| sander-tleap | sander-ours | angle | 5.1895 | 5.1895 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | dihed | 0.0705 | 0.0705 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw | -0.0034 | -0.0034 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 5.6814 | 5.6814 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | vdw14 | 0.0405 | 0.0405 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -8.2573 | -8.2573 | 0 | 0.000201 | True |
| gmx | numpy-ir | bond | 100.2241364 | 100.2243053 | 0.000169 | 0.00011 | False |
| gmx | numpy-ir | angle | 21.71297264 | 21.71297483 | 2.19e-06 | 0.000102 | True |
| gmx | numpy-ir | proper | 0.2951387465 | 0.2951387194 | 2.71e-08 | 0.0001 | True |
| gmx | numpy-ir | improper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj | -0.01425454766 | -0.01425454276 | 4.89e-09 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 23.77185631 | 23.77185701 | 7e-07 | 0.000102 | True |
| gmx | numpy-ir | lj14 | 0.1693737656 | 0.169373514 | 2.52e-07 | 0.0001 | True |
| gmx | numpy-ir | coulomb14 | -34.5496788 | -34.54968119 | 2.39e-06 | 0.000103 | True |
| gmx | numpy-ir | total | 111.6095505 | 111.6097137 | 0.000163 | 0.000111 | False |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 2.33e-08 | 0.001 | True |
| openmm-top | openmm-ir | 7.46e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 2  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 21.40412709 | 21.40412709 | 0 | 1.21e-06 | True |
| numpy-ir | reference-itself | angle | 13.58284109 | 13.58284109 | 0 | 1.14e-06 | True |
| numpy-ir | reference-itself | proper | 0.3103369131 | 0.3103369131 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | improper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | -0.01437304707 | -0.01437304707 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb | 23.52402371 | 23.52402371 | 0 | 1.24e-06 | True |
| numpy-ir | reference-itself | lj14 | 0.5829905996 | 0.5829905996 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -34.10783546 | -34.10783546 | 0 | 1.34e-06 | True |
| numpy-ir | reference-itself | total | 25.2821109 | 25.2821109 | 0 | 1.25e-06 | True |
| openmm-ir | numpy-ir | bond | 21.40412709 | 21.40412709 | 0 | 1.21e-06 | True |
| openmm-ir | numpy-ir | angle | 13.58284109 | 13.58284109 | 3.55e-15 | 1.14e-06 | True |
| openmm-ir | numpy-ir | proper | 0.3103369131 | 0.3103369131 | 2.78e-16 | 1e-06 | True |
| openmm-ir | numpy-ir | improper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | -0.01437304707 | -0.01437304707 | 1.73e-18 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb | 23.52402371 | 23.52402371 | 0 | 1.24e-06 | True |
| openmm-ir | numpy-ir | lj14 | 0.5829905996 | 0.5829905996 | 2.22e-16 | 1.01e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -34.10783546 | -34.10783546 | 0 | 1.34e-06 | True |
| openmm-ir | numpy-ir | total | 25.2821109 | 25.2821109 | 3.55e-15 | 1.34e-06 | True |
| openmm-xml | numpy-ir | bond | 21.40412709 | 21.40412709 | 0 | 1.21e-06 | True |
| openmm-xml | numpy-ir | angle | 13.58284109 | 13.58284109 | 3.55e-15 | 1.14e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 0.3103369131 | 0.3103369131 | 2.78e-16 | 1e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | -10.01519432 | -10.0151942 | 1.25e-07 | 1.1e-06 | True |
| openmm-xml | numpy-ir | total | 25.28211077 | 25.2821109 | 1.25e-07 | 1.34e-06 | True |
| openmm-prmtop | numpy-ir | bond | 21.40412709 | 21.40412709 | 0 | 0.000102 | True |
| openmm-prmtop | numpy-ir | angle | 13.58284109 | 13.58284109 | 1.24e-14 | 0.000101 | True |
| openmm-prmtop | numpy-ir | torsion_all | 0.3103369131 | 0.3103369131 | 2.78e-16 | 0.0001 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | -10.01519432 | -10.0151942 | 1.29e-07 | 0.000101 | True |
| openmm-prmtop | numpy-ir | total | 25.28211077 | 25.2821109 | 1.29e-07 | 0.000103 | True |
| openmm-top | numpy-ir | bond | 21.40412709 | 21.40412709 | 2.13e-14 | 0.000102 | True |
| openmm-top | numpy-ir | angle | 13.58284099 | 13.58284109 | 1.02e-07 | 0.000101 | True |
| openmm-top | numpy-ir | torsion_all | 0.3103369131 | 0.3103369131 | 6.97e-13 | 0.0001 | True |
| openmm-top | numpy-ir | nonbonded_all | -10.01519434 | -10.0151942 | 1.41e-07 | 0.000101 | True |
| openmm-top | numpy-ir | total | 25.28211066 | 25.2821109 | 2.43e-07 | 0.000103 | True |
| sander-ours | numpy-ir | bond | 5.1157 | 5.115709151 | 9.15e-06 | 0.000201 | True |
| sander-ours | numpy-ir | angle | 3.2464 | 3.246376935 | 2.31e-05 | 0.0002 | True |
| sander-ours | numpy-ir | dihed | 0.0742 | 0.07417230236 | 2.77e-05 | 0.0002 | True |
| sander-ours | numpy-ir | vdw | -0.0034 | -0.003435240695 | 3.52e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 5.6222 | 5.622376603 | 0.000177 | 0.000481 | True |
| sander-ours | numpy-ir | vdw14 | 0.1393 | 0.1393380974 | 3.81e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -8.1517 | -8.151968321 | 0.000268 | 0.000608 | True |
| sander-tleap | sander-ours | bond | 5.1157 | 5.1157 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | angle | 3.2464 | 3.2464 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | dihed | 0.0742 | 0.0742 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw | -0.0034 | -0.0034 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 5.6222 | 5.6222 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | vdw14 | 0.1393 | 0.1393 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -8.1517 | -8.1517 | 0 | 0.000201 | True |
| gmx | numpy-ir | bond | 21.40408897 | 21.40412709 | 3.81e-05 | 0.000102 | True |
| gmx | numpy-ir | angle | 13.58284664 | 13.58284109 | 5.55e-06 | 0.000101 | True |
| gmx | numpy-ir | proper | 0.3103371561 | 0.3103369131 | 2.43e-07 | 0.0001 | True |
| gmx | numpy-ir | improper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj | -0.01437304728 | -0.01437304707 | 2.1e-10 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 23.52402306 | 23.52402371 | 6.53e-07 | 0.000102 | True |
| gmx | numpy-ir | lj14 | 0.5829878449 | 0.5829905996 | 2.75e-06 | 0.0001 | True |
| gmx | numpy-ir | coulomb14 | -34.10782623 | -34.10783546 | 9.22e-06 | 0.000103 | True |
| gmx | numpy-ir | total | 25.28208351 | 25.2821109 | 2.74e-05 | 0.000103 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 2.24e-08 | 0.001 | True |
| openmm-top | openmm-ir | 7.81e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 3  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 24.65416149 | 24.65416149 | 0 | 1.25e-06 | True |
| numpy-ir | reference-itself | angle | 7.750801911 | 7.750801911 | 0 | 1.08e-06 | True |
| numpy-ir | reference-itself | proper | 0.166676559 | 0.166676559 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | improper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | -0.01706189404 | -0.01706189404 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb | 24.44033409 | 24.44033409 | 0 | 1.24e-06 | True |
| numpy-ir | reference-itself | lj14 | 1.186174426 | 1.186174426 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -35.21801192 | -35.21801192 | 0 | 1.35e-06 | True |
| numpy-ir | reference-itself | total | 22.96307467 | 22.96307467 | 0 | 1.23e-06 | True |
| openmm-ir | numpy-ir | bond | 24.65416149 | 24.65416149 | 7.11e-15 | 1.25e-06 | True |
| openmm-ir | numpy-ir | angle | 7.750801911 | 7.750801911 | 7.11e-15 | 1.08e-06 | True |
| openmm-ir | numpy-ir | proper | 0.166676559 | 0.166676559 | 3.05e-16 | 1e-06 | True |
| openmm-ir | numpy-ir | improper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | -0.01706189404 | -0.01706189404 | 3.47e-18 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb | 24.44033409 | 24.44033409 | 0 | 1.24e-06 | True |
| openmm-ir | numpy-ir | lj14 | 1.186174426 | 1.186174426 | 0 | 1.01e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -35.21801192 | -35.21801192 | 0 | 1.35e-06 | True |
| openmm-ir | numpy-ir | total | 22.96307467 | 22.96307467 | 0 | 1.35e-06 | True |
| openmm-xml | numpy-ir | bond | 24.65416149 | 24.65416149 | 7.11e-15 | 1.25e-06 | True |
| openmm-xml | numpy-ir | angle | 7.750801911 | 7.750801911 | 7.11e-15 | 1.08e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 0.166676559 | 0.166676559 | 3.05e-16 | 1e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | -9.60856542 | -9.608565293 | 1.28e-07 | 1.1e-06 | True |
| openmm-xml | numpy-ir | total | 22.96307454 | 22.96307467 | 1.28e-07 | 1.35e-06 | True |
| openmm-prmtop | numpy-ir | bond | 24.65416149 | 24.65416149 | 3.55e-15 | 0.000102 | True |
| openmm-prmtop | numpy-ir | angle | 7.750801911 | 7.750801911 | 1.07e-14 | 0.000101 | True |
| openmm-prmtop | numpy-ir | torsion_all | 0.166676559 | 0.166676559 | 3.05e-16 | 0.0001 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | -9.608565425 | -9.608565293 | 1.32e-07 | 0.000101 | True |
| openmm-prmtop | numpy-ir | total | 22.96307453 | 22.96307467 | 1.32e-07 | 0.000104 | True |
| openmm-top | numpy-ir | bond | 24.65416149 | 24.65416149 | 1.78e-14 | 0.000102 | True |
| openmm-top | numpy-ir | angle | 7.750801931 | 7.750801911 | 2.05e-08 | 0.000101 | True |
| openmm-top | numpy-ir | torsion_all | 0.166676559 | 0.166676559 | 6.07e-13 | 0.0001 | True |
| openmm-top | numpy-ir | nonbonded_all | -9.608565438 | -9.608565293 | 1.45e-07 | 0.000101 | True |
| openmm-top | numpy-ir | total | 22.96307454 | 22.96307467 | 1.24e-07 | 0.000104 | True |
| sander-ours | numpy-ir | bond | 5.8925 | 5.892486015 | 1.4e-05 | 0.000201 | True |
| sander-ours | numpy-ir | angle | 1.8525 | 1.852486116 | 1.39e-05 | 0.0002 | True |
| sander-ours | numpy-ir | dihed | 0.0398 | 0.03983665368 | 3.67e-05 | 0.0002 | True |
| sander-ours | numpy-ir | vdw | -0.0041 | -0.004077890544 | 2.21e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 5.8412 | 5.841380041 | 0.00018 | 0.000492 | True |
| sander-ours | numpy-ir | vdw14 | 0.2835 | 0.2835024918 | 2.49e-06 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -8.417 | -8.417306863 | 0.000307 | 0.000621 | True |
| sander-tleap | sander-ours | bond | 5.8925 | 5.8925 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | angle | 1.8525 | 1.8525 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | dihed | 0.0398 | 0.0398 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw | -0.0041 | -0.0041 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 5.8412 | 5.8412 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | vdw14 | 0.2835 | 0.2835 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -8.417 | -8.417 | 0 | 0.000201 | True |
| gmx | numpy-ir | bond | 24.6541729 | 24.65416149 | 1.14e-05 | 0.000102 | True |
| gmx | numpy-ir | angle | 7.750789642 | 7.750801911 | 1.23e-05 | 0.000101 | True |
| gmx | numpy-ir | proper | 0.1666763872 | 0.166676559 | 1.72e-07 | 0.0001 | True |
| gmx | numpy-ir | improper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj | -0.0170618929 | -0.01706189404 | 1.14e-09 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 24.44033051 | 24.44033409 | 3.59e-06 | 0.000102 | True |
| gmx | numpy-ir | lj14 | 1.186174393 | 1.186174426 | 3.3e-08 | 0.0001 | True |
| gmx | numpy-ir | coulomb14 | -35.21800995 | -35.21801192 | 1.97e-06 | 0.000104 | True |
| gmx | numpy-ir | total | 22.96307182 | 22.96307467 | 2.84e-06 | 0.000102 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 3.75e-08 | 0.001 | True |
| openmm-top | openmm-ir | 8.1e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

## benzene
### config 0  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 0.1126418815 | 0.1126418815 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | angle | 0.0113812304 | 0.0113812304 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | proper | 3.366673755e-05 | 3.366673755e-05 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | improper | 3.399055073e-06 | 3.399055073e-06 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | -0.7432418659 | -0.7432418659 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | coulomb | 10.62810423 | 10.62810423 | 0 | 1.11e-06 | True |
| numpy-ir | reference-itself | lj14 | 12.12007097 | 12.12007097 | 0 | 1.12e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -0.6926374943 | -0.6926374943 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | total | 21.43635602 | 21.43635602 | 0 | 1.21e-06 | True |
| openmm-ir | numpy-ir | bond | 0.1126418815 | 0.1126418815 | 1.11e-16 | 1e-06 | True |
| openmm-ir | numpy-ir | angle | 0.0113812304 | 0.0113812304 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | proper | 3.366673755e-05 | 3.366673755e-05 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | improper | 3.399055073e-06 | 3.399055073e-06 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | -0.7432418659 | -0.7432418659 | 0 | 1.01e-06 | True |
| openmm-ir | numpy-ir | coulomb | 10.62810423 | 10.62810423 | 3.55e-15 | 1.11e-06 | True |
| openmm-ir | numpy-ir | lj14 | 12.12007097 | 12.12007097 | 3.55e-15 | 1.12e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -0.6926374943 | -0.6926374943 | 0 | 1.01e-06 | True |
| openmm-ir | numpy-ir | total | 21.43635602 | 21.43635602 | 7.11e-15 | 1.12e-06 | True |
| openmm-xml | numpy-ir | bond | 0.1126418815 | 0.1126418815 | 1.11e-16 | 1e-06 | True |
| openmm-xml | numpy-ir | angle | 0.0113812304 | 0.0113812304 | 0 | 1e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 3.706579262e-05 | 3.706579262e-05 | 0 | 1e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | 21.31229596 | 21.31229584 | 1.18e-07 | 1.21e-06 | True |
| openmm-xml | numpy-ir | total | 21.43635614 | 21.43635602 | 1.18e-07 | 1.12e-06 | True |
| openmm-prmtop | numpy-ir | bond | 0.1126418815 | 0.1126418815 | 1.11e-16 | 0.0001 | True |
| openmm-prmtop | numpy-ir | angle | 0.0113812304 | 0.0113812304 | 0 | 0.0001 | True |
| openmm-prmtop | numpy-ir | torsion_all | 3.706579262e-05 | 3.706579262e-05 | 0 | 0.0001 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | 21.31229595 | 21.31229584 | 1.14e-07 | 0.000102 | True |
| openmm-prmtop | numpy-ir | total | 21.43635613 | 21.43635602 | 1.14e-07 | 0.000101 | True |
| openmm-top | numpy-ir | bond | 0.1126418815 | 0.1126418815 | 1.19e-15 | 0.0001 | True |
| openmm-top | numpy-ir | angle | 0.01138123904 | 0.0113812304 | 8.63e-09 | 0.0001 | True |
| openmm-top | numpy-ir | torsion_all | 3.706579897e-05 | 3.706579262e-05 | 6.35e-12 | 0.0001 | True |
| openmm-top | numpy-ir | nonbonded_all | 21.31229596 | 21.31229584 | 1.16e-07 | 0.000102 | True |
| openmm-top | numpy-ir | total | 21.43635614 | 21.43635602 | 1.25e-07 | 0.000101 | True |
| sander-ours | numpy-ir | bond | 0.0269 | 0.02692205581 | 2.21e-05 | 0.0002 | True |
| sander-ours | numpy-ir | angle | 0.0027 | 0.002720179351 | 2.02e-05 | 0.0002 | True |
| sander-ours | numpy-ir | dihed | 0 | 8.858937051e-06 | 8.86e-06 | 0.0002 | True |
| sander-ours | numpy-ir | vdw | -0.1776 | -0.1776390693 | 3.91e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 2.5401 | 2.540177876 | 7.79e-05 | 0.000327 | True |
| sander-ours | numpy-ir | vdw14 | 2.8968 | 2.896766483 | 3.35e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -0.1655 | -0.1655443342 | 4.43e-05 | 0.000208 | True |
| sander-tleap | sander-ours | bond | 0.0269 | 0.0269 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | angle | 0.0027 | 0.0027 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | dihed | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw | -0.1776 | -0.1776 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 2.5401 | 2.5401 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | 2.8968 | 2.8968 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -0.1655 | -0.1655 | 0 | 0.0002 | True |
| gmx | numpy-ir | bond | 0.1126409471 | 0.1126418815 | 9.34e-07 | 0.0001 | True |
| gmx | numpy-ir | angle | 0.01138130669 | 0.0113812304 | 7.63e-08 | 0.0001 | True |
| gmx | numpy-ir | proper | 3.7968988e-05 | 3.366673755e-05 | 4.3e-06 | 0.0001 | True |
| gmx | numpy-ir | improper | 2.743244e-06 | 3.399055073e-06 | 6.56e-07 | 0.0001 | True |
| gmx | numpy-ir | lj | -0.7432420254 | -0.7432418659 | 1.6e-07 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 10.62810802 | 10.62810423 | 3.79e-06 | 0.000101 | True |
| gmx | numpy-ir | lj14 | 12.12007427 | 12.12007097 | 3.31e-06 | 0.000101 | True |
| gmx | numpy-ir | coulomb14 | -0.6926355362 | -0.6926374943 | 1.96e-06 | 0.0001 | True |
| gmx | numpy-ir | total | 21.43636703 | 21.43635602 | 1.1e-05 | 0.000102 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 1.18e-08 | 0.001 | True |
| openmm-top | openmm-ir | 1.92e-07 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 1  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 119.1823692 | 119.1823692 | 0 | 2.19e-06 | True |
| numpy-ir | reference-itself | angle | 19.41553459 | 19.41553459 | 0 | 1.19e-06 | True |
| numpy-ir | reference-itself | proper | 11.55933001 | 11.55933001 | 0 | 1.12e-06 | True |
| numpy-ir | reference-itself | improper | 1.357188137 | 1.357188137 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | lj | -0.8178862794 | -0.8178862794 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | coulomb | 10.77626748 | 10.77626748 | 0 | 1.11e-06 | True |
| numpy-ir | reference-itself | lj14 | 18.5092036 | 18.5092036 | 0 | 1.19e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -0.7174716825 | -0.7174716825 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | total | 179.2645351 | 179.2645351 | 0 | 2.79e-06 | True |
| openmm-ir | numpy-ir | bond | 119.1823692 | 119.1823692 | 0 | 2.19e-06 | True |
| openmm-ir | numpy-ir | angle | 19.41553459 | 19.41553459 | 0 | 1.19e-06 | True |
| openmm-ir | numpy-ir | proper | 11.55933001 | 11.55933001 | 2.49e-14 | 1.12e-06 | True |
| openmm-ir | numpy-ir | improper | 1.357188137 | 1.357188137 | 6.44e-15 | 1.01e-06 | True |
| openmm-ir | numpy-ir | lj | -0.8178862794 | -0.8178862794 | 0 | 1.01e-06 | True |
| openmm-ir | numpy-ir | coulomb | 10.77626748 | 10.77626748 | 0 | 1.11e-06 | True |
| openmm-ir | numpy-ir | lj14 | 18.5092036 | 18.5092036 | 0 | 1.19e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -0.7174716825 | -0.7174716825 | 2.66e-15 | 1.01e-06 | True |
| openmm-ir | numpy-ir | total | 179.2645351 | 179.2645351 | 2.84e-14 | 2.19e-06 | True |
| openmm-xml | numpy-ir | bond | 119.1823692 | 119.1823692 | 0 | 2.19e-06 | True |
| openmm-xml | numpy-ir | angle | 19.41553459 | 19.41553459 | 0 | 1.19e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 12.91651815 | 12.91651815 | 3.02e-14 | 1.13e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | 27.75011323 | 27.75011311 | 1.19e-07 | 1.28e-06 | True |
| openmm-xml | numpy-ir | total | 179.2645352 | 179.2645351 | 1.19e-07 | 2.19e-06 | True |
| openmm-prmtop | numpy-ir | bond | 119.1823692 | 119.1823692 | 0 | 0.000112 | True |
| openmm-prmtop | numpy-ir | angle | 19.41553459 | 19.41553459 | 0 | 0.000102 | True |
| openmm-prmtop | numpy-ir | torsion_all | 12.91651815 | 12.91651815 | 3.38e-14 | 0.000101 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | 27.75011323 | 27.75011311 | 1.14e-07 | 0.000103 | True |
| openmm-prmtop | numpy-ir | total | 179.2645352 | 179.2645351 | 1.14e-07 | 0.000112 | True |
| openmm-top | numpy-ir | bond | 119.1823692 | 119.1823692 | 2.84e-14 | 0.000112 | True |
| openmm-top | numpy-ir | angle | 19.41553458 | 19.41553459 | 4.69e-09 | 0.000102 | True |
| openmm-top | numpy-ir | torsion_all | 12.91651815 | 12.91651815 | 2.45e-09 | 0.000101 | True |
| openmm-top | numpy-ir | nonbonded_all | 27.75011323 | 27.75011311 | 1.17e-07 | 0.000103 | True |
| openmm-top | numpy-ir | total | 179.2645352 | 179.2645351 | 1.1e-07 | 0.000112 | True |
| sander-ours | numpy-ir | bond | 28.4853 | 28.4852699 | 3.01e-05 | 0.000203 | True |
| sander-ours | numpy-ir | angle | 4.6404 | 4.640424137 | 2.41e-05 | 0.0002 | True |
| sander-ours | numpy-ir | dihed | 3.0871 | 3.087121929 | 2.19e-05 | 0.0002 | True |
| sander-ours | numpy-ir | vdw | -0.1955 | -0.1954795123 | 2.05e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 2.5755 | 2.575589741 | 8.97e-05 | 0.000329 | True |
| sander-ours | numpy-ir | vdw14 | 4.4238 | 4.423805832 | 5.83e-06 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -0.1715 | -0.1714798476 | 2.02e-05 | 0.000209 | True |
| sander-tleap | sander-ours | bond | 28.4853 | 28.4853 | 0 | 0.000203 | True |
| sander-tleap | sander-ours | angle | 4.6404 | 4.6404 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | dihed | 3.0871 | 3.0871 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw | -0.1955 | -0.1955 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 2.5755 | 2.5755 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | 4.4238 | 4.4238 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -0.1715 | -0.1715 | 0 | 0.0002 | True |
| gmx | numpy-ir | bond | 119.1824417 | 119.1823692 | 7.25e-05 | 0.000112 | True |
| gmx | numpy-ir | angle | 19.41554832 | 19.41553459 | 1.37e-05 | 0.000102 | True |
| gmx | numpy-ir | proper | 11.55932713 | 11.55933001 | 2.89e-06 | 0.000101 | True |
| gmx | numpy-ir | improper | 1.357187986 | 1.357188137 | 1.5e-07 | 0.0001 | True |
| gmx | numpy-ir | lj | -0.8178867102 | -0.8178862794 | 4.31e-07 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 10.77626705 | 10.77626748 | 4.24e-07 | 0.000101 | True |
| gmx | numpy-ir | lj14 | 18.50921822 | 18.5092036 | 1.46e-05 | 0.000102 | True |
| gmx | numpy-ir | coulomb14 | -0.7174773216 | -0.7174716825 | 5.64e-06 | 0.0001 | True |
| gmx | numpy-ir | total | 179.2646179 | 179.2645351 | 8.28e-05 | 0.000118 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 1.75e-08 | 0.001 | True |
| openmm-top | openmm-ir | 1.35e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 2  passed=False

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 61.89138275 | 61.89138275 | 0 | 1.62e-06 | True |
| numpy-ir | reference-itself | angle | 28.57377621 | 28.57377621 | 0 | 1.29e-06 | True |
| numpy-ir | reference-itself | proper | 1.90664744 | 1.90664744 | 0 | 1.02e-06 | True |
| numpy-ir | reference-itself | improper | 0.2296190583 | 0.2296190583 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | -0.7434275408 | -0.7434275408 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | coulomb | 10.95522143 | 10.95522143 | 0 | 1.11e-06 | True |
| numpy-ir | reference-itself | lj14 | 10.08723169 | 10.08723169 | 0 | 1.1e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -0.6674225375 | -0.6674225375 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | total | 112.2330285 | 112.2330285 | 0 | 2.12e-06 | True |
| openmm-ir | numpy-ir | bond | 61.89138275 | 61.89138275 | 6.39e-14 | 1.62e-06 | True |
| openmm-ir | numpy-ir | angle | 28.57377621 | 28.57377621 | 1.78e-14 | 1.29e-06 | True |
| openmm-ir | numpy-ir | proper | 1.90664744 | 1.90664744 | 1.33e-15 | 1.02e-06 | True |
| openmm-ir | numpy-ir | improper | 0.2296190583 | 0.2296190583 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | -0.7434275408 | -0.7434275408 | 2.22e-16 | 1.01e-06 | True |
| openmm-ir | numpy-ir | coulomb | 10.95522143 | 10.95522143 | 0 | 1.11e-06 | True |
| openmm-ir | numpy-ir | lj14 | 10.08723169 | 10.08723169 | 1.78e-15 | 1.1e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -0.6674225375 | -0.6674225375 | 8.88e-16 | 1.01e-06 | True |
| openmm-ir | numpy-ir | total | 112.2330285 | 112.2330285 | 8.53e-14 | 1.62e-06 | True |
| openmm-xml | numpy-ir | bond | 61.89138275 | 61.89138275 | 6.39e-14 | 1.62e-06 | True |
| openmm-xml | numpy-ir | angle | 28.57377621 | 28.57377621 | 1.78e-14 | 1.29e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 2.136266499 | 2.136266499 | 1.33e-15 | 1.02e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | 19.63160316 | 19.63160304 | 1.22e-07 | 1.2e-06 | True |
| openmm-xml | numpy-ir | total | 112.2330286 | 112.2330285 | 1.22e-07 | 1.62e-06 | True |
| openmm-prmtop | numpy-ir | bond | 61.89138275 | 61.89138275 | 6.39e-14 | 0.000106 | True |
| openmm-prmtop | numpy-ir | angle | 28.57377621 | 28.57377621 | 1.42e-14 | 0.000103 | True |
| openmm-prmtop | numpy-ir | torsion_all | 2.136266499 | 2.136266499 | 1.33e-15 | 0.0001 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | 19.63160316 | 19.63160304 | 1.18e-07 | 0.000102 | True |
| openmm-prmtop | numpy-ir | total | 112.2330286 | 112.2330285 | 1.18e-07 | 0.000106 | True |
| openmm-top | numpy-ir | bond | 61.89138275 | 61.89138275 | 1.42e-14 | 0.000106 | True |
| openmm-top | numpy-ir | angle | 28.57377622 | 28.57377621 | 6.01e-09 | 0.000103 | True |
| openmm-top | numpy-ir | torsion_all | 2.136266498 | 2.136266499 | 5.68e-10 | 0.0001 | True |
| openmm-top | numpy-ir | nonbonded_all | 19.63160316 | 19.63160304 | 1.2e-07 | 0.000102 | True |
| openmm-top | numpy-ir | total | 112.2330286 | 112.2330285 | 1.26e-07 | 0.000106 | True |
| sander-ours | numpy-ir | bond | 14.7924 | 14.7923955 | 4.5e-06 | 0.000201 | True |
| sander-ours | numpy-ir | angle | 6.8293 | 6.829296418 | 3.58e-06 | 0.000201 | True |
| sander-ours | numpy-ir | dihed | 0.5106 | 0.5105799471 | 2.01e-05 | 0.0002 | True |
| sander-ours | numpy-ir | vdw | -0.1777 | -0.1776834466 | 1.66e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 2.6183 | 2.618360761 | 6.08e-05 | 0.000331 | True |
| sander-ours | numpy-ir | vdw14 | 2.4109 | 2.410906236 | 6.24e-06 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -0.1595 | -0.1595178149 | 1.78e-05 | 0.000208 | True |
| sander-tleap | sander-ours | bond | 14.7924 | 14.7924 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | angle | 6.8293 | 6.8293 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | dihed | 0.5106 | 0.5106 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw | -0.1777 | -0.1777 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 2.6183 | 2.6183 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | 2.4109 | 2.4109 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -0.1595 | -0.1595 | 0 | 0.0002 | True |
| gmx | numpy-ir | bond | 61.89121628 | 61.89138275 | 0.000166 | 0.000106 | False |
| gmx | numpy-ir | angle | 28.57376099 | 28.57377621 | 1.52e-05 | 0.000103 | True |
| gmx | numpy-ir | proper | 1.906648993 | 1.90664744 | 1.55e-06 | 0.0001 | True |
| gmx | numpy-ir | improper | 0.2296183109 | 0.2296190583 | 7.47e-07 | 0.0001 | True |
| gmx | numpy-ir | lj | -0.7434278727 | -0.7434275408 | 3.32e-07 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 10.95521736 | 10.95522143 | 4.06e-06 | 0.000101 | True |
| gmx | numpy-ir | lj14 | 10.08724213 | 10.08723169 | 1.04e-05 | 0.000101 | True |
| gmx | numpy-ir | coulomb14 | -0.6674222946 | -0.6674225375 | 2.43e-07 | 0.0001 | True |
| gmx | numpy-ir | total | 112.2328568 | 112.2330285 | 0.000172 | 0.000111 | False |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 1.36e-08 | 0.001 | True |
| openmm-top | openmm-ir | 6.94e-07 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 3  passed=False

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 80.82374198 | 80.82374198 | 0 | 1.81e-06 | True |
| numpy-ir | reference-itself | angle | 37.76224293 | 37.76224293 | 0 | 1.38e-06 | True |
| numpy-ir | reference-itself | proper | 5.331389591 | 5.331389591 | 0 | 1.05e-06 | True |
| numpy-ir | reference-itself | improper | 0.6711358119 | 0.6711358119 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | lj | -0.7553759441 | -0.7553759441 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | coulomb | 10.77591317 | 10.77591317 | 0 | 1.11e-06 | True |
| numpy-ir | reference-itself | lj14 | 12.78587264 | 12.78587264 | 0 | 1.13e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -0.6841122596 | -0.6841122596 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | total | 146.7108079 | 146.7108079 | 0 | 2.47e-06 | True |
| openmm-ir | numpy-ir | bond | 80.82374198 | 80.82374198 | 1.42e-14 | 1.81e-06 | True |
| openmm-ir | numpy-ir | angle | 37.76224293 | 37.76224293 | 7.11e-15 | 1.38e-06 | True |
| openmm-ir | numpy-ir | proper | 5.331389591 | 5.331389591 | 1.69e-14 | 1.05e-06 | True |
| openmm-ir | numpy-ir | improper | 0.6711358119 | 0.6711358119 | 9.99e-16 | 1.01e-06 | True |
| openmm-ir | numpy-ir | lj | -0.7553759441 | -0.7553759441 | 0 | 1.01e-06 | True |
| openmm-ir | numpy-ir | coulomb | 10.77591317 | 10.77591317 | 0 | 1.11e-06 | True |
| openmm-ir | numpy-ir | lj14 | 12.78587264 | 12.78587264 | 1.42e-14 | 1.13e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -0.6841122596 | -0.6841122596 | 3.55e-15 | 1.01e-06 | True |
| openmm-ir | numpy-ir | total | 146.7108079 | 146.7108079 | 2.84e-14 | 1.81e-06 | True |
| openmm-xml | numpy-ir | bond | 80.82374198 | 80.82374198 | 1.42e-14 | 1.81e-06 | True |
| openmm-xml | numpy-ir | angle | 37.76224293 | 37.76224293 | 7.11e-15 | 1.38e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 6.002525403 | 6.002525403 | 1.78e-14 | 1.06e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | 22.12229773 | 22.12229761 | 1.19e-07 | 1.22e-06 | True |
| openmm-xml | numpy-ir | total | 146.710808 | 146.7108079 | 1.19e-07 | 1.81e-06 | True |
| openmm-prmtop | numpy-ir | bond | 80.82374198 | 80.82374198 | 1.42e-14 | 0.000108 | True |
| openmm-prmtop | numpy-ir | angle | 37.76224293 | 37.76224293 | 0 | 0.000104 | True |
| openmm-prmtop | numpy-ir | torsion_all | 6.002525403 | 6.002525403 | 1.95e-14 | 0.000101 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | 22.12229772 | 22.12229761 | 1.15e-07 | 0.000102 | True |
| openmm-prmtop | numpy-ir | total | 146.710808 | 146.7108079 | 1.15e-07 | 0.000108 | True |
| openmm-top | numpy-ir | bond | 80.82374198 | 80.82374198 | 2.84e-14 | 0.000108 | True |
| openmm-top | numpy-ir | angle | 37.76224293 | 37.76224293 | 7.09e-09 | 0.000104 | True |
| openmm-top | numpy-ir | torsion_all | 6.002525403 | 6.002525403 | 3.85e-10 | 0.000101 | True |
| openmm-top | numpy-ir | nonbonded_all | 22.12229772 | 22.12229761 | 1.18e-07 | 0.000102 | True |
| openmm-top | numpy-ir | total | 146.710808 | 146.7108079 | 1.24e-07 | 0.000108 | True |
| sander-ours | numpy-ir | bond | 19.3173 | 19.31733795 | 3.79e-05 | 0.000202 | True |
| sander-ours | numpy-ir | angle | 9.0254 | 9.025392669 | 7.33e-06 | 0.000201 | True |
| sander-ours | numpy-ir | dihed | 1.4346 | 1.434638003 | 3.8e-05 | 0.0002 | True |
| sander-ours | numpy-ir | vdw | -0.1805 | -0.1805391836 | 3.92e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 2.5754 | 2.575505061 | 0.000105 | 0.000329 | True |
| sander-ours | numpy-ir | vdw14 | 3.0559 | 3.055896901 | 3.1e-06 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -0.1635 | -0.1635067542 | 6.75e-06 | 0.000208 | True |
| sander-tleap | sander-ours | bond | 19.3173 | 19.3173 | 0 | 0.000202 | True |
| sander-tleap | sander-ours | angle | 9.0254 | 9.0254 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | dihed | 1.4346 | 1.4346 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw | -0.1805 | -0.1805 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 2.5754 | 2.5754 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | 3.0559 | 3.0559 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -0.1635 | -0.1635 | 0 | 0.0002 | True |
| gmx | numpy-ir | bond | 80.82362366 | 80.82374198 | 0.000118 | 0.000108 | False |
| gmx | numpy-ir | angle | 37.76218414 | 37.76224293 | 5.88e-05 | 0.000104 | True |
| gmx | numpy-ir | proper | 5.331389427 | 5.331389591 | 1.64e-07 | 0.000101 | True |
| gmx | numpy-ir | improper | 0.6711349487 | 0.6711358119 | 8.63e-07 | 0.0001 | True |
| gmx | numpy-ir | lj | -0.7553762197 | -0.7553759441 | 2.76e-07 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 10.77591515 | 10.77591317 | 1.97e-06 | 0.000101 | True |
| gmx | numpy-ir | lj14 | 12.78587437 | 12.78587264 | 1.73e-06 | 0.000101 | True |
| gmx | numpy-ir | coulomb14 | -0.6841130257 | -0.6841122596 | 7.66e-07 | 0.0001 | True |
| gmx | numpy-ir | total | 146.7106323 | 146.7108079 | 0.000176 | 0.000115 | False |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 1.48e-08 | 0.001 | True |
| openmm-top | openmm-ir | 7.96e-07 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

## nma
### config 0  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 0.9501068381 | 0.9501068381 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | angle | 1.634761665 | 1.634761665 | 0 | 1.02e-06 | True |
| numpy-ir | reference-itself | proper | 104.2050736 | 104.2050736 | 0 | 2.04e-06 | True |
| numpy-ir | reference-itself | improper | 0.8172111871 | 0.8172111871 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | lj | -0.509404987 | -0.509404987 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | coulomb | -1.240306379 | -1.240306379 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | lj14 | 4.869462934 | 4.869462934 | 0 | 1.05e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -158.7006254 | -158.7006254 | 0 | 2.59e-06 | True |
| numpy-ir | reference-itself | total | -47.97372051 | -47.97372051 | 0 | 1.48e-06 | True |
| openmm-ir | numpy-ir | bond | 0.9501068381 | 0.9501068381 | 2.66e-15 | 1.01e-06 | True |
| openmm-ir | numpy-ir | angle | 1.634761665 | 1.634761665 | 5.77e-15 | 1.02e-06 | True |
| openmm-ir | numpy-ir | proper | 104.2050736 | 104.2050736 | 0 | 2.04e-06 | True |
| openmm-ir | numpy-ir | improper | 0.8172111871 | 0.8172111871 | 0 | 1.01e-06 | True |
| openmm-ir | numpy-ir | lj | -0.509404987 | -0.509404987 | 3.33e-16 | 1.01e-06 | True |
| openmm-ir | numpy-ir | coulomb | -1.240306379 | -1.240306379 | 0 | 1.01e-06 | True |
| openmm-ir | numpy-ir | lj14 | 4.869462934 | 4.869462934 | 0 | 1.05e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -158.7006254 | -158.7006254 | 2.84e-14 | 2.59e-06 | True |
| openmm-ir | numpy-ir | total | -47.97372051 | -47.97372051 | 2.13e-14 | 2.59e-06 | True |
| openmm-xml | numpy-ir | bond | 0.9501068381 | 0.9501068381 | 2.66e-15 | 1.01e-06 | True |
| openmm-xml | numpy-ir | angle | 1.634761665 | 1.634761665 | 5.77e-15 | 1.02e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 105.0222848 | 105.0222848 | 0 | 2.05e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | -155.5808757 | -155.5808738 | 1.89e-06 | 2.56e-06 | True |
| openmm-xml | numpy-ir | total | -47.9737224 | -47.97372051 | 1.89e-06 | 2.59e-06 | True |
| openmm-prmtop | numpy-ir | bond | 0.9501068381 | 0.9501068381 | 2.44e-15 | 0.0001 | True |
| openmm-prmtop | numpy-ir | angle | 1.634761665 | 1.634761665 | 9.33e-15 | 0.0001 | True |
| openmm-prmtop | numpy-ir | torsion_all | 105.0222848 | 105.0222848 | 1.42e-14 | 0.000111 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | -155.5808757 | -155.5808738 | 1.9e-06 | 0.000116 | True |
| openmm-prmtop | numpy-ir | total | -47.97372241 | -47.97372051 | 1.9e-06 | 0.000116 | True |
| openmm-top | numpy-ir | bond | 0.9501068381 | 0.9501068381 | 7.44e-15 | 0.0001 | True |
| openmm-top | numpy-ir | angle | 1.634761669 | 1.634761665 | 4.11e-09 | 0.0001 | True |
| openmm-top | numpy-ir | torsion_all | 105.0222848 | 105.0222848 | 4.89e-09 | 0.000111 | True |
| openmm-top | numpy-ir | nonbonded_all | -155.5808757 | -155.5808738 | 1.91e-06 | 0.000116 | True |
| openmm-top | numpy-ir | total | -47.97372241 | -47.97372051 | 1.9e-06 | 0.000116 | True |
| sander-ours | numpy-ir | bond | 0.2271 | 0.2270809843 | 1.9e-05 | 0.0002 | True |
| sander-ours | numpy-ir | angle | 0.3907 | 0.3907174152 | 1.74e-05 | 0.0002 | True |
| sander-ours | numpy-ir | dihed | 25.1009 | 25.10092849 | 2.85e-05 | 0.000203 | True |
| sander-ours | numpy-ir | vdw | -0.1218 | -0.1217507139 | 4.93e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel | -0.2964 | -0.2964403392 | 4.03e-05 | 0.000215 | True |
| sander-ours | numpy-ir | vdw14 | 1.1638 | 1.163829573 | 2.96e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -37.929 | -37.9303598 | 0.00136 | 0.0021 | True |
| sander-tleap | sander-ours | bond | 0.2271 | 0.2271 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | angle | 0.3907 | 0.3907 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | dihed | 25.1009 | 25.1009 | 0 | 0.000203 | True |
| sander-tleap | sander-ours | vdw | -0.1218 | -0.1218 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | -0.2964 | -0.2964 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | 1.1638 | 1.1638 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -37.929 | -37.929 | 0 | 0.000204 | True |
| gmx | numpy-ir | bond | 0.9501005411 | 0.9501068381 | 6.3e-06 | 0.0001 | True |
| gmx | numpy-ir | angle | 1.634764194 | 1.634761665 | 2.53e-06 | 0.0001 | True |
| gmx | numpy-ir | proper | 104.2050858 | 104.2050736 | 1.21e-05 | 0.00011 | True |
| gmx | numpy-ir | improper | 0.8172117472 | 0.8172111871 | 5.6e-07 | 0.0001 | True |
| gmx | numpy-ir | lj | -0.5094050169 | -0.509404987 | 2.99e-08 | 0.0001 | True |
| gmx | numpy-ir | coulomb | -1.240307808 | -1.240306379 | 1.43e-06 | 0.0001 | True |
| gmx | numpy-ir | lj14 | 4.869464397 | 4.869462934 | 1.46e-06 | 0.0001 | True |
| gmx | numpy-ir | coulomb14 | -158.7006378 | -158.7006254 | 1.24e-05 | 0.000116 | True |
| gmx | numpy-ir | total | -47.97372437 | -47.97372051 | 3.86e-06 | 0.000105 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 1.9e-07 | 0.001 | True |
| openmm-top | openmm-ir | 7.76e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 1  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 66.27263928 | 66.27263928 | 0 | 1.66e-06 | True |
| numpy-ir | reference-itself | angle | 18.35968129 | 18.35968129 | 0 | 1.18e-06 | True |
| numpy-ir | reference-itself | proper | 103.4230996 | 103.4230996 | 0 | 2.03e-06 | True |
| numpy-ir | reference-itself | improper | 0.3765435898 | 0.3765435898 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | -0.4676485325 | -0.4676485325 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb | -1.511585096 | -1.511585096 | 0 | 1.02e-06 | True |
| numpy-ir | reference-itself | lj14 | 2.941760682 | 2.941760682 | 0 | 1.03e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -154.9775136 | -154.9775136 | 0 | 2.55e-06 | True |
| numpy-ir | reference-itself | total | 34.41697727 | 34.41697727 | 0 | 1.34e-06 | True |
| openmm-ir | numpy-ir | bond | 66.27263928 | 66.27263928 | 0 | 1.66e-06 | True |
| openmm-ir | numpy-ir | angle | 18.35968129 | 18.35968129 | 1.07e-14 | 1.18e-06 | True |
| openmm-ir | numpy-ir | proper | 103.4230996 | 103.4230996 | 0 | 2.03e-06 | True |
| openmm-ir | numpy-ir | improper | 0.3765435898 | 0.3765435898 | 7.77e-15 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | -0.4676485325 | -0.4676485325 | 2.22e-16 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb | -1.511585096 | -1.511585096 | 0 | 1.02e-06 | True |
| openmm-ir | numpy-ir | lj14 | 2.941760682 | 2.941760682 | 4.44e-16 | 1.03e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -154.9775136 | -154.9775136 | 2.84e-14 | 2.55e-06 | True |
| openmm-ir | numpy-ir | total | 34.41697727 | 34.41697727 | 7.11e-15 | 2.55e-06 | True |
| openmm-xml | numpy-ir | bond | 66.27263928 | 66.27263928 | 0 | 1.66e-06 | True |
| openmm-xml | numpy-ir | angle | 18.35968129 | 18.35968129 | 1.07e-14 | 1.18e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 103.7996432 | 103.7996432 | 0 | 2.04e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | -154.0149884 | -154.0149865 | 1.85e-06 | 2.54e-06 | True |
| openmm-xml | numpy-ir | total | 34.41697542 | 34.41697727 | 1.85e-06 | 2.55e-06 | True |
| openmm-prmtop | numpy-ir | bond | 66.27263928 | 66.27263928 | 1.42e-14 | 0.000107 | True |
| openmm-prmtop | numpy-ir | angle | 18.35968129 | 18.35968129 | 0 | 0.000102 | True |
| openmm-prmtop | numpy-ir | torsion_all | 103.7996432 | 103.7996432 | 1.42e-14 | 0.00011 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | -154.0149884 | -154.0149865 | 1.86e-06 | 0.000115 | True |
| openmm-prmtop | numpy-ir | total | 34.41697541 | 34.41697727 | 1.86e-06 | 0.000115 | True |
| openmm-top | numpy-ir | bond | 66.27263928 | 66.27263928 | 1.42e-14 | 0.000107 | True |
| openmm-top | numpy-ir | angle | 18.35968129 | 18.35968129 | 2.99e-09 | 0.000102 | True |
| openmm-top | numpy-ir | torsion_all | 103.7996432 | 103.7996432 | 9.09e-09 | 0.00011 | True |
| openmm-top | numpy-ir | nonbonded_all | -154.0149884 | -154.0149865 | 1.87e-06 | 0.000115 | True |
| openmm-top | numpy-ir | total | 34.41697541 | 34.41697727 | 1.86e-06 | 0.000115 | True |
| sander-ours | numpy-ir | bond | 15.8395 | 15.83954094 | 4.09e-05 | 0.000202 | True |
| sander-ours | numpy-ir | angle | 4.3881 | 4.388069142 | 3.09e-05 | 0.0002 | True |
| sander-ours | numpy-ir | dihed | 24.8087 | 24.80871014 | 1.01e-05 | 0.000202 | True |
| sander-ours | numpy-ir | vdw | -0.1118 | -0.1117706818 | 2.93e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel | -0.3613 | -0.3612775087 | 2.25e-05 | 0.000218 | True |
| sander-ours | numpy-ir | vdw14 | 0.7031 | 0.7030976773 | 2.32e-06 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -37.0392 | -37.04051471 | 0.00131 | 0.00205 | True |
| sander-tleap | sander-ours | bond | 15.8395 | 15.8395 | 0 | 0.000202 | True |
| sander-tleap | sander-ours | angle | 4.3881 | 4.3881 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | dihed | 24.8087 | 24.8087 | 0 | 0.000202 | True |
| sander-tleap | sander-ours | vdw | -0.1118 | -0.1118 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | -0.3613 | -0.3613 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | 0.7031 | 0.7031 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -37.0392 | -37.0392 | 0 | 0.000204 | True |
| gmx | numpy-ir | bond | 66.27272034 | 66.27263928 | 8.11e-05 | 0.000107 | True |
| gmx | numpy-ir | angle | 18.35967827 | 18.35968129 | 3.02e-06 | 0.000102 | True |
| gmx | numpy-ir | proper | 103.423111 | 103.4230996 | 1.13e-05 | 0.00011 | True |
| gmx | numpy-ir | improper | 0.3765461147 | 0.3765435898 | 2.52e-06 | 0.0001 | True |
| gmx | numpy-ir | lj | -0.4676497281 | -0.4676485325 | 1.2e-06 | 0.0001 | True |
| gmx | numpy-ir | coulomb | -1.511584282 | -1.511585096 | 8.14e-07 | 0.0001 | True |
| gmx | numpy-ir | lj14 | 2.941751242 | 2.941760682 | 9.44e-06 | 0.0001 | True |
| gmx | numpy-ir | coulomb14 | -154.9775085 | -154.9775136 | 5.02e-06 | 0.000115 | True |
| gmx | numpy-ir | total | 34.41706467 | 34.41697727 | 8.74e-05 | 0.000103 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 1.44e-07 | 0.001 | True |
| openmm-top | openmm-ir | 7.32e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 2  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 56.89593871 | 56.89593871 | 0 | 1.57e-06 | True |
| numpy-ir | reference-itself | angle | 30.09037899 | 30.09037899 | 0 | 1.3e-06 | True |
| numpy-ir | reference-itself | proper | 106.190717 | 106.190717 | 0 | 2.06e-06 | True |
| numpy-ir | reference-itself | improper | 0.9684522793 | 0.9684522793 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | lj | -0.448652211 | -0.448652211 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb | -1.699268773 | -1.699268773 | 0 | 1.02e-06 | True |
| numpy-ir | reference-itself | lj14 | 5.668838474 | 5.668838474 | 0 | 1.06e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -158.6497147 | -158.6497147 | 0 | 2.59e-06 | True |
| numpy-ir | reference-itself | total | 39.01668973 | 39.01668973 | 0 | 1.39e-06 | True |
| openmm-ir | numpy-ir | bond | 56.89593871 | 56.89593871 | 4.97e-14 | 1.57e-06 | True |
| openmm-ir | numpy-ir | angle | 30.09037899 | 30.09037899 | 4.62e-14 | 1.3e-06 | True |
| openmm-ir | numpy-ir | proper | 106.190717 | 106.190717 | 1.42e-14 | 2.06e-06 | True |
| openmm-ir | numpy-ir | improper | 0.9684522793 | 0.9684522793 | 4.88e-15 | 1.01e-06 | True |
| openmm-ir | numpy-ir | lj | -0.448652211 | -0.448652211 | 1.67e-16 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb | -1.699268773 | -1.699268773 | 3.55e-15 | 1.02e-06 | True |
| openmm-ir | numpy-ir | lj14 | 5.668838474 | 5.668838474 | 2.66e-15 | 1.06e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -158.6497147 | -158.6497147 | 2.84e-14 | 2.59e-06 | True |
| openmm-ir | numpy-ir | total | 39.01668973 | 39.01668973 | 0 | 2.59e-06 | True |
| openmm-xml | numpy-ir | bond | 56.89593871 | 56.89593871 | 4.97e-14 | 1.57e-06 | True |
| openmm-xml | numpy-ir | angle | 30.09037899 | 30.09037899 | 4.62e-14 | 1.3e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 107.1591692 | 107.1591692 | 1.42e-14 | 2.07e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | -155.1287991 | -155.1287972 | 1.9e-06 | 2.55e-06 | True |
| openmm-xml | numpy-ir | total | 39.01668783 | 39.01668973 | 1.9e-06 | 2.59e-06 | True |
| openmm-prmtop | numpy-ir | bond | 56.89593871 | 56.89593871 | 4.97e-14 | 0.000106 | True |
| openmm-prmtop | numpy-ir | angle | 30.09037899 | 30.09037899 | 5.68e-14 | 0.000103 | True |
| openmm-prmtop | numpy-ir | torsion_all | 107.1591692 | 107.1591692 | 0 | 0.000111 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | -155.1287991 | -155.1287972 | 1.91e-06 | 0.000116 | True |
| openmm-prmtop | numpy-ir | total | 39.01668782 | 39.01668973 | 1.91e-06 | 0.000116 | True |
| openmm-top | numpy-ir | bond | 56.89593871 | 56.89593871 | 5.68e-14 | 0.000106 | True |
| openmm-top | numpy-ir | angle | 30.09037909 | 30.09037899 | 1.04e-07 | 0.000103 | True |
| openmm-top | numpy-ir | torsion_all | 107.1591692 | 107.1591692 | 2.13e-09 | 0.000111 | True |
| openmm-top | numpy-ir | nonbonded_all | -155.1287991 | -155.1287972 | 1.92e-06 | 0.000116 | True |
| openmm-top | numpy-ir | total | 39.01668792 | 39.01668973 | 1.81e-06 | 0.000116 | True |
| sander-ours | numpy-ir | bond | 13.5985 | 13.59845571 | 4.43e-05 | 0.000201 | True |
| sander-ours | numpy-ir | angle | 7.1918 | 7.191773181 | 2.68e-05 | 0.000201 | True |
| sander-ours | numpy-ir | dihed | 25.6117 | 25.61165613 | 4.39e-05 | 0.000203 | True |
| sander-ours | numpy-ir | vdw | -0.1072 | -0.107230452 | 3.05e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel | -0.4061 | -0.4061349839 | 3.5e-05 | 0.00022 | True |
| sander-ours | numpy-ir | vdw14 | 1.3549 | 1.354884912 | 1.51e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -37.9169 | -37.91819185 | 0.00129 | 0.0021 | True |
| sander-tleap | sander-ours | bond | 13.5985 | 13.5985 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | angle | 7.1918 | 7.1918 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | dihed | 25.6117 | 25.6117 | 0 | 0.000203 | True |
| sander-tleap | sander-ours | vdw | -0.1072 | -0.1072 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | -0.4061 | -0.4061 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | 1.3549 | 1.3549 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -37.9169 | -37.9169 | 0 | 0.000204 | True |
| gmx | numpy-ir | bond | 56.89601898 | 56.89593871 | 8.03e-05 | 0.000106 | True |
| gmx | numpy-ir | angle | 30.09037971 | 30.09037899 | 7.26e-07 | 0.000103 | True |
| gmx | numpy-ir | proper | 106.190712 | 106.190717 | 4.99e-06 | 0.000111 | True |
| gmx | numpy-ir | improper | 0.9684525132 | 0.9684522793 | 2.34e-07 | 0.0001 | True |
| gmx | numpy-ir | lj | -0.4486527443 | -0.448652211 | 5.33e-07 | 0.0001 | True |
| gmx | numpy-ir | coulomb | -1.699271202 | -1.699268773 | 2.43e-06 | 0.0001 | True |
| gmx | numpy-ir | lj14 | 5.668827057 | 5.668838474 | 1.14e-05 | 0.000101 | True |
| gmx | numpy-ir | coulomb14 | -158.649704 | -158.6497147 | 1.07e-05 | 0.000116 | True |
| gmx | numpy-ir | total | 39.01676178 | 39.01668973 | 7.21e-05 | 0.000104 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 2.13e-07 | 0.001 | True |
| openmm-top | openmm-ir | 8.18e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 3  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 57.70295688 | 57.70295688 | 0 | 1.58e-06 | True |
| numpy-ir | reference-itself | angle | 19.24804174 | 19.24804174 | 0 | 1.19e-06 | True |
| numpy-ir | reference-itself | proper | 104.28802 | 104.28802 | 0 | 2.04e-06 | True |
| numpy-ir | reference-itself | improper | 2.423971464 | 2.423971464 | 0 | 1.02e-06 | True |
| numpy-ir | reference-itself | lj | -0.5264760139 | -0.5264760139 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | coulomb | -1.411356065 | -1.411356065 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | lj14 | 4.586631092 | 4.586631092 | 0 | 1.05e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -156.1945843 | -156.1945843 | 0 | 2.56e-06 | True |
| numpy-ir | reference-itself | total | 30.11720479 | 30.11720479 | 0 | 1.3e-06 | True |
| openmm-ir | numpy-ir | bond | 57.70295688 | 57.70295688 | 7.11e-14 | 1.58e-06 | True |
| openmm-ir | numpy-ir | angle | 19.24804174 | 19.24804174 | 2.84e-14 | 1.19e-06 | True |
| openmm-ir | numpy-ir | proper | 104.28802 | 104.28802 | 2.84e-14 | 2.04e-06 | True |
| openmm-ir | numpy-ir | improper | 2.423971464 | 2.423971464 | 3.55e-15 | 1.02e-06 | True |
| openmm-ir | numpy-ir | lj | -0.5264760139 | -0.5264760139 | 2.78e-15 | 1.01e-06 | True |
| openmm-ir | numpy-ir | coulomb | -1.411356065 | -1.411356065 | 0 | 1.01e-06 | True |
| openmm-ir | numpy-ir | lj14 | 4.586631092 | 4.586631092 | 0 | 1.05e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -156.1945843 | -156.1945843 | 0 | 2.56e-06 | True |
| openmm-ir | numpy-ir | total | 30.11720479 | 30.11720479 | 2.13e-14 | 2.56e-06 | True |
| openmm-xml | numpy-ir | bond | 57.70295688 | 57.70295688 | 7.11e-14 | 1.58e-06 | True |
| openmm-xml | numpy-ir | angle | 19.24804174 | 19.24804174 | 2.84e-14 | 1.19e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 106.7119915 | 106.7119915 | 2.84e-14 | 2.07e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | -153.5457872 | -153.5457853 | 1.87e-06 | 2.54e-06 | True |
| openmm-xml | numpy-ir | total | 30.11720293 | 30.11720479 | 1.87e-06 | 2.56e-06 | True |
| openmm-prmtop | numpy-ir | bond | 57.70295688 | 57.70295688 | 7.11e-14 | 0.000106 | True |
| openmm-prmtop | numpy-ir | angle | 19.24804174 | 19.24804174 | 2.13e-14 | 0.000102 | True |
| openmm-prmtop | numpy-ir | torsion_all | 106.7119915 | 106.7119915 | 0 | 0.000111 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | -153.5457872 | -153.5457853 | 1.87e-06 | 0.000115 | True |
| openmm-prmtop | numpy-ir | total | 30.11720292 | 30.11720479 | 1.87e-06 | 0.000116 | True |
| openmm-top | numpy-ir | bond | 57.70295688 | 57.70295688 | 1.99e-13 | 0.000106 | True |
| openmm-top | numpy-ir | angle | 19.24804182 | 19.24804174 | 8.64e-08 | 0.000102 | True |
| openmm-top | numpy-ir | torsion_all | 106.7119915 | 106.7119915 | 1.54e-08 | 0.000111 | True |
| openmm-top | numpy-ir | nonbonded_all | -153.5457872 | -153.5457853 | 1.89e-06 | 0.000115 | True |
| openmm-top | numpy-ir | total | 30.11720301 | 30.11720479 | 1.78e-06 | 0.000116 | True |
| sander-ours | numpy-ir | bond | 13.7913 | 13.79133769 | 3.77e-05 | 0.000201 | True |
| sander-ours | numpy-ir | angle | 4.6004 | 4.600392384 | 7.62e-06 | 0.0002 | True |
| sander-ours | numpy-ir | dihed | 25.5048 | 25.50477808 | 2.19e-05 | 0.000203 | True |
| sander-ours | numpy-ir | vdw | -0.1258 | -0.1258307873 | 3.08e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel | -0.3373 | -0.3373221952 | 2.22e-05 | 0.000217 | True |
| sander-ours | numpy-ir | vdw14 | 1.0962 | 1.09623114 | 3.11e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -37.3301 | -37.33140161 | 0.0013 | 0.00207 | True |
| sander-tleap | sander-ours | bond | 13.7913 | 13.7913 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | angle | 4.6004 | 4.6004 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | dihed | 25.5048 | 25.5048 | 0 | 0.000203 | True |
| sander-tleap | sander-ours | vdw | -0.1258 | -0.1258 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | -0.3373 | -0.3373 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | 1.0962 | 1.0962 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -37.3301 | -37.3301 | 0 | 0.000204 | True |
| gmx | numpy-ir | bond | 57.70298004 | 57.70295688 | 2.32e-05 | 0.000106 | True |
| gmx | numpy-ir | angle | 19.24802399 | 19.24804174 | 1.77e-05 | 0.000102 | True |
| gmx | numpy-ir | proper | 104.2880173 | 104.28802 | 2.75e-06 | 0.00011 | True |
| gmx | numpy-ir | improper | 2.423959732 | 2.423971464 | 1.17e-05 | 0.0001 | True |
| gmx | numpy-ir | lj | -0.5264770985 | -0.5264760139 | 1.08e-06 | 0.0001 | True |
| gmx | numpy-ir | coulomb | -1.411353111 | -1.411356065 | 2.95e-06 | 0.0001 | True |
| gmx | numpy-ir | lj14 | 4.586625576 | 4.586631092 | 5.52e-06 | 0.0001 | True |
| gmx | numpy-ir | coulomb14 | -156.1945953 | -156.1945843 | 1.1e-05 | 0.000116 | True |
| gmx | numpy-ir | total | 30.11717606 | 30.11720479 | 2.87e-05 | 0.000103 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 1.69e-07 | 0.001 | True |
| openmm-top | openmm-ir | 8.26e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

## cyclopropane
### config 0  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 0.1698281495 | 0.1698281495 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | angle | 0.062543369 | 0.062543369 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | proper | 26.97107427 | 26.97107427 | 0 | 1.27e-06 | True |
| numpy-ir | reference-itself | improper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj14 | 0.02599070105 | 0.02599070105 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb14 | 23.32701892 | 23.32701892 | 0 | 1.23e-06 | True |
| numpy-ir | reference-itself | total | 50.55645541 | 50.55645541 | 0 | 1.51e-06 | True |
| openmm-ir | numpy-ir | bond | 0.1698281495 | 0.1698281495 | 3.05e-16 | 1e-06 | True |
| openmm-ir | numpy-ir | angle | 0.062543369 | 0.062543369 | 8.33e-16 | 1e-06 | True |
| openmm-ir | numpy-ir | proper | 26.97107427 | 26.97107427 | 3.55e-15 | 1.27e-06 | True |
| openmm-ir | numpy-ir | improper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj14 | 0.02599070105 | 0.02599070105 | 1.18e-16 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | 23.32701892 | 23.32701892 | 0 | 1.23e-06 | True |
| openmm-ir | numpy-ir | total | 50.55645541 | 50.55645541 | 7.11e-15 | 1.27e-06 | True |
| openmm-xml | numpy-ir | bond | 0.1698281495 | 0.1698281495 | 3.05e-16 | 1e-06 | True |
| openmm-xml | numpy-ir | angle | 0.062543369 | 0.062543369 | 8.33e-16 | 1e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 26.97107427 | 26.97107427 | 3.55e-15 | 1.27e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | 23.3530099 | 23.35300962 | 2.76e-07 | 1.23e-06 | True |
| openmm-xml | numpy-ir | total | 50.55645568 | 50.55645541 | 2.76e-07 | 1.27e-06 | True |
| openmm-prmtop | numpy-ir | bond | 0.1698281495 | 0.1698281495 | 3.05e-16 | 0.0001 | True |
| openmm-prmtop | numpy-ir | angle | 0.062543369 | 0.062543369 | 8.33e-16 | 0.0001 | True |
| openmm-prmtop | numpy-ir | torsion_all | 26.97107427 | 26.97107427 | 3.55e-15 | 0.000103 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | 23.3530099 | 23.35300962 | 2.76e-07 | 0.000102 | True |
| openmm-prmtop | numpy-ir | total | 50.55645568 | 50.55645541 | 2.76e-07 | 0.000103 | True |
| openmm-top | numpy-ir | bond | 0.1698281495 | 0.1698281495 | 6.69e-15 | 0.0001 | True |
| openmm-top | numpy-ir | angle | 0.06254337336 | 0.062543369 | 4.37e-09 | 0.0001 | True |
| openmm-top | numpy-ir | torsion_all | 26.97107427 | 26.97107427 | 1.66e-10 | 0.000103 | True |
| openmm-top | numpy-ir | nonbonded_all | 23.35300992 | 23.35300962 | 2.98e-07 | 0.000102 | True |
| openmm-top | numpy-ir | total | 50.55645571 | 50.55645541 | 3.02e-07 | 0.000103 | True |
| sander-ours | numpy-ir | bond | 0.0406 | 0.04058990188 | 1.01e-05 | 0.0002 | True |
| sander-ours | numpy-ir | angle | 0.0149 | 0.01494822395 | 4.82e-05 | 0.0002 | True |
| sander-ours | numpy-ir | dihed | 6.4462 | 6.446241459 | 4.15e-05 | 0.000201 | True |
| sander-ours | numpy-ir | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | vdw14 | 0.0062 | 0.006211926637 | 1.19e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | 5.5751 | 5.575291329 | 0.000191 | 0.000479 | True |
| sander-tleap | sander-ours | bond | 0.0406 | 0.0406 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | angle | 0.0149 | 0.0149 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | dihed | 6.4462 | 6.4462 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | 0.0062 | 0.0062 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | 5.5751 | 5.5751 | 0 | 0.000201 | True |
| gmx | numpy-ir | bond | 0.1698301882 | 0.1698281495 | 2.04e-06 | 0.0001 | True |
| gmx | numpy-ir | angle | 0.06254481524 | 0.062543369 | 1.45e-06 | 0.0001 | True |
| gmx | numpy-ir | proper | 26.97107315 | 26.97107427 | 1.11e-06 | 0.000103 | True |
| gmx | numpy-ir | improper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj14 | 0.02599062771 | 0.02599070105 | 7.33e-08 | 0.0001 | True |
| gmx | numpy-ir | coulomb14 | 23.32702065 | 23.32701892 | 1.72e-06 | 0.000102 | True |
| gmx | numpy-ir | total | 50.55645752 | 50.55645541 | 2.11e-06 | 0.000105 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 2.84e-14 | 0.001 | True |
| openmm-top | openmm-ir | 1.8e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 1  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 62.14399434 | 62.14399434 | 0 | 1.62e-06 | True |
| numpy-ir | reference-itself | angle | 16.54110787 | 16.54110787 | 0 | 1.17e-06 | True |
| numpy-ir | reference-itself | proper | 26.39153906 | 26.39153906 | 0 | 1.26e-06 | True |
| numpy-ir | reference-itself | improper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj14 | 0.2320824432 | 0.2320824432 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb14 | 23.76271771 | 23.76271771 | 0 | 1.24e-06 | True |
| numpy-ir | reference-itself | total | 129.0714414 | 129.0714414 | 0 | 2.29e-06 | True |
| openmm-ir | numpy-ir | bond | 62.14399434 | 62.14399434 | 0 | 1.62e-06 | True |
| openmm-ir | numpy-ir | angle | 16.54110787 | 16.54110787 | 1.42e-14 | 1.17e-06 | True |
| openmm-ir | numpy-ir | proper | 26.39153906 | 26.39153906 | 0 | 1.26e-06 | True |
| openmm-ir | numpy-ir | improper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj14 | 0.2320824432 | 0.2320824432 | 1.67e-16 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | 23.76271771 | 23.76271771 | 7.11e-15 | 1.24e-06 | True |
| openmm-ir | numpy-ir | total | 129.0714414 | 129.0714414 | 2.84e-14 | 1.62e-06 | True |
| openmm-xml | numpy-ir | bond | 62.14399434 | 62.14399434 | 0 | 1.62e-06 | True |
| openmm-xml | numpy-ir | angle | 16.54110787 | 16.54110787 | 1.42e-14 | 1.17e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 26.39153906 | 26.39153906 | 0 | 1.26e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | 23.99480043 | 23.99480015 | 2.81e-07 | 1.24e-06 | True |
| openmm-xml | numpy-ir | total | 129.0714417 | 129.0714414 | 2.81e-07 | 1.62e-06 | True |
| openmm-prmtop | numpy-ir | bond | 62.14399434 | 62.14399434 | 7.11e-15 | 0.000106 | True |
| openmm-prmtop | numpy-ir | angle | 16.54110787 | 16.54110787 | 1.42e-14 | 0.000102 | True |
| openmm-prmtop | numpy-ir | torsion_all | 26.39153906 | 26.39153906 | 0 | 0.000103 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | 23.99480043 | 23.99480015 | 2.81e-07 | 0.000102 | True |
| openmm-prmtop | numpy-ir | total | 129.0714417 | 129.0714414 | 2.81e-07 | 0.000106 | True |
| openmm-top | numpy-ir | bond | 62.14399434 | 62.14399434 | 9.95e-14 | 0.000106 | True |
| openmm-top | numpy-ir | angle | 16.54110786 | 16.54110787 | 5.25e-09 | 0.000102 | True |
| openmm-top | numpy-ir | torsion_all | 26.39153906 | 26.39153906 | 1.62e-10 | 0.000103 | True |
| openmm-top | numpy-ir | nonbonded_all | 23.99480045 | 23.99480015 | 3.03e-07 | 0.000102 | True |
| openmm-top | numpy-ir | total | 129.0714417 | 129.0714414 | 2.98e-07 | 0.000106 | True |
| sander-ours | numpy-ir | bond | 14.8528 | 14.85277111 | 2.89e-05 | 0.000201 | True |
| sander-ours | numpy-ir | angle | 3.9534 | 3.953419663 | 1.97e-05 | 0.0002 | True |
| sander-ours | numpy-ir | dihed | 6.3077 | 6.30772922 | 2.92e-05 | 0.000201 | True |
| sander-ours | numpy-ir | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | vdw14 | 0.0555 | 0.05546903519 | 3.1e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | 5.6792 | 5.679425838 | 0.000226 | 0.000484 | True |
| sander-tleap | sander-ours | bond | 14.8528 | 14.8528 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | angle | 3.9534 | 3.9534 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | dihed | 6.3077 | 6.3077 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | 0.0555 | 0.0555 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | 5.6792 | 5.6792 | 0 | 0.000201 | True |
| gmx | numpy-ir | bond | 62.14395523 | 62.14399434 | 3.91e-05 | 0.000106 | True |
| gmx | numpy-ir | angle | 16.54111099 | 16.54110787 | 3.12e-06 | 0.000102 | True |
| gmx | numpy-ir | proper | 26.39154053 | 26.39153906 | 1.47e-06 | 0.000103 | True |
| gmx | numpy-ir | improper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj14 | 0.2320820689 | 0.2320824432 | 3.74e-07 | 0.0001 | True |
| gmx | numpy-ir | coulomb14 | 23.76272011 | 23.76271771 | 2.4e-06 | 0.000102 | True |
| gmx | numpy-ir | total | 129.0714111 | 129.0714414 | 3.03e-05 | 0.000113 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 9.09e-13 | 0.001 | True |
| openmm-top | openmm-ir | 1.88e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 2  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 34.51630139 | 34.51630139 | 0 | 1.35e-06 | True |
| numpy-ir | reference-itself | angle | 24.68405195 | 24.68405195 | 0 | 1.25e-06 | True |
| numpy-ir | reference-itself | proper | 26.42852355 | 26.42852355 | 0 | 1.26e-06 | True |
| numpy-ir | reference-itself | improper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj14 | 0.214469467 | 0.214469467 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb14 | 23.48557173 | 23.48557173 | 0 | 1.23e-06 | True |
| numpy-ir | reference-itself | total | 109.3289181 | 109.3289181 | 0 | 2.09e-06 | True |
| openmm-ir | numpy-ir | bond | 34.51630139 | 34.51630139 | 1.42e-14 | 1.35e-06 | True |
| openmm-ir | numpy-ir | angle | 24.68405195 | 24.68405195 | 2.13e-14 | 1.25e-06 | True |
| openmm-ir | numpy-ir | proper | 26.42852355 | 26.42852355 | 3.55e-15 | 1.26e-06 | True |
| openmm-ir | numpy-ir | improper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj14 | 0.214469467 | 0.214469467 | 5e-16 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | 23.48557173 | 23.48557173 | 3.55e-15 | 1.23e-06 | True |
| openmm-ir | numpy-ir | total | 109.3289181 | 109.3289181 | 2.84e-14 | 1.35e-06 | True |
| openmm-xml | numpy-ir | bond | 34.51630139 | 34.51630139 | 1.42e-14 | 1.35e-06 | True |
| openmm-xml | numpy-ir | angle | 24.68405195 | 24.68405195 | 2.13e-14 | 1.25e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 26.42852355 | 26.42852355 | 3.55e-15 | 1.26e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | 23.70004148 | 23.7000412 | 2.78e-07 | 1.24e-06 | True |
| openmm-xml | numpy-ir | total | 109.3289184 | 109.3289181 | 2.78e-07 | 1.35e-06 | True |
| openmm-prmtop | numpy-ir | bond | 34.51630139 | 34.51630139 | 7.11e-15 | 0.000103 | True |
| openmm-prmtop | numpy-ir | angle | 24.68405195 | 24.68405195 | 2.13e-14 | 0.000102 | True |
| openmm-prmtop | numpy-ir | torsion_all | 26.42852355 | 26.42852355 | 3.55e-15 | 0.000103 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | 23.70004148 | 23.7000412 | 2.78e-07 | 0.000102 | True |
| openmm-prmtop | numpy-ir | total | 109.3289184 | 109.3289181 | 2.78e-07 | 0.000103 | True |
| openmm-top | numpy-ir | bond | 34.51630139 | 34.51630139 | 1.42e-14 | 0.000103 | True |
| openmm-top | numpy-ir | angle | 24.68405196 | 24.68405195 | 3.22e-09 | 0.000102 | True |
| openmm-top | numpy-ir | torsion_all | 26.42852355 | 26.42852355 | 1.62e-10 | 0.000103 | True |
| openmm-top | numpy-ir | nonbonded_all | 23.7000415 | 23.7000412 | 3e-07 | 0.000102 | True |
| openmm-top | numpy-ir | total | 109.3289184 | 109.3289181 | 3.03e-07 | 0.000103 | True |
| sander-ours | numpy-ir | bond | 8.2496 | 8.249594023 | 5.98e-06 | 0.000201 | True |
| sander-ours | numpy-ir | angle | 5.8996 | 5.899630008 | 3e-05 | 0.000201 | True |
| sander-ours | numpy-ir | dihed | 6.3166 | 6.316568727 | 3.13e-05 | 0.000201 | True |
| sander-ours | numpy-ir | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | vdw14 | 0.0513 | 0.05125943283 | 4.06e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | 5.613 | 5.61318636 | 0.000186 | 0.000481 | True |
| sander-tleap | sander-ours | bond | 8.2496 | 8.2496 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | angle | 5.8996 | 5.8996 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | dihed | 6.3166 | 6.3166 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | 0.0513 | 0.0513 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | 5.613 | 5.613 | 0 | 0.000201 | True |
| gmx | numpy-ir | bond | 34.51631165 | 34.51630139 | 1.03e-05 | 0.000103 | True |
| gmx | numpy-ir | angle | 24.68403435 | 24.68405195 | 1.76e-05 | 0.000102 | True |
| gmx | numpy-ir | proper | 26.42852592 | 26.42852355 | 2.37e-06 | 0.000103 | True |
| gmx | numpy-ir | improper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj14 | 0.2144681662 | 0.214469467 | 1.3e-06 | 0.0001 | True |
| gmx | numpy-ir | coulomb14 | 23.48557663 | 23.48557173 | 4.9e-06 | 0.000102 | True |
| gmx | numpy-ir | total | 109.3289185 | 109.3289181 | 3.6e-07 | 0.000111 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 9.09e-13 | 0.001 | True |
| openmm-top | openmm-ir | 1.99e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 3  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 36.90279134 | 36.90279134 | 0 | 1.37e-06 | True |
| numpy-ir | reference-itself | angle | 16.05337891 | 16.05337891 | 0 | 1.16e-06 | True |
| numpy-ir | reference-itself | proper | 26.6779418 | 26.6779418 | 0 | 1.27e-06 | True |
| numpy-ir | reference-itself | improper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj14 | 0.009424188908 | 0.009424188908 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb14 | 23.07031128 | 23.07031128 | 0 | 1.23e-06 | True |
| numpy-ir | reference-itself | total | 102.7138475 | 102.7138475 | 0 | 2.03e-06 | True |
| openmm-ir | numpy-ir | bond | 36.90279134 | 36.90279134 | 7.11e-15 | 1.37e-06 | True |
| openmm-ir | numpy-ir | angle | 16.05337891 | 16.05337891 | 2.13e-14 | 1.16e-06 | True |
| openmm-ir | numpy-ir | proper | 26.6779418 | 26.6779418 | 3.55e-15 | 1.27e-06 | True |
| openmm-ir | numpy-ir | improper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj14 | 0.009424188908 | 0.009424188908 | 5.69e-16 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | 23.07031128 | 23.07031128 | 3.55e-15 | 1.23e-06 | True |
| openmm-ir | numpy-ir | total | 102.7138475 | 102.7138475 | 2.84e-14 | 1.37e-06 | True |
| openmm-xml | numpy-ir | bond | 36.90279134 | 36.90279134 | 7.11e-15 | 1.37e-06 | True |
| openmm-xml | numpy-ir | angle | 16.05337891 | 16.05337891 | 2.13e-14 | 1.16e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 26.6779418 | 26.6779418 | 3.55e-15 | 1.27e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | 23.07973574 | 23.07973547 | 2.73e-07 | 1.23e-06 | True |
| openmm-xml | numpy-ir | total | 102.7138478 | 102.7138475 | 2.73e-07 | 1.37e-06 | True |
| openmm-prmtop | numpy-ir | bond | 36.90279134 | 36.90279134 | 1.42e-14 | 0.000104 | True |
| openmm-prmtop | numpy-ir | angle | 16.05337891 | 16.05337891 | 2.13e-14 | 0.000102 | True |
| openmm-prmtop | numpy-ir | torsion_all | 26.6779418 | 26.6779418 | 3.55e-15 | 0.000103 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | 23.07973574 | 23.07973547 | 2.73e-07 | 0.000102 | True |
| openmm-prmtop | numpy-ir | total | 102.7138478 | 102.7138475 | 2.73e-07 | 0.000104 | True |
| openmm-top | numpy-ir | bond | 36.90279134 | 36.90279134 | 6.39e-14 | 0.000104 | True |
| openmm-top | numpy-ir | angle | 16.05337892 | 16.05337891 | 6.51e-09 | 0.000102 | True |
| openmm-top | numpy-ir | torsion_all | 26.6779418 | 26.6779418 | 1.64e-10 | 0.000103 | True |
| openmm-top | numpy-ir | nonbonded_all | 23.07973577 | 23.07973547 | 2.94e-07 | 0.000102 | True |
| openmm-top | numpy-ir | total | 102.7138478 | 102.7138475 | 3.01e-07 | 0.000104 | True |
| sander-ours | numpy-ir | bond | 8.82 | 8.81997881 | 2.12e-05 | 0.000201 | True |
| sander-ours | numpy-ir | angle | 3.8368 | 3.836849644 | 4.96e-05 | 0.0002 | True |
| sander-ours | numpy-ir | dihed | 6.3762 | 6.376181119 | 1.89e-05 | 0.000201 | True |
| sander-ours | numpy-ir | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | vdw14 | 0.0023 | 0.002252435207 | 4.76e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | 5.5137 | 5.513936731 | 0.000237 | 0.000476 | True |
| sander-tleap | sander-ours | bond | 8.82 | 8.82 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | angle | 3.8368 | 3.8368 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | dihed | 6.3762 | 6.3762 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | 0.0023 | 0.0023 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | 5.5137 | 5.5137 | 0 | 0.000201 | True |
| gmx | numpy-ir | bond | 36.9027977 | 36.90279134 | 6.36e-06 | 0.000104 | True |
| gmx | numpy-ir | angle | 16.05339241 | 16.05337891 | 1.35e-05 | 0.000102 | True |
| gmx | numpy-ir | proper | 26.67794228 | 26.6779418 | 4.72e-07 | 0.000103 | True |
| gmx | numpy-ir | improper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj14 | 0.00942396 | 0.009424188908 | 2.29e-07 | 0.0001 | True |
| gmx | numpy-ir | coulomb14 | 23.07031441 | 23.07031128 | 3.12e-06 | 0.000102 | True |
| gmx | numpy-ir | total | 102.7138672 | 102.7138475 | 1.97e-05 | 0.00011 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 4.55e-13 | 0.001 | True |
| openmm-top | openmm-ir | 1.75e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

## cyclobutane
### config 0  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 0.8492151726 | 0.8492151726 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | angle | 1.479556376 | 1.479556376 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | proper | 27.1562796 | 27.1562796 | 0 | 1.27e-06 | True |
| numpy-ir | reference-itself | improper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | -0.3284134243 | -0.3284134243 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb | 7.225466808 | 7.225466808 | 0 | 1.07e-06 | True |
| numpy-ir | reference-itself | lj14 | 3.62935762 | 3.62935762 | 0 | 1.04e-06 | True |
| numpy-ir | reference-itself | coulomb14 | 0.5983036016 | 0.5983036016 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | total | 40.60976575 | 40.60976575 | 0 | 1.41e-06 | True |
| openmm-ir | numpy-ir | bond | 0.8492151726 | 0.8492151726 | 6.66e-16 | 1.01e-06 | True |
| openmm-ir | numpy-ir | angle | 1.479556376 | 1.479556376 | 1.55e-15 | 1.01e-06 | True |
| openmm-ir | numpy-ir | proper | 27.1562796 | 27.1562796 | 0 | 1.27e-06 | True |
| openmm-ir | numpy-ir | improper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | -0.3284134243 | -0.3284134243 | 1.11e-16 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb | 7.225466808 | 7.225466808 | 0 | 1.07e-06 | True |
| openmm-ir | numpy-ir | lj14 | 3.62935762 | 3.62935762 | 1.33e-15 | 1.04e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | 0.5983036016 | 0.5983036016 | 2.89e-15 | 1.01e-06 | True |
| openmm-ir | numpy-ir | total | 40.60976575 | 40.60976575 | 0 | 1.27e-06 | True |
| openmm-xml | numpy-ir | bond | 0.8492151726 | 0.8492151726 | 6.66e-16 | 1.01e-06 | True |
| openmm-xml | numpy-ir | angle | 1.479556376 | 1.479556376 | 1.55e-15 | 1.01e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 27.1562796 | 27.1562796 | 0 | 1.27e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | 11.1247147 | 11.12471461 | 9.26e-08 | 1.11e-06 | True |
| openmm-xml | numpy-ir | total | 40.60976584 | 40.60976575 | 9.26e-08 | 1.27e-06 | True |
| openmm-prmtop | numpy-ir | bond | 0.8492151726 | 0.8492151726 | 5.55e-16 | 0.0001 | True |
| openmm-prmtop | numpy-ir | angle | 1.479556376 | 1.479556376 | 2e-15 | 0.0001 | True |
| openmm-prmtop | numpy-ir | torsion_all | 27.1562796 | 27.1562796 | 3.55e-15 | 0.000103 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | 11.12471469 | 11.12471461 | 8.12e-08 | 0.000101 | True |
| openmm-prmtop | numpy-ir | total | 40.60976583 | 40.60976575 | 8.12e-08 | 0.000103 | True |
| openmm-top | numpy-ir | bond | 0.8492151726 | 0.8492151726 | 3.19e-14 | 0.0001 | True |
| openmm-top | numpy-ir | angle | 1.47955638 | 1.479556376 | 3.88e-09 | 0.0001 | True |
| openmm-top | numpy-ir | torsion_all | 27.1562796 | 27.1562796 | 1.67e-10 | 0.000103 | True |
| openmm-top | numpy-ir | nonbonded_all | 11.12471468 | 11.12471461 | 7.95e-08 | 0.000101 | True |
| openmm-top | numpy-ir | total | 40.60976583 | 40.60976575 | 8.32e-08 | 0.000103 | True |
| sander-ours | numpy-ir | bond | 0.203 | 0.2029672975 | 3.27e-05 | 0.0002 | True |
| sander-ours | numpy-ir | angle | 0.3536 | 0.3536224608 | 2.25e-05 | 0.0002 | True |
| sander-ours | numpy-ir | dihed | 6.4905 | 6.490506596 | 6.6e-06 | 0.000201 | True |
| sander-ours | numpy-ir | vdw | -0.0785 | -0.07849269224 | 7.31e-06 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 1.7269 | 1.726928013 | 2.8e-05 | 0.000286 | True |
| sander-ours | numpy-ir | vdw14 | 0.8674 | 0.8674372897 | 3.73e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | 0.143 | 0.1429979927 | 2.01e-06 | 0.000207 | True |
| sander-tleap | sander-ours | bond | 0.203 | 0.203 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | angle | 0.3536 | 0.3536 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | dihed | 6.4905 | 6.4905 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | vdw | -0.0785 | -0.0785 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 1.7269 | 1.7269 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | 0.8674 | 0.8674 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | 0.143 | 0.143 | 0 | 0.0002 | True |
| gmx | numpy-ir | bond | 0.8492215872 | 0.8492151726 | 6.41e-06 | 0.0001 | True |
| gmx | numpy-ir | angle | 1.479562044 | 1.479556376 | 5.67e-06 | 0.0001 | True |
| gmx | numpy-ir | proper | 27.15628433 | 27.1562796 | 4.74e-06 | 0.000103 | True |
| gmx | numpy-ir | improper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj | -0.3284133673 | -0.3284134243 | 5.7e-08 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 7.225466728 | 7.225466808 | 7.96e-08 | 0.000101 | True |
| gmx | numpy-ir | lj14 | 3.629358053 | 3.62935762 | 4.33e-07 | 0.0001 | True |
| gmx | numpy-ir | coulomb14 | 0.598300755 | 0.5983036016 | 2.85e-06 | 0.0001 | True |
| gmx | numpy-ir | total | 40.60977936 | 40.60976575 | 1.36e-05 | 0.000104 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 5.23e-08 | 0.001 | True |
| openmm-top | openmm-ir | 3.38e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 1  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 115.2557076 | 115.2557076 | 0 | 2.15e-06 | True |
| numpy-ir | reference-itself | angle | 75.38676741 | 75.38676741 | 0 | 1.75e-06 | True |
| numpy-ir | reference-itself | proper | 25.92207054 | 25.92207054 | 0 | 1.26e-06 | True |
| numpy-ir | reference-itself | improper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | -0.2949938265 | -0.2949938265 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb | 7.251141814 | 7.251141814 | 0 | 1.07e-06 | True |
| numpy-ir | reference-itself | lj14 | 4.241364517 | 4.241364517 | 0 | 1.04e-06 | True |
| numpy-ir | reference-itself | coulomb14 | 0.6083074285 | 0.6083074285 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | total | 228.3703655 | 228.3703655 | 0 | 3.28e-06 | True |
| openmm-ir | numpy-ir | bond | 115.2557076 | 115.2557076 | 2.84e-14 | 2.15e-06 | True |
| openmm-ir | numpy-ir | angle | 75.38676741 | 75.38676741 | 5.68e-14 | 1.75e-06 | True |
| openmm-ir | numpy-ir | proper | 25.92207054 | 25.92207054 | 3.55e-15 | 1.26e-06 | True |
| openmm-ir | numpy-ir | improper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | -0.2949938265 | -0.2949938265 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb | 7.251141814 | 7.251141814 | 0 | 1.07e-06 | True |
| openmm-ir | numpy-ir | lj14 | 4.241364517 | 4.241364517 | 3.55e-15 | 1.04e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | 0.6083074285 | 0.6083074285 | 1.44e-15 | 1.01e-06 | True |
| openmm-ir | numpy-ir | total | 228.3703655 | 228.3703655 | 8.53e-14 | 2.15e-06 | True |
| openmm-xml | numpy-ir | bond | 115.2557076 | 115.2557076 | 2.84e-14 | 2.15e-06 | True |
| openmm-xml | numpy-ir | angle | 75.38676741 | 75.38676741 | 5.68e-14 | 1.75e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 25.92207054 | 25.92207054 | 3.55e-15 | 1.26e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | 11.80582003 | 11.80581993 | 9.3e-08 | 1.12e-06 | True |
| openmm-xml | numpy-ir | total | 228.3703656 | 228.3703655 | 9.3e-08 | 2.15e-06 | True |
| openmm-prmtop | numpy-ir | bond | 115.2557076 | 115.2557076 | 2.84e-14 | 0.000112 | True |
| openmm-prmtop | numpy-ir | angle | 75.38676741 | 75.38676741 | 2.84e-14 | 0.000108 | True |
| openmm-prmtop | numpy-ir | torsion_all | 25.92207054 | 25.92207054 | 0 | 0.000103 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | 11.80582001 | 11.80581993 | 8.16e-08 | 0.000101 | True |
| openmm-prmtop | numpy-ir | total | 228.3703655 | 228.3703655 | 8.16e-08 | 0.000112 | True |
| openmm-top | numpy-ir | bond | 115.2557076 | 115.2557076 | 2.84e-14 | 0.000112 | True |
| openmm-top | numpy-ir | angle | 75.38676741 | 75.38676741 | 3.82e-09 | 0.000108 | True |
| openmm-top | numpy-ir | torsion_all | 25.92207054 | 25.92207054 | 1.59e-10 | 0.000103 | True |
| openmm-top | numpy-ir | nonbonded_all | 11.80582001 | 11.80581993 | 7.94e-08 | 0.000101 | True |
| openmm-top | numpy-ir | total | 228.3703655 | 228.3703655 | 7.54e-08 | 0.000112 | True |
| sander-ours | numpy-ir | bond | 27.5468 | 27.54677524 | 2.48e-05 | 0.000203 | True |
| sander-ours | numpy-ir | angle | 18.0179 | 18.01786984 | 3.02e-05 | 0.000202 | True |
| sander-ours | numpy-ir | dihed | 6.1955 | 6.195523551 | 2.36e-05 | 0.000201 | True |
| sander-ours | numpy-ir | vdw | -0.0705 | -0.07050521666 | 5.22e-06 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 1.733 | 1.733064487 | 6.45e-05 | 0.000287 | True |
| sander-ours | numpy-ir | vdw14 | 1.0137 | 1.013710449 | 1.04e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | 0.1454 | 0.1453889647 | 1.1e-05 | 0.000207 | True |
| sander-tleap | sander-ours | bond | 27.5468 | 27.5468 | 0 | 0.000203 | True |
| sander-tleap | sander-ours | angle | 18.0179 | 18.0179 | 0 | 0.000202 | True |
| sander-tleap | sander-ours | dihed | 6.1955 | 6.1955 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | vdw | -0.0705 | -0.0705 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 1.733 | 1.733 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | 1.0137 | 1.0137 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | 0.1454 | 0.1454 | 0 | 0.0002 | True |
| gmx | numpy-ir | bond | 115.2557297 | 115.2557076 | 2.21e-05 | 0.000112 | True |
| gmx | numpy-ir | angle | 75.38674164 | 75.38676741 | 2.58e-05 | 0.000108 | True |
| gmx | numpy-ir | proper | 25.92206955 | 25.92207054 | 9.86e-07 | 0.000103 | True |
| gmx | numpy-ir | improper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj | -0.2949940562 | -0.2949938265 | 2.3e-07 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 7.251141548 | 7.251141814 | 2.66e-07 | 0.000101 | True |
| gmx | numpy-ir | lj14 | 4.241368294 | 4.241364517 | 3.78e-06 | 0.0001 | True |
| gmx | numpy-ir | coulomb14 | 0.6083080173 | 0.6083074285 | 5.89e-07 | 0.0001 | True |
| gmx | numpy-ir | total | 228.3703613 | 228.3703655 | 4.14e-06 | 0.000123 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 6.25e-08 | 0.001 | True |
| openmm-top | openmm-ir | 3.44e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 2  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 119.9050121 | 119.9050121 | 0 | 2.2e-06 | True |
| numpy-ir | reference-itself | angle | 33.37048999 | 33.37048999 | 0 | 1.33e-06 | True |
| numpy-ir | reference-itself | proper | 27.2783065 | 27.2783065 | 0 | 1.27e-06 | True |
| numpy-ir | reference-itself | improper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | -0.3330607831 | -0.3330607831 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb | 7.130106598 | 7.130106598 | 0 | 1.07e-06 | True |
| numpy-ir | reference-itself | lj14 | 4.082253954 | 4.082253954 | 0 | 1.04e-06 | True |
| numpy-ir | reference-itself | coulomb14 | 0.2621358553 | 0.2621358553 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | total | 191.6952442 | 191.6952442 | 0 | 2.92e-06 | True |
| openmm-ir | numpy-ir | bond | 119.9050121 | 119.9050121 | 1.42e-14 | 2.2e-06 | True |
| openmm-ir | numpy-ir | angle | 33.37048999 | 33.37048999 | 1.42e-14 | 1.33e-06 | True |
| openmm-ir | numpy-ir | proper | 27.2783065 | 27.2783065 | 0 | 1.27e-06 | True |
| openmm-ir | numpy-ir | improper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | -0.3330607831 | -0.3330607831 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb | 7.130106598 | 7.130106598 | 0 | 1.07e-06 | True |
| openmm-ir | numpy-ir | lj14 | 4.082253954 | 4.082253954 | 0 | 1.04e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | 0.2621358553 | 0.2621358553 | 2.11e-15 | 1e-06 | True |
| openmm-ir | numpy-ir | total | 191.6952442 | 191.6952442 | 0 | 2.2e-06 | True |
| openmm-xml | numpy-ir | bond | 119.9050121 | 119.9050121 | 1.42e-14 | 2.2e-06 | True |
| openmm-xml | numpy-ir | angle | 33.37048999 | 33.37048999 | 1.42e-14 | 1.33e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 27.2783065 | 27.2783065 | 0 | 1.27e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | 11.14143571 | 11.14143562 | 8.75e-08 | 1.11e-06 | True |
| openmm-xml | numpy-ir | total | 191.6952443 | 191.6952442 | 8.75e-08 | 2.2e-06 | True |
| openmm-prmtop | numpy-ir | bond | 119.9050121 | 119.9050121 | 1.42e-14 | 0.000112 | True |
| openmm-prmtop | numpy-ir | angle | 33.37048999 | 33.37048999 | 2.13e-14 | 0.000103 | True |
| openmm-prmtop | numpy-ir | torsion_all | 27.2783065 | 27.2783065 | 3.55e-15 | 0.000103 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | 11.1414357 | 11.14143562 | 7.55e-08 | 0.000101 | True |
| openmm-prmtop | numpy-ir | total | 191.6952443 | 191.6952442 | 7.55e-08 | 0.000112 | True |
| openmm-top | numpy-ir | bond | 119.9050121 | 119.9050121 | 2.84e-14 | 0.000112 | True |
| openmm-top | numpy-ir | angle | 33.3704899 | 33.37048999 | 9.47e-08 | 0.000103 | True |
| openmm-top | numpy-ir | torsion_all | 27.2783065 | 27.2783065 | 1.68e-10 | 0.000103 | True |
| openmm-top | numpy-ir | nonbonded_all | 11.1414357 | 11.14143562 | 7.43e-08 | 0.000101 | True |
| openmm-top | numpy-ir | total | 191.6952442 | 191.6952442 | 2.05e-08 | 0.000112 | True |
| sander-ours | numpy-ir | bond | 28.658 | 28.65798568 | 1.43e-05 | 0.000203 | True |
| sander-ours | numpy-ir | angle | 7.9757 | 7.975738526 | 3.85e-05 | 0.000201 | True |
| sander-ours | numpy-ir | dihed | 6.5197 | 6.519671727 | 2.83e-05 | 0.000201 | True |
| sander-ours | numpy-ir | vdw | -0.0796 | -0.07960343764 | 3.44e-06 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 1.7041 | 1.704136376 | 3.64e-05 | 0.000285 | True |
| sander-ours | numpy-ir | vdw14 | 0.9757 | 0.9756821114 | 1.79e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | 0.0626 | 0.06265197307 | 5.2e-05 | 0.000203 | True |
| sander-tleap | sander-ours | bond | 28.658 | 28.658 | 0 | 0.000203 | True |
| sander-tleap | sander-ours | angle | 7.9757 | 7.9757 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | dihed | 6.5197 | 6.5197 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | vdw | -0.0796 | -0.0796 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 1.7041 | 1.7041 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | 0.9757 | 0.9757 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | 0.0626 | 0.0626 | 0 | 0.0002 | True |
| gmx | numpy-ir | bond | 119.9050522 | 119.9050121 | 4.01e-05 | 0.000112 | True |
| gmx | numpy-ir | angle | 33.37050629 | 33.37048999 | 1.63e-05 | 0.000103 | True |
| gmx | numpy-ir | proper | 27.27830887 | 27.2783065 | 2.36e-06 | 0.000103 | True |
| gmx | numpy-ir | improper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj | -0.3330608606 | -0.3330607831 | 7.76e-08 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 7.130106926 | 7.130106598 | 3.28e-07 | 0.000101 | True |
| gmx | numpy-ir | lj14 | 4.082254887 | 4.082253954 | 9.32e-07 | 0.0001 | True |
| gmx | numpy-ir | coulomb14 | 0.2621349692 | 0.2621358553 | 8.86e-07 | 0.0001 | True |
| gmx | numpy-ir | total | 191.6953125 | 191.6952442 | 6.83e-05 | 0.000119 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 5.78e-08 | 0.001 | True |
| openmm-top | openmm-ir | 3.5e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 3  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 127.7463512 | 127.7463512 | 0 | 2.28e-06 | True |
| numpy-ir | reference-itself | angle | 49.87253576 | 49.87253576 | 0 | 1.5e-06 | True |
| numpy-ir | reference-itself | proper | 28.11874407 | 28.11874407 | 0 | 1.28e-06 | True |
| numpy-ir | reference-itself | improper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | -0.3603461807 | -0.3603461807 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb | 7.197534313 | 7.197534313 | 0 | 1.07e-06 | True |
| numpy-ir | reference-itself | lj14 | 3.774944577 | 3.774944577 | 0 | 1.04e-06 | True |
| numpy-ir | reference-itself | coulomb14 | 0.6417521681 | 0.6417521681 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | total | 216.9915159 | 216.9915159 | 0 | 3.17e-06 | True |
| openmm-ir | numpy-ir | bond | 127.7463512 | 127.7463512 | 5.68e-14 | 2.28e-06 | True |
| openmm-ir | numpy-ir | angle | 49.87253576 | 49.87253576 | 0 | 1.5e-06 | True |
| openmm-ir | numpy-ir | proper | 28.11874407 | 28.11874407 | 3.55e-15 | 1.28e-06 | True |
| openmm-ir | numpy-ir | improper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | -0.3603461807 | -0.3603461807 | 5.55e-17 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb | 7.197534313 | 7.197534313 | 0 | 1.07e-06 | True |
| openmm-ir | numpy-ir | lj14 | 3.774944577 | 3.774944577 | 5.33e-15 | 1.04e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | 0.6417521681 | 0.6417521681 | 1.22e-15 | 1.01e-06 | True |
| openmm-ir | numpy-ir | total | 216.9915159 | 216.9915159 | 5.68e-14 | 2.28e-06 | True |
| openmm-xml | numpy-ir | bond | 127.7463512 | 127.7463512 | 5.68e-14 | 2.28e-06 | True |
| openmm-xml | numpy-ir | angle | 49.87253576 | 49.87253576 | 0 | 1.5e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 28.11874407 | 28.11874407 | 3.55e-15 | 1.28e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | 11.25388497 | 11.25388488 | 9.28e-08 | 1.11e-06 | True |
| openmm-xml | numpy-ir | total | 216.991516 | 216.9915159 | 9.28e-08 | 2.28e-06 | True |
| openmm-prmtop | numpy-ir | bond | 127.7463512 | 127.7463512 | 5.68e-14 | 0.000113 | True |
| openmm-prmtop | numpy-ir | angle | 49.87253576 | 49.87253576 | 1.42e-14 | 0.000105 | True |
| openmm-prmtop | numpy-ir | torsion_all | 28.11874407 | 28.11874407 | 0 | 0.000103 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | 11.25388496 | 11.25388488 | 8.14e-08 | 0.000101 | True |
| openmm-prmtop | numpy-ir | total | 216.991516 | 216.9915159 | 8.14e-08 | 0.000113 | True |
| openmm-top | numpy-ir | bond | 127.7463512 | 127.7463512 | 0 | 0.000113 | True |
| openmm-top | numpy-ir | angle | 49.87253572 | 49.87253576 | 3.32e-08 | 0.000105 | True |
| openmm-top | numpy-ir | torsion_all | 28.11874407 | 28.11874407 | 1.73e-10 | 0.000103 | True |
| openmm-top | numpy-ir | nonbonded_all | 11.25388496 | 11.25388488 | 7.96e-08 | 0.000101 | True |
| openmm-top | numpy-ir | total | 216.9915159 | 216.9915159 | 4.62e-08 | 0.000113 | True |
| sander-ours | numpy-ir | bond | 30.5321 | 30.5321107 | 1.07e-05 | 0.000203 | True |
| sander-ours | numpy-ir | angle | 11.9198 | 11.91982212 | 2.21e-05 | 0.000201 | True |
| sander-ours | numpy-ir | dihed | 6.7205 | 6.720541126 | 4.11e-05 | 0.000201 | True |
| sander-ours | numpy-ir | vdw | -0.0861 | -0.08612480419 | 2.48e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 1.7202 | 1.720251987 | 5.2e-05 | 0.000286 | True |
| sander-ours | numpy-ir | vdw14 | 0.9022 | 0.9022334076 | 3.34e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | 0.1534 | 0.1533824493 | 1.76e-05 | 0.000208 | True |
| sander-tleap | sander-ours | bond | 30.5321 | 30.5321 | 0 | 0.000203 | True |
| sander-tleap | sander-ours | angle | 11.9198 | 11.9198 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | dihed | 6.7205 | 6.7205 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | vdw | -0.0861 | -0.0861 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 1.7202 | 1.7202 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | 0.9022 | 0.9022 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | 0.1534 | 0.1534 | 0 | 0.0002 | True |
| gmx | numpy-ir | bond | 127.746315 | 127.7463512 | 3.62e-05 | 0.000113 | True |
| gmx | numpy-ir | angle | 49.87258148 | 49.87253576 | 4.57e-05 | 0.000105 | True |
| gmx | numpy-ir | proper | 28.11874199 | 28.11874407 | 2.08e-06 | 0.000103 | True |
| gmx | numpy-ir | improper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj | -0.3603461385 | -0.3603461807 | 4.23e-08 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 7.197534561 | 7.197534313 | 2.49e-07 | 0.000101 | True |
| gmx | numpy-ir | lj14 | 3.774942636 | 3.774944577 | 1.94e-06 | 0.0001 | True |
| gmx | numpy-ir | coulomb14 | 0.6417523026 | 0.6417521681 | 1.35e-07 | 0.0001 | True |
| gmx | numpy-ir | total | 216.9915314 | 216.9915159 | 1.55e-05 | 0.000122 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 6.22e-08 | 0.001 | True |
| openmm-top | openmm-ir | 3.54e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

## norbornane
### config 0  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 1.431229985 | 1.431229985 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | angle | 19.87787454 | 19.87787454 | 0 | 1.2e-06 | True |
| numpy-ir | reference-itself | proper | 39.15169399 | 39.15169399 | 0 | 1.39e-06 | True |
| numpy-ir | reference-itself | improper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | -1.734386613 | -1.734386613 | 0 | 1.02e-06 | True |
| numpy-ir | reference-itself | coulomb | 21.84776847 | 21.84776847 | 0 | 1.22e-06 | True |
| numpy-ir | reference-itself | lj14 | 12.66723714 | 12.66723714 | 0 | 1.13e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -27.00675596 | -27.00675596 | 0 | 1.27e-06 | True |
| numpy-ir | reference-itself | total | 66.23466156 | 66.23466156 | 0 | 1.66e-06 | True |
| openmm-ir | numpy-ir | bond | 1.431229985 | 1.431229985 | 2.66e-14 | 1.01e-06 | True |
| openmm-ir | numpy-ir | angle | 19.87787454 | 19.87787454 | 3.55e-15 | 1.2e-06 | True |
| openmm-ir | numpy-ir | proper | 39.15169399 | 39.15169399 | 0 | 1.39e-06 | True |
| openmm-ir | numpy-ir | improper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | -1.734386613 | -1.734386613 | 1.78e-15 | 1.02e-06 | True |
| openmm-ir | numpy-ir | coulomb | 21.84776847 | 21.84776847 | 0 | 1.22e-06 | True |
| openmm-ir | numpy-ir | lj14 | 12.66723714 | 12.66723714 | 5.33e-15 | 1.13e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -27.00675596 | -27.00675596 | 3.55e-15 | 1.27e-06 | True |
| openmm-ir | numpy-ir | total | 66.23466156 | 66.23466156 | 2.84e-14 | 1.39e-06 | True |
| openmm-xml | numpy-ir | bond | 1.431229985 | 1.431229985 | 2.66e-14 | 1.01e-06 | True |
| openmm-xml | numpy-ir | angle | 19.87787454 | 19.87787454 | 3.55e-15 | 1.2e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 39.15169399 | 39.15169399 | 0 | 1.39e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | 5.773862986 | 5.773863047 | 6.11e-08 | 1.06e-06 | True |
| openmm-xml | numpy-ir | total | 66.2346615 | 66.23466156 | 6.11e-08 | 1.39e-06 | True |
| openmm-prmtop | numpy-ir | bond | 1.431229985 | 1.431229985 | 2.64e-14 | 0.0001 | True |
| openmm-prmtop | numpy-ir | angle | 19.87787454 | 19.87787454 | 3.55e-15 | 0.000102 | True |
| openmm-prmtop | numpy-ir | torsion_all | 39.15169399 | 39.15169399 | 2.13e-14 | 0.000104 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | 5.77386295 | 5.773863047 | 9.67e-08 | 0.000101 | True |
| openmm-prmtop | numpy-ir | total | 66.23466146 | 66.23466156 | 9.67e-08 | 0.000104 | True |
| openmm-top | numpy-ir | bond | 1.431229985 | 1.431229985 | 2.73e-14 | 0.0001 | True |
| openmm-top | numpy-ir | angle | 19.87787457 | 19.87787454 | 2.82e-08 | 0.000102 | True |
| openmm-top | numpy-ir | torsion_all | 39.15169399 | 39.15169399 | 0 | 0.000104 | True |
| openmm-top | numpy-ir | nonbonded_all | 5.773863009 | 5.773863047 | 3.8e-08 | 0.000101 | True |
| openmm-top | numpy-ir | total | 66.23466155 | 66.23466156 | 9.73e-09 | 0.000104 | True |
| sander-ours | numpy-ir | bond | 0.3421 | 0.3420721761 | 2.78e-05 | 0.0002 | True |
| sander-ours | numpy-ir | angle | 4.7509 | 4.750926038 | 2.6e-05 | 0.0002 | True |
| sander-ours | numpy-ir | dihed | 9.3575 | 9.357479442 | 2.06e-05 | 0.000201 | True |
| sander-ours | numpy-ir | vdw | -0.4145 | -0.4145283492 | 2.83e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 5.2216 | 5.221741987 | 0.000142 | 0.000461 | True |
| sander-ours | numpy-ir | vdw14 | 3.0275 | 3.027542338 | 4.23e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -6.4545 | -6.454769588 | 0.00027 | 0.000523 | True |
| sander-tleap | sander-ours | bond | 0.3421 | 0.3421 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | angle | 4.7509 | 4.7509 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | dihed | 9.3575 | 9.3575 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | vdw | -0.4145 | -0.4145 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 5.2216 | 5.2216 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | vdw14 | 3.0275 | 3.0275 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -6.4545 | -6.4545 | 0 | 0.000201 | True |
| gmx | numpy-ir | bond | 1.431214571 | 1.431229985 | 1.54e-05 | 0.0001 | True |
| gmx | numpy-ir | angle | 19.87790108 | 19.87787454 | 2.65e-05 | 0.000102 | True |
| gmx | numpy-ir | proper | 39.15169144 | 39.15169399 | 2.55e-06 | 0.000104 | True |
| gmx | numpy-ir | improper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj | -1.73439002 | -1.734386613 | 3.41e-06 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 21.84776878 | 21.84776847 | 3.09e-07 | 0.000102 | True |
| gmx | numpy-ir | lj14 | 12.66723442 | 12.66723714 | 2.72e-06 | 0.000101 | True |
| gmx | numpy-ir | coulomb14 | -27.00675774 | -27.00675596 | 1.78e-06 | 0.000103 | True |
| gmx | numpy-ir | total | 66.23466492 | 66.23466156 | 3.36e-06 | 0.000107 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 7.34e-08 | 0.001 | True |
| openmm-top | openmm-ir | 1.52e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 1  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 88.75652229 | 88.75652229 | 0 | 1.89e-06 | True |
| numpy-ir | reference-itself | angle | 72.52331226 | 72.52331226 | 0 | 1.73e-06 | True |
| numpy-ir | reference-itself | proper | 39.84011166 | 39.84011166 | 0 | 1.4e-06 | True |
| numpy-ir | reference-itself | improper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | -1.655154283 | -1.655154283 | 0 | 1.02e-06 | True |
| numpy-ir | reference-itself | coulomb | 21.78738535 | 21.78738535 | 0 | 1.22e-06 | True |
| numpy-ir | reference-itself | lj14 | 13.24585313 | 13.24585313 | 0 | 1.13e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -26.89885219 | -26.89885219 | 0 | 1.27e-06 | True |
| numpy-ir | reference-itself | total | 207.5991782 | 207.5991782 | 0 | 3.08e-06 | True |
| openmm-ir | numpy-ir | bond | 88.75652229 | 88.75652229 | 0 | 1.89e-06 | True |
| openmm-ir | numpy-ir | angle | 72.52331226 | 72.52331226 | 2.84e-14 | 1.73e-06 | True |
| openmm-ir | numpy-ir | proper | 39.84011166 | 39.84011166 | 0 | 1.4e-06 | True |
| openmm-ir | numpy-ir | improper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | -1.655154283 | -1.655154283 | 2.22e-16 | 1.02e-06 | True |
| openmm-ir | numpy-ir | coulomb | 21.78738535 | 21.78738535 | 0 | 1.22e-06 | True |
| openmm-ir | numpy-ir | lj14 | 13.24585313 | 13.24585313 | 8.88e-15 | 1.13e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -26.89885219 | -26.89885219 | 0 | 1.27e-06 | True |
| openmm-ir | numpy-ir | total | 207.5991782 | 207.5991782 | 2.84e-14 | 1.89e-06 | True |
| openmm-xml | numpy-ir | bond | 88.75652229 | 88.75652229 | 0 | 1.89e-06 | True |
| openmm-xml | numpy-ir | angle | 72.52331226 | 72.52331226 | 2.84e-14 | 1.73e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 39.84011166 | 39.84011166 | 0 | 1.4e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | 6.479231953 | 6.479232013 | 6.05e-08 | 1.06e-06 | True |
| openmm-xml | numpy-ir | total | 207.5991782 | 207.5991782 | 6.05e-08 | 1.89e-06 | True |
| openmm-prmtop | numpy-ir | bond | 88.75652229 | 88.75652229 | 0 | 0.000109 | True |
| openmm-prmtop | numpy-ir | angle | 72.52331226 | 72.52331226 | 2.84e-14 | 0.000107 | True |
| openmm-prmtop | numpy-ir | torsion_all | 39.84011166 | 39.84011166 | 7.11e-15 | 0.000104 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | 6.479231918 | 6.479232013 | 9.59e-08 | 0.000101 | True |
| openmm-prmtop | numpy-ir | total | 207.5991781 | 207.5991782 | 9.59e-08 | 0.000109 | True |
| openmm-top | numpy-ir | bond | 88.75652229 | 88.75652229 | 7.11e-14 | 0.000109 | True |
| openmm-top | numpy-ir | angle | 72.5233123 | 72.52331226 | 4.13e-08 | 0.000107 | True |
| openmm-top | numpy-ir | torsion_all | 39.84011166 | 39.84011166 | 7.11e-15 | 0.000104 | True |
| openmm-top | numpy-ir | nonbonded_all | 6.479231976 | 6.479232013 | 3.79e-08 | 0.000101 | True |
| openmm-top | numpy-ir | total | 207.5991782 | 207.5991782 | 3.38e-09 | 0.000109 | True |
| sander-ours | numpy-ir | bond | 21.2133 | 21.21331795 | 1.79e-05 | 0.000202 | True |
| sander-ours | numpy-ir | angle | 17.3335 | 17.33348763 | 1.24e-05 | 0.000202 | True |
| sander-ours | numpy-ir | dihed | 9.522 | 9.522015216 | 1.52e-05 | 0.000201 | True |
| sander-ours | numpy-ir | vdw | -0.3956 | -0.3955913678 | 8.63e-06 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 5.2071 | 5.207310075 | 0.00021 | 0.00046 | True |
| sander-ours | numpy-ir | vdw14 | 3.1658 | 3.165834879 | 3.49e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -6.4288 | -6.42897997 | 0.00018 | 0.000521 | True |
| sander-tleap | sander-ours | bond | 21.2133 | 21.2133 | 0 | 0.000202 | True |
| sander-tleap | sander-ours | angle | 17.3335 | 17.3335 | 0 | 0.000202 | True |
| sander-tleap | sander-ours | dihed | 9.522 | 9.522 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | vdw | -0.3956 | -0.3956 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 5.2071 | 5.2071 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | vdw14 | 3.1658 | 3.1658 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -6.4288 | -6.4288 | 0 | 0.000201 | True |
| gmx | numpy-ir | bond | 88.75649261 | 88.75652229 | 2.97e-05 | 0.000109 | True |
| gmx | numpy-ir | angle | 72.52331543 | 72.52331226 | 3.17e-06 | 0.000107 | True |
| gmx | numpy-ir | proper | 39.84010696 | 39.84011166 | 4.7e-06 | 0.000104 | True |
| gmx | numpy-ir | improper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj | -1.65515554 | -1.655154283 | 1.26e-06 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 21.78738403 | 21.78738535 | 1.32e-06 | 0.000102 | True |
| gmx | numpy-ir | lj14 | 13.24584961 | 13.24585313 | 3.52e-06 | 0.000101 | True |
| gmx | numpy-ir | coulomb14 | -26.89886093 | -26.89885219 | 8.74e-06 | 0.000103 | True |
| gmx | numpy-ir | total | 207.5991364 | 207.5991782 | 4.19e-05 | 0.000121 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 7.89e-08 | 0.001 | True |
| openmm-top | openmm-ir | 1.57e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 2  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 78.59119922 | 78.59119922 | 0 | 1.79e-06 | True |
| numpy-ir | reference-itself | angle | 72.0636049 | 72.0636049 | 0 | 1.72e-06 | True |
| numpy-ir | reference-itself | proper | 38.30560271 | 38.30560271 | 0 | 1.38e-06 | True |
| numpy-ir | reference-itself | improper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | -1.944750647 | -1.944750647 | 0 | 1.02e-06 | True |
| numpy-ir | reference-itself | coulomb | 21.6010553 | 21.6010553 | 0 | 1.22e-06 | True |
| numpy-ir | reference-itself | lj14 | 12.04085823 | 12.04085823 | 0 | 1.12e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -26.72913645 | -26.72913645 | 0 | 1.27e-06 | True |
| numpy-ir | reference-itself | total | 193.9284333 | 193.9284333 | 0 | 2.94e-06 | True |
| openmm-ir | numpy-ir | bond | 78.59119922 | 78.59119922 | 1.42e-14 | 1.79e-06 | True |
| openmm-ir | numpy-ir | angle | 72.0636049 | 72.0636049 | 1.42e-14 | 1.72e-06 | True |
| openmm-ir | numpy-ir | proper | 38.30560271 | 38.30560271 | 1.42e-14 | 1.38e-06 | True |
| openmm-ir | numpy-ir | improper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | -1.944750647 | -1.944750647 | 2.22e-16 | 1.02e-06 | True |
| openmm-ir | numpy-ir | coulomb | 21.6010553 | 21.6010553 | 0 | 1.22e-06 | True |
| openmm-ir | numpy-ir | lj14 | 12.04085823 | 12.04085823 | 1.78e-15 | 1.12e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -26.72913645 | -26.72913645 | 1.42e-14 | 1.27e-06 | True |
| openmm-ir | numpy-ir | total | 193.9284333 | 193.9284333 | 5.68e-14 | 1.79e-06 | True |
| openmm-xml | numpy-ir | bond | 78.59119922 | 78.59119922 | 1.42e-14 | 1.79e-06 | True |
| openmm-xml | numpy-ir | angle | 72.0636049 | 72.0636049 | 1.42e-14 | 1.72e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 38.30560271 | 38.30560271 | 1.42e-14 | 1.38e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | 4.968026376 | 4.968026437 | 6.07e-08 | 1.05e-06 | True |
| openmm-xml | numpy-ir | total | 193.9284332 | 193.9284333 | 6.07e-08 | 1.79e-06 | True |
| openmm-prmtop | numpy-ir | bond | 78.59119922 | 78.59119922 | 2.84e-14 | 0.000108 | True |
| openmm-prmtop | numpy-ir | angle | 72.0636049 | 72.0636049 | 0 | 0.000107 | True |
| openmm-prmtop | numpy-ir | torsion_all | 38.30560271 | 38.30560271 | 3.55e-14 | 0.000104 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | 4.968026342 | 4.968026437 | 9.47e-08 | 0.0001 | True |
| openmm-prmtop | numpy-ir | total | 193.9284332 | 193.9284333 | 9.47e-08 | 0.000108 | True |
| openmm-top | numpy-ir | bond | 78.59119922 | 78.59119922 | 2.84e-14 | 0.000108 | True |
| openmm-top | numpy-ir | angle | 72.06360494 | 72.0636049 | 4.21e-08 | 0.000107 | True |
| openmm-top | numpy-ir | torsion_all | 38.30560271 | 38.30560271 | 2.84e-14 | 0.000104 | True |
| openmm-top | numpy-ir | nonbonded_all | 4.968026399 | 4.968026437 | 3.74e-08 | 0.0001 | True |
| openmm-top | numpy-ir | total | 193.9284333 | 193.9284333 | 4.68e-09 | 0.000108 | True |
| sander-ours | numpy-ir | bond | 18.7837 | 18.78374742 | 4.74e-05 | 0.000202 | True |
| sander-ours | numpy-ir | angle | 17.2236 | 17.22361494 | 1.49e-05 | 0.000202 | True |
| sander-ours | numpy-ir | dihed | 9.1553 | 9.155258774 | 4.12e-05 | 0.000201 | True |
| sander-ours | numpy-ir | vdw | -0.4648 | -0.46480656 | 6.56e-06 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 5.1626 | 5.162776124 | 0.000176 | 0.000458 | True |
| sander-ours | numpy-ir | vdw14 | 2.8778 | 2.877834184 | 3.42e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -6.3882 | -6.388416932 | 0.000217 | 0.000519 | True |
| sander-tleap | sander-ours | bond | 18.7837 | 18.7837 | 0 | 0.000202 | True |
| sander-tleap | sander-ours | angle | 17.2236 | 17.2236 | 0 | 0.000202 | True |
| sander-tleap | sander-ours | dihed | 9.1553 | 9.1553 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | vdw | -0.4648 | -0.4648 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 5.1626 | 5.1626 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | vdw14 | 2.8778 | 2.8778 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -6.3882 | -6.3882 | 0 | 0.000201 | True |
| gmx | numpy-ir | bond | 78.59111023 | 78.59119922 | 8.9e-05 | 0.000108 | True |
| gmx | numpy-ir | angle | 72.06359863 | 72.0636049 | 6.27e-06 | 0.000107 | True |
| gmx | numpy-ir | proper | 38.30560303 | 38.30560271 | 3.16e-07 | 0.000104 | True |
| gmx | numpy-ir | improper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj | -1.944751859 | -1.944750647 | 1.21e-06 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 21.60105515 | 21.6010553 | 1.59e-07 | 0.000102 | True |
| gmx | numpy-ir | lj14 | 12.04085159 | 12.04085823 | 6.63e-06 | 0.000101 | True |
| gmx | numpy-ir | coulomb14 | -26.72912788 | -26.72913645 | 8.56e-06 | 0.000103 | True |
| gmx | numpy-ir | total | 193.9283447 | 193.9284333 | 8.85e-05 | 0.000119 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 6.98e-08 | 0.001 | True |
| openmm-top | openmm-ir | 1.82e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 3  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 77.96302343 | 77.96302343 | 0 | 1.78e-06 | True |
| numpy-ir | reference-itself | angle | 64.74494503 | 64.74494503 | 0 | 1.65e-06 | True |
| numpy-ir | reference-itself | proper | 39.06591784 | 39.06591784 | 0 | 1.39e-06 | True |
| numpy-ir | reference-itself | improper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | -1.374448556 | -1.374448556 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | coulomb | 22.07590215 | 22.07590215 | 0 | 1.22e-06 | True |
| numpy-ir | reference-itself | lj14 | 13.23482912 | 13.23482912 | 0 | 1.13e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -26.91805424 | -26.91805424 | 0 | 1.27e-06 | True |
| numpy-ir | reference-itself | total | 188.7921148 | 188.7921148 | 0 | 2.89e-06 | True |
| openmm-ir | numpy-ir | bond | 77.96302343 | 77.96302343 | 0 | 1.78e-06 | True |
| openmm-ir | numpy-ir | angle | 64.74494503 | 64.74494503 | 5.68e-14 | 1.65e-06 | True |
| openmm-ir | numpy-ir | proper | 39.06591784 | 39.06591784 | 0 | 1.39e-06 | True |
| openmm-ir | numpy-ir | improper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | -1.374448556 | -1.374448556 | 0 | 1.01e-06 | True |
| openmm-ir | numpy-ir | coulomb | 22.07590215 | 22.07590215 | 0 | 1.22e-06 | True |
| openmm-ir | numpy-ir | lj14 | 13.23482912 | 13.23482912 | 1.78e-15 | 1.13e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -26.91805424 | -26.91805424 | 0 | 1.27e-06 | True |
| openmm-ir | numpy-ir | total | 188.7921148 | 188.7921148 | 5.68e-14 | 1.78e-06 | True |
| openmm-xml | numpy-ir | bond | 77.96302343 | 77.96302343 | 0 | 1.78e-06 | True |
| openmm-xml | numpy-ir | angle | 64.74494503 | 64.74494503 | 5.68e-14 | 1.65e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 39.06591784 | 39.06591784 | 0 | 1.39e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | 7.018228405 | 7.018228463 | 5.73e-08 | 1.07e-06 | True |
| openmm-xml | numpy-ir | total | 188.7921147 | 188.7921148 | 5.73e-08 | 1.78e-06 | True |
| openmm-prmtop | numpy-ir | bond | 77.96302343 | 77.96302343 | 1.42e-14 | 0.000108 | True |
| openmm-prmtop | numpy-ir | angle | 64.74494503 | 64.74494503 | 5.68e-14 | 0.000106 | True |
| openmm-prmtop | numpy-ir | torsion_all | 39.06591784 | 39.06591784 | 0 | 0.000104 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | 7.018228369 | 7.018228463 | 9.34e-08 | 0.000101 | True |
| openmm-prmtop | numpy-ir | total | 188.7921147 | 188.7921148 | 9.34e-08 | 0.000108 | True |
| openmm-top | numpy-ir | bond | 77.96302343 | 77.96302343 | 8.53e-14 | 0.000108 | True |
| openmm-top | numpy-ir | angle | 64.74494507 | 64.74494503 | 4.41e-08 | 0.000106 | True |
| openmm-top | numpy-ir | torsion_all | 39.06591784 | 39.06591784 | 0 | 0.000104 | True |
| openmm-top | numpy-ir | nonbonded_all | 7.018228428 | 7.018228463 | 3.45e-08 | 0.000101 | True |
| openmm-top | numpy-ir | total | 188.7921148 | 188.7921148 | 9.54e-09 | 0.000108 | True |
| sander-ours | numpy-ir | bond | 18.6336 | 18.63360981 | 9.81e-06 | 0.000202 | True |
| sander-ours | numpy-ir | angle | 15.4744 | 15.47441325 | 1.32e-05 | 0.000202 | True |
| sander-ours | numpy-ir | dihed | 9.337 | 9.33697845 | 2.15e-05 | 0.000201 | True |
| sander-ours | numpy-ir | vdw | -0.3285 | -0.328501089 | 1.09e-06 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 5.2761 | 5.276267243 | 0.000167 | 0.000464 | True |
| sander-ours | numpy-ir | vdw14 | 3.1632 | 3.163200075 | 7.53e-08 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -6.4333 | -6.43356937 | 0.000269 | 0.000522 | True |
| sander-tleap | sander-ours | bond | 18.6336 | 18.6336 | 0 | 0.000202 | True |
| sander-tleap | sander-ours | angle | 15.4744 | 15.4744 | 0 | 0.000202 | True |
| sander-tleap | sander-ours | dihed | 9.337 | 9.337 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | vdw | -0.3285 | -0.3285 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 5.2761 | 5.2761 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | vdw14 | 3.1632 | 3.1632 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -6.4333 | -6.4333 | 0 | 0.000201 | True |
| gmx | numpy-ir | bond | 77.96299744 | 77.96302343 | 2.6e-05 | 0.000108 | True |
| gmx | numpy-ir | angle | 64.74491119 | 64.74494503 | 3.38e-05 | 0.000106 | True |
| gmx | numpy-ir | proper | 39.0659256 | 39.06591784 | 7.76e-06 | 0.000104 | True |
| gmx | numpy-ir | improper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj | -1.374451876 | -1.374448556 | 3.32e-06 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 22.07590294 | 22.07590215 | 7.92e-07 | 0.000102 | True |
| gmx | numpy-ir | lj14 | 13.23481941 | 13.23482912 | 9.7e-06 | 0.000101 | True |
| gmx | numpy-ir | coulomb14 | -26.91804886 | -26.91805424 | 5.38e-06 | 0.000103 | True |
| gmx | numpy-ir | total | 188.7920532 | 188.7921148 | 6.15e-05 | 0.000119 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 7.67e-08 | 0.001 | True |
| openmm-top | openmm-ir | 1.73e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

## acetate
### config 0  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 9.441754753 | 9.441754753 | 0 | 1.09e-06 | True |
| numpy-ir | reference-itself | angle | 12.58729995 | 12.58729995 | 0 | 1.13e-06 | True |
| numpy-ir | reference-itself | proper | 192.8405594 | 192.8405594 | 0 | 2.93e-06 | True |
| numpy-ir | reference-itself | improper | 1.751809337e-06 | 1.751809337e-06 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj14 | 1.638880082 | 1.638880082 | 0 | 1.02e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -14.84882028 | -14.84882028 | 0 | 1.15e-06 | True |
| numpy-ir | reference-itself | total | 201.6596757 | 201.6596757 | 0 | 3.02e-06 | True |
| openmm-ir | numpy-ir | bond | 9.441754753 | 9.441754753 | 0 | 1.09e-06 | True |
| openmm-ir | numpy-ir | angle | 12.58729995 | 12.58729995 | 0 | 1.13e-06 | True |
| openmm-ir | numpy-ir | proper | 192.8405594 | 192.8405594 | 0 | 2.93e-06 | True |
| openmm-ir | numpy-ir | improper | 1.751809337e-06 | 1.751809337e-06 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj14 | 1.638880082 | 1.638880082 | 8.88e-16 | 1.02e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -14.84882028 | -14.84882028 | 0 | 1.15e-06 | True |
| openmm-ir | numpy-ir | total | 201.6596757 | 201.6596757 | 0 | 2.93e-06 | True |
| openmm-xml | numpy-ir | bond | 9.441754753 | 9.441754753 | 0 | 1.09e-06 | True |
| openmm-xml | numpy-ir | angle | 12.58729995 | 12.58729995 | 0 | 1.13e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 192.8405612 | 192.8405612 | 0 | 2.93e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | -13.20994038 | -13.2099402 | 1.76e-07 | 1.13e-06 | True |
| openmm-xml | numpy-ir | total | 201.6596755 | 201.6596757 | 1.76e-07 | 2.93e-06 | True |
| openmm-prmtop | numpy-ir | bond | 9.441754753 | 9.441754753 | 3.55e-15 | 0.000101 | True |
| openmm-prmtop | numpy-ir | angle | 12.58729995 | 12.58729995 | 0 | 0.000101 | True |
| openmm-prmtop | numpy-ir | torsion_all | 192.8405612 | 192.8405612 | 0 | 0.000119 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | -13.20994039 | -13.2099402 | 1.89e-07 | 0.000101 | True |
| openmm-prmtop | numpy-ir | total | 201.6596755 | 201.6596757 | 1.89e-07 | 0.000119 | True |
| openmm-top | numpy-ir | bond | 9.441754753 | 9.441754753 | 5.86e-14 | 0.000101 | True |
| openmm-top | numpy-ir | angle | 12.58730004 | 12.58729995 | 8.34e-08 | 0.000101 | True |
| openmm-top | numpy-ir | torsion_all | 192.8405612 | 192.8405612 | 1.35e-11 | 0.000119 | True |
| openmm-top | numpy-ir | nonbonded_all | -13.20994036 | -13.2099402 | 1.6e-07 | 0.000101 | True |
| openmm-top | numpy-ir | total | 201.6596756 | 201.6596757 | 7.61e-08 | 0.000119 | True |
| sander-ours | numpy-ir | bond | 2.2566 | 2.256633545 | 3.35e-05 | 0.0002 | True |
| sander-ours | numpy-ir | angle | 3.0084 | 3.008436891 | 3.69e-05 | 0.0002 | True |
| sander-ours | numpy-ir | dihed | 46.09 | 46.09000028 | 2.79e-07 | 0.000205 | True |
| sander-ours | numpy-ir | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | vdw14 | 0.3917 | 0.3917017405 | 1.74e-06 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -3.5488 | -3.548953222 | 0.000153 | 0.000377 | True |
| sander-tleap | sander-ours | bond | 2.2566 | 2.2566 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | angle | 3.0084 | 3.0084 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | dihed | 46.09 | 46.09 | 0 | 0.000205 | True |
| sander-tleap | sander-ours | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | 0.3917 | 0.3917 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -3.5488 | -3.5488 | 0 | 0.0002 | True |
| gmx | numpy-ir | bond | 9.441771507 | 9.441754753 | 1.68e-05 | 0.000101 | True |
| gmx | numpy-ir | angle | 12.58731747 | 12.58729995 | 1.75e-05 | 0.000101 | True |
| gmx | numpy-ir | proper | 192.8405762 | 192.8405594 | 1.68e-05 | 0.000119 | True |
| gmx | numpy-ir | improper | 1.645946e-06 | 1.751809337e-06 | 1.06e-07 | 0.0001 | True |
| gmx | numpy-ir | lj | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj14 | 1.638878703 | 1.638880082 | 1.38e-06 | 0.0001 | True |
| gmx | numpy-ir | coulomb14 | -14.84881878 | -14.84882028 | 1.5e-06 | 0.000101 | True |
| gmx | numpy-ir | total | 201.659729 | 201.6596757 | 5.33e-05 | 0.00012 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 1.51e-07 | 0.001 | True |
| openmm-top | openmm-ir | 6.96e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 1  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 37.72496229 | 37.72496229 | 0 | 1.38e-06 | True |
| numpy-ir | reference-itself | angle | 20.44942325 | 20.44942325 | 0 | 1.2e-06 | True |
| numpy-ir | reference-itself | proper | 193.4507104 | 193.4507104 | 0 | 2.93e-06 | True |
| numpy-ir | reference-itself | improper | 0.5090147767 | 0.5090147767 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | lj | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj14 | 0.9996293506 | 0.9996293506 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -14.5382572 | -14.5382572 | 0 | 1.15e-06 | True |
| numpy-ir | reference-itself | total | 238.5954829 | 238.5954829 | 0 | 3.39e-06 | True |
| openmm-ir | numpy-ir | bond | 37.72496229 | 37.72496229 | 1.42e-14 | 1.38e-06 | True |
| openmm-ir | numpy-ir | angle | 20.44942325 | 20.44942325 | 2.84e-14 | 1.2e-06 | True |
| openmm-ir | numpy-ir | proper | 193.4507104 | 193.4507104 | 0 | 2.93e-06 | True |
| openmm-ir | numpy-ir | improper | 0.5090147767 | 0.5090147767 | 2e-15 | 1.01e-06 | True |
| openmm-ir | numpy-ir | lj | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj14 | 0.9996293506 | 0.9996293506 | 4.44e-16 | 1.01e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -14.5382572 | -14.5382572 | 0 | 1.15e-06 | True |
| openmm-ir | numpy-ir | total | 238.5954829 | 238.5954829 | 0 | 2.93e-06 | True |
| openmm-xml | numpy-ir | bond | 37.72496229 | 37.72496229 | 1.42e-14 | 1.38e-06 | True |
| openmm-xml | numpy-ir | angle | 20.44942325 | 20.44942325 | 2.84e-14 | 1.2e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 193.9597252 | 193.9597252 | 0 | 2.94e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | -13.53862802 | -13.53862785 | 1.72e-07 | 1.14e-06 | True |
| openmm-xml | numpy-ir | total | 238.5954827 | 238.5954829 | 1.72e-07 | 2.93e-06 | True |
| openmm-prmtop | numpy-ir | bond | 37.72496229 | 37.72496229 | 2.13e-14 | 0.000104 | True |
| openmm-prmtop | numpy-ir | angle | 20.44942325 | 20.44942325 | 2.13e-14 | 0.000102 | True |
| openmm-prmtop | numpy-ir | torsion_all | 193.9597252 | 193.9597252 | 0 | 0.000119 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | -13.53862804 | -13.53862785 | 1.84e-07 | 0.000101 | True |
| openmm-prmtop | numpy-ir | total | 238.5954827 | 238.5954829 | 1.84e-07 | 0.000119 | True |
| openmm-top | numpy-ir | bond | 37.72496229 | 37.72496229 | 1.21e-13 | 0.000104 | True |
| openmm-top | numpy-ir | angle | 20.44942331 | 20.44942325 | 6.18e-08 | 0.000102 | True |
| openmm-top | numpy-ir | torsion_all | 193.9597252 | 193.9597252 | 4.62e-09 | 0.000119 | True |
| openmm-top | numpy-ir | nonbonded_all | -13.53862801 | -13.53862785 | 1.56e-07 | 0.000101 | True |
| openmm-top | numpy-ir | total | 238.5954828 | 238.5954829 | 9.92e-08 | 0.000119 | True |
| sander-ours | numpy-ir | bond | 9.0165 | 9.016482384 | 1.76e-05 | 0.000201 | True |
| sander-ours | numpy-ir | angle | 4.8875 | 4.887529457 | 2.95e-05 | 0.0002 | True |
| sander-ours | numpy-ir | dihed | 46.3575 | 46.35748691 | 1.31e-05 | 0.000205 | True |
| sander-ours | numpy-ir | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | vdw14 | 0.2389 | 0.2389171488 | 1.71e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -3.4746 | -3.474726865 | 0.000127 | 0.000374 | True |
| sander-tleap | sander-ours | bond | 9.0165 | 9.0165 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | angle | 4.8875 | 4.8875 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | dihed | 46.3575 | 46.3575 | 0 | 0.000205 | True |
| sander-tleap | sander-ours | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | 0.2389 | 0.2389 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -3.4746 | -3.4746 | 0 | 0.0002 | True |
| gmx | numpy-ir | bond | 37.72489929 | 37.72496229 | 6.3e-05 | 0.000104 | True |
| gmx | numpy-ir | angle | 20.4494381 | 20.44942325 | 1.48e-05 | 0.000102 | True |
| gmx | numpy-ir | proper | 193.4507294 | 193.4507104 | 1.89e-05 | 0.000119 | True |
| gmx | numpy-ir | improper | 0.5090144277 | 0.5090147767 | 3.49e-07 | 0.0001 | True |
| gmx | numpy-ir | lj | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj14 | 0.9996300936 | 0.9996293506 | 7.43e-07 | 0.0001 | True |
| gmx | numpy-ir | coulomb14 | -14.5382576 | -14.5382572 | 3.96e-07 | 0.000101 | True |
| gmx | numpy-ir | total | 238.595459 | 238.5954829 | 2.39e-05 | 0.000124 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 1.41e-07 | 0.001 | True |
| openmm-top | openmm-ir | 7.06e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 2  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 61.14815051 | 61.14815051 | 0 | 1.61e-06 | True |
| numpy-ir | reference-itself | angle | 81.13997439 | 81.13997439 | 0 | 1.81e-06 | True |
| numpy-ir | reference-itself | proper | 203.9757186 | 203.9757186 | 0 | 3.04e-06 | True |
| numpy-ir | reference-itself | improper | 0.05472210537 | 0.05472210537 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj14 | 3.645523707 | 3.645523707 | 0 | 1.04e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -15.02097654 | -15.02097654 | 0 | 1.15e-06 | True |
| numpy-ir | reference-itself | total | 334.9431127 | 334.9431127 | 0 | 4.35e-06 | True |
| openmm-ir | numpy-ir | bond | 61.14815051 | 61.14815051 | 0 | 1.61e-06 | True |
| openmm-ir | numpy-ir | angle | 81.13997439 | 81.13997439 | 1.42e-14 | 1.81e-06 | True |
| openmm-ir | numpy-ir | proper | 203.9757186 | 203.9757186 | 0 | 3.04e-06 | True |
| openmm-ir | numpy-ir | improper | 0.05472210537 | 0.05472210537 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj14 | 3.645523707 | 3.645523707 | 4.44e-16 | 1.04e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -15.02097654 | -15.02097654 | 0 | 1.15e-06 | True |
| openmm-ir | numpy-ir | total | 334.9431127 | 334.9431127 | 0 | 3.04e-06 | True |
| openmm-xml | numpy-ir | bond | 61.14815051 | 61.14815051 | 0 | 1.61e-06 | True |
| openmm-xml | numpy-ir | angle | 81.13997439 | 81.13997439 | 1.42e-14 | 1.81e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 204.0304407 | 204.0304407 | 0 | 3.04e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | -11.37545301 | -11.37545283 | 1.78e-07 | 1.11e-06 | True |
| openmm-xml | numpy-ir | total | 334.9431126 | 334.9431127 | 1.78e-07 | 3.04e-06 | True |
| openmm-prmtop | numpy-ir | bond | 61.14815051 | 61.14815051 | 0 | 0.000106 | True |
| openmm-prmtop | numpy-ir | angle | 81.13997439 | 81.13997439 | 1.42e-14 | 0.000108 | True |
| openmm-prmtop | numpy-ir | torsion_all | 204.0304407 | 204.0304407 | 0 | 0.00012 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | -11.37545302 | -11.37545283 | 1.93e-07 | 0.000101 | True |
| openmm-prmtop | numpy-ir | total | 334.9431126 | 334.9431127 | 1.93e-07 | 0.00012 | True |
| openmm-top | numpy-ir | bond | 61.14815051 | 61.14815051 | 2.84e-14 | 0.000106 | True |
| openmm-top | numpy-ir | angle | 81.1399746 | 81.13997439 | 2.03e-07 | 0.000108 | True |
| openmm-top | numpy-ir | torsion_all | 204.0304407 | 204.0304407 | 7.11e-09 | 0.00012 | True |
| openmm-top | numpy-ir | nonbonded_all | -11.37545299 | -11.37545283 | 1.61e-07 | 0.000101 | True |
| openmm-top | numpy-ir | total | 334.9431128 | 334.9431127 | 4.97e-08 | 0.00012 | True |
| sander-ours | numpy-ir | bond | 14.6148 | 14.61475873 | 4.13e-05 | 0.000201 | True |
| sander-ours | numpy-ir | angle | 19.3929 | 19.39291931 | 1.93e-05 | 0.000202 | True |
| sander-ours | numpy-ir | dihed | 48.7644 | 48.76444567 | 4.57e-05 | 0.000205 | True |
| sander-ours | numpy-ir | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | vdw14 | 0.8713 | 0.8713010772 | 1.08e-06 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -3.59 | -3.590099555 | 9.96e-05 | 0.00038 | True |
| sander-tleap | sander-ours | bond | 14.6148 | 14.6148 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | angle | 19.3929 | 19.3929 | 0 | 0.000202 | True |
| sander-tleap | sander-ours | dihed | 48.7644 | 48.7644 | 0 | 0.000205 | True |
| sander-tleap | sander-ours | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | 0.8713 | 0.8713 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -3.59 | -3.59 | 0 | 0.0002 | True |
| gmx | numpy-ir | bond | 61.14809036 | 61.14815051 | 6.02e-05 | 0.000106 | True |
| gmx | numpy-ir | angle | 81.14002991 | 81.13997439 | 5.55e-05 | 0.000108 | True |
| gmx | numpy-ir | proper | 203.9757233 | 203.9757186 | 4.7e-06 | 0.00012 | True |
| gmx | numpy-ir | improper | 0.05472168326 | 0.05472210537 | 4.22e-07 | 0.0001 | True |
| gmx | numpy-ir | lj | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj14 | 3.645513535 | 3.645523707 | 1.02e-05 | 0.0001 | True |
| gmx | numpy-ir | coulomb14 | -15.02097416 | -15.02097654 | 2.38e-06 | 0.000102 | True |
| gmx | numpy-ir | total | 334.9431152 | 334.9431127 | 2.49e-06 | 0.000133 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 2.3e-07 | 0.001 | True |
| openmm-top | openmm-ir | 6.61e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 3  passed=False

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 79.6765043 | 79.6765043 | 0 | 1.8e-06 | True |
| numpy-ir | reference-itself | angle | 31.37476311 | 31.37476311 | 0 | 1.31e-06 | True |
| numpy-ir | reference-itself | proper | 191.3070689 | 191.3070689 | 0 | 2.91e-06 | True |
| numpy-ir | reference-itself | improper | 0.03547588222 | 0.03547588222 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj14 | 1.732409674 | 1.732409674 | 0 | 1.02e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -14.87343537 | -14.87343537 | 0 | 1.15e-06 | True |
| numpy-ir | reference-itself | total | 289.2527865 | 289.2527865 | 0 | 3.89e-06 | True |
| openmm-ir | numpy-ir | bond | 79.6765043 | 79.6765043 | 1.42e-14 | 1.8e-06 | True |
| openmm-ir | numpy-ir | angle | 31.37476311 | 31.37476311 | 3.55e-15 | 1.31e-06 | True |
| openmm-ir | numpy-ir | proper | 191.3070689 | 191.3070689 | 0 | 2.91e-06 | True |
| openmm-ir | numpy-ir | improper | 0.03547588222 | 0.03547588222 | 5.07e-16 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj14 | 1.732409674 | 1.732409674 | 1.11e-15 | 1.02e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -14.87343537 | -14.87343537 | 0 | 1.15e-06 | True |
| openmm-ir | numpy-ir | total | 289.2527865 | 289.2527865 | 0 | 2.91e-06 | True |
| openmm-xml | numpy-ir | bond | 79.6765043 | 79.6765043 | 1.42e-14 | 1.8e-06 | True |
| openmm-xml | numpy-ir | angle | 31.37476311 | 31.37476311 | 3.55e-15 | 1.31e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 191.3425448 | 191.3425448 | 0 | 2.91e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | -13.14102587 | -13.1410257 | 1.76e-07 | 1.13e-06 | True |
| openmm-xml | numpy-ir | total | 289.2527863 | 289.2527865 | 1.76e-07 | 2.91e-06 | True |
| openmm-prmtop | numpy-ir | bond | 79.6765043 | 79.6765043 | 2.84e-14 | 0.000108 | True |
| openmm-prmtop | numpy-ir | angle | 31.37476311 | 31.37476311 | 0 | 0.000103 | True |
| openmm-prmtop | numpy-ir | torsion_all | 191.3425448 | 191.3425448 | 0 | 0.000119 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | -13.14102589 | -13.1410257 | 1.9e-07 | 0.000101 | True |
| openmm-prmtop | numpy-ir | total | 289.2527863 | 289.2527865 | 1.9e-07 | 0.000119 | True |
| openmm-top | numpy-ir | bond | 79.6765043 | 79.6765043 | 9.95e-14 | 0.000108 | True |
| openmm-top | numpy-ir | angle | 31.37476324 | 31.37476311 | 1.34e-07 | 0.000103 | True |
| openmm-top | numpy-ir | torsion_all | 191.3425448 | 191.3425448 | 1.26e-09 | 0.000119 | True |
| openmm-top | numpy-ir | nonbonded_all | -13.14102586 | -13.1410257 | 1.6e-07 | 0.000101 | True |
| openmm-top | numpy-ir | total | 289.2527865 | 289.2527865 | 2.49e-08 | 0.000119 | True |
| sander-ours | numpy-ir | bond | 19.0431 | 19.04314156 | 4.16e-05 | 0.000202 | True |
| sander-ours | numpy-ir | angle | 7.4987 | 7.498748353 | 4.84e-05 | 0.000201 | True |
| sander-ours | numpy-ir | dihed | 45.732 | 45.73196577 | 3.42e-05 | 0.000205 | True |
| sander-ours | numpy-ir | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | vdw14 | 0.4141 | 0.4140558493 | 4.42e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -3.5547 | -3.55483637 | 0.000136 | 0.000378 | True |
| sander-tleap | sander-ours | bond | 19.0431 | 19.0431 | 0 | 0.000202 | True |
| sander-tleap | sander-ours | angle | 7.4987 | 7.4987 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | dihed | 45.732 | 45.732 | 0 | 0.000205 | True |
| sander-tleap | sander-ours | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | 0.4141 | 0.4141 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -3.5547 | -3.5547 | 0 | 0.0002 | True |
| gmx | numpy-ir | bond | 79.67641449 | 79.6765043 | 8.98e-05 | 0.000108 | True |
| gmx | numpy-ir | angle | 31.37478638 | 31.37476311 | 2.33e-05 | 0.000103 | True |
| gmx | numpy-ir | proper | 191.3070526 | 191.3070689 | 1.63e-05 | 0.000119 | True |
| gmx | numpy-ir | improper | 0.03547590599 | 0.03547588222 | 2.38e-08 | 0.0001 | True |
| gmx | numpy-ir | lj | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj14 | 1.732406497 | 1.732409674 | 3.18e-06 | 0.0001 | True |
| gmx | numpy-ir | coulomb14 | -14.87343407 | -14.87343537 | 1.3e-06 | 0.000101 | True |
| gmx | numpy-ir | total | 289.252655 | 289.2527865 | 0.000131 | 0.000129 | False |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 1.52e-07 | 0.001 | True |
| openmm-top | openmm-ir | 6.72e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

## methylammonium
### config 0  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 0.09632231627 | 0.09632231627 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | angle | 0.4434451887 | 0.4434451887 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | proper | 1.040343704e-05 | 1.040343704e-05 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | improper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj14 | -0.04580930017 | -0.04580930017 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb14 | 227.1018033 | 227.1018033 | 0 | 3.27e-06 | True |
| numpy-ir | reference-itself | total | 227.5957719 | 227.5957719 | 0 | 3.28e-06 | True |
| openmm-ir | numpy-ir | bond | 0.09632231627 | 0.09632231627 | 1.11e-15 | 1e-06 | True |
| openmm-ir | numpy-ir | angle | 0.4434451887 | 0.4434451887 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | proper | 1.040343704e-05 | 1.040343704e-05 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | improper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj14 | -0.04580930017 | -0.04580930017 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | 227.1018033 | 227.1018033 | 0 | 3.27e-06 | True |
| openmm-ir | numpy-ir | total | 227.5957719 | 227.5957719 | 0 | 3.27e-06 | True |
| openmm-xml | numpy-ir | bond | 0.09632231627 | 0.09632231627 | 1.11e-15 | 1e-06 | True |
| openmm-xml | numpy-ir | angle | 0.4434451887 | 0.4434451887 | 0 | 1e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 1.040343704e-05 | 1.040343704e-05 | 0 | 1e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | 227.0559967 | 227.055994 | 2.69e-06 | 3.27e-06 | True |
| openmm-xml | numpy-ir | total | 227.5957746 | 227.5957719 | 2.69e-06 | 3.27e-06 | True |
| openmm-prmtop | numpy-ir | bond | 0.09632231627 | 0.09632231627 | 1.1e-15 | 0.0001 | True |
| openmm-prmtop | numpy-ir | angle | 0.4434451887 | 0.4434451887 | 0 | 0.0001 | True |
| openmm-prmtop | numpy-ir | torsion_all | 1.040343704e-05 | 1.040343704e-05 | 0 | 0.0001 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | 227.0559967 | 227.055994 | 2.69e-06 | 0.000123 | True |
| openmm-prmtop | numpy-ir | total | 227.5957746 | 227.5957719 | 2.69e-06 | 0.000123 | True |
| openmm-top | numpy-ir | bond | 0.09632231627 | 0.09632231627 | 1.11e-15 | 0.0001 | True |
| openmm-top | numpy-ir | angle | 0.4434451991 | 0.4434451887 | 1.04e-08 | 0.0001 | True |
| openmm-top | numpy-ir | torsion_all | 1.040343704e-05 | 1.040343704e-05 | 0 | 0.0001 | True |
| openmm-top | numpy-ir | nonbonded_all | 227.0559967 | 227.055994 | 2.68e-06 | 0.000123 | True |
| openmm-top | numpy-ir | total | 227.5957746 | 227.5957719 | 2.69e-06 | 0.000123 | True |
| sander-ours | numpy-ir | bond | 0.023 | 0.02302158611 | 2.16e-05 | 0.0002 | True |
| sander-ours | numpy-ir | angle | 0.106 | 0.1059859438 | 1.41e-05 | 0.0002 | True |
| sander-ours | numpy-ir | dihed | 0 | 2.486481127e-06 | 2.49e-06 | 0.0002 | True |
| sander-ours | numpy-ir | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | vdw14 | -0.0109 | -0.01094868551 | 4.87e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | 54.2768 | 54.27863368 | 0.00183 | 0.00291 | True |
| sander-tleap | sander-ours | bond | 0.023 | 0.023 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | angle | 0.106 | 0.106 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | dihed | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | -0.0109 | -0.0109 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | 54.2768 | 54.2768 | 0 | 0.000205 | True |
| gmx | numpy-ir | bond | 0.0963236168 | 0.09632231627 | 1.3e-06 | 0.0001 | True |
| gmx | numpy-ir | angle | 0.4434462786 | 0.4434451887 | 1.09e-06 | 0.0001 | True |
| gmx | numpy-ir | proper | 1.0438293e-05 | 1.040343704e-05 | 3.49e-08 | 0.0001 | True |
| gmx | numpy-ir | improper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj14 | -0.04580929503 | -0.04580930017 | 5.14e-09 | 0.0001 | True |
| gmx | numpy-ir | coulomb14 | 227.1017914 | 227.1018033 | 1.19e-05 | 0.000123 | True |
| gmx | numpy-ir | total | 227.5957642 | 227.5957719 | 7.75e-06 | 0.000123 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 2.67e-11 | 0.001 | True |
| openmm-top | openmm-ir | 6.19e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 1  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 138.7391139 | 138.7391139 | 0 | 2.39e-06 | True |
| numpy-ir | reference-itself | angle | 27.01909036 | 27.01909036 | 0 | 1.27e-06 | True |
| numpy-ir | reference-itself | proper | 0.1712388679 | 0.1712388679 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | improper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj14 | -0.05065713365 | -0.05065713365 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb14 | 230.4050256 | 230.4050256 | 0 | 3.3e-06 | True |
| numpy-ir | reference-itself | total | 396.2838116 | 396.2838116 | 0 | 4.96e-06 | True |
| openmm-ir | numpy-ir | bond | 138.7391139 | 138.7391139 | 8.53e-14 | 2.39e-06 | True |
| openmm-ir | numpy-ir | angle | 27.01909036 | 27.01909036 | 7.11e-15 | 1.27e-06 | True |
| openmm-ir | numpy-ir | proper | 0.1712388679 | 0.1712388679 | 1.94e-16 | 1e-06 | True |
| openmm-ir | numpy-ir | improper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj14 | -0.05065713365 | -0.05065713365 | 1.39e-17 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | 230.4050256 | 230.4050256 | 2.84e-14 | 3.3e-06 | True |
| openmm-ir | numpy-ir | total | 396.2838116 | 396.2838116 | 5.68e-14 | 3.3e-06 | True |
| openmm-xml | numpy-ir | bond | 138.7391139 | 138.7391139 | 8.53e-14 | 2.39e-06 | True |
| openmm-xml | numpy-ir | angle | 27.01909036 | 27.01909036 | 7.11e-15 | 1.27e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 0.1712388679 | 0.1712388679 | 1.94e-16 | 1e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | 230.3543712 | 230.3543685 | 2.73e-06 | 3.3e-06 | True |
| openmm-xml | numpy-ir | total | 396.2838144 | 396.2838116 | 2.73e-06 | 3.3e-06 | True |
| openmm-prmtop | numpy-ir | bond | 138.7391139 | 138.7391139 | 1.14e-13 | 0.000114 | True |
| openmm-prmtop | numpy-ir | angle | 27.01909036 | 27.01909036 | 7.11e-15 | 0.000103 | True |
| openmm-prmtop | numpy-ir | torsion_all | 0.1712388679 | 0.1712388679 | 1.94e-16 | 0.0001 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | 230.3543712 | 230.3543685 | 2.73e-06 | 0.000123 | True |
| openmm-prmtop | numpy-ir | total | 396.2838144 | 396.2838116 | 2.73e-06 | 0.000123 | True |
| openmm-top | numpy-ir | bond | 138.7391139 | 138.7391139 | 1.14e-13 | 0.000114 | True |
| openmm-top | numpy-ir | angle | 27.01909045 | 27.01909036 | 9.19e-08 | 0.000103 | True |
| openmm-top | numpy-ir | torsion_all | 0.1712388679 | 0.1712388679 | 1.94e-16 | 0.0001 | True |
| openmm-top | numpy-ir | nonbonded_all | 230.3543712 | 230.3543685 | 2.72e-06 | 0.000123 | True |
| openmm-top | numpy-ir | total | 396.2838145 | 396.2838116 | 2.81e-06 | 0.000123 | True |
| sander-ours | numpy-ir | bond | 33.1594 | 33.15944406 | 4.41e-05 | 0.000203 | True |
| sander-ours | numpy-ir | angle | 6.4577 | 6.457717581 | 1.76e-05 | 0.000201 | True |
| sander-ours | numpy-ir | dihed | 0.0409 | 0.04092707169 | 2.71e-05 | 0.0002 | True |
| sander-ours | numpy-ir | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | vdw14 | -0.0121 | -0.01210734552 | 7.35e-06 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | 55.0662 | 55.06812276 | 0.00192 | 0.00295 | True |
| sander-tleap | sander-ours | bond | 33.1594 | 33.1594 | 0 | 0.000203 | True |
| sander-tleap | sander-ours | angle | 6.4577 | 6.4577 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | dihed | 0.0409 | 0.0409 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | -0.0121 | -0.0121 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | 55.0662 | 55.0662 | 0 | 0.000206 | True |
| gmx | numpy-ir | bond | 138.73909 | 138.7391139 | 2.4e-05 | 0.000114 | True |
| gmx | numpy-ir | angle | 27.019104 | 27.01909036 | 1.36e-05 | 0.000103 | True |
| gmx | numpy-ir | proper | 0.1712389141 | 0.1712388679 | 4.62e-08 | 0.0001 | True |
| gmx | numpy-ir | improper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj14 | -0.05065713078 | -0.05065713365 | 2.87e-09 | 0.0001 | True |
| gmx | numpy-ir | coulomb14 | 230.4050293 | 230.4050256 | 3.69e-06 | 0.000123 | True |
| gmx | numpy-ir | total | 396.2838135 | 396.2838116 | 1.83e-06 | 0.00014 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 2.91e-11 | 0.001 | True |
| openmm-top | openmm-ir | 5.84e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 2  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 36.7595768 | 36.7595768 | 0 | 1.37e-06 | True |
| numpy-ir | reference-itself | angle | 15.47669894 | 15.47669894 | 0 | 1.15e-06 | True |
| numpy-ir | reference-itself | proper | 0.1801345844 | 0.1801345844 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | improper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj14 | -0.04819571467 | -0.04819571467 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb14 | 228.9277561 | 228.9277561 | 0 | 3.29e-06 | True |
| numpy-ir | reference-itself | total | 281.2959707 | 281.2959707 | 0 | 3.81e-06 | True |
| openmm-ir | numpy-ir | bond | 36.7595768 | 36.7595768 | 7.11e-15 | 1.37e-06 | True |
| openmm-ir | numpy-ir | angle | 15.47669894 | 15.47669894 | 1.78e-14 | 1.15e-06 | True |
| openmm-ir | numpy-ir | proper | 0.1801345844 | 0.1801345844 | 1.53e-15 | 1e-06 | True |
| openmm-ir | numpy-ir | improper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj14 | -0.04819571467 | -0.04819571467 | 6.94e-18 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | 228.9277561 | 228.9277561 | 2.84e-14 | 3.29e-06 | True |
| openmm-ir | numpy-ir | total | 281.2959707 | 281.2959707 | 5.68e-14 | 3.29e-06 | True |
| openmm-xml | numpy-ir | bond | 36.7595768 | 36.7595768 | 7.11e-15 | 1.37e-06 | True |
| openmm-xml | numpy-ir | angle | 15.47669894 | 15.47669894 | 1.78e-14 | 1.15e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 0.1801345844 | 0.1801345844 | 1.53e-15 | 1e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | 228.8795631 | 228.8795604 | 2.71e-06 | 3.29e-06 | True |
| openmm-xml | numpy-ir | total | 281.2959734 | 281.2959707 | 2.71e-06 | 3.29e-06 | True |
| openmm-prmtop | numpy-ir | bond | 36.7595768 | 36.7595768 | 7.11e-15 | 0.000104 | True |
| openmm-prmtop | numpy-ir | angle | 15.47669894 | 15.47669894 | 1.78e-14 | 0.000102 | True |
| openmm-prmtop | numpy-ir | torsion_all | 0.1801345844 | 0.1801345844 | 1.53e-15 | 0.0001 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | 228.8795631 | 228.8795604 | 2.71e-06 | 0.000123 | True |
| openmm-prmtop | numpy-ir | total | 281.2959734 | 281.2959707 | 2.71e-06 | 0.000123 | True |
| openmm-top | numpy-ir | bond | 36.7595768 | 36.7595768 | 0 | 0.000104 | True |
| openmm-top | numpy-ir | angle | 15.47669892 | 15.47669894 | 2.32e-08 | 0.000102 | True |
| openmm-top | numpy-ir | torsion_all | 0.1801345844 | 0.1801345844 | 1.53e-15 | 0.0001 | True |
| openmm-top | numpy-ir | nonbonded_all | 228.8795631 | 228.8795604 | 2.7e-06 | 0.000123 | True |
| openmm-top | numpy-ir | total | 281.2959734 | 281.2959707 | 2.68e-06 | 0.000123 | True |
| sander-ours | numpy-ir | bond | 8.7858 | 8.785749714 | 5.03e-05 | 0.000201 | True |
| sander-ours | numpy-ir | angle | 3.699 | 3.699019824 | 1.98e-05 | 0.0002 | True |
| sander-ours | numpy-ir | dihed | 0.0431 | 0.04305319896 | 4.68e-05 | 0.0002 | True |
| sander-ours | numpy-ir | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | vdw14 | -0.0115 | -0.01151905226 | 1.91e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | 54.7132 | 54.71504687 | 0.00185 | 0.00294 | True |
| sander-tleap | sander-ours | bond | 8.7858 | 8.7858 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | angle | 3.699 | 3.699 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | dihed | 0.0431 | 0.0431 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | -0.0115 | -0.0115 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | 54.7132 | 54.7132 | 0 | 0.000205 | True |
| gmx | numpy-ir | bond | 36.75956345 | 36.7595768 | 1.34e-05 | 0.000104 | True |
| gmx | numpy-ir | angle | 15.47669983 | 15.47669894 | 8.84e-07 | 0.000102 | True |
| gmx | numpy-ir | proper | 0.1801347137 | 0.1801345844 | 1.29e-07 | 0.0001 | True |
| gmx | numpy-ir | improper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj14 | -0.04819569364 | -0.04819571467 | 2.1e-08 | 0.0001 | True |
| gmx | numpy-ir | coulomb14 | 228.9277344 | 228.9277561 | 2.17e-05 | 0.000123 | True |
| gmx | numpy-ir | total | 281.295929 | 281.2959707 | 4.18e-05 | 0.000128 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 2.82e-11 | 0.001 | True |
| openmm-top | openmm-ir | 5.99e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 3  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 13.78383517 | 13.78383517 | 0 | 1.14e-06 | True |
| numpy-ir | reference-itself | angle | 15.62546729 | 15.62546729 | 0 | 1.16e-06 | True |
| numpy-ir | reference-itself | proper | 0.2046293527 | 0.2046293527 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | improper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj14 | -0.04540661059 | -0.04540661059 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb14 | 226.5016797 | 226.5016797 | 0 | 3.27e-06 | True |
| numpy-ir | reference-itself | total | 256.0702049 | 256.0702049 | 0 | 3.56e-06 | True |
| openmm-ir | numpy-ir | bond | 13.78383517 | 13.78383517 | 4.09e-14 | 1.14e-06 | True |
| openmm-ir | numpy-ir | angle | 15.62546729 | 15.62546729 | 2.13e-14 | 1.16e-06 | True |
| openmm-ir | numpy-ir | proper | 0.2046293527 | 0.2046293527 | 2.78e-16 | 1e-06 | True |
| openmm-ir | numpy-ir | improper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj14 | -0.04540661059 | -0.04540661059 | 6.94e-18 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | 226.5016797 | 226.5016797 | 0 | 3.27e-06 | True |
| openmm-ir | numpy-ir | total | 256.0702049 | 256.0702049 | 5.68e-14 | 3.27e-06 | True |
| openmm-xml | numpy-ir | bond | 13.78383517 | 13.78383517 | 4.09e-14 | 1.14e-06 | True |
| openmm-xml | numpy-ir | angle | 15.62546729 | 15.62546729 | 2.13e-14 | 1.16e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 0.2046293527 | 0.2046293527 | 2.78e-16 | 1e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | 226.4562758 | 226.4562731 | 2.68e-06 | 3.26e-06 | True |
| openmm-xml | numpy-ir | total | 256.0702076 | 256.0702049 | 2.68e-06 | 3.27e-06 | True |
| openmm-prmtop | numpy-ir | bond | 13.78383517 | 13.78383517 | 4.26e-14 | 0.000101 | True |
| openmm-prmtop | numpy-ir | angle | 15.62546729 | 15.62546729 | 2.13e-14 | 0.000102 | True |
| openmm-prmtop | numpy-ir | torsion_all | 0.2046293527 | 0.2046293527 | 2.78e-16 | 0.0001 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | 226.4562758 | 226.4562731 | 2.68e-06 | 0.000123 | True |
| openmm-prmtop | numpy-ir | total | 256.0702076 | 256.0702049 | 2.68e-06 | 0.000123 | True |
| openmm-top | numpy-ir | bond | 13.78383517 | 13.78383517 | 4.09e-14 | 0.000101 | True |
| openmm-top | numpy-ir | angle | 15.6254673 | 15.62546729 | 9.92e-09 | 0.000102 | True |
| openmm-top | numpy-ir | torsion_all | 0.2046293527 | 0.2046293527 | 2.78e-16 | 0.0001 | True |
| openmm-top | numpy-ir | nonbonded_all | 226.4562758 | 226.4562731 | 2.67e-06 | 0.000123 | True |
| openmm-top | numpy-ir | total | 256.0702076 | 256.0702049 | 2.68e-06 | 0.000123 | True |
| sander-ours | numpy-ir | bond | 3.2944 | 3.294415673 | 1.57e-05 | 0.0002 | True |
| sander-ours | numpy-ir | angle | 3.7346 | 3.734576312 | 2.37e-05 | 0.0002 | True |
| sander-ours | numpy-ir | dihed | 0.0489 | 0.04890758907 | 7.59e-06 | 0.0002 | True |
| sander-ours | numpy-ir | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | vdw14 | -0.0109 | -0.01085244039 | 4.76e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | 54.1333 | 54.13520069 | 0.0019 | 0.00291 | True |
| sander-tleap | sander-ours | bond | 3.2944 | 3.2944 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | angle | 3.7346 | 3.7346 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | dihed | 0.0489 | 0.0489 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | -0.0109 | -0.0109 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | 54.1333 | 54.1333 | 0 | 0.000205 | True |
| gmx | numpy-ir | bond | 13.78384781 | 13.78383517 | 1.26e-05 | 0.000101 | True |
| gmx | numpy-ir | angle | 15.62547493 | 15.62546729 | 7.64e-06 | 0.000102 | True |
| gmx | numpy-ir | proper | 0.2046293616 | 0.2046293527 | 8.98e-09 | 0.0001 | True |
| gmx | numpy-ir | improper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj14 | -0.0454065986 | -0.04540661059 | 1.2e-08 | 0.0001 | True |
| gmx | numpy-ir | coulomb14 | 226.5016479 | 226.5016797 | 3.17e-05 | 0.000123 | True |
| gmx | numpy-ir | total | 256.0701904 | 256.0702049 | 1.45e-05 | 0.000126 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 2.8e-11 | 0.001 | True |
| openmm-top | openmm-ir | 6.09e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

## chlorobenzene
### config 0  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 1.166230953 | 1.166230953 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | angle | 0.09933909937 | 0.09933909937 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | proper | 0.0001286021831 | 0.0001286021831 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | improper | 1.303640113e-05 | 1.303640113e-05 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | -1.156288482 | -1.156288482 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | coulomb | 5.839891223 | 5.839891223 | 0 | 1.06e-06 | True |
| numpy-ir | reference-itself | lj14 | 12.26611398 | 12.26611398 | 0 | 1.12e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -3.554685334 | -3.554685334 | 0 | 1.04e-06 | True |
| numpy-ir | reference-itself | total | 14.66074308 | 14.66074308 | 0 | 1.15e-06 | True |
| openmm-ir | numpy-ir | bond | 1.166230953 | 1.166230953 | 6.66e-16 | 1.01e-06 | True |
| openmm-ir | numpy-ir | angle | 0.09933909937 | 0.09933909937 | 2.15e-15 | 1e-06 | True |
| openmm-ir | numpy-ir | proper | 0.0001286021831 | 0.0001286021831 | 1.68e-15 | 1e-06 | True |
| openmm-ir | numpy-ir | improper | 1.303640113e-05 | 1.303640113e-05 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | -1.156288482 | -1.156288482 | 0 | 1.01e-06 | True |
| openmm-ir | numpy-ir | coulomb | 5.839891223 | 5.839891223 | 0 | 1.06e-06 | True |
| openmm-ir | numpy-ir | lj14 | 12.26611398 | 12.26611398 | 0 | 1.12e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -3.554685334 | -3.554685334 | 1.42e-14 | 1.04e-06 | True |
| openmm-ir | numpy-ir | total | 14.66074308 | 14.66074308 | 1.42e-14 | 1.12e-06 | True |
| openmm-xml | numpy-ir | bond | 1.166230953 | 1.166230953 | 6.66e-16 | 1.01e-06 | True |
| openmm-xml | numpy-ir | angle | 0.09933909937 | 0.09933909937 | 2.15e-15 | 1e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 0.0001416385842 | 0.0001416385842 | 1.68e-15 | 1e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | 13.39503141 | 13.39503139 | 2.7e-08 | 1.13e-06 | True |
| openmm-xml | numpy-ir | total | 14.6607431 | 14.66074308 | 2.7e-08 | 1.12e-06 | True |
| openmm-prmtop | numpy-ir | bond | 1.166230953 | 1.166230953 | 6.66e-16 | 0.0001 | True |
| openmm-prmtop | numpy-ir | angle | 0.09933909937 | 0.09933909937 | 2.15e-15 | 0.0001 | True |
| openmm-prmtop | numpy-ir | torsion_all | 0.0001416385842 | 0.0001416385842 | 1.68e-15 | 0.0001 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | 13.39503141 | 13.39503139 | 2.76e-08 | 0.000101 | True |
| openmm-prmtop | numpy-ir | total | 14.66074311 | 14.66074308 | 2.76e-08 | 0.000101 | True |
| openmm-top | numpy-ir | bond | 1.166230953 | 1.166230953 | 1.35e-14 | 0.0001 | True |
| openmm-top | numpy-ir | angle | 0.09933911358 | 0.09933909937 | 1.42e-08 | 0.0001 | True |
| openmm-top | numpy-ir | torsion_all | 0.0001416385733 | 0.0001416385842 | 1.09e-11 | 0.0001 | True |
| openmm-top | numpy-ir | nonbonded_all | 13.39503139 | 13.39503139 | 1.3e-10 | 0.000101 | True |
| openmm-top | numpy-ir | total | 14.66074309 | 14.66074308 | 1.41e-08 | 0.000101 | True |
| sander-ours | numpy-ir | bond | 0.2787 | 0.2787358875 | 3.59e-05 | 0.0002 | True |
| sander-ours | numpy-ir | angle | 0.0237 | 0.02374261457 | 4.26e-05 | 0.0002 | True |
| sander-ours | numpy-ir | dihed | 0 | 3.385243409e-05 | 3.39e-05 | 0.0002 | True |
| sander-ours | numpy-ir | vdw | -0.2764 | -0.2763595799 | 4.04e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 1.3957 | 1.395767501 | 6.75e-05 | 0.00027 | True |
| sander-ours | numpy-ir | vdw14 | 2.9317 | 2.931671601 | 2.84e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -0.8496 | -0.849590185 | 9.81e-06 | 0.000242 | True |
| sander-tleap | sander-ours | bond | 0.2787 | 0.2787 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | angle | 0.0237 | 0.0237 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | dihed | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw | -0.2764 | -0.2764 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 1.3957 | 1.3957 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | 2.9317 | 2.9317 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -0.8496 | -0.8496 | 0 | 0.0002 | True |
| gmx | numpy-ir | bond | 1.166232347 | 1.166230953 | 1.39e-06 | 0.0001 | True |
| gmx | numpy-ir | angle | 0.0993386358 | 0.09933909937 | 4.64e-07 | 0.0001 | True |
| gmx | numpy-ir | proper | 0.000127467327 | 0.0001286021831 | 1.13e-06 | 0.0001 | True |
| gmx | numpy-ir | improper | 1.3167572e-05 | 1.303640113e-05 | 1.31e-07 | 0.0001 | True |
| gmx | numpy-ir | lj | -1.156288624 | -1.156288482 | 1.42e-07 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 5.839890957 | 5.839891223 | 2.66e-07 | 0.000101 | True |
| gmx | numpy-ir | lj14 | 12.266119 | 12.26611398 | 5.02e-06 | 0.000101 | True |
| gmx | numpy-ir | coulomb14 | -3.554683685 | -3.554685334 | 1.65e-06 | 0.0001 | True |
| gmx | numpy-ir | total | 14.66074944 | 14.66074308 | 6.36e-06 | 0.000101 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 7.3e-08 | 0.001 | True |
| openmm-top | openmm-ir | 1.92e-07 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 1  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 42.01667352 | 42.01667352 | 0 | 1.42e-06 | True |
| numpy-ir | reference-itself | angle | 19.86722601 | 19.86722601 | 0 | 1.2e-06 | True |
| numpy-ir | reference-itself | proper | 18.57472493 | 18.57472493 | 0 | 1.19e-06 | True |
| numpy-ir | reference-itself | improper | 2.569099414 | 2.569099414 | 0 | 1.03e-06 | True |
| numpy-ir | reference-itself | lj | -1.202163624 | -1.202163624 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | coulomb | 5.927143804 | 5.927143804 | 0 | 1.06e-06 | True |
| numpy-ir | reference-itself | lj14 | 13.83416165 | 13.83416165 | 0 | 1.14e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -4.261994596 | -4.261994596 | 0 | 1.04e-06 | True |
| numpy-ir | reference-itself | total | 97.32487111 | 97.32487111 | 0 | 1.97e-06 | True |
| openmm-ir | numpy-ir | bond | 42.01667352 | 42.01667352 | 1.42e-14 | 1.42e-06 | True |
| openmm-ir | numpy-ir | angle | 19.86722601 | 19.86722601 | 3.55e-15 | 1.2e-06 | True |
| openmm-ir | numpy-ir | proper | 18.57472493 | 18.57472493 | 2.84e-14 | 1.19e-06 | True |
| openmm-ir | numpy-ir | improper | 2.569099414 | 2.569099414 | 8.88e-16 | 1.03e-06 | True |
| openmm-ir | numpy-ir | lj | -1.202163624 | -1.202163624 | 0 | 1.01e-06 | True |
| openmm-ir | numpy-ir | coulomb | 5.927143804 | 5.927143804 | 0 | 1.06e-06 | True |
| openmm-ir | numpy-ir | lj14 | 13.83416165 | 13.83416165 | 1.78e-15 | 1.14e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -4.261994596 | -4.261994596 | 5.33e-15 | 1.04e-06 | True |
| openmm-ir | numpy-ir | total | 97.32487111 | 97.32487111 | 4.26e-14 | 1.42e-06 | True |
| openmm-xml | numpy-ir | bond | 42.01667352 | 42.01667352 | 1.42e-14 | 1.42e-06 | True |
| openmm-xml | numpy-ir | angle | 19.86722601 | 19.86722601 | 3.55e-15 | 1.2e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 21.14382435 | 21.14382435 | 2.84e-14 | 1.21e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | 14.29714725 | 14.29714723 | 1.97e-08 | 1.14e-06 | True |
| openmm-xml | numpy-ir | total | 97.32487113 | 97.32487111 | 1.97e-08 | 1.42e-06 | True |
| openmm-prmtop | numpy-ir | bond | 42.01667352 | 42.01667352 | 2.13e-14 | 0.000104 | True |
| openmm-prmtop | numpy-ir | angle | 19.86722601 | 19.86722601 | 0 | 0.000102 | True |
| openmm-prmtop | numpy-ir | torsion_all | 21.14382435 | 21.14382435 | 2.84e-14 | 0.000102 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | 14.29714725 | 14.29714723 | 2.11e-08 | 0.000101 | True |
| openmm-prmtop | numpy-ir | total | 97.32487113 | 97.32487111 | 2.11e-08 | 0.000104 | True |
| openmm-top | numpy-ir | bond | 42.01667352 | 42.01667352 | 3.55e-14 | 0.000104 | True |
| openmm-top | numpy-ir | angle | 19.86722601 | 19.86722601 | 6.06e-09 | 0.000102 | True |
| openmm-top | numpy-ir | torsion_all | 21.14382435 | 21.14382435 | 7.16e-09 | 0.000102 | True |
| openmm-top | numpy-ir | nonbonded_all | 14.29714723 | 14.29714723 | 7.73e-09 | 0.000101 | True |
| openmm-top | numpy-ir | total | 97.3248711 | 97.32487111 | 6.64e-09 | 0.000104 | True |
| sander-ours | numpy-ir | bond | 10.0422 | 10.04222598 | 2.6e-05 | 0.000201 | True |
| sander-ours | numpy-ir | angle | 4.7484 | 4.748380978 | 1.9e-05 | 0.0002 | True |
| sander-ours | numpy-ir | dihed | 5.0535 | 5.053495303 | 4.7e-06 | 0.000201 | True |
| sander-ours | numpy-ir | vdw | -0.2873 | -0.287324002 | 2.4e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 1.4166 | 1.416621368 | 2.14e-05 | 0.000271 | True |
| sander-ours | numpy-ir | vdw14 | 3.3064 | 3.306443989 | 4.4e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -1.0186 | -1.018641156 | 4.12e-05 | 0.000251 | True |
| sander-tleap | sander-ours | bond | 10.0422 | 10.0422 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | angle | 4.7484 | 4.7484 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | dihed | 5.0535 | 5.0535 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | vdw | -0.2873 | -0.2873 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 1.4166 | 1.4166 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | 3.3064 | 3.3064 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -1.0186 | -1.0186 | 0 | 0.0002 | True |
| gmx | numpy-ir | bond | 42.01672745 | 42.01667352 | 5.39e-05 | 0.000104 | True |
| gmx | numpy-ir | angle | 19.867239 | 19.86722601 | 1.3e-05 | 0.000102 | True |
| gmx | numpy-ir | proper | 18.57472801 | 18.57472493 | 3.08e-06 | 0.000102 | True |
| gmx | numpy-ir | improper | 2.569097757 | 2.569099414 | 1.66e-06 | 0.0001 | True |
| gmx | numpy-ir | lj | -1.202163935 | -1.202163624 | 3.1e-07 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 5.927143574 | 5.927143804 | 2.31e-07 | 0.000101 | True |
| gmx | numpy-ir | lj14 | 13.83414364 | 13.83416165 | 1.8e-05 | 0.000101 | True |
| gmx | numpy-ir | coulomb14 | -4.261990547 | -4.261994596 | 4.05e-06 | 0.0001 | True |
| gmx | numpy-ir | total | 97.32492828 | 97.32487111 | 5.72e-05 | 0.00011 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 9.24e-08 | 0.001 | True |
| openmm-top | openmm-ir | 1.09e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 2  passed=False

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 107.3909043 | 107.3909043 | 0 | 2.07e-06 | True |
| numpy-ir | reference-itself | angle | 16.40529929 | 16.40529929 | 0 | 1.16e-06 | True |
| numpy-ir | reference-itself | proper | 8.484664367 | 8.484664367 | 0 | 1.08e-06 | True |
| numpy-ir | reference-itself | improper | 1.125782576 | 1.125782576 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | lj | -1.162571945 | -1.162571945 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | coulomb | 6.401359456 | 6.401359456 | 0 | 1.06e-06 | True |
| numpy-ir | reference-itself | lj14 | 9.366010817 | 9.366010817 | 0 | 1.09e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -3.593168603 | -3.593168603 | 0 | 1.04e-06 | True |
| numpy-ir | reference-itself | total | 144.4182802 | 144.4182802 | 0 | 2.44e-06 | True |
| openmm-ir | numpy-ir | bond | 107.3909043 | 107.3909043 | 7.11e-14 | 2.07e-06 | True |
| openmm-ir | numpy-ir | angle | 16.40529929 | 16.40529929 | 3.55e-15 | 1.16e-06 | True |
| openmm-ir | numpy-ir | proper | 8.484664367 | 8.484664367 | 0 | 1.08e-06 | True |
| openmm-ir | numpy-ir | improper | 1.125782576 | 1.125782576 | 4.44e-16 | 1.01e-06 | True |
| openmm-ir | numpy-ir | lj | -1.162571945 | -1.162571945 | 2.22e-16 | 1.01e-06 | True |
| openmm-ir | numpy-ir | coulomb | 6.401359456 | 6.401359456 | 0 | 1.06e-06 | True |
| openmm-ir | numpy-ir | lj14 | 9.366010817 | 9.366010817 | 1.78e-15 | 1.09e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -3.593168603 | -3.593168603 | 5.33e-15 | 1.04e-06 | True |
| openmm-ir | numpy-ir | total | 144.4182802 | 144.4182802 | 5.68e-14 | 2.07e-06 | True |
| openmm-xml | numpy-ir | bond | 107.3909043 | 107.3909043 | 7.11e-14 | 2.07e-06 | True |
| openmm-xml | numpy-ir | angle | 16.40529929 | 16.40529929 | 3.55e-15 | 1.16e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 9.610446943 | 9.610446943 | 1.78e-15 | 1.1e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | 11.01162976 | 11.01162973 | 3.32e-08 | 1.11e-06 | True |
| openmm-xml | numpy-ir | total | 144.4182803 | 144.4182802 | 3.32e-08 | 2.07e-06 | True |
| openmm-prmtop | numpy-ir | bond | 107.3909043 | 107.3909043 | 8.53e-14 | 0.000111 | True |
| openmm-prmtop | numpy-ir | angle | 16.40529929 | 16.40529929 | 3.55e-15 | 0.000102 | True |
| openmm-prmtop | numpy-ir | torsion_all | 9.610446943 | 9.610446943 | 1.78e-15 | 0.000101 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | 11.01162976 | 11.01162973 | 3.35e-08 | 0.000101 | True |
| openmm-prmtop | numpy-ir | total | 144.4182803 | 144.4182802 | 3.35e-08 | 0.000111 | True |
| openmm-top | numpy-ir | bond | 107.3909043 | 107.3909043 | 2.27e-13 | 0.000111 | True |
| openmm-top | numpy-ir | angle | 16.4052993 | 16.40529929 | 4.04e-09 | 0.000102 | True |
| openmm-top | numpy-ir | torsion_all | 9.610446938 | 9.610446943 | 4.5e-09 | 0.000101 | True |
| openmm-top | numpy-ir | nonbonded_all | 11.01162973 | 11.01162973 | 5.62e-09 | 0.000101 | True |
| openmm-top | numpy-ir | total | 144.4182802 | 144.4182802 | 5.15e-09 | 0.000111 | True |
| sander-ours | numpy-ir | bond | 25.667 | 25.66704213 | 4.21e-05 | 0.000203 | True |
| sander-ours | numpy-ir | angle | 3.921 | 3.920960634 | 3.94e-05 | 0.0002 | True |
| sander-ours | numpy-ir | dihed | 2.297 | 2.296951946 | 4.81e-05 | 0.0002 | True |
| sander-ours | numpy-ir | vdw | -0.2779 | -0.2778613635 | 3.86e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 1.5299 | 1.529961629 | 6.16e-05 | 0.000276 | True |
| sander-ours | numpy-ir | vdw14 | 2.2385 | 2.23853031 | 3.03e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -0.8588 | -0.8587879069 | 1.21e-05 | 0.000243 | True |
| sander-tleap | sander-ours | bond | 25.667 | 25.667 | 0 | 0.000203 | True |
| sander-tleap | sander-ours | angle | 3.921 | 3.921 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | dihed | 2.297 | 2.297 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw | -0.2779 | -0.2779 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 1.5299 | 1.5299 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | 2.2385 | 2.2385 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -0.8588 | -0.8588 | 0 | 0.0002 | True |
| gmx | numpy-ir | bond | 107.3907928 | 107.3909043 | 0.000111 | 0.000111 | False |
| gmx | numpy-ir | angle | 16.40527534 | 16.40529929 | 2.39e-05 | 0.000102 | True |
| gmx | numpy-ir | proper | 8.484657288 | 8.484664367 | 7.08e-06 | 0.000101 | True |
| gmx | numpy-ir | improper | 1.125781894 | 1.125782576 | 6.82e-07 | 0.0001 | True |
| gmx | numpy-ir | lj | -1.162572026 | -1.162571945 | 8.13e-08 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 6.401360512 | 6.401359456 | 1.06e-06 | 0.000101 | True |
| gmx | numpy-ir | lj14 | 9.366008759 | 9.366010817 | 2.06e-06 | 0.000101 | True |
| gmx | numpy-ir | coulomb14 | -3.593166351 | -3.593168603 | 2.25e-06 | 0.0001 | True |
| gmx | numpy-ir | total | 144.4181519 | 144.4182802 | 0.000128 | 0.000114 | False |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 6.74e-08 | 0.001 | True |
| openmm-top | openmm-ir | 1e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 3  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 20.81102674 | 20.81102674 | 0 | 1.21e-06 | True |
| numpy-ir | reference-itself | angle | 27.51123927 | 27.51123927 | 0 | 1.28e-06 | True |
| numpy-ir | reference-itself | proper | 1.56998816 | 1.56998816 | 0 | 1.02e-06 | True |
| numpy-ir | reference-itself | improper | 0.189222842 | 0.189222842 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | -1.160602053 | -1.160602053 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | coulomb | 6.04935329 | 6.04935329 | 0 | 1.06e-06 | True |
| numpy-ir | reference-itself | lj14 | 12.43837605 | 12.43837605 | 0 | 1.12e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -3.324239133 | -3.324239133 | 0 | 1.03e-06 | True |
| numpy-ir | reference-itself | total | 64.08436517 | 64.08436517 | 0 | 1.64e-06 | True |
| openmm-ir | numpy-ir | bond | 20.81102674 | 20.81102674 | 0 | 1.21e-06 | True |
| openmm-ir | numpy-ir | angle | 27.51123927 | 27.51123927 | 7.11e-15 | 1.28e-06 | True |
| openmm-ir | numpy-ir | proper | 1.56998816 | 1.56998816 | 1.18e-14 | 1.02e-06 | True |
| openmm-ir | numpy-ir | improper | 0.189222842 | 0.189222842 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | -1.160602053 | -1.160602053 | 0 | 1.01e-06 | True |
| openmm-ir | numpy-ir | coulomb | 6.04935329 | 6.04935329 | 1.78e-15 | 1.06e-06 | True |
| openmm-ir | numpy-ir | lj14 | 12.43837605 | 12.43837605 | 7.11e-15 | 1.12e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -3.324239133 | -3.324239133 | 0 | 1.03e-06 | True |
| openmm-ir | numpy-ir | total | 64.08436517 | 64.08436517 | 0 | 1.28e-06 | True |
| openmm-xml | numpy-ir | bond | 20.81102674 | 20.81102674 | 0 | 1.21e-06 | True |
| openmm-xml | numpy-ir | angle | 27.51123927 | 27.51123927 | 7.11e-15 | 1.28e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 1.759211002 | 1.759211002 | 1.18e-14 | 1.02e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | 14.00288819 | 14.00288816 | 3.23e-08 | 1.14e-06 | True |
| openmm-xml | numpy-ir | total | 64.0843652 | 64.08436517 | 3.23e-08 | 1.28e-06 | True |
| openmm-prmtop | numpy-ir | bond | 20.81102674 | 20.81102674 | 7.11e-15 | 0.000102 | True |
| openmm-prmtop | numpy-ir | angle | 27.51123927 | 27.51123927 | 1.07e-14 | 0.000103 | True |
| openmm-prmtop | numpy-ir | torsion_all | 1.759211002 | 1.759211002 | 1.15e-14 | 0.0001 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | 14.00288819 | 14.00288816 | 3.27e-08 | 0.000101 | True |
| openmm-prmtop | numpy-ir | total | 64.0843652 | 64.08436517 | 3.27e-08 | 0.000103 | True |
| openmm-top | numpy-ir | bond | 20.81102674 | 20.81102674 | 3.55e-15 | 0.000102 | True |
| openmm-top | numpy-ir | angle | 27.51123928 | 27.51123927 | 1.13e-08 | 0.000103 | True |
| openmm-top | numpy-ir | torsion_all | 1.759211001 | 1.759211002 | 5.96e-10 | 0.0001 | True |
| openmm-top | numpy-ir | nonbonded_all | 14.00288816 | 14.00288816 | 4.98e-09 | 0.000101 | True |
| openmm-top | numpy-ir | total | 64.08436518 | 64.08436517 | 1.57e-08 | 0.000103 | True |
| sander-ours | numpy-ir | bond | 4.974 | 4.973954766 | 4.52e-05 | 0.0002 | True |
| sander-ours | numpy-ir | angle | 6.5753 | 6.575343993 | 4.4e-05 | 0.000201 | True |
| sander-ours | numpy-ir | dihed | 0.4205 | 0.4204615205 | 3.85e-05 | 0.0002 | True |
| sander-ours | numpy-ir | vdw | -0.2774 | -0.2773905481 | 9.45e-06 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 1.4458 | 1.445830136 | 3.01e-05 | 0.000272 | True |
| sander-ours | numpy-ir | vdw14 | 2.9728 | 2.972843225 | 4.32e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -0.7945 | -0.7945122211 | 1.22e-05 | 0.00024 | True |
| sander-tleap | sander-ours | bond | 4.974 | 4.974 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | angle | 6.5753 | 6.5753 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | dihed | 0.4205 | 0.4205 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw | -0.2774 | -0.2774 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 1.4458 | 1.4458 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | 2.9728 | 2.9728 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -0.7945 | -0.7945 | 0 | 0.0002 | True |
| gmx | numpy-ir | bond | 20.81105232 | 20.81102674 | 2.56e-05 | 0.000102 | True |
| gmx | numpy-ir | angle | 27.51125145 | 27.51123927 | 1.22e-05 | 0.000103 | True |
| gmx | numpy-ir | proper | 1.569988847 | 1.56998816 | 6.87e-07 | 0.0001 | True |
| gmx | numpy-ir | improper | 0.1892226636 | 0.189222842 | 1.78e-07 | 0.0001 | True |
| gmx | numpy-ir | lj | -1.160602212 | -1.160602053 | 1.59e-07 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 6.049351215 | 6.04935329 | 2.07e-06 | 0.000101 | True |
| gmx | numpy-ir | lj14 | 12.43838024 | 12.43837605 | 4.19e-06 | 0.000101 | True |
| gmx | numpy-ir | coulomb14 | -3.324239731 | -3.324239133 | 5.98e-07 | 0.0001 | True |
| gmx | numpy-ir | total | 64.08440399 | 64.08436517 | 3.88e-05 | 0.000106 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 6.65e-08 | 0.001 | True |
| openmm-top | openmm-ir | 6.09e-07 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

## methyl_benzoate
### config 0  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 0.2238031184 | 0.2238031184 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | angle | 3.175003666 | 3.175003666 | 0 | 1.03e-06 | True |
| numpy-ir | reference-itself | proper | 13.39880152 | 13.39880152 | 0 | 1.13e-06 | True |
| numpy-ir | reference-itself | improper | 1.513171003e-05 | 1.513171003e-05 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | 1.619221082 | 1.619221082 | 0 | 1.02e-06 | True |
| numpy-ir | reference-itself | coulomb | -93.35445106 | -93.35445106 | 0 | 1.93e-06 | True |
| numpy-ir | reference-itself | lj14 | 24.78902596 | 24.78902596 | 0 | 1.25e-06 | True |
| numpy-ir | reference-itself | coulomb14 | 85.17636724 | 85.17636724 | 0 | 1.85e-06 | True |
| numpy-ir | reference-itself | total | 35.02778665 | 35.02778665 | 0 | 1.35e-06 | True |
| openmm-ir | numpy-ir | bond | 0.2238031184 | 0.2238031184 | 4.66e-15 | 1e-06 | True |
| openmm-ir | numpy-ir | angle | 3.175003666 | 3.175003666 | 6.22e-15 | 1.03e-06 | True |
| openmm-ir | numpy-ir | proper | 13.39880152 | 13.39880152 | 0 | 1.13e-06 | True |
| openmm-ir | numpy-ir | improper | 1.513171003e-05 | 1.513171003e-05 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | 1.619221082 | 1.619221082 | 1.78e-15 | 1.02e-06 | True |
| openmm-ir | numpy-ir | coulomb | -93.35445106 | -93.35445106 | 0 | 1.93e-06 | True |
| openmm-ir | numpy-ir | lj14 | 24.78902596 | 24.78902596 | 1.07e-14 | 1.25e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | 85.17636724 | 85.17636724 | 0 | 1.85e-06 | True |
| openmm-ir | numpy-ir | total | 35.02778665 | 35.02778665 | 7.11e-15 | 1.93e-06 | True |
| openmm-xml | numpy-ir | bond | 0.2238031184 | 0.2238031184 | 4.66e-15 | 1e-06 | True |
| openmm-xml | numpy-ir | angle | 3.175003666 | 3.175003666 | 6.22e-15 | 1.03e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 13.39881665 | 13.39881665 | 0 | 1.13e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | 18.23016312 | 18.23016322 | 9.68e-08 | 1.18e-06 | True |
| openmm-xml | numpy-ir | total | 35.02778655 | 35.02778665 | 9.68e-08 | 1.93e-06 | True |
| openmm-prmtop | numpy-ir | bond | 0.2238031184 | 0.2238031184 | 4.69e-15 | 0.0001 | True |
| openmm-prmtop | numpy-ir | angle | 3.175003666 | 3.175003666 | 1.47e-14 | 0.0001 | True |
| openmm-prmtop | numpy-ir | torsion_all | 13.39881665 | 13.39881665 | 1.78e-15 | 0.000101 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | 18.23016309 | 18.23016322 | 1.21e-07 | 0.000102 | True |
| openmm-prmtop | numpy-ir | total | 35.02778653 | 35.02778665 | 1.21e-07 | 0.000109 | True |
| openmm-top | numpy-ir | bond | 0.2238031184 | 0.2238031184 | 1.36e-15 | 0.0001 | True |
| openmm-top | numpy-ir | angle | 3.1750037 | 3.175003666 | 3.39e-08 | 0.0001 | True |
| openmm-top | numpy-ir | torsion_all | 13.39881665 | 13.39881665 | 8.85e-12 | 0.000101 | True |
| openmm-top | numpy-ir | nonbonded_all | 18.23016312 | 18.23016322 | 9.96e-08 | 0.000102 | True |
| openmm-top | numpy-ir | total | 35.02778658 | 35.02778665 | 6.56e-08 | 0.000109 | True |
| sander-ours | numpy-ir | bond | 0.0535 | 0.05349022907 | 9.77e-06 | 0.0002 | True |
| sander-ours | numpy-ir | angle | 0.7588 | 0.7588440884 | 4.41e-05 | 0.0002 | True |
| sander-ours | numpy-ir | dihed | 3.2024 | 3.202394036 | 5.96e-06 | 0.0002 | True |
| sander-ours | numpy-ir | vdw | 0.387 | 0.3870031267 | 3.13e-06 | 0.0002 | True |
| sander-ours | numpy-ir | eel | -22.3115 | -22.3122493 | 0.000749 | 0.00132 | True |
| sander-ours | numpy-ir | vdw14 | 5.9247 | 5.924719397 | 1.94e-05 | 0.000201 | True |
| sander-ours | numpy-ir | eel14 | 20.3569 | 20.35764035 | 0.00074 | 0.00122 | True |
| sander-tleap | sander-ours | bond | 0.0535 | 0.0535 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | angle | 0.7588 | 0.7588 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | dihed | 3.2024 | 3.2024 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw | 0.387 | 0.387 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | -22.3115 | -22.3115 | 0 | 0.000202 | True |
| sander-tleap | sander-ours | vdw14 | 5.9247 | 5.9247 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | eel14 | 20.3569 | 20.3569 | 0 | 0.000202 | True |
| gmx | numpy-ir | bond | 0.2238037884 | 0.2238031184 | 6.7e-07 | 0.0001 | True |
| gmx | numpy-ir | angle | 3.175035477 | 3.175003666 | 3.18e-05 | 0.0001 | True |
| gmx | numpy-ir | proper | 13.39880276 | 13.39880152 | 1.24e-06 | 0.000101 | True |
| gmx | numpy-ir | improper | 1.5087842e-05 | 1.513171003e-05 | 4.39e-08 | 0.0001 | True |
| gmx | numpy-ir | lj | 1.619218826 | 1.619221082 | 2.26e-06 | 0.0001 | True |
| gmx | numpy-ir | coulomb | -93.35446167 | -93.35445106 | 1.06e-05 | 0.000109 | True |
| gmx | numpy-ir | lj14 | 24.78902054 | 24.78902596 | 5.42e-06 | 0.000102 | True |
| gmx | numpy-ir | coulomb14 | 85.17636871 | 85.17636724 | 1.47e-06 | 0.000109 | True |
| gmx | numpy-ir | total | 35.02780151 | 35.02778665 | 1.49e-05 | 0.000104 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 2.26e-07 | 0.001 | True |
| openmm-top | openmm-ir | 5.98e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 1  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 104.3147563 | 104.3147563 | 0 | 2.04e-06 | True |
| numpy-ir | reference-itself | angle | 34.21740396 | 34.21740396 | 0 | 1.34e-06 | True |
| numpy-ir | reference-itself | proper | 26.41726215 | 26.41726215 | 0 | 1.26e-06 | True |
| numpy-ir | reference-itself | improper | 1.262131846 | 1.262131846 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | lj | 0.6941367119 | 0.6941367119 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | coulomb | -91.85164462 | -91.85164462 | 0 | 1.92e-06 | True |
| numpy-ir | reference-itself | lj14 | 30.49273032 | 30.49273032 | 0 | 1.3e-06 | True |
| numpy-ir | reference-itself | coulomb14 | 84.0012077 | 84.0012077 | 0 | 1.84e-06 | True |
| numpy-ir | reference-itself | total | 189.5479844 | 189.5479844 | 0 | 2.9e-06 | True |
| openmm-ir | numpy-ir | bond | 104.3147563 | 104.3147563 | 2.84e-14 | 2.04e-06 | True |
| openmm-ir | numpy-ir | angle | 34.21740396 | 34.21740396 | 7.11e-15 | 1.34e-06 | True |
| openmm-ir | numpy-ir | proper | 26.41726215 | 26.41726215 | 0 | 1.26e-06 | True |
| openmm-ir | numpy-ir | improper | 1.262131846 | 1.262131846 | 2.22e-16 | 1.01e-06 | True |
| openmm-ir | numpy-ir | lj | 0.6941367119 | 0.6941367119 | 0 | 1.01e-06 | True |
| openmm-ir | numpy-ir | coulomb | -91.85164462 | -91.85164462 | 0 | 1.92e-06 | True |
| openmm-ir | numpy-ir | lj14 | 30.49273032 | 30.49273032 | 0 | 1.3e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | 84.0012077 | 84.0012077 | 2.84e-14 | 1.84e-06 | True |
| openmm-ir | numpy-ir | total | 189.5479844 | 189.5479844 | 8.53e-14 | 2.04e-06 | True |
| openmm-xml | numpy-ir | bond | 104.3147563 | 104.3147563 | 2.84e-14 | 2.04e-06 | True |
| openmm-xml | numpy-ir | angle | 34.21740396 | 34.21740396 | 7.11e-15 | 1.34e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 27.67939399 | 27.67939399 | 0 | 1.28e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | 23.33643002 | 23.33643012 | 9.29e-08 | 1.23e-06 | True |
| openmm-xml | numpy-ir | total | 189.5479843 | 189.5479844 | 9.29e-08 | 2.04e-06 | True |
| openmm-prmtop | numpy-ir | bond | 104.3147563 | 104.3147563 | 0 | 0.00011 | True |
| openmm-prmtop | numpy-ir | angle | 34.21740396 | 34.21740396 | 2.84e-14 | 0.000103 | True |
| openmm-prmtop | numpy-ir | torsion_all | 27.67939399 | 27.67939399 | 0 | 0.000103 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | 23.33643 | 23.33643012 | 1.19e-07 | 0.000102 | True |
| openmm-prmtop | numpy-ir | total | 189.5479842 | 189.5479844 | 1.19e-07 | 0.00011 | True |
| openmm-top | numpy-ir | bond | 104.3147563 | 104.3147563 | 8.53e-14 | 0.00011 | True |
| openmm-top | numpy-ir | angle | 34.21740401 | 34.21740396 | 5.14e-08 | 0.000103 | True |
| openmm-top | numpy-ir | torsion_all | 27.67939399 | 27.67939399 | 5.29e-09 | 0.000103 | True |
| openmm-top | numpy-ir | nonbonded_all | 23.33643002 | 23.33643012 | 9.56e-08 | 0.000102 | True |
| openmm-top | numpy-ir | total | 189.5479843 | 189.5479844 | 4.95e-08 | 0.00011 | True |
| sander-ours | numpy-ir | bond | 24.9318 | 24.93182512 | 2.51e-05 | 0.000202 | True |
| sander-ours | numpy-ir | angle | 8.1782 | 8.178155822 | 4.42e-05 | 0.000201 | True |
| sander-ours | numpy-ir | dihed | 6.6155 | 6.615533938 | 3.39e-05 | 0.000201 | True |
| sander-ours | numpy-ir | vdw | 0.1659 | 0.1659026558 | 2.66e-06 | 0.0002 | True |
| sander-ours | numpy-ir | eel | -21.9523 | -21.95306994 | 0.00077 | 0.0013 | True |
| sander-ours | numpy-ir | vdw14 | 7.2879 | 7.287937457 | 3.75e-05 | 0.000201 | True |
| sander-ours | numpy-ir | eel14 | 20.0761 | 20.07677048 | 0.00067 | 0.0012 | True |
| sander-tleap | sander-ours | bond | 24.9318 | 24.9318 | 0 | 0.000202 | True |
| sander-tleap | sander-ours | angle | 8.1782 | 8.1782 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | dihed | 6.6155 | 6.6155 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | vdw | 0.1659 | 0.1659 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | -21.9523 | -21.9523 | 0 | 0.000202 | True |
| sander-tleap | sander-ours | vdw14 | 7.2879 | 7.2879 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | eel14 | 20.0761 | 20.0761 | 0 | 0.000202 | True |
| gmx | numpy-ir | bond | 104.3146744 | 104.3147563 | 8.19e-05 | 0.00011 | True |
| gmx | numpy-ir | angle | 34.21741867 | 34.21740396 | 1.47e-05 | 0.000103 | True |
| gmx | numpy-ir | proper | 26.41725922 | 26.41726215 | 2.93e-06 | 0.000103 | True |
| gmx | numpy-ir | improper | 1.262130976 | 1.262131846 | 8.71e-07 | 0.0001 | True |
| gmx | numpy-ir | lj | 0.6941315532 | 0.6941367119 | 5.16e-06 | 0.0001 | True |
| gmx | numpy-ir | coulomb | -91.85163879 | -91.85164462 | 5.83e-06 | 0.000109 | True |
| gmx | numpy-ir | lj14 | 30.49274445 | 30.49273032 | 1.41e-05 | 0.000103 | True |
| gmx | numpy-ir | coulomb14 | 84.00119781 | 84.0012077 | 9.89e-06 | 0.000108 | True |
| gmx | numpy-ir | total | 189.5478973 | 189.5479844 | 8.7e-05 | 0.000119 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 3.25e-07 | 0.001 | True |
| openmm-top | openmm-ir | 5.7e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 2  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 118.2685698 | 118.2685698 | 0 | 2.18e-06 | True |
| numpy-ir | reference-itself | angle | 57.52590708 | 57.52590708 | 0 | 1.58e-06 | True |
| numpy-ir | reference-itself | proper | 28.62501753 | 28.62501753 | 0 | 1.29e-06 | True |
| numpy-ir | reference-itself | improper | 2.056022689 | 2.056022689 | 0 | 1.02e-06 | True |
| numpy-ir | reference-itself | lj | 3.349488378 | 3.349488378 | 0 | 1.03e-06 | True |
| numpy-ir | reference-itself | coulomb | -95.17644843 | -95.17644843 | 0 | 1.95e-06 | True |
| numpy-ir | reference-itself | lj14 | 25.10370149 | 25.10370149 | 0 | 1.25e-06 | True |
| numpy-ir | reference-itself | coulomb14 | 87.54689716 | 87.54689716 | 0 | 1.88e-06 | True |
| numpy-ir | reference-itself | total | 227.2991557 | 227.2991557 | 0 | 3.27e-06 | True |
| openmm-ir | numpy-ir | bond | 118.2685698 | 118.2685698 | 5.68e-14 | 2.18e-06 | True |
| openmm-ir | numpy-ir | angle | 57.52590708 | 57.52590708 | 7.11e-15 | 1.58e-06 | True |
| openmm-ir | numpy-ir | proper | 28.62501753 | 28.62501753 | 0 | 1.29e-06 | True |
| openmm-ir | numpy-ir | improper | 2.056022689 | 2.056022689 | 2.22e-15 | 1.02e-06 | True |
| openmm-ir | numpy-ir | lj | 3.349488378 | 3.349488378 | 2.66e-15 | 1.03e-06 | True |
| openmm-ir | numpy-ir | coulomb | -95.17644843 | -95.17644843 | 2.84e-14 | 1.95e-06 | True |
| openmm-ir | numpy-ir | lj14 | 25.10370149 | 25.10370149 | 3.55e-15 | 1.25e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | 87.54689716 | 87.54689716 | 1.42e-14 | 1.88e-06 | True |
| openmm-ir | numpy-ir | total | 227.2991557 | 227.2991557 | 2.84e-14 | 2.18e-06 | True |
| openmm-xml | numpy-ir | bond | 118.2685698 | 118.2685698 | 5.68e-14 | 2.18e-06 | True |
| openmm-xml | numpy-ir | angle | 57.52590708 | 57.52590708 | 7.11e-15 | 1.58e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 30.68104022 | 30.68104022 | 3.55e-15 | 1.31e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | 20.82363851 | 20.8236386 | 9.03e-08 | 1.21e-06 | True |
| openmm-xml | numpy-ir | total | 227.2991556 | 227.2991557 | 9.03e-08 | 2.18e-06 | True |
| openmm-prmtop | numpy-ir | bond | 118.2685698 | 118.2685698 | 1.42e-14 | 0.000112 | True |
| openmm-prmtop | numpy-ir | angle | 57.52590708 | 57.52590708 | 7.11e-15 | 0.000106 | True |
| openmm-prmtop | numpy-ir | torsion_all | 30.68104022 | 30.68104022 | 3.55e-15 | 0.000103 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | 20.82363849 | 20.8236386 | 1.13e-07 | 0.000102 | True |
| openmm-prmtop | numpy-ir | total | 227.2991556 | 227.2991557 | 1.13e-07 | 0.000112 | True |
| openmm-top | numpy-ir | bond | 118.2685698 | 118.2685698 | 1.85e-13 | 0.000112 | True |
| openmm-top | numpy-ir | angle | 57.52590705 | 57.52590708 | 3.01e-08 | 0.000106 | True |
| openmm-top | numpy-ir | torsion_all | 30.68104021 | 30.68104022 | 5.57e-09 | 0.000103 | True |
| openmm-top | numpy-ir | nonbonded_all | 20.82363851 | 20.8236386 | 9.48e-08 | 0.000102 | True |
| openmm-top | numpy-ir | total | 227.2991556 | 227.2991557 | 1.3e-07 | 0.000112 | True |
| sander-ours | numpy-ir | bond | 28.2669 | 28.26686659 | 3.34e-05 | 0.000203 | True |
| sander-ours | numpy-ir | angle | 13.749 | 13.74902177 | 2.18e-05 | 0.000201 | True |
| sander-ours | numpy-ir | dihed | 7.3329 | 7.332944602 | 4.46e-05 | 0.000201 | True |
| sander-ours | numpy-ir | vdw | 0.8005 | 0.8005469355 | 4.69e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel | -22.7469 | -22.74771712 | 0.000817 | 0.00134 | True |
| sander-ours | numpy-ir | vdw14 | 5.9999 | 5.999928655 | 2.87e-05 | 0.000201 | True |
| sander-ours | numpy-ir | eel14 | 20.9235 | 20.9242106 | 0.000711 | 0.00125 | True |
| sander-tleap | sander-ours | bond | 28.2669 | 28.2669 | 0 | 0.000203 | True |
| sander-tleap | sander-ours | angle | 13.749 | 13.749 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | dihed | 7.3329 | 7.3329 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | vdw | 0.8005 | 0.8005 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | -22.7469 | -22.7469 | 0 | 0.000202 | True |
| sander-tleap | sander-ours | vdw14 | 5.9999 | 5.9999 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | eel14 | 20.9235 | 20.9235 | 0 | 0.000202 | True |
| gmx | numpy-ir | bond | 118.2685623 | 118.2685698 | 7.52e-06 | 0.000112 | True |
| gmx | numpy-ir | angle | 57.52585602 | 57.52590708 | 5.11e-05 | 0.000106 | True |
| gmx | numpy-ir | proper | 28.62501526 | 28.62501753 | 2.27e-06 | 0.000103 | True |
| gmx | numpy-ir | improper | 2.056020021 | 2.056022689 | 2.67e-06 | 0.0001 | True |
| gmx | numpy-ir | lj | 3.349482775 | 3.349488378 | 5.6e-06 | 0.0001 | True |
| gmx | numpy-ir | coulomb | -95.17643738 | -95.17644843 | 1.11e-05 | 0.00011 | True |
| gmx | numpy-ir | lj14 | 25.10369682 | 25.10370149 | 4.67e-06 | 0.000103 | True |
| gmx | numpy-ir | coulomb14 | 87.54689789 | 87.54689716 | 7.24e-07 | 0.000109 | True |
| gmx | numpy-ir | total | 227.2990875 | 227.2991557 | 6.82e-05 | 0.000123 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 1.48e-07 | 0.001 | True |
| openmm-top | openmm-ir | 6.13e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 3  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 169.0501323 | 169.0501323 | 0 | 2.69e-06 | True |
| numpy-ir | reference-itself | angle | 51.0448828 | 51.0448828 | 0 | 1.51e-06 | True |
| numpy-ir | reference-itself | proper | 21.91409867 | 21.91409867 | 0 | 1.22e-06 | True |
| numpy-ir | reference-itself | improper | 1.144592205 | 1.144592205 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | lj | 1.107404929 | 1.107404929 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | coulomb | -92.63573562 | -92.63573562 | 0 | 1.93e-06 | True |
| numpy-ir | reference-itself | lj14 | 38.38404921 | 38.38404921 | 0 | 1.38e-06 | True |
| numpy-ir | reference-itself | coulomb14 | 82.95307417 | 82.95307417 | 0 | 1.83e-06 | True |
| numpy-ir | reference-itself | total | 272.9624987 | 272.9624987 | 0 | 3.73e-06 | True |
| openmm-ir | numpy-ir | bond | 169.0501323 | 169.0501323 | 5.68e-14 | 2.69e-06 | True |
| openmm-ir | numpy-ir | angle | 51.0448828 | 51.0448828 | 1.42e-14 | 1.51e-06 | True |
| openmm-ir | numpy-ir | proper | 21.91409867 | 21.91409867 | 3.2e-14 | 1.22e-06 | True |
| openmm-ir | numpy-ir | improper | 1.144592205 | 1.144592205 | 4.44e-15 | 1.01e-06 | True |
| openmm-ir | numpy-ir | lj | 1.107404929 | 1.107404929 | 3.55e-15 | 1.01e-06 | True |
| openmm-ir | numpy-ir | coulomb | -92.63573562 | -92.63573562 | 0 | 1.93e-06 | True |
| openmm-ir | numpy-ir | lj14 | 38.38404921 | 38.38404921 | 0 | 1.38e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | 82.95307417 | 82.95307417 | 4.26e-14 | 1.83e-06 | True |
| openmm-ir | numpy-ir | total | 272.9624987 | 272.9624987 | 1.14e-13 | 2.69e-06 | True |
| openmm-xml | numpy-ir | bond | 169.0501323 | 169.0501323 | 5.68e-14 | 2.69e-06 | True |
| openmm-xml | numpy-ir | angle | 51.0448828 | 51.0448828 | 1.42e-14 | 1.51e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 23.05869088 | 23.05869088 | 3.55e-14 | 1.23e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | 29.80879258 | 29.80879269 | 1.15e-07 | 1.3e-06 | True |
| openmm-xml | numpy-ir | total | 272.9624986 | 272.9624987 | 1.15e-07 | 2.69e-06 | True |
| openmm-prmtop | numpy-ir | bond | 169.0501323 | 169.0501323 | 2.84e-14 | 0.000117 | True |
| openmm-prmtop | numpy-ir | angle | 51.0448828 | 51.0448828 | 4.26e-14 | 0.000105 | True |
| openmm-prmtop | numpy-ir | torsion_all | 23.05869088 | 23.05869088 | 3.2e-14 | 0.000102 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | 29.80879254 | 29.80879269 | 1.5e-07 | 0.000103 | True |
| openmm-prmtop | numpy-ir | total | 272.9624985 | 272.9624987 | 1.5e-07 | 0.000117 | True |
| openmm-top | numpy-ir | bond | 169.0501323 | 169.0501323 | 8.53e-14 | 0.000117 | True |
| openmm-top | numpy-ir | angle | 51.04488272 | 51.0448828 | 8.48e-08 | 0.000105 | True |
| openmm-top | numpy-ir | torsion_all | 23.05869088 | 23.05869088 | 1.72e-09 | 0.000102 | True |
| openmm-top | numpy-ir | nonbonded_all | 29.80879257 | 29.80879269 | 1.19e-07 | 0.000103 | True |
| openmm-top | numpy-ir | total | 272.9624985 | 272.9624987 | 2.05e-07 | 0.000117 | True |
| sander-ours | numpy-ir | bond | 40.4039 | 40.40395131 | 5.13e-05 | 0.000204 | True |
| sander-ours | numpy-ir | angle | 12.2 | 12.20001979 | 1.98e-05 | 0.000201 | True |
| sander-ours | numpy-ir | dihed | 5.5112 | 5.511159387 | 4.06e-05 | 0.000201 | True |
| sander-ours | numpy-ir | vdw | 0.2647 | 0.2646761302 | 2.39e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel | -22.1397 | -22.14047218 | 0.000772 | 0.00131 | True |
| sander-ours | numpy-ir | vdw14 | 9.174 | 9.174007937 | 7.94e-06 | 0.000201 | True |
| sander-ours | numpy-ir | eel14 | 19.8256 | 19.82626056 | 0.000661 | 0.00119 | True |
| sander-tleap | sander-ours | bond | 40.4039 | 40.4039 | 0 | 0.000204 | True |
| sander-tleap | sander-ours | angle | 12.2 | 12.2 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | dihed | 5.5112 | 5.5112 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | vdw | 0.2647 | 0.2647 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | -22.1397 | -22.1397 | 0 | 0.000202 | True |
| sander-tleap | sander-ours | vdw14 | 9.174 | 9.174 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | eel14 | 19.8256 | 19.8256 | 0 | 0.000202 | True |
| gmx | numpy-ir | bond | 169.0500946 | 169.0501323 | 3.77e-05 | 0.000117 | True |
| gmx | numpy-ir | angle | 51.04486847 | 51.0448828 | 1.43e-05 | 0.000105 | True |
| gmx | numpy-ir | proper | 21.91410446 | 21.91409867 | 5.79e-06 | 0.000102 | True |
| gmx | numpy-ir | improper | 1.144592762 | 1.144592205 | 5.57e-07 | 0.0001 | True |
| gmx | numpy-ir | lj | 1.107403755 | 1.107404929 | 1.17e-06 | 0.0001 | True |
| gmx | numpy-ir | coulomb | -92.63574219 | -92.63573562 | 6.57e-06 | 0.000109 | True |
| gmx | numpy-ir | lj14 | 38.3840065 | 38.38404921 | 4.27e-05 | 0.000104 | True |
| gmx | numpy-ir | coulomb14 | 82.95305634 | 82.95307417 | 1.78e-05 | 0.000108 | True |
| gmx | numpy-ir | total | 272.9624023 | 272.9624987 | 9.63e-05 | 0.000127 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 5.84e-07 | 0.001 | True |
| openmm-top | openmm-ir | 5.89e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

## dmso
### config 0  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 2.191600575 | 2.191600575 | 0 | 1.02e-06 | True |
| numpy-ir | reference-itself | angle | 1.538698594 | 1.538698594 | 0 | 1.02e-06 | True |
| numpy-ir | reference-itself | proper | 0.3864881364 | 0.3864881364 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | improper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | -0.3403176606 | -0.3403176606 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb | 25.807557 | 25.807557 | 0 | 1.26e-06 | True |
| numpy-ir | reference-itself | lj14 | -0.3718384765 | -0.3718384765 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -118.7572208 | -118.7572208 | 0 | 2.19e-06 | True |
| numpy-ir | reference-itself | total | -89.5450326 | -89.5450326 | 0 | 1.9e-06 | True |
| openmm-ir | numpy-ir | bond | 2.191600575 | 2.191600575 | 1.78e-15 | 1.02e-06 | True |
| openmm-ir | numpy-ir | angle | 1.538698594 | 1.538698594 | 3.11e-15 | 1.02e-06 | True |
| openmm-ir | numpy-ir | proper | 0.3864881364 | 0.3864881364 | 5.55e-16 | 1e-06 | True |
| openmm-ir | numpy-ir | improper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | -0.3403176606 | -0.3403176606 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb | 25.807557 | 25.807557 | 0 | 1.26e-06 | True |
| openmm-ir | numpy-ir | lj14 | -0.3718384765 | -0.3718384765 | 2.78e-16 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -118.7572208 | -118.7572208 | 0 | 2.19e-06 | True |
| openmm-ir | numpy-ir | total | -89.5450326 | -89.5450326 | 0 | 2.19e-06 | True |
| openmm-xml | numpy-ir | bond | 2.191600575 | 2.191600575 | 1.78e-15 | 1.02e-06 | True |
| openmm-xml | numpy-ir | angle | 1.538698594 | 1.538698594 | 3.11e-15 | 1.02e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 0.3864881364 | 0.3864881364 | 5.55e-16 | 1e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | -93.66182101 | -93.66181991 | 1.1e-06 | 1.94e-06 | True |
| openmm-xml | numpy-ir | total | -89.5450337 | -89.5450326 | 1.1e-06 | 2.19e-06 | True |
| openmm-prmtop | numpy-ir | bond | 2.191600575 | 2.191600575 | 1.78e-15 | 0.0001 | True |
| openmm-prmtop | numpy-ir | angle | 1.538698594 | 1.538698594 | 3.55e-15 | 0.0001 | True |
| openmm-prmtop | numpy-ir | torsion_all | 0.3864881364 | 0.3864881364 | 5.55e-16 | 0.0001 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | -93.66182101 | -93.66181991 | 1.1e-06 | 0.000109 | True |
| openmm-prmtop | numpy-ir | total | -89.5450337 | -89.5450326 | 1.1e-06 | 0.000112 | True |
| openmm-top | numpy-ir | bond | 2.191600575 | 2.191600575 | 4.26e-14 | 0.0001 | True |
| openmm-top | numpy-ir | angle | 1.538698537 | 1.538698594 | 5.77e-08 | 0.0001 | True |
| openmm-top | numpy-ir | torsion_all | 0.3864881364 | 0.3864881364 | 5e-16 | 0.0001 | True |
| openmm-top | numpy-ir | nonbonded_all | -93.661821 | -93.66181991 | 1.1e-06 | 0.000109 | True |
| openmm-top | numpy-ir | total | -89.54503376 | -89.5450326 | 1.16e-06 | 0.000112 | True |
| sander-ours | numpy-ir | bond | 0.5238 | 0.5238051087 | 5.11e-06 | 0.0002 | True |
| sander-ours | numpy-ir | angle | 0.3678 | 0.3677577902 | 4.22e-05 | 0.0002 | True |
| sander-ours | numpy-ir | dihed | 0.0924 | 0.09237288156 | 2.71e-05 | 0.0002 | True |
| sander-ours | numpy-ir | vdw | -0.0813 | -0.081337873 | 3.79e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 6.1679 | 6.168154159 | 0.000254 | 0.000508 | True |
| sander-ours | numpy-ir | vdw14 | -0.0889 | -0.08887152879 | 2.85e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -28.3827 | -28.38365697 | 0.000957 | 0.00162 | True |
| sander-tleap | sander-ours | bond | 0.5238 | 0.5238 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | angle | 0.3678 | 0.3678 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | dihed | 0.0924 | 0.0924 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw | -0.0813 | -0.0813 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 6.1679 | 6.1679 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | vdw14 | -0.0889 | -0.0889 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -28.3827 | -28.3827 | 0 | 0.000203 | True |
| gmx | numpy-ir | bond | 2.191618204 | 2.191600575 | 1.76e-05 | 0.0001 | True |
| gmx | numpy-ir | angle | 1.538696408 | 1.538698594 | 2.19e-06 | 0.0001 | True |
| gmx | numpy-ir | proper | 0.3864878118 | 0.3864881364 | 3.25e-07 | 0.0001 | True |
| gmx | numpy-ir | improper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj | -0.3403177261 | -0.3403176606 | 6.55e-08 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 25.80755997 | 25.807557 | 2.97e-06 | 0.000103 | True |
| gmx | numpy-ir | lj14 | -0.37184003 | -0.3718384765 | 1.55e-06 | 0.0001 | True |
| gmx | numpy-ir | coulomb14 | -118.7572327 | -118.7572208 | 1.19e-05 | 0.000112 | True |
| gmx | numpy-ir | total | -89.54502869 | -89.5450326 | 3.91e-06 | 0.000109 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 2.06e-08 | 0.001 | True |
| openmm-top | openmm-ir | 7.35e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 1  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 51.77425323 | 51.77425323 | 0 | 1.52e-06 | True |
| numpy-ir | reference-itself | angle | 14.91846254 | 14.91846254 | 0 | 1.15e-06 | True |
| numpy-ir | reference-itself | proper | 0.6329380944 | 0.6329380944 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | improper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | -0.3495564117 | -0.3495564117 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb | 25.54736947 | 25.54736947 | 0 | 1.26e-06 | True |
| numpy-ir | reference-itself | lj14 | -0.2019880524 | -0.2019880524 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -118.1707214 | -118.1707214 | 0 | 2.18e-06 | True |
| numpy-ir | reference-itself | total | -25.84924254 | -25.84924254 | 0 | 1.26e-06 | True |
| openmm-ir | numpy-ir | bond | 51.77425323 | 51.77425323 | 0 | 1.52e-06 | True |
| openmm-ir | numpy-ir | angle | 14.91846254 | 14.91846254 | 8.88e-15 | 1.15e-06 | True |
| openmm-ir | numpy-ir | proper | 0.6329380944 | 0.6329380944 | 2.22e-16 | 1.01e-06 | True |
| openmm-ir | numpy-ir | improper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | -0.3495564117 | -0.3495564117 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb | 25.54736947 | 25.54736947 | 0 | 1.26e-06 | True |
| openmm-ir | numpy-ir | lj14 | -0.2019880524 | -0.2019880524 | 1.17e-15 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -118.1707214 | -118.1707214 | 1.42e-14 | 2.18e-06 | True |
| openmm-ir | numpy-ir | total | -25.84924254 | -25.84924254 | 2.49e-14 | 2.18e-06 | True |
| openmm-xml | numpy-ir | bond | 51.77425323 | 51.77425323 | 0 | 1.52e-06 | True |
| openmm-xml | numpy-ir | angle | 14.91846254 | 14.91846254 | 8.88e-15 | 1.15e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 0.6329380944 | 0.6329380944 | 2.22e-16 | 1.01e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | -93.1748975 | -93.1748964 | 1.1e-06 | 1.93e-06 | True |
| openmm-xml | numpy-ir | total | -25.84924364 | -25.84924254 | 1.1e-06 | 2.18e-06 | True |
| openmm-prmtop | numpy-ir | bond | 51.77425323 | 51.77425323 | 0 | 0.000105 | True |
| openmm-prmtop | numpy-ir | angle | 14.91846254 | 14.91846254 | 7.11e-15 | 0.000101 | True |
| openmm-prmtop | numpy-ir | torsion_all | 0.6329380944 | 0.6329380944 | 2.22e-16 | 0.0001 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | -93.1748975 | -93.1748964 | 1.1e-06 | 0.000109 | True |
| openmm-prmtop | numpy-ir | total | -25.84924364 | -25.84924254 | 1.1e-06 | 0.000112 | True |
| openmm-top | numpy-ir | bond | 51.77425323 | 51.77425323 | 6.39e-14 | 0.000105 | True |
| openmm-top | numpy-ir | angle | 14.91846244 | 14.91846254 | 9.23e-08 | 0.000101 | True |
| openmm-top | numpy-ir | torsion_all | 0.6329380944 | 0.6329380944 | 2.22e-16 | 0.0001 | True |
| openmm-top | numpy-ir | nonbonded_all | -93.1748975 | -93.1748964 | 1.09e-06 | 0.000109 | True |
| openmm-top | numpy-ir | total | -25.84924373 | -25.84924254 | 1.19e-06 | 0.000112 | True |
| sander-ours | numpy-ir | bond | 12.3743 | 12.37434351 | 4.35e-05 | 0.000201 | True |
| sander-ours | numpy-ir | angle | 3.5656 | 3.565598121 | 1.88e-06 | 0.0002 | True |
| sander-ours | numpy-ir | dihed | 0.1513 | 0.1512758352 | 2.42e-05 | 0.0002 | True |
| sander-ours | numpy-ir | vdw | -0.0835 | -0.08354598751 | 4.6e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 6.1058 | 6.105967846 | 0.000168 | 0.000505 | True |
| sander-ours | numpy-ir | vdw14 | -0.0483 | -0.04827630316 | 2.37e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -28.2425 | -28.24348026 | 0.00098 | 0.00161 | True |
| sander-tleap | sander-ours | bond | 12.3743 | 12.3743 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | angle | 3.5656 | 3.5656 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | dihed | 0.1513 | 0.1513 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw | -0.0835 | -0.0835 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 6.1058 | 6.1058 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | vdw14 | -0.0483 | -0.0483 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -28.2425 | -28.2425 | 0 | 0.000203 | True |
| gmx | numpy-ir | bond | 51.77426529 | 51.77425323 | 1.21e-05 | 0.000105 | True |
| gmx | numpy-ir | angle | 14.91847706 | 14.91846254 | 1.45e-05 | 0.000101 | True |
| gmx | numpy-ir | proper | 0.6329388022 | 0.6329380944 | 7.08e-07 | 0.0001 | True |
| gmx | numpy-ir | improper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj | -0.3495565057 | -0.3495564117 | 9.4e-08 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 25.54737091 | 25.54736947 | 1.44e-06 | 0.000103 | True |
| gmx | numpy-ir | lj14 | -0.2019871771 | -0.2019880524 | 8.75e-07 | 0.0001 | True |
| gmx | numpy-ir | coulomb14 | -118.1707382 | -118.1707214 | 1.68e-05 | 0.000112 | True |
| gmx | numpy-ir | total | -25.84922791 | -25.84924254 | 1.46e-05 | 0.000103 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 2e-08 | 0.001 | True |
| openmm-top | openmm-ir | 6.92e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 2  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 50.85031062 | 50.85031062 | 0 | 1.51e-06 | True |
| numpy-ir | reference-itself | angle | 20.55717215 | 20.55717215 | 0 | 1.21e-06 | True |
| numpy-ir | reference-itself | proper | 0.2098537599 | 0.2098537599 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | improper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | -0.3519711231 | -0.3519711231 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb | 25.12416132 | 25.12416132 | 0 | 1.25e-06 | True |
| numpy-ir | reference-itself | lj14 | -0.5122862072 | -0.5122862072 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -118.5857362 | -118.5857362 | 0 | 2.19e-06 | True |
| numpy-ir | reference-itself | total | -22.70849564 | -22.70849564 | 0 | 1.23e-06 | True |
| openmm-ir | numpy-ir | bond | 50.85031062 | 50.85031062 | 4.26e-14 | 1.51e-06 | True |
| openmm-ir | numpy-ir | angle | 20.55717215 | 20.55717215 | 3.2e-14 | 1.21e-06 | True |
| openmm-ir | numpy-ir | proper | 0.2098537599 | 0.2098537599 | 1.94e-16 | 1e-06 | True |
| openmm-ir | numpy-ir | improper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | -0.3519711231 | -0.3519711231 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb | 25.12416132 | 25.12416132 | 0 | 1.25e-06 | True |
| openmm-ir | numpy-ir | lj14 | -0.5122862072 | -0.5122862072 | 1.11e-16 | 1.01e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -118.5857362 | -118.5857362 | 0 | 2.19e-06 | True |
| openmm-ir | numpy-ir | total | -22.70849564 | -22.70849564 | 1.07e-14 | 2.19e-06 | True |
| openmm-xml | numpy-ir | bond | 50.85031062 | 50.85031062 | 4.26e-14 | 1.51e-06 | True |
| openmm-xml | numpy-ir | angle | 20.55717215 | 20.55717215 | 3.2e-14 | 1.21e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 0.2098537599 | 0.2098537599 | 1.94e-16 | 1e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | -94.32583327 | -94.32583217 | 1.11e-06 | 1.94e-06 | True |
| openmm-xml | numpy-ir | total | -22.70849674 | -22.70849564 | 1.11e-06 | 2.19e-06 | True |
| openmm-prmtop | numpy-ir | bond | 50.85031062 | 50.85031062 | 4.26e-14 | 0.000105 | True |
| openmm-prmtop | numpy-ir | angle | 20.55717215 | 20.55717215 | 3.55e-14 | 0.000102 | True |
| openmm-prmtop | numpy-ir | torsion_all | 0.2098537599 | 0.2098537599 | 1.94e-16 | 0.0001 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | -94.32583327 | -94.32583217 | 1.11e-06 | 0.000109 | True |
| openmm-prmtop | numpy-ir | total | -22.70849674 | -22.70849564 | 1.11e-06 | 0.000112 | True |
| openmm-top | numpy-ir | bond | 50.85031062 | 50.85031062 | 4.97e-14 | 0.000105 | True |
| openmm-top | numpy-ir | angle | 20.55717202 | 20.55717215 | 1.23e-07 | 0.000102 | True |
| openmm-top | numpy-ir | torsion_all | 0.2098537599 | 0.2098537599 | 1.94e-16 | 0.0001 | True |
| openmm-top | numpy-ir | nonbonded_all | -94.32583327 | -94.32583217 | 1.1e-06 | 0.000109 | True |
| openmm-top | numpy-ir | total | -22.70849686 | -22.70849564 | 1.23e-06 | 0.000112 | True |
| sander-ours | numpy-ir | bond | 12.1535 | 12.15351592 | 1.59e-05 | 0.000201 | True |
| sander-ours | numpy-ir | angle | 4.9133 | 4.913282062 | 1.79e-05 | 0.0002 | True |
| sander-ours | numpy-ir | dihed | 0.0502 | 0.05015625235 | 4.37e-05 | 0.0002 | True |
| sander-ours | numpy-ir | vdw | -0.0841 | -0.08412311737 | 2.31e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 6.0046 | 6.00481867 | 0.000219 | 0.0005 | True |
| sander-ours | numpy-ir | vdw14 | -0.1224 | -0.1224393421 | 3.93e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -28.3417 | -28.34267116 | 0.000971 | 0.00162 | True |
| sander-tleap | sander-ours | bond | 12.1535 | 12.1535 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | angle | 4.9133 | 4.9133 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | dihed | 0.0502 | 0.0502 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw | -0.0841 | -0.0841 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 6.0046 | 6.0046 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | vdw14 | -0.1224 | -0.1224 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -28.3417 | -28.3417 | 0 | 0.000203 | True |
| gmx | numpy-ir | bond | 50.85039902 | 50.85031062 | 8.84e-05 | 0.000105 | True |
| gmx | numpy-ir | angle | 20.55717468 | 20.55717215 | 2.54e-06 | 0.000102 | True |
| gmx | numpy-ir | proper | 0.2098535746 | 0.2098537599 | 1.85e-07 | 0.0001 | True |
| gmx | numpy-ir | improper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj | -0.35197106 | -0.3519711231 | 6.3e-08 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 25.12416458 | 25.12416132 | 3.27e-06 | 0.000103 | True |
| gmx | numpy-ir | lj14 | -0.512286067 | -0.5122862072 | 1.4e-07 | 0.0001 | True |
| gmx | numpy-ir | coulomb14 | -118.5857468 | -118.5857362 | 1.06e-05 | 0.000112 | True |
| gmx | numpy-ir | total | -22.70840454 | -22.70849564 | 9.11e-05 | 0.000102 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 2.23e-08 | 0.001 | True |
| openmm-top | openmm-ir | 7.12e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 3  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 61.17628609 | 61.17628609 | 0 | 1.61e-06 | True |
| numpy-ir | reference-itself | angle | 56.94028593 | 56.94028593 | 0 | 1.57e-06 | True |
| numpy-ir | reference-itself | proper | 0.4833224749 | 0.4833224749 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | improper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | -0.3647664252 | -0.3647664252 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb | 25.27517364 | 25.27517364 | 0 | 1.25e-06 | True |
| numpy-ir | reference-itself | lj14 | -0.6057023539 | -0.6057023539 | 0 | 1.01e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -116.4305032 | -116.4305032 | 0 | 2.16e-06 | True |
| numpy-ir | reference-itself | total | 26.47409621 | 26.47409621 | 0 | 1.26e-06 | True |
| openmm-ir | numpy-ir | bond | 61.17628609 | 61.17628609 | 5.68e-14 | 1.61e-06 | True |
| openmm-ir | numpy-ir | angle | 56.94028593 | 56.94028593 | 0 | 1.57e-06 | True |
| openmm-ir | numpy-ir | proper | 0.4833224749 | 0.4833224749 | 3.33e-16 | 1e-06 | True |
| openmm-ir | numpy-ir | improper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | -0.3647664252 | -0.3647664252 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb | 25.27517364 | 25.27517364 | 0 | 1.25e-06 | True |
| openmm-ir | numpy-ir | lj14 | -0.6057023539 | -0.6057023539 | 5.55e-16 | 1.01e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -116.4305032 | -116.4305032 | 1.42e-14 | 2.16e-06 | True |
| openmm-ir | numpy-ir | total | 26.47409621 | 26.47409621 | 4.26e-14 | 2.16e-06 | True |
| openmm-xml | numpy-ir | bond | 61.17628609 | 61.17628609 | 5.68e-14 | 1.61e-06 | True |
| openmm-xml | numpy-ir | angle | 56.94028593 | 56.94028593 | 0 | 1.57e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 0.4833224749 | 0.4833224749 | 3.33e-16 | 1e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | -92.12579937 | -92.12579829 | 1.08e-06 | 1.92e-06 | True |
| openmm-xml | numpy-ir | total | 26.47409513 | 26.47409621 | 1.08e-06 | 2.16e-06 | True |
| openmm-prmtop | numpy-ir | bond | 61.17628609 | 61.17628609 | 5.68e-14 | 0.000106 | True |
| openmm-prmtop | numpy-ir | angle | 56.94028593 | 56.94028593 | 7.11e-15 | 0.000106 | True |
| openmm-prmtop | numpy-ir | torsion_all | 0.4833224749 | 0.4833224749 | 3.33e-16 | 0.0001 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | -92.12579937 | -92.12579829 | 1.08e-06 | 0.000109 | True |
| openmm-prmtop | numpy-ir | total | 26.47409513 | 26.47409621 | 1.08e-06 | 0.000112 | True |
| openmm-top | numpy-ir | bond | 61.17628609 | 61.17628609 | 7.11e-14 | 0.000106 | True |
| openmm-top | numpy-ir | angle | 56.94028576 | 56.94028593 | 1.76e-07 | 0.000106 | True |
| openmm-top | numpy-ir | torsion_all | 0.4833224749 | 0.4833224749 | 3.33e-16 | 0.0001 | True |
| openmm-top | numpy-ir | nonbonded_all | -92.12579937 | -92.12579829 | 1.08e-06 | 0.000109 | True |
| openmm-top | numpy-ir | total | 26.47409496 | 26.47409621 | 1.25e-06 | 0.000112 | True |
| sander-ours | numpy-ir | bond | 14.6215 | 14.62148329 | 1.67e-05 | 0.000201 | True |
| sander-ours | numpy-ir | angle | 13.6091 | 13.60905496 | 4.5e-05 | 0.000201 | True |
| sander-ours | numpy-ir | dihed | 0.1155 | 0.1155168439 | 1.68e-05 | 0.0002 | True |
| sander-ours | numpy-ir | vdw | -0.0872 | -0.08718126798 | 1.87e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 6.0407 | 6.040911483 | 0.000211 | 0.000502 | True |
| sander-ours | numpy-ir | vdw14 | -0.1448 | -0.144766337 | 3.37e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -27.8266 | -27.82755812 | 0.000958 | 0.00159 | True |
| sander-tleap | sander-ours | bond | 14.6215 | 14.6215 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | angle | 13.6091 | 13.6091 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | dihed | 0.1155 | 0.1155 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw | -0.0872 | -0.0872 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 6.0407 | 6.0407 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | vdw14 | -0.1448 | -0.1448 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -27.8266 | -27.8266 | 0 | 0.000203 | True |
| gmx | numpy-ir | bond | 61.17627335 | 61.17628609 | 1.27e-05 | 0.000106 | True |
| gmx | numpy-ir | angle | 56.94029617 | 56.94028593 | 1.02e-05 | 0.000106 | True |
| gmx | numpy-ir | proper | 0.4833222032 | 0.4833224749 | 2.72e-07 | 0.0001 | True |
| gmx | numpy-ir | improper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj | -0.3647663593 | -0.3647664252 | 6.59e-08 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 25.27517509 | 25.27517364 | 1.45e-06 | 0.000103 | True |
| gmx | numpy-ir | lj14 | -0.6057023406 | -0.6057023539 | 1.33e-08 | 0.0001 | True |
| gmx | numpy-ir | coulomb14 | -116.4305191 | -116.4305032 | 1.6e-05 | 0.000112 | True |
| gmx | numpy-ir | total | 26.47407532 | 26.47409621 | 2.09e-05 | 0.000103 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 1.6e-08 | 0.001 | True |
| openmm-top | openmm-ir | 7.91e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

## acetonitrile
### config 0  passed=False

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 0.04587591683 | 0.04587591683 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | angle | 0.4745081711 | 0.4745081711 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | proper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | improper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj14 | -0.2825473015 | -0.2825473015 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -29.01939985 | -29.01939985 | 0 | 1.29e-06 | True |
| numpy-ir | reference-itself | total | -28.78156307 | -28.78156307 | 0 | 1.29e-06 | True |
| openmm-ir | numpy-ir | bond | 0.04587591683 | 0.04587591683 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | angle | 0.4745081711 | 0.4745081711 | 8.52e-12 | 1e-06 | True |
| openmm-ir | numpy-ir | proper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | improper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj14 | -0.2825473015 | -0.2825473015 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -29.01939985 | -29.01939985 | 7.11e-15 | 1.29e-06 | True |
| openmm-ir | numpy-ir | total | -28.78156307 | -28.78156307 | 8.53e-12 | 1.29e-06 | True |
| openmm-xml | numpy-ir | bond | 0.04587591683 | 0.04587591683 | 0 | 1e-06 | True |
| openmm-xml | numpy-ir | angle | 0.4745081711 | 0.4745081711 | 8.52e-12 | 1e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 0 | 0 | 0 | 1e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | -29.3019475 | -29.30194716 | 3.43e-07 | 1.29e-06 | True |
| openmm-xml | numpy-ir | total | -28.78156341 | -28.78156307 | 3.43e-07 | 1.29e-06 | True |
| openmm-prmtop | numpy-ir | bond | 0.04587591683 | 0.04587591683 | 0 | 0.0001 | True |
| openmm-prmtop | numpy-ir | angle | 0.4745081711 | 0.4745081711 | 8.52e-12 | 0.0001 | True |
| openmm-prmtop | numpy-ir | torsion_all | 0 | 0 | 0 | 0.0001 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | -29.3019475 | -29.30194716 | 3.41e-07 | 0.000103 | True |
| openmm-prmtop | numpy-ir | total | -28.78156341 | -28.78156307 | 3.41e-07 | 0.000103 | True |
| openmm-top | numpy-ir | bond | 0.04587591683 | 0.04587591683 | 2.65e-15 | 0.0001 | True |
| openmm-top | numpy-ir | angle | 0.4745081678 | 0.4745081711 | 3.3e-09 | 0.0001 | True |
| openmm-top | numpy-ir | torsion_all | 0 | 0 | 0 | 0.0001 | True |
| openmm-top | numpy-ir | nonbonded_all | -29.3019475 | -29.30194716 | 3.42e-07 | 0.000103 | True |
| openmm-top | numpy-ir | total | -28.78156341 | -28.78156307 | 3.46e-07 | 0.000103 | True |
| sander-ours | numpy-ir | bond | 0.011 | 0.01096460727 | 3.54e-05 | 0.0002 | True |
| sander-ours | numpy-ir | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | vdw14 | -0.0675 | -0.06753042578 | 3.04e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -6.9356 | -6.935803024 | 0.000203 | 0.000547 | True |
| sander-tleap | sander-ours | bond | 0.011 | 0.011 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | angle | 0.1088 | 0.1088 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | dihed | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | -0.0675 | -0.0675 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -6.9356 | -6.9356 | 0 | 0.000201 | True |
| gmx | numpy-ir | bond | 0.04587664083 | 0.04587591683 | 7.24e-07 | 0.0001 | True |
| gmx | numpy-ir | angle | 0.4733843803 | 0.4745081711 | 0.00112 | 0.0001 | False |
| gmx | numpy-ir | proper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | improper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj14 | -0.2825472355 | -0.2825473015 | 6.6e-08 | 0.0001 | True |
| gmx | numpy-ir | coulomb14 | -29.01939774 | -29.01939985 | 2.12e-06 | 0.000103 | True |
| gmx | numpy-ir | total | -28.78268433 | -28.78156307 | 0.00112 | 0.000103 | False |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 5.29e-08 | 0.001 | True |
| openmm-top | openmm-ir | 4.89e-06 | 0.001 | True |
SKIPPED sander-ours:angle/dihed: near-linear angle(s) [(0, 1, 2)] (> 175.0 deg); engine-specific sander-vs-numpy difference (decision 0004.2)
GROMACS precision=mixed; strict file tolerance

### config 1  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 5.270075831 | 5.270075831 | 0 | 1.05e-06 | True |
| numpy-ir | reference-itself | angle | 3.600624038 | 3.600624038 | 0 | 1.04e-06 | True |
| numpy-ir | reference-itself | proper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | improper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj14 | -0.2703468752 | -0.2703468752 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -29.20949354 | -29.20949354 | 0 | 1.29e-06 | True |
| numpy-ir | reference-itself | total | -20.60914055 | -20.60914055 | 0 | 1.21e-06 | True |
| openmm-ir | numpy-ir | bond | 5.270075831 | 5.270075831 | 0 | 1.05e-06 | True |
| openmm-ir | numpy-ir | angle | 3.600624038 | 3.600624038 | 3.24e-14 | 1.04e-06 | True |
| openmm-ir | numpy-ir | proper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | improper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj14 | -0.2703468752 | -0.2703468752 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -29.20949354 | -29.20949354 | 7.11e-15 | 1.29e-06 | True |
| openmm-ir | numpy-ir | total | -20.60914055 | -20.60914055 | 2.84e-14 | 1.29e-06 | True |
| openmm-xml | numpy-ir | bond | 5.270075831 | 5.270075831 | 0 | 1.05e-06 | True |
| openmm-xml | numpy-ir | angle | 3.600624038 | 3.600624038 | 3.24e-14 | 1.04e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 0 | 0 | 0 | 1e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | -29.47984077 | -29.47984042 | 3.46e-07 | 1.29e-06 | True |
| openmm-xml | numpy-ir | total | -20.6091409 | -20.60914055 | 3.46e-07 | 1.29e-06 | True |
| openmm-prmtop | numpy-ir | bond | 5.270075831 | 5.270075831 | 8.88e-16 | 0.000101 | True |
| openmm-prmtop | numpy-ir | angle | 3.600624038 | 3.600624038 | 3.2e-14 | 0.0001 | True |
| openmm-prmtop | numpy-ir | torsion_all | 0 | 0 | 0 | 0.0001 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | -29.47984076 | -29.47984042 | 3.43e-07 | 0.000103 | True |
| openmm-prmtop | numpy-ir | total | -20.60914089 | -20.60914055 | 3.43e-07 | 0.000103 | True |
| openmm-top | numpy-ir | bond | 5.270075831 | 5.270075831 | 9.77e-15 | 0.000101 | True |
| openmm-top | numpy-ir | angle | 3.600624065 | 3.600624038 | 2.68e-08 | 0.0001 | True |
| openmm-top | numpy-ir | torsion_all | 0 | 0 | 0 | 0.0001 | True |
| openmm-top | numpy-ir | nonbonded_all | -29.47984076 | -29.47984042 | 3.45e-07 | 0.000103 | True |
| openmm-top | numpy-ir | total | -20.60914087 | -20.60914055 | 3.18e-07 | 0.000103 | True |
| sander-ours | numpy-ir | bond | 1.2596 | 1.259578354 | 2.16e-05 | 0.0002 | True |
| sander-ours | numpy-ir | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | vdw14 | -0.0646 | -0.06461445391 | 1.45e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -6.981 | -6.981236507 | 0.000237 | 0.000549 | True |
| sander-tleap | sander-ours | bond | 1.2596 | 1.2596 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | angle | 0.8606 | 0.8606 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | dihed | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | -0.0646 | -0.0646 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -6.981 | -6.981 | 0 | 0.000201 | True |
| gmx | numpy-ir | bond | 5.270064354 | 5.270075831 | 1.15e-05 | 0.000101 | True |
| gmx | numpy-ir | angle | 3.600635767 | 3.600624038 | 1.17e-05 | 0.0001 | True |
| gmx | numpy-ir | proper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | improper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj14 | -0.27034688 | -0.2703468752 | 4.81e-09 | 0.0001 | True |
| gmx | numpy-ir | coulomb14 | -29.20949173 | -29.20949354 | 1.82e-06 | 0.000103 | True |
| gmx | numpy-ir | total | -20.60913849 | -20.60914055 | 2.06e-06 | 0.000102 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 5.57e-08 | 0.001 | True |
| openmm-top | openmm-ir | 4.7e-06 | 0.001 | True |
SKIPPED sander-ours:angle/dihed: near-linear angle(s) [(0, 1, 2)] (> 175.0 deg); engine-specific sander-vs-numpy difference (decision 0004.2)
GROMACS precision=mixed; strict file tolerance

### config 2  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 48.05368473 | 48.05368473 | 0 | 1.48e-06 | True |
| numpy-ir | reference-itself | angle | 34.30677353 | 34.30677353 | 0 | 1.34e-06 | True |
| numpy-ir | reference-itself | proper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | improper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj14 | -0.2752731221 | -0.2752731221 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -28.95685035 | -28.95685035 | 0 | 1.29e-06 | True |
| numpy-ir | reference-itself | total | 53.12833479 | 53.12833479 | 0 | 1.53e-06 | True |
| openmm-ir | numpy-ir | bond | 48.05368473 | 48.05368473 | 6.39e-14 | 1.48e-06 | True |
| openmm-ir | numpy-ir | angle | 34.30677353 | 34.30677353 | 4.97e-14 | 1.34e-06 | True |
| openmm-ir | numpy-ir | proper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | improper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj14 | -0.2752731221 | -0.2752731221 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -28.95685035 | -28.95685035 | 7.11e-15 | 1.29e-06 | True |
| openmm-ir | numpy-ir | total | 53.12833479 | 53.12833479 | 1.07e-13 | 1.48e-06 | True |
| openmm-xml | numpy-ir | bond | 48.05368473 | 48.05368473 | 6.39e-14 | 1.48e-06 | True |
| openmm-xml | numpy-ir | angle | 34.30677353 | 34.30677353 | 4.97e-14 | 1.34e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 0 | 0 | 0 | 1e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | -29.23212381 | -29.23212347 | 3.43e-07 | 1.29e-06 | True |
| openmm-xml | numpy-ir | total | 53.12833445 | 53.12833479 | 3.43e-07 | 1.48e-06 | True |
| openmm-prmtop | numpy-ir | bond | 48.05368473 | 48.05368473 | 6.39e-14 | 0.000105 | True |
| openmm-prmtop | numpy-ir | angle | 34.30677353 | 34.30677353 | 4.97e-14 | 0.000103 | True |
| openmm-prmtop | numpy-ir | torsion_all | 0 | 0 | 0 | 0.0001 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | -29.23212381 | -29.23212347 | 3.4e-07 | 0.000103 | True |
| openmm-prmtop | numpy-ir | total | 53.12833445 | 53.12833479 | 3.4e-07 | 0.000105 | True |
| openmm-top | numpy-ir | bond | 48.05368473 | 48.05368473 | 4.26e-14 | 0.000105 | True |
| openmm-top | numpy-ir | angle | 34.30677353 | 34.30677353 | 3e-09 | 0.000103 | True |
| openmm-top | numpy-ir | torsion_all | 0 | 0 | 0 | 0.0001 | True |
| openmm-top | numpy-ir | nonbonded_all | -29.23212381 | -29.23212347 | 3.42e-07 | 0.000103 | True |
| openmm-top | numpy-ir | total | 53.12833445 | 53.12833479 | 3.39e-07 | 0.000105 | True |
| sander-ours | numpy-ir | bond | 11.4851 | 11.48510629 | 6.29e-06 | 0.000201 | True |
| sander-ours | numpy-ir | angle | 8.1995 | 8.199515662 | 1.57e-05 | 0.000201 | True |
| sander-ours | numpy-ir | dihed | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | vdw14 | -0.0658 | -0.06579185519 | 8.14e-06 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -6.9206 | -6.920853334 | 0.000253 | 0.000546 | True |
| sander-tleap | sander-ours | bond | 11.4851 | 11.4851 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | angle | 8.1995 | 8.1995 | 0 | 0.000201 | True |
| sander-tleap | sander-ours | dihed | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | -0.0658 | -0.0658 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -6.9206 | -6.9206 | 0 | 0.000201 | True |
| gmx | numpy-ir | bond | 48.05369568 | 48.05368473 | 1.09e-05 | 0.000105 | True |
| gmx | numpy-ir | angle | 34.30678177 | 34.30677353 | 8.24e-06 | 0.000103 | True |
| gmx | numpy-ir | proper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | improper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj14 | -0.2752732038 | -0.2752731221 | 8.17e-08 | 0.0001 | True |
| gmx | numpy-ir | coulomb14 | -28.95684433 | -28.95685035 | 6.02e-06 | 0.000103 | True |
| gmx | numpy-ir | total | 53.12835312 | 53.12833479 | 1.83e-05 | 0.000105 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 5.29e-08 | 0.001 | True |
| openmm-top | openmm-ir | 5.11e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance

### config 3  passed=True

| leg | reference | category | E_leg | E_ref | delta | tol | passed |
|---|---|---|---|---|---|---|---|
| numpy-ir | reference-itself | bond | 7.599469008 | 7.599469008 | 0 | 1.08e-06 | True |
| numpy-ir | reference-itself | angle | 12.96587611 | 12.96587611 | 0 | 1.13e-06 | True |
| numpy-ir | reference-itself | proper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | improper | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb | 0 | 0 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | lj14 | -0.2948437013 | -0.2948437013 | 0 | 1e-06 | True |
| numpy-ir | reference-itself | coulomb14 | -28.48778091 | -28.48778091 | 0 | 1.28e-06 | True |
| numpy-ir | reference-itself | total | -8.21727949 | -8.21727949 | 0 | 1.08e-06 | True |
| openmm-ir | numpy-ir | bond | 7.599469008 | 7.599469008 | 0 | 1.08e-06 | True |
| openmm-ir | numpy-ir | angle | 12.96587611 | 12.96587611 | 3.91e-14 | 1.13e-06 | True |
| openmm-ir | numpy-ir | proper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | improper | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb | 0 | 0 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | lj14 | -0.2948437013 | -0.2948437013 | 0 | 1e-06 | True |
| openmm-ir | numpy-ir | coulomb14 | -28.48778091 | -28.48778091 | 7.11e-15 | 1.28e-06 | True |
| openmm-ir | numpy-ir | total | -8.21727949 | -8.21727949 | 3.2e-14 | 1.28e-06 | True |
| openmm-xml | numpy-ir | bond | 7.599469008 | 7.599469008 | 0 | 1.08e-06 | True |
| openmm-xml | numpy-ir | angle | 12.96587611 | 12.96587611 | 3.91e-14 | 1.13e-06 | True |
| openmm-xml | numpy-ir | torsion_all | 0 | 0 | 0 | 1e-06 | True |
| openmm-xml | numpy-ir | nonbonded_all | -28.78262495 | -28.78262461 | 3.37e-07 | 1.29e-06 | True |
| openmm-xml | numpy-ir | total | -8.217279827 | -8.21727949 | 3.37e-07 | 1.28e-06 | True |
| openmm-prmtop | numpy-ir | bond | 7.599469008 | 7.599469008 | 8.88e-16 | 0.000101 | True |
| openmm-prmtop | numpy-ir | angle | 12.96587611 | 12.96587611 | 4.09e-14 | 0.000101 | True |
| openmm-prmtop | numpy-ir | torsion_all | 0 | 0 | 0 | 0.0001 | True |
| openmm-prmtop | numpy-ir | nonbonded_all | -28.78262495 | -28.78262461 | 3.35e-07 | 0.000103 | True |
| openmm-prmtop | numpy-ir | total | -8.217279825 | -8.21727949 | 3.35e-07 | 0.000103 | True |
| openmm-top | numpy-ir | bond | 7.599469008 | 7.599469008 | 4.44e-15 | 0.000101 | True |
| openmm-top | numpy-ir | angle | 12.96587609 | 12.96587611 | 2.22e-08 | 0.000101 | True |
| openmm-top | numpy-ir | torsion_all | 0 | 0 | 0 | 0.0001 | True |
| openmm-top | numpy-ir | nonbonded_all | -28.78262495 | -28.78262461 | 3.36e-07 | 0.000103 | True |
| openmm-top | numpy-ir | total | -8.217279848 | -8.21727949 | 3.58e-07 | 0.000103 | True |
| sander-ours | numpy-ir | bond | 1.8163 | 1.816316685 | 1.67e-05 | 0.0002 | True |
| sander-ours | numpy-ir | angle | 3.0989 | 3.098918765 | 1.88e-05 | 0.0002 | True |
| sander-ours | numpy-ir | dihed | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-ours | numpy-ir | vdw14 | -0.0705 | -0.07046933587 | 3.07e-05 | 0.0002 | True |
| sander-ours | numpy-ir | eel14 | -6.8085 | -6.808743047 | 0.000243 | 0.00054 | True |
| sander-tleap | sander-ours | bond | 1.8163 | 1.8163 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | angle | 3.0989 | 3.0989 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | dihed | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel | 0 | 0 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | vdw14 | -0.0705 | -0.0705 | 0 | 0.0002 | True |
| sander-tleap | sander-ours | eel14 | -6.8085 | -6.8085 | 0 | 0.000201 | True |
| gmx | numpy-ir | bond | 7.599452496 | 7.599469008 | 1.65e-05 | 0.000101 | True |
| gmx | numpy-ir | angle | 12.96587753 | 12.96587611 | 1.42e-06 | 0.000101 | True |
| gmx | numpy-ir | proper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | improper | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | coulomb | 0 | 0 | 0 | 0.0001 | True |
| gmx | numpy-ir | lj14 | -0.2948437333 | -0.2948437013 | 3.2e-08 | 0.0001 | True |
| gmx | numpy-ir | coulomb14 | -28.48778152 | -28.48778091 | 6.14e-07 | 0.000103 | True |
| gmx | numpy-ir | total | -8.217294693 | -8.21727949 | 1.52e-05 | 0.000101 | True |

| leg | reference | max|dF| | tol | passed |
|---|---|---|---|---|
| openmm-xml | openmm-ir | 0 | 1e-06 | True |
| openmm-prmtop | openmm-ir | 4.65e-08 | 0.001 | True |
| openmm-top | openmm-ir | 5.14e-06 | 0.001 | True |
GROMACS precision=mixed; strict file tolerance
