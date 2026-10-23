# Copper (fcc) phonon frequencies at X and L (experiment)

These are room-temperature neutron values. The numbers come from Savrasov & Savrasov [3, Table I, "exp." rows], which takes them from the Landolt-Börnstein compilation [4]. Which measurement [4] used is **not confirmed**. The numbers are most likely the room-temperature data of Svensson, Brockhouse & Rowe [1], the standard Cu dataset that [5] also compares against. I could not read the original tables in [1] or [2] (APS paywall; McMaster thesis copies are behind a bot check). The values are therefore **secondary-source** values.

1 THz = 33.356 cm⁻¹ = 4.1357 meV.

| Point | Branch | ν (THz) | cm⁻¹ | meV | Source |
|---|---|---|---|---|---|
| X | TA (doubly degenerate) | 5.13 | 171.1 | 21.22 | [3] (exp.), citing [4] |
| X | LA | 7.25 | 241.8 | 29.98 | [3] (exp.), citing [4] |
| L | TA (doubly degenerate) | 3.42 | 114.1 | 14.14 | [3] (exp.), citing [4] |
| L | LA | 7.30 | 243.5 | 30.19 | [3] (exp.), citing [4] |

[3] gives no error bars. The cm⁻¹ and meV columns are my own conversions.

**Temperature and lattice constant:** The abstract of [1] says "room temperature". Xu et al. [5] describe the neutron data of [1] as "taken at 296 K". The lattice constant used in [1] is **not confirmed**. a = 3.615 Å is the usual room-temperature value, but I did not trace it to a source.

## Check against the nearest-neighbour model

The nearest-neighbour central-spring model of `examples/crystals_3d/copper.py` has K fitted to the measured LA(X) = 7.25 THz (output of the example):

| Point | Model (THz) | Measured (THz) | Model − measured |
|---|---|---|---|
| X TA | 5.13 | 5.13 | 0 (model: LA/√2) |
| X LA | 7.25 (fit) | 7.25 | 0 |
| L TA | 3.62 | 3.42 | +6 % |
| L LA | 7.25 | 7.30 | −0.7 % |

- In this model LA(L) = LA(X); the measured values differ by 0.7 %. The largest error of the model is TA(L), about 6 % too high. In [1] this is attributed to forces beyond the nearest neighbour (abstract: weaker forces out to at least sixth neighbours).

## Discrepancies

- The measured LA(L) of 7.30 THz is the likely origin of "about 7.3 THz" being quoted for LA(X). The measured LA(X) is 7.25 THz.
- Nicklow et al. [2] measured Cu at 49 K and 298 K. Their abstract says that frequencies above ~3 THz shift with temperature "often 1 % or less". The 49 K values at X and L should therefore lie up to ~1 % above the room-temperature values. A 49 K LA(X) close to 7.3 THz is plausible, but **not confirmed**.
- Nilsson & Rolandson [6] measured Cu at 80 K, mostly at off-symmetry points. I could not read their values.

## Not confirmed

- All values in the table: they were taken from [3], which cites the compilation [4], and not from the tables of [1].
- The error bars and lattice constant of [1].
- Values at 49 K [2] and 80 K [6].
- Frequencies at K and W: none of the sources I could read tabulate them. [1] measured along [0ζ1] (X–W), according to its abstract, but I did not see the values.

## Citations

1. E. C. Svensson, B. N. Brockhouse, J. M. Rowe, "Crystal dynamics of copper", Phys. Rev. **155**, 619 (1967). doi:10.1103/PhysRev.155.619. *(abstract only read)*
2. R. M. Nicklow, G. Gilat, H. G. Smith, L. J. Raubenheimer, M. K. Wilkinson, "Phonon frequencies in copper at 49 and 298°K", Phys. Rev. **164**, 922 (1967). doi:10.1103/PhysRev.164.922. *(abstract only read)*
3. S. Y. Savrasov, D. Y. Savrasov, "Electron-phonon interactions and related physical properties of metals from linear-response theory", Phys. Rev. B **54**, 16487 (1996), Table I. arXiv:cond-mat/9604124. *(secondary for the experimental values; cites [4])*
4. Landolt-Börnstein, New Series, Group III, Vol. 13a (Springer, Berlin, 1981). *(compilation, not read)*
5. R. Xu, H. Hong, P. Zschack, T.-C. Chiang, "Direct mapping of phonon dispersion relations in copper by momentum-resolved x-ray calorimetry", Phys. Rev. Lett. **101**, 085504 (2008). doi:10.1103/PhysRevLett.101.085504. *(source for 296 K; data only shown in a figure)*
6. G. Nilsson, S. Rolandson, "Lattice dynamics of copper at 80 K", Phys. Rev. B **7**, 2393 (1973). doi:10.1103/PhysRevB.7.2393. *(abstract only read)*
