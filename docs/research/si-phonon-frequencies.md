# Silicon phonon frequencies at Γ, X, L (experiment)

Neutron values of Nilsson & Nelin [1] (with Dolling [2]). The numbers come from a table in Savrasov [3, Table I], which cites [2] and Nilsson & Nelin. All four are checked against the rounded values in the ioffe database [4]. I could not read the original table in [1] (paywalled), so the values are **secondary-source** values.

| Point | Branch | ν (THz) | Original | Source |
|---|---|---|---|---|
| Γ | LO = TO (Γ25') | 15.53 | 15.53 THz | [1,2] via [3]; [4]: 15.5 |
| X | TA (X3) | 4.49 | 4.49 THz | [1,2] via [3]; [4]: 4.5 |
| X | LA = LO (X1) | 12.32 | 12.32 THz | [1,2] via [3]; [4]: 12.3 |
| X | TO (X4) | 13.90 | 13.90 THz | [1,2] via [3]; [4]: 13.9 |
| L | TA (L3) | 3.43 | 3.43 THz | [1,2] via [3]; [4]: 3.45 |
| L | LA (L2') | 11.35 | 11.35 THz | [1,2] via [3]; [4]: 11.3 |
| L | LO (L1) | 12.60 | 12.60 THz | [1,2] via [3]; [4]: 12.6 |
| L | TO (L3') | 14.68 | 14.68 THz | [1,2] via [3]; [4]: 14.7 |
| Γ (Raman) | LO = TO | 15.60 | 520.2 ± 0.5 cm⁻¹ at 300 K | [5] (abstract) |

Error bars from [1] are **not confirmed**. In diamond, LA and LO are degenerate at X.

**Temperature and lattice constant:** room temperature. [4] gives T = 300 K for the Γ value. The exact temperature in [1] is **not confirmed**: their earlier Ge work was at 80 K, but the Si data are commonly quoted as room temperature. The lattice constant is a = 5.431 Å at 300 K [4], and 5.4309 Å in the IFJ phonon database [6] (citing Yim & Paff 1974). The value used in [1] is **not confirmed**.

## Discrepancies

- Γ optical: 15.53 THz (neutron [1,2] via [3]) vs 15.60 THz (Raman, 520.2 cm⁻¹ [5]). The 0.4 % gap is within neutron resolution. Raman is the more precise number for Γ. Use 15.53 THz if all points should come from one neutron dataset (recommended for a consistent fit). A value of 15.3 THz was not found in any source checked, so it is **not confirmed**.
- L TA: 3.43 THz [3] vs 3.45 THz [4]. [4] is rounded or comes from another source; use 3.43.
- [3] cites Nilsson & Nelin as "Phys. Rev. B 8, 5711 (1973)". That reference is wrong: PRB 8, 5711 is a paper on AlAs/AlP. The correct reference is [1].
- Kulda et al. [7] measured LA and LO frequencies and eigenvectors along the three main directions (abstract). I could not access their table or temperature, so their values are **not confirmed** and not used here.

## References

1. G. Nilsson and G. Nelin, "Study of the homology between silicon and germanium by thermal-neutron spectrometry", Phys. Rev. B **6**, 3777 (1972). doi:10.1103/PhysRevB.6.3777
2. G. Dolling, "Lattice vibrations in crystals with the diamond structure", in *Inelastic Scattering of Neutrons in Solids and Liquids*, Vol. II, p. 37 (IAEA, Vienna, 1963).
3. S. Y. Savrasov, "Linear-response theory and lattice dynamics: a muffin-tin-orbital approach", Phys. Rev. B **54**, 16470 (1996), Table I. arXiv:cond-mat/9604125. *(secondary for the experimental values)*
4. Ioffe Institute, NSM semiconductor database, Si: <https://www.ioffe.ru/SVA/NSM/Semicond/Si/mechanic.html> (cites Dolling 1963; Tubino, Piseri, Zerbi, J. Chem. Phys. **56**, 1022 (1972)). *(secondary)*
5. J. H. Parker Jr., D. W. Feldman, M. Ashkin, "Raman scattering by silicon and germanium", Phys. Rev. **155**, 712 (1967). doi:10.1103/PhysRev.155.712
6. IFJ phonon dispersion database, Si: <https://wolf.ifj.edu.pl/phonondb/si.htm>. *(secondary)*
7. J. Kulda, D. Strauch, P. Pavone, Y. Ishii, "Inelastic-neutron-scattering study of phonon eigenvectors and frequencies in Si", Phys. Rev. B **50**, 13347 (1994). doi:10.1103/PhysRevB.50.13347
