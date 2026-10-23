# Bilayer graphene (AB): interlayer modes at Γ

Measured Raman frequencies of the two rigid-layer modes of Bernal bilayer graphene (2LG), with the interlayer force constants per unit area derived from them. 1 THz = 33.356 cm⁻¹.

## Frequencies

| Mode | cm⁻¹ | THz | Source |
|---|---|---|---|
| Shear (C mode, E_g) | ~31 | 0.93 | Tan 2012 (abstract and text: "Pos(C)∼31 cm⁻¹ for BLG") |
| Shear, from the same 2LG spectrum | 30 | 0.90 | Lui 2014 (peak list "30, 81, 113, 174, 198 cm⁻¹") |
| Layer breathing (LBM, ZO′), fundamental | 81 | 2.43 | Lui 2014 (linewidth 4.5 cm⁻¹, suspended 2LG, seen only above ~600 K) |
| LBM extrapolated from 2ZO′ overtone to q = 0 | 80 ± 2 | 2.40 ± 0.06 | Lui & Heinz 2013 |
| Bulk graphite shear (E₂g) | ~43.5 | 1.30 | Tan 2012 |

No primary measurement gives 32 cm⁻¹ for the bilayer shear mode: Tan 2012 says ~31 cm⁻¹.

## Force constants per unit area

Linear-chain model with mass per area μ = 7.6×10⁻²⁷ kg/Å² (single-layer graphene, Tan 2012 and Wu 2015). Bilayer: ω = √(2α/μ). Graphite (N → ∞): ω = 2√(α/μ).

| Constant | Value (N/m³) | Source | Check: bilayer ω from √(2α/μ) |
|---|---|---|---|
| Shear α∥ | 12.8×10¹⁸ | Tan 2012 | 30.8 cm⁻¹ |
| Breathing α⊥ (nearest layer) | 106×10¹⁸ | Wu 2015, fitted to LBMs of twisted 4LG, not to 2LG | 88.7 cm⁻¹ |
| Breathing α⊥ (2nd-nearest layer model), plus β⊥ ≈ 9.3×10¹⁸ | 116×10¹⁸ (DFT average 115.6×10¹⁸) | Wu 2015 | 92.6 cm⁻¹ (α only) |
| Breathing α⊥ needed to match 81 cm⁻¹ | ≈ 88×10¹⁸ | my own calculation, not printed in any source | 81 cm⁻¹ |

Lui & Heinz 2013 and Lui 2014 give no force constant, only the bulk frequency ω₀ of their chain model (ω_bilayer = ω₀/√2).

## Geometry

- Graphite, room temperature: a = 2.464 Å, c = 6.711 Å, so the interlayer distance is c/2 = 3.356 Å (Trucano & Chen 1975). I took these values from the AMCSD and AIST crystal-structure databases, which cite that paper. I did not read the paper itself.
- Lui & Heinz 2013 use "c ~ 0.34 nm" for the interlayer spacing.
- I found no measured interlayer distance for the bilayer itself. The usual 3.35 Å is the graphite value.
- 2.46 Å for the in-plane lattice constant is the rounded graphite value.

## Discrepancies

- **Breathing mode, 81 vs. ~89 cm⁻¹.** The bilayer measurement (80–81 cm⁻¹) implies α⊥ ≈ 88×10¹⁸ N/m³. Wu 2015 fits α⊥ = 106×10¹⁸ N/m³ to thicker twisted samples, and that value gives ~89 cm⁻¹ for 2LG. A chain model with only nearest-layer springs does not fit 2LG and thick samples with one α⊥. Wu 2015 adds a second-neighbour spring to fix this.
- **Bulk value in Lui's two papers.** Lui & Heinz 2013 fit ω₀ = 132.3 cm⁻¹ from the 2ZO′ overtone of graphite (half of 264.5 cm⁻¹). This predicts 93.5 cm⁻¹ for 2LG, against their own extrapolated 80 ± 2 cm⁻¹. Lui 2014 instead takes ω₀ = 114.6 cm⁻¹ from the 81 cm⁻¹ bilayer peak. Neutron scattering gives ~128 cm⁻¹ for the graphite B₂g mode (as cited by Wu 2015, not checked).
- **Shear mode.** Measurements agree on 30–31 cm⁻¹.
- **Example to compare against.** For the shear mode, use α∥ = 12.8×10¹⁸ N/m³ against ~31 cm⁻¹. For the breathing mode, use 81 cm⁻¹ (directly measured), and note that the published α⊥ (106×10¹⁸) predicts ~89 cm⁻¹.

## Not confirmed

- A bilayer shear frequency of 32 cm⁻¹: the primary source says ~31 cm⁻¹.
- A force constant for the breathing mode fitted directly to the bilayer: none of the papers read prints one.
- A measured interlayer distance for the bilayer itself.
- The Trucano & Chen lattice constants: taken from database entries that cite the paper, not from the paper.
- The 89 cm⁻¹ LBM value that appears in some secondary sources: not traced to a primary measurement.

## Citations

- P. H. Tan, W. P. Han, W. J. Zhao, Z. H. Wu, K. Chang, H. Wang, Y. F. Wang, N. Bonini, N. Marzari, G. Savini, A. Lombardo, A. C. Ferrari, "The shear mode of multilayer graphene", *Nature Materials* **11**, 294–300 (2012). doi:10.1038/nmat3245. arXiv:1106.1146. (Authors as on arXiv v1; the published author list was not checked.)
- C. H. Lui, T. F. Heinz, "Layer breathing modes in few-layer graphene" (arXiv title; published title not checked), *Phys. Rev. B* **87**, 121404(R) (2013). doi:10.1103/PhysRevB.87.121404. arXiv:1210.0960.
- C. H. Lui, Z. Ye, C. Keiser, X. Xiao, R. He, "Temperature-activated layer-breathing vibrations in few-layer graphene", *Nano Letters* **14**, 4615 (2014). doi:10.1021/nl501678j. arXiv:1405.6735.
- J.-B. Wu, Z.-X. Hu, X. Zhang, W.-P. Han, Y. Lu, W. Shi, X.-F. Qiao, M. Ijäs, S. Milana, W. Ji, A. C. Ferrari, P.-H. Tan, "Interface coupling in twisted multilayer graphene by resonant Raman spectroscopy of layer breathing modes", *ACS Nano* **9**, 7440–7449 (2015). doi:10.1021/acsnano.5b02502. arXiv:1505.00095.
- P. Trucano, R. Chen, "Structure of graphite by neutron diffraction", *Nature* **258**, 136–137 (1975). doi:10.1038/258136a0.
