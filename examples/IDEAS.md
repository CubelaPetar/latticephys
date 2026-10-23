# Planned examples

Sorted by folder. **Needs** names the library feature from `TODO.md` that must come first; no entry means the library can do it already. **Test** is the quantitative check against theory.

## Molecules

- **Unequal springs (e.g. acetylene H–C≡C–H)** via `force_constants`; stretch modes in cm⁻¹ against experiment.
- **CO₂ bending mode, H₂O (three modes).** Needs: angle-bending forces (central springs give no restoring force for bending a linear molecule).

## 1D chains

- **Spring-dimerized chain** (equal masses, K1 ≠ K2): gap at the zone edge, phonon analogue of the SSH model. Test: analytic two-branch ω(q).
- **Second-neighbour springs** `(0, 0, 2, K2)`: the dispersion changes shape; negative K2 makes the chain unstable.
- **Backward waves:** a wave packet on the optical branch, where ω falls with q, so crests and packet move in opposite directions. Test: group velocity dω/dq < 0.
- **Integrator accuracy:** energy error against time step, velocity Verlet vs. explicit Euler.
- **Fermi–Pasta–Ulam–Tsingou chain:** energy returns to the first mode instead of spreading. Needs: `simulate()` with a force function.

## 2D lattices

Animations here should colour atoms by displacement: `animate(..., color_scale=...)`, as in `wave_2d.py`.

- **Refraction and reflection at an interface:** heavier half of a triangular supercell, plane wave packet at an angle; total internal reflection above the critical angle. Test: Snell's law sin θ₁ / c₁ = sin θ₂ / c₂.
- **Rayleigh edge waves:** a kick at the edge of a wide strip travels along the edge, slower than c_T. Test: lowest strip branch at large q against the Rayleigh equation (2D value not yet known to us, compute it).
- **Phononic crystal and waveguide:** a periodic pattern of heavy atoms opens a gap; removing one row guides a wave at a gap frequency.
- **Thermal motion:** all modes with random amplitude and phase for temperature T; slow long-wave drifts (Mermin–Wagner in real time). Test: time-averaged ⟨u²⟩ equals `mean_square_displacement`. Needs: `thermal_state()`.
- **Phonon focusing:** anisotropic wave fronts on the square lattice with second neighbours. Test: front shape against the group-velocity surface.
- **Isotope defect in 2D:** localized mode above the band, like `chains/mass_defect.py`.
- **Chiral phonons at K in hBN:** atoms move on circles, opposite sense at K and K′. Needs: travelling-wave start in `bloch_mode`.
- **Topological edge modes (twisted kagome):** zero modes on one edge of a strip only (Kane & Lubensky). Geometry still to look up.
- **Strained lattices:** pre-tension gives transverse stiffness to otherwise floppy lattices. Needs: extended spring model.
- **Anderson localization with mass disorder:** high-frequency waves stay near the source. Needs: sparse matrices (2D localization lengths need large systems).

## Textbook exercises

From G. Czycholl, *Theoretische Festkörperphysik* (2016), copy in `references/` (not in git). The library features they need are still to be checked against the exercise texts.

- **Chapter 2 exercises**, especially 2.4.
- **Chapter 4 exercises**, especially 4.8 and 4.9.
- **NaCl dispersion against Abb. 4.6:** compare `crystals_3d/nacl.py` with the book's curves.
