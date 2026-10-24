# Models

Phonons are computed from harmonic springs between atoms. What they can and cannot describe:

| Model | Bond parameter | Notes |
|---|---|---|
| Central spring | K | Force only along the bond. Rotationally invariant. |
| Born–von Kármán | (K_L, K_T) | Adds stiffness K_T against displacements perpendicular to the bond. Not rotationally invariant: a free molecule or ribbon gets no zero rotation modes. |

## Limits

- **Floppy lattices.** Central springs to the nearest neighbours alone leave zero-frequency branches in the square, honeycomb, kagome and Lieb lattices, rocksalt (TA along Γ–X) and bcc (along Γ–N). Second-neighbour springs or K_T fix this.
- **No bending forces.** Without angle terms the graphene ZA branch is linear instead of ∝ q², and the flat TA and falling optical branches of Si are missed ([`silicon.py`](../../examples/crystals_3d/silicon.py)).
- **No long-range Coulomb forces.** Ionic crystals (NaCl, hBN) show no LO–TO splitting at Γ.
- **Harmonic only.** No thermal expansion, no phonon–phonon scattering.

Spring constants in the examples are fitted to a few measured frequencies or are illustrative; each example says which. Planned extensions are in [`TODO.md`](../../TODO.md).
