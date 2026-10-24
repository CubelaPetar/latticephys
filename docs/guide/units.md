# Units

Inputs: lengths in Å, masses in amu, energies in eV, spring constants in eV/Å², wave vectors in 1/Å.

Time then comes out in units of √(amu Å²/eV) ≈ 10.18 fs (`units.TIME_FS`), and angular frequencies ω in its inverse. Convert only at the output, with the factors in `latticephys.units`:

| Factor | Gives |
|---|---|
| `OMEGA_TO_THZ` | ν = ω / 2π in THz |
| `OMEGA_TO_MEV` | ħω in meV |
| `OMEGA_TO_INV_CM` | ν / c in cm⁻¹ |
| `VELOCITY_TO_KM_S` | group and sound velocities in km/s |

Frequencies are always ω, never ω²: zero modes are exactly 0, and an unstable mode shows up as a negative ω.
