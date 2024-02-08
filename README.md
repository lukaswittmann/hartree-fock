A simple implementation of the Hartree-Fock algorithm. STO-nG basis sets are used, which are generated from Slater exponents given in the input file.
Some parts like the input reader and straightforward functions were taken from [marvinfriede](https://github.com/marvinfriede)'s [hartree-fock implementation](https://github.com/marvinfriede/hartree-fock).
Currently, the following algorithms are implemented:

- RHF and UHF
- RMP2
- Numerical gradient
- Charge density
- Mulliken population analysis

This program was written in the context of teaching the Quantum Chemistry II course myself.
