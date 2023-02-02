# QuDiPy (Quantum Dots in Python)

This respository aims to be a general purpose tool for simulating electron dynamics in quantum dots with an emphasis on quantum information processing. We adopt an effective mass formalism which enables simulations across a variety of material systems and heterostructures. The code base is under construction by the [Coherent Spintronics Group](https://research.iqc.uwaterloo.ca/csg/) (CSG) at the University of Waterloo. The two codebase maintainers are Bohdan Khromets (@hromecB) and Zach D. Merino (@zmerino). Reach out to either of us if you have questions regarding this project.

Version 1.0 is currently under development. When completed, we aim to have the following capabilities:
- Integration with [nextnano](https://www.nextnano.de/) for the simulation of realistic nanoscale devices
- Real space (orbital) simulations:
  - Determination of electric potentials and fields in realistic quantum-dot based devices
  - Time evolution simulations of electron transfers
- Estimation of the Stark shift in silicon heterostructures
- Calculation of the many-electron spectra of general quantum dot networks using a [modified LCHO-CI approach](https://arxiv.org/pdf/2012.10512.pdf). These energy spectra can be mapped onto:
  -    Heisenberg Hamiltonian to determine the pairwise exchange interactions of electrons in the dot network, _or_
  -    Hubbard Hamiltonian to determine charge configuration and charge stability regions
- General $n$-spin system simulation (unitary _or_ non-unitary within Lindblad formalism) under time-varying ESR field, $g$-factors and exchange couplings in an effective spin Hamiltonian
- Simulation of charge stability diagrams using either a constant interaction model or a Hubbard Hamiltonian model. In addition, this module can fit charge stability data to extract capacitances for the constant interaciton model
- Tools for finding optimal control voltage pulses to realize desired quantum algorithms:
  -   constrained [effective shape engineering](https://uwspace.uwaterloo.ca/handle/10012/17823) with 1-to-1 mapping between effective and voltage parameters
  -   GRAPE
  -   constant-adiabaticity pulses [for electron shuttling](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.102.125406)

You can see the general progress in the module by going to the tutorials folder and following along the Jupyter notebooks.  Currently, some of the data for the tutorials are hosted elsewhere.  If you would like these available to you, please email either bohdan.khromets@uwaterloo.ca or zmerino@uwaterloo.ca (subject line: QuDiPy).  We emphasize that this is work in progress.

For development guidelines, see https://github.com/mainCSG/QuDiPy/blob/master/Development%20Guidelines.md.
