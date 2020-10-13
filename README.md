# matlab-lammps-phasediagram

Code for the article "Decoding the Physical Principles of Two-Component Biomolecular Phase Separation".

All data in the article are generated using the molecular-dynamics simulation platform LAMMPS (https://lammps.sandia.gov/). Lammps input files are produced with the MatLab codes in this repository. 

Running order: 
Units_nano.m -> InitialState_Valence_MediumSystem.m -> In_Valence_MediumSystem_Anneal.m -> In_Valence_MediumSystem_Relax.m -> In_Valence_MediumSystem_Record.m.
