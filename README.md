# Nuclear Deexcitation Simulator (NucDeEx)
NUcDeEx is a nuclear deexcitation simulator for neutrino interactions and nucleon decays shared in [GitHub](https://github.com/SeishoAbe/NucDeEx).  
It is composed of two parts: [TALYS v1.96](https://tendl.web.psi.ch/tendl_2019/talys.html) and a kinematics simulator based on [ROOT](https://root.cern/).  
Note that this simulator reads files containing branching ratios calculated with TALYS.  
Therefore, the simulator codes themselves are independent of the TALYS code.

## Contact 
Seisho Abe (Kamioka Obs., ICRR, the University of Tokyo): seisho@icrr.u-tokyo.ac.jp

## Citing NucDeEx
**If you use NucDeEx, please always cite [this paper](https://link.aps.org/doi/10.1103/PhysRevD.109.036009) and [TALYS's paper](https://doi.org/10.1140/epja/s10050-023-01034-3).**
```
@article{PhysRevD.109.036009,
  title = {Nuclear deexcitation simulator for neutrino interactions and nucleon decays of $^{12}\mathrm{C}$ and $^{16}\mathrm{O}$ based on TALYS},
  author = {Abe, Seisho},
  journal = {Phys. Rev. D},
  volume = {109},
  issue = {3},
  pages = {036009},
  numpages = {13},
  year = {2024},
  month = {Feb},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevD.109.036009},
  url = {https://link.aps.org/doi/10.1103/PhysRevD.109.036009}
}
```
```
@Article{Koning2023,
  author={Koning, Arjan and Hilaire, Stephane and Goriely, Stephane},
  title={TALYS: modeling of nuclear reactions},
  journal={The European Physical Journal A},
  year={2023},
  month={Jun},
  day={14},
  volume={59},
  number={6},
  pages={131},
  issn={1434-601X},
  doi={10.1140/epja/s10050-023-01034-3},
  url={https://doi.org/10.1140/epja/s10050-023-01034-3}
}
```
In addition, please consider citing the following references, explaining [the predecessor deexcitation simulator to NucDeEx](https://dx.doi.org/10.1088/1742-6596/2156/1/012189) 
and [application to atmospheric neutrino analysis at KamLAND](https://link.aps.org/doi/10.1103/PhysRevD.107.072006).

## Compile

This simulator depends on **ROOT**.
At least it works with the following gcc versions and ROOT versions.
- gcc 7.5.0 and ROOT v6.18.04 (ubuntu 18.04)
- gcc 8.5.0 and ROOT v5.34.38 (redhut 8.5.0-4)

To build NucDeEx, after preparing general root enviroment, just type:
```
source setup.sh
make
```

## Run Simulation
You need to declare `NucDeExDeexcitation` object and call a function:
```
NucDeExDeexcitation::DoDeex(Zt,Nt,Z,N,shell,Ex,Pinit)
```
- Zt(Nt): Z(N) number of the target nucleus
- Z(N): Z(N) number of the residual nucleus
- shell: shell flag for the hole state
  - 0: Automatically determined according to the excitation energy.
  - 1: s1/2-hole state
  - 2: p3/2-hole state
  - 3: p1/2-hole state
- Ex: Excitation energy in MeV 
- Pinit: 3D momentum (TVector3) of the residual nucleus

An example code for simulation is provided:
```
./bin/simulation (TARGET) 2 1 1 (SEED)
```
- argv[1]: Target nucleus. 11C, 11B, 15O, 15N are available now.
- argv[2]: Flag for level density model. "2 (Back-shifted Fermi gas)" is recommended
- argv[3]: Flag for the optical model. "1" is recommended.
- argv[4]: Flag for J^\pi specification. Set "1".
- argv[5]: Random seed (optional). If this is not specified, seed "1" is set. <be>

## TALYS information & calculation
- The branching ratios are calculated using [TALYS v1.96](https://tendl.web.psi.ch/tendl_2019/talys.html).
- If you want to use TALYS, please follow the documents in tar and install it. 
- Then, you can use the following scripts and codes 
```
 run_talys.sh
 ./bin/plot_decay
```
- `./bin/plot_decay` is an analyzer on TALYS output.
- The files created by this executable are used in the simulator (`./bin/simulation`).

## Interfaces and build scripts for INCLXX and Geant4
NucDeEx also provides intefaces and build scripts for use in [INCL++](https://irfu.cea.fr/dphn/Spallation/incl.html) and [Geant4](https://geant4.web.cern.ch/).
### INCL++
- Interfaces: `include/G4INCLNucDeExInterface.hh`, `src/G4INCLNucDeExInterface.cc`
- Build scrits: `./inclxx`
### Geant4
- Interfaces: `include/G4NucDeExInterface.hh`, `src/G4NucDeExInterface.cc`
- Build scrits: `./geant4`

## Directory

### Codes
- `./main`: main codes for TALYS calculation and event simulation
- `./src` & `./include`: TALYS output analyzer and simulator codes
- `./bin`: Executables
- `./lib`: Library directory
- `./obj`: Object file directory
- `./macro`: ROOT macro files to check the output. This is a private submodule.

### Tables, output, etc.
- `./output`
	- Files converted from TALYS' output
	- These contain branching ratios.
- `./sim_out`
	- Event simulation output 
- `./tables`
  - `./nucleus`
    - List of nuclei considered in TALYS calculation
  - `./separation_energy`
	- Separation energy tables obtained from TALYS
  - `./energy_distribution`
	- Energy distribution file input into TALYS
  - `./sf`
	- Benhar SF used in event simulation
- `./data`
  - Reference data (other predictions, experimental data) are stored. A private submodule.
- `./inclxx`
  - Build scripts for use in INCLXX.
- `./geant4`
  - Build scripts for use in Geant4.
