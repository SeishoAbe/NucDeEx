# Nuclear Deexcitation Simulator (NucDeEx)
NUcDeEx is a nuclear deexcitation simulator for neutrino interactions and nucleon decays shared in [GitHub](https://github.com/SeishoAbe/NucDeEx).  
It is composed of two parts: [TALYS v1.96](https://tendl.web.psi.ch/tendl_2019/talys.html) and a kinematics simulator based on [ROOT](https://root.cern/).  
Note that this simulator reads files containing branching ratios calculated with TALYS.  
Therefore, the simulator codes themselves are independent of the TALYS code.

## Contact 
Seisho Abe (ICRR, the University of Tokyo): seisho@icrr.u-tokyo.ac.jp

## Citing NucDeEx
**If you use NucDeEx, please always cite this paper [URL](https://arxiv.org/abs/2310.10394).**
```
@misc{abe2023nuclear,
      title={Nuclear deexcitation simulator for neutrino interactions and nucleon decays of $^{12}\text{C}$ and $^{16}\text{O}$ based on TALYS}, 
      author={Seisho Abe},
      year={2023},
      eprint={2310.10394},
      archivePrefix={arXiv},
      primaryClass={hep-ph}
}
```
It would be better to also cite [TALYS' paper](https://doi.org/10.1016/j.nds.2012.11.002).  
In addition, please consider citing the following references, explaining the predecessor deexcitation simulator to NucDeEx.
```
@article{Abe_2021,
doi = {10.1088/1742-6596/2156/1/012189},
url = {https://dx.doi.org/10.1088/1742-6596/2156/1/012189},
year = {2021},
month = {dec},
publisher = {IOP Publishing},
volume = {2156},
number = {1},
pages = {012189},
author = {Seisho Abe and the KamLAND Collaboration},
title = {Nuclear de-excitation associated with neutrino-carbon interactions},
journal = {Journal of Physics: Conference Series},
}
```
```
@article{PhysRevD.107.072006,
  title = {First measurement of the strange axial coupling constant using neutral-current quasielastic interactions of atmospheric neutrinos at KamLAND},
  author = {Abe, S. and others},
  collaboration = {KamLAND Collaboration},
  journal = {Phys. Rev. D},
  volume = {107},
  issue = {7},
  pages = {072006},
  numpages = {17},
  year = {2023},
  month = {Apr},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevD.107.072006},
  url = {https://link.aps.org/doi/10.1103/PhysRevD.107.072006}
}

```


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

## Deexcitation Simulation
You need to declare "Deexcitation" object and call a function:
```
Deexcitation::DoDeex(Zt,Nt,Z,N,shell,Ex,Pinit)
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
  - Reference data (other predictions, experimental data) are stored.
