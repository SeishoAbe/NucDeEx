# Nuclear Deexcitation Simulator (NucDeEx)
This is a nuclear deexcitation simulator (NucDeEx) for neutrino interactions and nucleon decays.  
It is composed of two parts: [TALYS v1.96](https://tendl.web.psi.ch/tendl_2019/talys.html) and a kinematics simulator based on [ROOT](https://root.cern/).  
Note that this simulator reads files containing branching ratios calculated with TALYS.  
Therefore, the simulator codes themselves are independent of the TALYS code.


## Contact 
Seisho Abe (ICRR, the University of Tokyo): seisho@icrr.u-tokyo.ac.jp

## Citing NucDeEx
Please cite this on arXiv [URL](https://arxiv.org/abs/2310.10394).
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
You also should also cite [TALYS' paper](https://doi.org/10.1016/j.nds.2012.11.002).

## Compile

This simulator depends on ROOT.
At least it works with the following gcc versions and ROOT versions.
- gcc 7.5.0 and ROOT v6.18.04 (ubuntu 18.04)
- gcc 8.5.0 and ROOT v5.34.38 (redhut 8.5.0-4)

To build NucDeEx, just type:
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
- Ex: Excitation energy in MeV 
- Pinit: 3D momentum of the residual nucleus

An example code for simulation is provided:
```
./bin/simulation (TARGET) 2 1 
```
- argv[1]: Target nucleus. 11C, 11B, 15O, 15N is available as of June 23rd.
- argv[2]: Flag for level density model. "2 (Back-shifted Fermi gas)" is recommended
- argv[3]: Flag for the optical model. "1" is recommended. <br>


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
- `./src` & `include/`: TALYS output analyzer and simulator codes
- `./bin`: Executables
- `./lib`: Library directory
- `./obj`: Object file directory
- `./macro`: ROOT macro files to check the output

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
