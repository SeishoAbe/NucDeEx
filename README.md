# Nuclear Deexcitation Simulator
For neutrino interaction and proton decay

## Contact 
Seisho Abe: seisho@km.icrr.u-tokyo.ac.jp

## Simulation
### How to run the simulation
```
make 
./bin/simulation (TARGET) 2 1 
```
- argv[1]: Target nucleus. 11C, 11B, 15O, 15N is available as of June 23rd.
- argv[2]: Flag for level density model. "2 (Back-shifted Fermi gas)" is recommended
- argv[3]: Flag for the optical model. "1" is recommended. <br>

### Notes
- This simulation reads files that contain branching ratios, which are converted from TALYS output.
- Therefore, these simulator codes are independent of TALYS codes.

## TALYS information & calculation
- The branching ratios are calculated using [TALYS](https://tendl.web.psi.ch/tendl_2019/talys.html).
- If you want to use TALYS, please follow the documents in tar and install it.
- Then you can use the following scripts and codes 
```
 run_talys.sh
 ./bin/plot_decay
```
- `./bin/plot_decay` is an analyzer on TALYS output.
- The files created by this executable are used in the simulator (`./bin/simulation`).

## Directory

### Codes
- `main/`: main codes for TALYS calculation and event simulation
- `src/` & `include/`: TALYS output analyzer and simulator codes
- `bin/`: Executables
- `lib/`: Library directory
- `obj/`: Object file directory
- `macro/`: ROOT macro files to check the output

### Tables, output, etc.
- `output/`
	- Files converted from TALYS' output
	- These contain branching ratios.
- `sim_out/`
	- Event simulation output 
- `tables/`
  - `nucleus/`
    - List of nuclei considered in TALYS calculation
  - `separation_energy/`
	- Separation energy tables obtained from TALYS
  - `energy_distribution/`
	- Energy distribution file input into TALYS
  - `sf/`
	- Benhar SF used in event simulation
- `data/`
  - Reference data (other predictions, experimental data) are stored.
