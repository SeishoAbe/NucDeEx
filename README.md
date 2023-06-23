# Nuclear Deexcitation Simulator for Neutrino Interaction and Proton Decay

## Contact 
Seisho Abe: seisho@km.icrr.u-tokyo.ac.jp

## How to simulate 
```
make 
./bin/simulation (TARGET) 2 1 
```
argv[2]: Flag for level density model. "2 (Back-shifted Fermi gas)" is recommended
argv[3]: Flag for optical model. "1" is recomennded.


## Directory

## Code
- `main/`: main codes for TALYS calulation and event simulation
- `src/` & `incldue/`: TALYS output analyzer and simulator codes
- `bin/`: Executables
- `lib/`: Library directory
- `obj/`: Object file directory

## Tables, etc.
- `output/`: Files converted from TALYS' output
- `sim_out/`: Event simulation output 
- `tables/`: 
	- `nucleus/`: List of nucleus considered in TALYS calculation
  - `separation_energy/`: Separation energy tables obtained from TALYS
  - `energy_distribution/`: Energy distrubition file inputed into TALYS
	- `sf/`: Benhar SF used in simulation
