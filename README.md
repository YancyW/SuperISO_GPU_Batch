# SuperIso Version 5.0 (November 2025)
*By Farvah Nazila Mahmoudi (mahmoudi@in2p3.fr) and Siavash Neshatpour (siavash.neshatpour@gmail.com)*


## INTRODUCTION

The most recent version of this program can be obtained from:  
http://superiso.in2p3.fr

This C code (in C99 standard) provides the calculation of several flavour physics observables as well as the muon anomalous magnetic moment, for the Standard Model and a set of supersymmetric or 2HDM parameters which can be either calculated within this program using the provided interfaces with Isajet, Softsusy, SuSpect, SPheno, NMSSMTools or 2HDMC, or directly supplied by the user thanks to a SUSY Les Houches Accord (SLHA1 or SLHA2) file.

If you use this program to write a paper, please cite:

- F. Mahmoudi, arXiv:0710.2067 [hep-ph], Comput. Phys. Commun. 178 (2008) 745
- F. Mahmoudi, arXiv:0808.3144 [hep-ph], Comput. Phys. Commun. 180 (2009) 1579
- F. Mahmoudi, Comput. Phys. Commun. 180 (2009) 1718
- S. Neshatpour, F. Mahmoudi, arXiv: 2105.03428 [hep-ph], PoS TOOLS2020 (2021) 036

The program has been tested on Linux and Mac machines with gcc, icc and clang.

## Installation and Compilation

```bash
# Extract the archive
tar xzvf superiso_vX.X.tgz

# Navigate to the directory
cd superiso_vX.X

# In Makefile, first define your C compiler, and the path to Isajet, Softsusy, SPheno, 
# SuSpect, NMSSMTools, or 2HDMC, if you want to use them. 
# If not, you must comment them in the main programs.

# Compile with:
make
```

17 main programs are available:

| Program | Description |
|---------|-------------|
| `predict.c` | computation of SM or NP predictions with uncertainties |
| `scan_coeffs.c` | tool for scanning over Wilson coefficients |
| `sm.c` | calculation of a few observables in the SM |
| `sm_chi2.c` | calculation of the b->sll observables with chi2 evaluation in the SM |
| `cmssm.c` | calculation of the observables for a CMSSM point, linked to isajet, softsusy, spheno or suspect |
| `amsb.c` | calculation of the observables for an AMSB point, linked to isajet (>=7.80), softsusy, spheno or suspect |
| `hcamsb.c` | calculation of the observables for a hypercharged AMSB point, linked to isajet (>=7.80) |
| `mmamsb.c` | calculation of the observables for a mixed moduli AMSB point, linked to isajet (>=7.80) |
| `gmsb.c` | calculation of the observables for a GMSB point, linked to isajet (>=7.80), softsusy, spheno or suspect |
| `nuhm.c` | calculation of the observables for a NUHM point, linked to isajet (>=7.80), softsusy or suspect |
| `thdm.c` | calculation of the observables for a 2HDM type I-IV point, linked to 2HDMC |
| `slha.c` | scan of the user provided SLHA file and calculation of the observables |
| `slha_chi2.c` | calculation of the b->sll observables with chi2 evaluation for a given SLHA file |
| `cnmssm.c` | calculation of the observables for a semi constrained NMSSM point, linked to NMSSMTools |
| `ngmsb.c` | calculation of the observables for a NMSSM GMSB point, linked to NMSSMTools |
| `nnuhm.c` | calculation of the observables for a NMSSM NUHM point, linked to NMSSMTools |
| `modelindep_chi2.c` | calculation of the b->sll observables with chi2 evaluation for a given set of (C9e,C9mu,C10e,C10mu) |

They can be compiled with: `make name` (or `make name.c`), where "name" stands for predict, sm, cmssm, amsb, gmsb, nuhm, ...

They can be run with: `./name.x`

Note that an example SLHA file is also provided: `example.lha`

## Included Files

### Procedures in src/
- **alphas.c**: calculation of the strong coupling constant
- **bdtaunu.c**: calculation of the observables related to B-> D(*) tau nu
- **bkll.c**: calculation of the observables related to B-> K l+ l- with full or soft form factor approaches
- **bkstargamma.c**: calculation of the branching ratio of B-> K* gamma
- **bkstarll.c**: calculation of the observables related to B-> K* l+ l- with full or soft form factor approaches
- **bsgamma.c**: calculation of the inclusive branching ratio of B-> Xs gamma
- **bsll.c**: calculation of the observables related to B-> Xs l+ l-
- **bsll_extra.c**: extra functions for bsll.c
- **bsmumu.c**: calculation of the branching ratios of B(s,d)-> l+ l-
- **bsphill.c**: calculation of the observables related to Bs-> phi l+ l- with full or soft form factor approaches
- **btaunu.c**: calculation of the branching ratio of B-> tau nu
- **chi2.c**: calculation of the chi2
- **dmunu.c**: calculation of the branching ratios of D-> mu nu
- **dslnu.c**: calculation of the branching ratios of Ds-> tau nu and Ds-> mu nu
- **excluded_masses.c**: procedures which check the charge of the LSP, and the agreement with the collider lower mass limits (LEP and Tevatron)
- **fleshouches.c**: generation of FLHA output files
- **general.c**: some useful mathematical functions 
- **gmuon.c**: calculation of the muon anomalous magnetic moment
- **interpreter.c**: interpreter of observable names
- **isajet.c**: procedures calling Isajet
- **isospin.c**: calculation of the isospin asymmetry in B-> K* gamma
- **kaon.c**: calculation of the observables related to leptonic and semileptonic decays
- **kmunu.c**: calculation of the observables related to K-> mu nu
- **lambdablambdall.c**: calculation of the observables related to Lambda_b-> Lambda l+ l-
- **leshouches.c**: procedures reading files in the SLHA format
- **nmssmtools.c**: procedures calling NMSSMTools
- **quark_masses.c**: calculation of the quark masses
- **spheno.c**: procedures calling SPheno
- **softsusy.c**: procedures calling Softsusy
- **suspect.c**: procedures calling SuSpect
- **wilson.c**: calculation of the Wilson coefficients and RGE runnings
- **zwidths.c**: calculation of partial decay withs of the Z boson into SUSY particles
- **2hdmc.c**: procedures calling 2HDMC

### Headers in src
- **include.h**: definitions and prototypes

### Input files in chi2_input
- **myobs.in**: (long) list of b->s observables
- **exp_values.in**: experimental values with uncertainties
- **exp_corr.in**: experimental correlation matrix
- **nuisance.in**: central values and errors of the nuisance parameters
- **nuisance_corr.in**: correlation matrix of the nuisance parameters
- **ckm.in**: list of CKM matrix elements (considered as observables)
- **wilsons.in**: list of Wilson coefficients (considered as observables)
- **obs_kaons.in**: list of kaon observables
- **obs_BKstrmm_Pi_CMS.in**: list of B->K*mumu ang. obs. with CMS binning
- **obs_BKstrmm_Pi_LHCb.in**: list of B->K*mumu ang. obs with LHCb binning
- **mycoeffs.in**: NP Wilson coeffs. values to be used with `predict.c`
- **coeff_list.txt**: list of Wilson coeffs. to be scanned over with `scan_coeffs.c`

### Other
- **example.lha**: example SLHA file
- **Makefile(s)**
- **README.md**

## History

### v5.0 - November 2025
- Updated parameters, form factors, and other numerical values
- Addition of KL->pi ll decay observables
- Implementation of lepton flavor choice for KL->ll
- New routines to facilitate calculations and parameter scanning:
  - `predict.c`: For computing SM or NP predictions with uncertainties (run `./predict.x` for usage instructions)
  - `scan_coeffs.c`: For scanning over Wilson coefficients (run `./scan_coeffs.x` for usage instructions)
- Added several new angular observables for B->K* ll and Bs->phi ll (including CP-asymmetries)
- Multiple form factor options for B->K*, Bs->phi and B->K decays:
  - Default based on LCSR + Lattice results from GKvDV23 (arXiv:2305.06301)
  - Can be configured via variables in `leshouches.c`: `BK_FormFactor_choice`, `Bsphi_FormFactor_choice`, `BKstar_FormFactor_choice`
- Fixed ARFLAG "bug" which was not working on macOS (rcsU -> rcs in case former not found)
- Fixed missing beta_l in B->K*ll calculations (relevant only for decay to tau leptons)
- Fixed conjugation of Wilson coefficients (relevant for studying complex NP effects)
- Added non-physical observables (e.g. form factors) to be outputted

### v4.1 - 11/04/2020
- Improvement of the automatic error calculation
- Implementation of an interpreter to ease calls to observable calculations
- Bug fixed in CQ2 calculation in 2HDM, with no significant change in the results (thanks to S. Heinemeyer)
- Addition of Kaon observables
- Addition of Lambda_b-> Lambda l+l- observables
- Updated parameters, form factors, ...
- Speed optimisations

### v4.0 - 06/30/2018
- Addition of the automatic calculation of the theory errors and correlations (thanks to A. Arbey)
- Addition of routines for computing and manipulating the covariance matrices (thanks to A. Arbey)
- Addition of B->K* gamma branching ratio (thanks to S. Neshatpour)
- Addition of several implementations for B->K*ll (thanks to S. Neshatpour)
- Improved FLHA reader
- Speed optimisations
- Possibility to use OpenMP for the covariance matrix calculation
- Bug fixes

### v3.6 - 04/22/2017
- Addition of the chi2 calculation for b->sll observables, based on Hurth, Mahmoudi, Neshatpour, Nucl. Phys. B909 (2016) 737 [arXiv:1603.00865]
- No change in the results for the observables since v3.5
- Improved implementation of BR(B->D l nu) and addition of other observables
- Addition of B->D* l nu branching ratio, forward backward asymmetry, tau and D* polarizations
- Minor improvements in the structure of the code
- Bug fixes

### v3.5 - 04/19/2016
- Addition of observables related to B->K* ell+ ell-, B->K ell+ ell- and Bs->phi ell+ ell- with full and soft form factor implementations (thanks to contributions from S. Neshatpour)
- Lepton flavour generic routines for most of the leptonic and semileptonic observables
- Addition of four main programs for b->sll observables
- Minor corrections and improvements
- Updated parameter values

### v3.4 - 10/27/2014
- Addition of more observables related to B->K* mu+ mu-
- Addition of three loop 2HDM corrections for b->s gamma
- Addition of three loop QCD corrections and NLO electroweak corrections for BR(Bs -> mu+ mu-)
- Addition of an FLHA input file program
- Minor corrections and improvements
- Updated parameter values (PDG 14, HFAG, FLAG)

### v3.3 - 06/16/2012
- Addition of observables related to B->K* mu+ mu-
- Addition of observables related to B->Xs l+ l-
- Addition of untagged BR(Bs -> mu+ mu-)
- Addition of more Cprime Wilson coefficients
- Minor corrections and improvements
- Updated parameter values

### v3.2 - 10/26/2011
- Addition of BR(Bd -> mu+ mu-)
- Minor corrections and improvements
- Updated parameter values

### v3.1 - 05/23/2011
- Improved interface for Higgsbounds
- Wolfenstein parameterization for the CKM matrix

### v3.0 - 02/23/2011
- NMSSM calculations included
- Interface with NMSSMTools

### v2.8 - 11/02/2010
- FLHA compliant
- Interface with HiggsBounds v2 and HDECAY
- Interface with SPheno
- Interface with SuSpect
- Addition of HC-AMSB and MM-AMSB models
- Improvement of the BR(Bs -> mu+ mu-) calculation
- Updated parameter values

### v2.7 - 11/16/2009
- Compatibility with the new Isajet 7.80
- Implementation of D-> mu nu
- Improvement in the CKM matrix implementation
- Updated parameter values

### v2.6 - 08/31/2009
- Implementation of the general 2HDM
- Interface with 2HDMC
- Addition of a main program for the SM
- Improvement of g-2 calculation

### v2.5 - 03/29/2009
- Replacement of "float" type by "double" to match the relic version conventions
- Minor improvements

### v2.4 - 02/15/2009
- Implementation of the branching ratios of Ds-> tau nu and Ds-> mu nu
- Bug fixed in the SLHA2 reader

### v2.3 - 08/25/2008
- Implementation of the branching ratio of Bs-> mu+ mu-
- Update of the numerical inputs (based on PDG 2008)
- Some minor corrections and improvements
- Change of the name of some routines

### v2.2 - 06/26/2008
- Implementation of B-> tau nu, B-> D tau nu and K-> mu nu observables

### v2.1 - 04/15/2008
- SLHA2 compliant

### v2.0 - 02/18/2008
- Full NNLO implementation for the branching ratio (b -> s gamma)
- Improvement of the calculation of the Wilson coefficients
- Updated parameters for the computation of the isospin asymmetry
- Implementation of the muon g-2 calculation

### v1.1 - 01/14/2008
- NUHM module has been added, as well as some minor corrections

### v1.0 - 09/17/2007
- First public release

## LICENSE

SuperIso Copyright (C) 2007-2020 F. Mahmoudi

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or any 
later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

See <http://www.gnu.org/licenses/>.