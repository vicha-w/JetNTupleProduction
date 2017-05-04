# Hi there
This GitHub repo is forked from the official CMS 2011 OpenData Validation by @tamshai et al, and is modified for my bachelor's thesis. Now my bachelor's thesis is finally finished, I will keep updating this README to better reflect the files in this repo.

All the best,

VW

---

# Simple (?) Jet Tuple production 2011

This project is a CMSSW module producing flat tuples from 2011A Jet data.

Source code was originally forked from the official CMS 2011 OpenData Validation by @tamshai et al:
https://github.com/tamshai/cms-opendata-validation/

The scripts in this repo has been modified in order to support both VM distributed on CERN OpenData portal, available on http://opendata.cern.ch/VM/CMS, as well as lxplus.

## Creating the working area

First setup your own git configuration, or alternatively use these dummy values (required by the command ```git cms-addpkg```): 

```bash
git config --global user.name 'Your Name'
git config --global user.email 'your@ema.il'
git config --global user.github 'username'
```

Next create the working area:
```bash
mkdir WorkingArea # make another directory
cd ./WorkingArea
cmsrel CMSSW_5_3_32 # checkout CMSSW 5_3_32
cd ./CMSSW_5_3_32/src

# initialise CMSSW environment
cmsenv
git cms-addpkg PhysicsTools/PatAlgos

# clone this repo
git clone https://github.com/vicha-w/JetNTupleProduction
cp JetNTupleProduction/jetProducer_cfi.py PhysicsTools/PatAlgos/python/producersLayer1/

# compile everything
scram b

# move to our working directory
cd JetNTupleProduction/AnalysisFW/python/

```

## Setting up additional files

With `JetNTupleProduction/AnalysisFW/python/` as the current folder, run the following commands:

1. Download index files. Here you may download any index files from CERN OpenData portal, such as these two index files : 
    
    ```bash
    # Primary dataset
    wget http://opendata.cern.ch/record/21/files/CMS_Run2011A_Jet_AOD_12Oct2013-v1_20000_file_index.txt

    # MC dataset
    wget http://opendata.cern.ch/record/1562/files/CMS_MonteCarlo2011_Summer11LegDR_QCD_Pt-80to120_TuneZ2_7TeV_pythia6_AODSIM_PU_S13_START53_LV6-v1_00000_file_index.txt 
    ```
    You can choose any index files from OpenData by substituing the above URLs with the URL of your index file of your choice.

2. Download JSON of good runs. This JSON file is required when you need to analyse data from 2011 primary dataset:

    ```bash
    wget http://opendata.cern.ch/record/1001/files/Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON.txt
    ```
    
3. **(Optional for VM)** Create links to the condition databases:

    ```bash
    ln -sf /cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA FT_53_LV5_AN1 
    ln -sf /cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA.db FT_53_LV5_AN1_RUNA.db
    
    ln -sf /cvmfs/cms-opendata-conddb.cern.ch/START53_LV6A1 START53_LV6A1
    ln -sf /cvmfs/cms-opendata-conddb.cern.ch/START53_LV6A1.db START53_LV6A1.db
    ```

4. Breakdown downloaded index files using provided custom script. Each downloaded index file can contain hundreds of links to actual file on the database. In case you want maximum number of events, you should break down index files into many files containing just one link to actual file on the database. The python script just for the task is `indexFragmenter.py`.

    ```bash
    python indexFragmenter.py <YOUR INDEX FILE HERE>
    ```
    
    By default, the script will extract data index file. If you want to extract MC index files, simply change the variable `isData` on line 6 in the script to `False`. The output index file will have a name ending in `<NUMBER>_fragment.txt`. The heading part, `CMS_MonteCarlo2011_Summer11LegDR_` for MC index files and `CMS_Run2011A_` for data index files are also removed.

## Run the program:
To create tuples from data or Mote Carlo simulations, modify the first few lines of `OpenDataCheckout.py` script.

```python
isMC = True # MC or data
customGlobalTag = 'START53_LV6::All'
customIndexFile = 'OpenIndex/DY/DYJetsToLL_M-10To50_TuneZ2_7TeV-pythia6_00000_1_fragment.txt'
customOutFileName = 'OpenTuple/DYJetsToLL_M-10To50_TuneZ2_7TeV-pythia6_00000_1_test_fragment.root'
customBTagDiscrim = 'combinedSecondaryVertexBJetTags'

numberOfEvents = 50
```

* `isMC`: Boolean value. `True` if the following dataset is MC. `False` if the following dataset is primary dataset.
* `customGlobalTag`: String. If unspecified, default GlobalTag is used according to `isMC` value. (`START53_LV6A1::All` if `isMC` is `True`, and `FT_53_LV5_AN1::All` if `isMC` is `False`.)
* `customIndexFile`: String. If unspecified, `CMS_Run2011A_Jet_AOD_12Oct2013-v1_20000_file_index.txt`.
* `customOutFileName`: String. If unspecified, `OpenDataTree_mc.root` or `OpenDataTree_data.root`.
* `customBTagDiscrim`: String - b-tagging discriminator scheme. **Required.**
* `numberOfEvents`: Integer. Enter `-1` to run all events from every file linked from the index file. (`indexFragmenter.py` can help here!)
 
After running the code, you can browse the tuples by initialising `TBrowser`, ROOT's default GUI tree browser, in ROOT:

```bash
root -l
root [] TBrowser bowser
```

## Simple list of tuple variables

* Jets reconstructed using the anti-kT algorithm with a parameter R = 0.5 (short. AK5).

```cpp
int     njet;                     // Number of AK5 jets
float   jet_pt[njet];             // Corrected transverse momentum
float   jet_eta[njet];            // Pseudorapidity
float   jet_phi[njet];            // Azimuthal angle
float   jet_btag[njet];           // b-tagging discriminator value
float   jet_E[njet];              // Energy
```

* Missing energy
```cpp
float   met_et;                   // Transverse energy
float   met_phi;                  // Azimuthal angle
// No pseudorapidity available from missing energy!
```

* Muons from tight muon criteria https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#The2011Data
```cpp
int     nmu;                      // Number of muons
float   muon_pt[njet];            // Corrected transverse momentum
float   muon_eta[njet];           // Pseudorapidity
float   muon_phi[njet];           // Azimuthal angle
float   muon_E[njet];             // Energy
float   muon_charge[njet];        // Charge
```

* Electrons from loose electron criteria (See https://arxiv.org/abs/1208.2671 ref 1 - 5)
```cpp
int     nele;                     // Number of electrons
float   electron_pt[njet];        // Corrected transverse momentum
float   electron_eta[njet];       // Pseudorapidity
float   electron_phi[njet];       // Azimuthal angle
float   electron_E[njet];         // Energy
float   electron_charge[njet];    // Charge
```

* All variables starting with `b_` are values originally available from database before cuts.

## Bash script to run cmsRun through multiple index files

In this repo, two bash scripts are added, named `loopOpenIndex.sh` and `loopOpenIndex.mc.sh`. These two scripts are designed to help you execute cmsRun through multiple fragmented index files (made by `indexFragmenter.py`) easily, without manually modifying `OpenDataCheckout.py` every single time. Before using the scripts, you may store the fragmented index files that came from the same original index file in the same folder.

## Condensing tuple file

In order to condense the tuple file into machine-learning-friendly format, you can "condense" the tuple file by using C++ script `condenseNTuple.cpp` with ROOT. To execute, simply call:

```bash
root -l -q "condenseNTuple.cpp(\"<YOUR TUPLE FILE HERE>\")"
```

Option `-l` suppresses splash screen, and `-q` quits ROOT after finished executing the script. These two options can be ideal when you have lots of tuple files to work with, by simply writing bash script such as the following:

```bash
for file in *.root;
    do root -l -q "condenseNTuple.cpp(\"$file\")";
    done
```

As per usual, you can freely modify every line of code here. Tips and comments are appreciated! ;-)