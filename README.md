# Hi there
This GitHub repo is forked from the official CMS 2011 OpenData Validation by @tamshai et al, and will be further modified for my bachelor's thesis. Now my bachelor's thesis is finally finished, I will keep updating this README to better reflect the files in this repo.

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

```
git config --global user.name 'Your Name'
git config --global user.email 'your@ema.il'
git config --global user.github 'username'
```

Next create the working area:
```
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
    
    ```
    # Primary dataset
    wget http://opendata.cern.ch/record/21/files/CMS_Run2011A_Jet_AOD_12Oct2013-v1_20000_file_index.txt

    # MC dataset
    wget http://opendata.cern.ch/record/1562/files/CMS_MonteCarlo2011_Summer11LegDR_QCD_Pt-80to120_TuneZ2_7TeV_pythia6_AODSIM_PU_S13_START53_LV6-v1_00000_file_index.txt 
    ```
    You can choose any index files from OpenData by substituing the above URLs with the URL of your index file of your choice.

2. Download JSON of good runs. This JSON file is required when you need to analyse data from 2011 primary dataset:

    ```
    wget http://opendata.cern.ch/record/1001/files/Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON.txt
    ```
    
3. **(Optional for VM)** Create links to the condition databases:

    ```
    ln -sf /cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA FT_53_LV5_AN1 
    ln -sf /cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA.db FT_53_LV5_AN1_RUNA.db
    
    ln -sf /cvmfs/cms-opendata-conddb.cern.ch/START53_LV6A1 START53_LV6A1
    ln -sf /cvmfs/cms-opendata-conddb.cern.ch/START53_LV6A1.db START53_LV6A1.db
    ```

4. Breakdown downloaded index files using provided custom script. Each downloaded index file can contain hundreds of links to actual file on the database. In case you want maximum number of events, you should break down index files into many files containing just one link to actual file on the database. The python script just for the task is `indexFragmenter.py`.

    ```
    python indexFragmenter.py <YOUR INDEX FILE HERE>
    ```
    
    The output index file will have a name ending in `<NUMBER>_fragment.txt`.

## Run the program:
To create tuples from data run the following command:

```
    cmsRun OpenDataTreeProducer_dataPAT_2011_cfg.py
```
    
This command creates tuples from Monte Carlo simulations:

```
    cmsRun OpenDataTreeProducer_mcPAT_2011_cfg.py
```

To create tuples from data or Mote Carlo simulations, modify `OpenDataCheckout.py` script.

```
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
 
After running the code, you can browse the tuples by opening the produced files in ROOT:

```
    root OpenDataTree_*
```
 
Finally, run this command in the ROOT command prompt:

```
    TBrowser t
```
 
## Tuple files
More info about the data and Monte Carlo tuples can be found in files ```tuple_info_data``` and  ```tuple_info_mc```.

## Tuple variables

* Properties of the event:

```cpp
    int             run;                // Run number
    float           lumi;               // Luminosity section
    long long       event;              // Event number
    float           ntrg;               // Number of triggers
    bool            triggers[ntrg];     // Trigger bits
    vector<string>  *triggernames;      // Trigger names
    float           prescales[ntrg];    // Trigger prescales
    float           met;                // Missing transverse energy
    float           sumet;              // Sum of transverse energy
    float           rho;                // Energy density
```


* Jets reconstructed using the anti-kT algorithm with a parameter R = 0.5 (short. AK5).

```cpp
    int     njet;           // Number of AK5 jets
    float   jet_pt[njet];   // Corrected transverse momentum
    float   jet_eta[njet];  // Pseudorapidity
    float   jet_phi[njet];  // Azimuthal angle
    float   jet_E[njet];    // Energy
```

* Other AK5 jet information

```cpp
    bool    jet_tightID[njet];  // Tight selection pass/fail
    float   jet_area[njet];     // Jet area in eta-phi plane
    float   jet_jes[njet];      // Jet energy correction
    int     jet_igen[njet];     // Index of the matching generated jet
```

* Composition values of the AK5 jets

```cpp
    float    chf[njet];      // Charged hadron energy fraction
    float    nhf[njet];      // Neutral hadron energy fraction
    float    phf[njet];      // Photon energy fraction
    float    elf[njet];      // Electron energy fraction
    float    muf[njet];      // Muon energy fraction
    float    hf_hf[njet];    // Forward calorimeter (HF) hadron energy fraction
    float    hf_phf[njet];   // HF photon energy fraction
    int      hf_hm[njet];    // HF hadron multiplicity
    int      hf_phm[njet];   // HF photon multiplicity
    int      chm[njet];      // Charged hadron multiplicity
    int      nhm[njet];      // Neutral hadron multiplicity
    int      phm[njet];      // Photon multiplicity
    int      elm[njet];      // Electron multiplicity
    int      mum[njet];      // Muon multiplicity
    float    beta[njet];     // Fraction of chf associated to the hard process
    float    bstar[njet];    // Fraction of chf associated to pile-up
```

* Jets reconstructed using the anti-kT algorithm with a parameter R = 0.7 (short. AK7)

```cpp
    int     njet_ak7;               // Number of jets
    float   jet_pt_ak7[njet_ak7];   // Transverse momentum
    float   jet_eta_ak7[njet_ak7];  // Pseudorapidity
    float   jet_phi_ak7[njet_ak7];  // Azimuthal angle
    float   jet_E_ak7[njet_ak7];    // Energy
    float   jet_area_ak7[njet_ak7]; // Jet area
    float   jet_jes_ak7[njet_ak7];  // Jet energy corection factor
    int     ak7_to_ak5[njet_ak7];   // Index of the corresponding AK5 jet 
```

* True properties of jets generated in the Monte Carlo simulation (only MC datasets)

```cpp
    int     ngen;           // Number of jets generated
    float   gen_pt[ngen];   // Transverse momentum
    float   gen_eta[ngen];  // Pseudorapidity
    float   gen_phi[ngen];  // Azimuthal angle
    float   gen_E[ngen];    // Energy

    float   pthat;          // Transverse momentum in the rest frame of the hard interaction
    float   mcweight;       // Monte Carlo weight of the event
```
