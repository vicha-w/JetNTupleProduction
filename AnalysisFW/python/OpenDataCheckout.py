
# Forked from 2011 OpenData Validation repo
# https://github.com/cms-smpj/SMPJ/tree/v1.0

## Custom switches

isMC = True # MC or data
customGlobalTag = 'START53_LV6::All'
customIndexFile = 'OpenIndex/DY/DYJetsToLL_M-10To50_TuneZ2_7TeV-pythia6_00000_1_fragment.txt'
customOutFileName = 'OpenTuple/DYJetsToLL_M-10To50_TuneZ2_7TeV-pythia6_00000_1_test_fragment.root'
customBTagDiscrim = 'combinedSecondaryVertexBJetTags'

numberOfEvents = 50

# True : when running in OpenData virtual machine
# False: when runing in lxplus 
runOnVM = False

## Import skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *
import FWCore.Utilities.FileUtils as FileUtils

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Index of data files
if customIndexFile != '' and type(customIndexFile) == str:
    files2011data = FileUtils.loadListFromFile(customIndexFile)
else:
    files2011data = FileUtils.loadListFromFile('CMS_Run2011A_Jet_AOD_12Oct2013-v1_20000_file_index.txt')
process.source.fileNames = cms.untracked.vstring(*files2011data)

# Read condition data from local sqlite database
if runOnVM:
    process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA.db')

if not isMC:
    # Read good luminosity sections from local JSON file
    import FWCore.PythonUtilities.LumiList as LumiList 
    goodJSON = './Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON.txt' 
    myLumis = LumiList.LumiList(filename = goodJSON).getCMSSWString().split(',') 
    process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange() 
    process.source.lumisToProcess.extend(myLumis)

# Specify custom GlobalTag with customGlobalTag variable
# Defaults from cms-opendata-validation:
# Global tag for 2011A data 'FT_53_LV5_AN1::All'
# Global tag for Summer11LegDR-PU_S13_START53_LV6-v1 'START53_LV6A1::All'

if customGlobalTag == '' or type(customGlobalTag) != str:
    if isMC: process.GlobalTag.globaltag = cms.string('START53_LV6A1::All')
    else: process.GlobalTag.globaltag = 'FT_53_LV5_AN1::All'
else: process.GlobalTag.globaltag = cms.string(customGlobalTag)

# Load PAT config
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff") # re-run tau discriminators (new version)
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load('RecoJets.Configuration.RecoJets_cff')
process.load('RecoJets.JetProducers.TrackJetParameters_cfi')
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')


# Import PAT tools
from PhysicsTools.PatAlgos.tools.pfTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.metTools import *
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.tools.muonTools import *
from PhysicsTools.PatAlgos.tools.electronTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector

# Select good vertices
process.goodOfflinePrimaryVertices = cms.EDFilter(
    "VertexSelector",
    filter = cms.bool(False),
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
    )


process.ak5PFJets.doAreaFastjet = True
process.ak7PFJets.doAreaFastjet = True
process.kt6PFJets.doRhoFastjet = True

if not isMC:
    removeMCMatchingPF2PAT( process, '' )
    runOnData(process)

# Choose PF met
addPfMET(process, 'PF')

# Adding non CHS jets to process
addJetCollection(process,cms.InputTag('ak5PFJets'),
                 'AK5', 'PFCorr',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5PF', cms.vstring(['L1FastJet','L2Relative','L3Absolute','L2L3Residual'])),
                 doType1MET   = True,
                 doL1Cleaning = True,
                 doL1Counters = False,
                 doJetID      = True,
                 jetIdLabel   = "ak5"
                 )

addJetCollection(process,cms.InputTag('ak7PFJets'),
                 'AK7', 'PFCorr',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK7PF', cms.vstring(['L1FastJet','L2Relative','L3Absolute','L2L3Residual'])),
                 doType1MET   = True,
                 doL1Cleaning = True,
                 doL1Counters = False,
                 doJetID      = True,
                 jetIdLabel   = "ak7"
                 )

# Tracking failure filter
from RecoMET.METFilters.trackingFailureFilter_cfi import trackingFailureFilter
process.trackingFailureFilter = trackingFailureFilter.clone()
process.trackingFailureFilter.VertexSource = cms.InputTag('goodOfflinePrimaryVertices')


################### EDAnalyzer ##############################
if not isMC:
    process.ak5ak7 = cms.EDAnalyzer('OpenDataTreeProducer',
        ## jet collections ###########################
        pfak7jets       = cms.InputTag('selectedPatJetsAK7PFCorr'),
        pfak5jets       = cms.InputTag('selectedPatJetsAK5PFCorr'),
        ## MET collection ####
        pfmet           = cms.InputTag('pfMET5'),
        ## set the conditions for good Vtx counting ##
        offlineVertices = cms.InputTag('goodOfflinePrimaryVertices'),
        goodVtxNdof     = cms.double(4), 
        goodVtxZ        = cms.double(24),
        ## rho #######################################
        srcPFRho        = cms.InputTag('kt6PFJets','rho'),
        ## preselection cuts #########################
        maxY            = cms.double(5.0), 
        minPFPtJets     = cms.double(30),
        minNPFJets      = cms.int32(2),
        minJJMass       = cms.double(-1),
        isMCarlo        = cms.untracked.bool(False),
        ## trigger ###################################
        printTriggerMenu = cms.untracked.bool(True),
        processName     = cms.string('HLT'),
        triggerNames    = cms.vstring(
                                    'HLT_Jet30', 'HLT_Jet60', 'HLT_Jet80', 'HLT_Jet110', 
                                    'HLT_Jet150', 'HLT_Jet190','HLT_Jet240','HLT_Jet370',
                                    ),
        triggerResults  = cms.InputTag("TriggerResults","","HLT"),
        triggerEvent    = cms.InputTag("hltTriggerSummaryAOD","","HLT"),

        muon            = cms.InputTag('cleanPatMuons'),
        electron        = cms.InputTag('cleanPatElectrons'),
        bTagDiscriminator = cms.string(customBTagDiscrim),

        minPtMuons      = cms.untracked.double(20),
        maxEtaMuons     = cms.untracked.double(2.4),
        globalMuon      = cms.untracked.bool(True),
        trackerMuon     = cms.untracked.bool(True),
        numValidHitsMuon= cms.untracked.int32(10),
        chi2OverNdofMuon= cms.untracked.double(10.),
        muonID          = cms.string('GlobalMuonPromptTight'),
        muonTIP         = cms.untracked.double(0.04),
        RMI             = cms.untracked.double(0.20),

        minPtElectrons  = cms.untracked.double(20),
        maxEtaElectrons = cms.untracked.double(2.5),
        electronID      = cms.string('eidRobustLoose'),
        electronTIP     = cms.untracked.double(0.04),
        REI             = cms.untracked.double(0.17)

    )
else:
    process.ak5ak7 = cms.EDAnalyzer('OpenDataTreeProducer',
        ## jet collections ###########################
        pfak7jets       = cms.InputTag('selectedPatJetsAK7PFCorr'),
        pfak5jets       = cms.InputTag('selectedPatJetsAK5PFCorr'),
        ## MET collection ####
        pfmet           = cms.InputTag('pfMET5'),
        ## set the conditions for good Vtx counting ##
        offlineVertices = cms.InputTag('goodOfflinePrimaryVertices'),
        goodVtxNdof     = cms.double(4), 
        goodVtxZ        = cms.double(24),
        ## rho #######################################
        srcPFRho        = cms.InputTag('kt6PFJets','rho'),
        ## preselection cuts #########################
        maxY            = cms.double(5.0), 
        minPFPtJets     = cms.double(30),
        minNPFJets      = cms.int32(2),
        minGenPt        = cms.untracked.double(30),
        minJJMass       = cms.double(-1),
        isMCarlo        = cms.untracked.bool(True),
        genjets         = cms.untracked.InputTag('ak5GenJets'),
        useGenInfo      = cms.untracked.bool(True),
        ## trigger ###################################
        printTriggerMenu = cms.untracked.bool(True),
        processName     = cms.string('HLT'),
        triggerNames    = cms.vstring(
                                    'HLT_Jet30', 'HLT_Jet60', 'HLT_Jet80', 'HLT_Jet110', 
                                    'HLT_Jet150','HLT_Jet190','HLT_Jet240','HLT_Jet370',
                                    ),
        triggerResults  = cms.InputTag("TriggerResults","","HLT"),
        triggerEvent    = cms.InputTag("hltTriggerSummaryAOD","","HLT"),

        muon            = cms.InputTag('cleanPatMuons'),
        electron        = cms.InputTag('cleanPatElectrons'),
        bTagDiscriminator = cms.string(customBTagDiscrim),

        minPtMuons      = cms.untracked.double(20),
        maxEtaMuons     = cms.untracked.double(2.4),
        globalMuon      = cms.untracked.bool(True),
        trackerMuon     = cms.untracked.bool(True),
        numValidHitsMuon= cms.untracked.int32(10),
        chi2OverNdofMuon= cms.untracked.double(10.),
        muonID          = cms.string('GlobalMuonPromptTight'),
        muonTIP         = cms.untracked.double(0.04),
        RMI             = cms.untracked.double(0.20),

        minPtElectrons  = cms.untracked.double(20),
        maxEtaElectrons = cms.untracked.double(2.5),
        electronID      = cms.string('eidRobustLoose'),
        electronTIP     = cms.untracked.double(0.04),
        REI             = cms.untracked.double(0.17)
    )

# HLT filter
process.hltFilter = cms.EDFilter('HLTHighLevel',
    TriggerResultsTag  = cms.InputTag('TriggerResults','','HLT'),
    HLTPaths           = cms.vstring('HLT_Jet*', 'HLT_DiJetAve*'),
    eventSetupPathsKey = cms.string(''),
    andOr              = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
    throw              = cms.bool(False)
)

# Run everything
process.p = cms.Path(
    process.goodOfflinePrimaryVertices*
    # process.hltFilter *
    process.trackingFailureFilter *
    process.patDefaultSequence *
    process.ak5ak7
)

# Processing time on VM (2011 laptop)
# - DATA: 50000 events / 4 hours
# - MC:   50000 events / 5 hours

# Change number of events here:
process.maxEvents.input = numberOfEvents

process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Output file
if customOutFileName == '' or type(customOutFileName) != str:
    if isMC: process.TFileService = cms.Service("TFileService", fileName = cms.string('OpenDataTree_mc.root'))
    else: process.TFileService = cms.Service("TFileService", fileName = cms.string('OpenDataTree_data.root'))
else: process.TFileService = cms.Service("TFileService", fileName = cms.string(customOutFileName))

# To suppress long output at the end of the job
process.options.wantSummary = False

del process.outpath
