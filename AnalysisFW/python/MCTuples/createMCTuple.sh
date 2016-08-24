

cwd=$(pwd)

# Environment
echo "CMSSW environment"
cd /afs/cern.ch/user/m/mhaapale/work/public/CMSSW_5_3_32
eval `scramv1 runtime -sh`
cd -

echo "EOS environment"
source /afs/cern.ch/project/eos/installation/cms/etc/setup.sh

# Mount EOS
echo "Mounting EOS"
mkdir -p eos
eosmount eos/

# Top-level directory
cd eos/cms/store/group/phys_smp/mhaapale/MCTuples/

# Executing script
root mergeMC.C

cd $cwd

echo "Success!"
