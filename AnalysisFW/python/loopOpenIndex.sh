# Before using this script, remember to initialise CMSSW environment first
# by executing cmsenv within CMSSW folder

for ((x=181; x<=200; x+=1)) do
    for file in DataOpenIndex/DoubleElectron_AOD_12Oct2013-v1_20000_${x}_fragment.txt
    do
        # specify that the index file $file is for data
        sed -i -e "7s@.*@isMC = False # MC or data@" OpenDataCheckout.py

        # set index file
        sed -i -e "9s@'.*'@'${file}'@" OpenDataCheckout.py

        # print out the index file being used
        sed -e "9q;d" OpenDataCheckout.py

        #modify the output file name
        sed -i -e "10s@'.*'@'${file}'@" OpenDataCheckout.py
        sed -i -e "10s@DataOpenIndex@DataOpenTuple@" OpenDataCheckout.py
        sed -i -e "10s@.txt@.root@" OpenDataCheckout.py
        
        # specify the number of events
        sed -i -e "13s@.*@numberOfEvents = -1@" OpenDataCheckout.py

        # print out the output file name
        sed -e "10q;d" OpenDataCheckout.py
        
        # cmsRun
        # Run cmsenv first.
        cmsRun OpenDataCheckout.py
    done
done
