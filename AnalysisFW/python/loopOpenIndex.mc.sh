# Before using this script, remember to initialise CMSSW environment first
# by executing cmsenv within CMSSW folder

for ((x=117; x<=119; x+=1)) do
    for folder in OpenIndex/WJets # YOUR INDEX FOLDER HERE
    do
        for file in $folder/*_${x}_fragment.txt
        do
            # specify that the index file $file is for MC
            sed -i -e "7s@.*@isMC = True # MC or data@" OpenDataCheckout.py

            # set index file
            sed -i -e "9s@'.*'@'${file}'@" OpenDataCheckout.py

            # print out the index file being used
            sed -e "9q;d" OpenDataCheckout.py

            # modify the output file name
            sed -i -e "10s@'.*'@'${file}'@" OpenDataCheckout.py
            sed -i -e "10s@$folder@OpenTuple@" OpenDataCheckout.py
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
done
