############### MUON TAG ###########################
#./bin/splitGenericJobsOnLxplus_byFiles.py -l data_DoubleMuon -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/muon/ntuple_presel/DoubleMuon/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/mu -f plots -c cfg/analysis_data_DoubleMuon.cfg -n 28 -q workday -s -v --htcondor


############### ELECTRON TAG ###########################
./bin/splitGenericJobsOnLxplus_byFiles.py -l data_DoubleEG -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/ele/ntuple_presel/DoubleEG/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/ele -f plots -c cfg/analysis_data_DoubleEG.cfg -n 27 -q workday -s -v --htcondor



