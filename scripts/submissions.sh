############### MUON TAG ###########################
./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_ggH -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/muon/ntuple_newVar/GluGluHToZG_M-125_13TeV_powheg_pythia8_CUETP8M1Up/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/mu/new -f plots -c cfg/analysis_mc_ggH_mu.cfg -n 1 -q workday -s -v --htcondor

#./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_ttH -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/muon/ntuple_newVar/ttHToZG_M125_13TeV_powheg_pythia8/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/mu/new -f plots -c cfg/analysis_mc_ttH_mu.cfg -n 1 -q workday -s -v --htcondor

#./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_WplusH -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/muon/ntuple_newVar/WplusH_HToZG_WToAll_M125_13TeV_powheg_pythia8/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/mu/new -f plots -c cfg/analysis_mc_W+H_mu.cfg -n 1 -q workday -s -v --htcondor

#./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_WminusH -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/muon/ntuple_newVar/WminusH_HToZG_WToAll_M125_13TeV_powheg_pythia8/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/mu/new -f plots -c cfg/analysis_mc_W-H_mu.cfg -n 1 -q workday -s -v --htcondor

#./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_VBF -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/muon/ntuple_newVar/VBFHToZG_M-125_13TeV_powheg_pythia8/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/mu/new -f plots -c cfg/analysis_mc_VBF_mu.cfg -n 1 -q workday -s -v --htcondor

#./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_ZH -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/muon/ntuple_newVar/ZH_HToZG_ZToAll_M-125_13TeV_powheg_pythia8/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/mu/new -f plots -c cfg/analysis_mc_ZH_mu.cfg -n 1 -q workday -s -v --htcondor

#./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_ZG -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i  /eos/user/f/fcetorel/HtoZg/muon/ntuple_newVar/ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/HtoZg_ntuples_v1/190708_152300/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/mu/new -f plots -c cfg/analysis_mc_ZG_mu.cfg -n 11 -q workday -s -v --htcondor

#./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_DY -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/muon/ntuple_newVar/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/HtoZg_ntuples_v1/190708_163732/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/mu/new -f plots -c cfg/analysis_mc_DY_mu.cfg -n 36 -q workday -s -v --htcondor


############### ELECTRON TAG ###########################
#./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_ggH -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/ele/ntuple_newVar/GluGluHToZG_M-125_13TeV_powheg_pythia8_CUETP8M1Up/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/ele/new -f plots -c cfg/analysis_mc_ggH_ele.cfg -n 1 -q workday -s -v --htcondor

#./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_ttH -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/ele/ntuple_newVar/ttHToZG_M125_13TeV_powheg_pythia8/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/ele/new -f plots -c cfg/analysis_mc_ttH_ele.cfg -n 1 -q workday -s -v --htcondor

#./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_WplusH -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/ele/ntuple_newVar/WplusH_HToZG_WToAll_M125_13TeV_powheg_pythia8/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/ele/new -f plots -c cfg/analysis_mc_W+H_ele.cfg -n 1 -q workday -s -v --htcondor

#./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_WminusH -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/ele/ntuple_newVar/WminusH_HToZG_WToAll_M125_13TeV_powheg_pythia8/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/ele/new -f plots -c cfg/analysis_mc_W-H_ele.cfg -n 1 -q workday -s -v --htcondor

#./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_VBF -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/ele/ntuple_newVar/VBFHToZG_M-125_13TeV_powheg_pythia8/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/ele/new -f plots -c cfg/analysis_mc_VBF_ele.cfg -n 1 -q workday -s -v --htcondor

#./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_ZH -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/ele/ntuple_newVar/ZH_HToZG_ZToAll_M-125_13TeV_powheg_pythia8/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/ele/new -f plots -c cfg/analysis_mc_ZH_ele.cfg -n 1 -q workday -s -v --htcondor

##./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_ZG -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/ele/ntuple_newVar/ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/HtoZg_ntuples_v2/190712_075959/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/ele/new -f plots -c cfg/analysis_mc_ZG_ele.cfg -n 11 -q workday -s -v --htcondor

##./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_DY -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/ele/ntuple_newVar/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/ele/new -f plots -c cfg/analysis_mc_DY_ele.cfg -n 36 -q workday -s -v --htcondor


