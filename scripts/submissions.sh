############### MUON TAG ###########################
#./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_ggH -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/muon/ntuple_presel/GluGluHToZG_M-125_13TeV_powheg_pythia8_CUETP8M1Up/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/mu -f plots -c cfg/analysis_mc_ggH_mu.cfg -n 1 -q workday -s -v --htcondor

#./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_ttH -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/muon/ntuple_presel/ttHToZG_M125_13TeV_powheg_pythia8_CUETP8M1Up/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/mu -f plots -c cfg/analysis_mc_ttH_mu.cfg -n 1 -q workday -s -v --htcondor

#./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_W+H -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/muon/ntuple_presel/WplusH_HToZG_WToAll_M125_13TeV_powheg_pythia8_CUETP8M1Up/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/mu -f plots -c cfg/analysis_mc_W+H_mu.cfg -n 1 -q workday -s -v --htcondor

#./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_W-H -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/muon/ntuple_presel/WminusH_HToZG_WToAll_M125_13TeV_powheg_pythia8_CUETP8M1Up/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/mu -f plots -c cfg/analysis_mc_W-H_mu.cfg -n 1 -q workday -s -v --htcondor

#./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_VBF -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/muon/ntuple_presel/VBFHToZG_M-125_13TeV_powheg_pythia8_CUETP8M1Up/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/mu -f plots -c cfg/analysis_mc_VBF_mu.cfg -n 1 -q workday -s -v --htcondor

#./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_ZH -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/muon/ntuple_presel/ZH_HToZG_ZToAll_M-125_13TeV_powheg_pythia8_CUETP8M1Up/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/mu -f plots -c cfg/analysis_mc_ZH_mu.cfg -n 1 -q workday -s -v --htcondor

#./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_ZG -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/muon/ntuple_presel/ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/mu -f plots -c cfg/analysis_mc_ZG_mu.cfg -n 11 -q workday -s -v --htcondor

#./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_DY -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/muon/ntuple_presel/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_CUETP8M1Up/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/mu -f plots -c cfg/analysis_mc_DY_mu.cfg -n 1 -q workday -s -v --htcondor


############### ELECTRON TAG ###########################
./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_ggH -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/ele/ntuple_presel/GluGluHToZG_M-125_13TeV_powheg_pythia8_CUETP8M1Up/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/ele -f plots -c cfg/analysis_mc_ggH_ele.cfg -n 1 -q workday -s -v --htcondor

./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_ttH -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/ele/ntuple_presel/ttHToZG_M125_13TeV_powheg_pythia8_CUETP8M1Up/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/ele -f plots -c cfg/analysis_mc_ttH_ele.cfg -n 1 -q workday -s -v --htcondor

./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_W+H -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/ele/ntuple_presel/WplusH_HToZG_WToAll_M125_13TeV_powheg_pythia8_CUETP8M1Up/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/ele -f plots -c cfg/analysis_mc_W+H_ele.cfg -n 1 -q workday -s -v --htcondor

./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_W-H -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/ele/ntuple_presel/WminusH_HToZG_WToAll_M125_13TeV_powheg_pythia8_CUETP8M1Up/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/ele -f plots -c cfg/analysis_mc_W-H_ele.cfg -n 1 -q workday -s -v --htcondor

./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_VBF -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/ele/ntuple_presel/VBFHToZG_M-125_13TeV_powheg_pythia8_CUETP8M1Up/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/ele -f plots -c cfg/analysis_mc_VBF_ele.cfg -n 1 -q workday -s -v --htcondor

./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_ZH -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/ele/ntuple_presel/ZH_HToZG_ZToAll_M-125_13TeV_powheg_pythia8_CUETP8M1Up/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/ele -f plots -c cfg/analysis_mc_ZH_ele.cfg -n 1 -q workday -s -v --htcondor

./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_ZG -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/ele/ntuple_presel/ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/ele -f plots -c cfg/analysis_mc_ZG_ele.cfg -n 11 -q workday -s -v --htcondor

./bin/splitGenericJobsOnLxplus_byFiles.py -l mc_DY -b /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis -e bin/hzg_analysis.exe -i /eos/user/f/fcetorel/HtoZg/ele/ntuple_presel/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_CUETP8M1Up/ -o /afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/ele -f plots -c cfg/analysis_mc_DY_ele.cfg -n 1 -q workday -s -v --htcondor


