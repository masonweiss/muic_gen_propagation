import ROOT
import array
import numpy as np

ROOT.gErrorIgnoreLevel = ROOT.kFatal

print("===========\nBEGIN PROCESSING MUIC SIMULATION RESULTS\nLAUNCHING ROOT WITH VERSION : " + ROOT.__version__ + "\n===========")
# fnames = ["1", "10", "100",
#           "1k", "10k", "100k"]
g4_fname = "muic_ntuple_t0.root"   # propagated gen muons root file
csv_fname = "all_unmatched.csv"    # source gen muon csv table

h = 12 # use 12 decimals in csv file
verbose = False  # print out info 

g4_file = ROOT.TFile.Open("in/"+g4_fname)
g4_tree = g4_file.Get("events")
num_propagated_events = g4_tree.GetEntries()

print(f"Processing file: in/{g4_fname} with {num_propagated_events} entries")

original_csvfile = open("in/"+csv_fname, "r")
original_csv = original_csvfile.readlines()[1:]
num_unmatched_events = len(original_csv)
new_csvfile = open("out/prop_"+csv_fname, "w")

print(f"Processing file: in/{csv_fname} with {num_unmatched_events} entries")

new_csvfile.write("file_idx,event_idx,gen_particle_idx,eta,phi,pt,vx,vy,vz,sum_e,sum_e_minus_pt,Q2,x,y,eta_before_nozzle,eta_after_nozzle,energy_after_nozzle\n")

unmatched_idx = 0
for prop_idx in range(num_propagated_events):
    entry = original_csv[prop_idx].strip()
    new_csvfile.write(entry)

    g4_tree.GetEntry(prop_idx)
    if len(g4_tree.gen_e_after_nozzle) == 1:  # gen muon goes thru nozzle
        energy_after_nozzle = float(g4_tree.gen_e_after_nozzle[0])
        eidx = int(g4_tree.gen_event_idx[0])
        fidx = int(g4_tree.gen_file_idx[0])
        deta = float(g4_tree.gen_pos_eta_scatter_nozzle[0])
        eta_after_nozzle = float(g4_tree.gen_pos_eta_after_nozzle[0])
        eta_before_nozzle = float(g4_tree.gen_pos_eta_before_nozzle[0])
        entry_list = entry.split(",")

        assert int(entry_list[0]) == fidx and int(entry_list[1]) == eidx, "Error matching muons from geant 4 to existing data table"

        new_csvfile.write(f",{round(eta_before_nozzle, h)},{round(eta_after_nozzle, h)},{round(energy_after_nozzle, h)}\n")
    else: # gen muon does not go thru nozzle, too low eta
        new_csvfile.write(f",0,0,0\n")


original_csvfile.close()
