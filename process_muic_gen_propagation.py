import ROOT
import array
import numpy as np

ROOT.gErrorIgnoreLevel = ROOT.kFatal

print("===========\nBEGIN PROCESSING MUIC SIMULATION RESULTS\nLAUNCHING ROOT WITH VERSION : " + ROOT.__version__ + "\n===========")
dd4hep_fnames = ["dd4hep_1000events/imcc_1", "dd4hep_1000events/imcc_10", 
                 "dd4hep_1000events/imcc_100", "dd4hep_1000events/imcc_1k", 
                 "dd4hep_4000events/imcc_10k", "dd4hep_4000events/imcc_100k"]  # dd4hep output file names (without .root)
g4_fname = "muicsimulation/muic_ntuple_t0.root"   # propagated gen muons root file
in_csv_fname = "gen_reco_match/all_unmatched.csv"    # source gen muon csv table
out_csv_fname = "prop_all_unmatched.csv"    # source gen muon csv table

h = 12 # use 12 decimals in csv file

make_csv = False
run_dd4hep_merge = True

toroid_etas = [[-7.2,-5.6,660],[-5.6,-4.0,140],[-4.0,-2.4,28]]  # (min eta, max eta, z)
np.random.seed(0)

if make_csv:
    g4_file = ROOT.TFile.Open("in/"+g4_fname)
    g4_tree = g4_file.Get("events")
    num_propagated_events = g4_tree.GetEntries()

    ####################################################
    #  SAVE MUIC SIMULATION OF UNMATCHED MUONS TO CSV  #
    ####################################################

    print(f"Processing file: in/{g4_fname} with {num_propagated_events} entries")

    original_csvfile = open("in/"+in_csv_fname, "r")
    original_csv = original_csvfile.readlines()[1:]
    num_unmatched_events = len(original_csv)

    new_csvfile = open("out/"+out_csv_fname, "w")

    print(f"Processing file: in/{in_csv_fname} with {num_unmatched_events} entries")

    new_csvfile.write("file_idx,event_idx,gen_particle_idx,eta,phi,pt,vx,vy,vz,sum_e,sum_e_minus_pz,Q2,x,y,eta_before_nozzle,eta_after_nozzle,energy_after_nozzle,smeared_e,smeared_pt,hit_nozzle\n")

    print(f"Writing file: out/{out_csv_fname} with {num_unmatched_events} entries")

    for prop_idx in range(num_propagated_events):
        entry = original_csv[prop_idx].strip()
        new_csvfile.write(entry)

        g4_tree.GetEntry(prop_idx)

        if len(g4_tree.gen_e_after_nozzle) == 1:  # gen muon hits nozzle
            energy_after_nozzle = float(g4_tree.gen_e_after_nozzle[0])
            eidx = int(g4_tree.gen_event_idx[0])
            fidx = int(g4_tree.gen_file_idx[0])
            deta = float(g4_tree.gen_pos_eta_scatter_nozzle[0])
            eta_after_nozzle = float(g4_tree.gen_pos_eta_after_nozzle[0])
            eta_before_nozzle = float(g4_tree.gen_pos_eta_before_nozzle[0])
            entry_list = entry.split(",")
            
            assert int(entry_list[0]) == fidx and int(entry_list[1]) == eidx, "Error matching muons from geant 4 to existing data table"

            #############################
            #  ASSIGN TOROID AND SMEAR  #
            #############################

            energy_after_nozzle  # GeV
            hit_z = 0

            # assign toroid
            for toroid in toroid_etas:
                if eta_after_nozzle >= toroid[0] and eta_after_nozzle < toroid[1]:
                    hit_z = toroid[2]
                    break
            
            if hit_z:
                MU_MASS = 0.10565835982561111 # extracted from the dd4hep results GeV

                hit_r = hit_z * (2*np.exp(eta_after_nozzle)) / (1-np.exp(2*eta_after_nozzle))
                unsmeared_e = energy_after_nozzle # unsmeared energy in GeV

                assert unsmeared_e > MU_MASS, "Simulated total energy of muon is less than mass"
                unsmeared_p = np.sqrt(pow(unsmeared_e,2) - pow(MU_MASS,2)) / 1000 # unsmeared momentum in TeV

                p_unc = 0.016 * (hit_r - 1.05) + 0.017  # % uncertainty in momentum in TeV
                smeared_p = np.random.normal(loc=unsmeared_p,scale=p_unc*unsmeared_p) * 1000  # smear with gaussian, std dev = absolute uncertainty, units GeV
                smeared_e = np.sqrt(pow(smeared_p,2) + pow(MU_MASS,2)) # smeared energy in GeV 
                smeared_pt = smeared_p / np.cosh(eta_after_nozzle)
            
            filled_with_geant4 = True

            new_csvfile.write(f",{round(eta_before_nozzle, h)},{round(eta_after_nozzle, h)},{round(energy_after_nozzle, h)},{round(smeared_e, h)},{round(smeared_pt, h)},1\n")

        else: # gen muon does not go thru nozzle, too low eta
            new_csvfile.write(f",0,0,0,0,0,0\n")
        
    original_csvfile.close()
    new_csvfile.close()
    g4_file.Close()

################################################
#  MERGE GEANT4 SIMULATION RESULTS WITH DD4HEP #
################################################

if run_dd4hep_merge:
    print("===========\nBEGIN MERGING MUICSIMULATION RESULTS WITH DD4HEP RESULTS\n===========")
    prop_csv_f = open("out/"+out_csv_fname, "r")
    prop_csv = prop_csv_f.readlines()[1:]
    prop_csv_f.close()

    curr_csvline = 0

    for dd4hep_fname in dd4hep_fnames:
        num_unmatched = 0
        num_hits = 0

        dd4hep_file = ROOT.TFile.Open("in/"+dd4hep_fname+".root")
        output_file = ROOT.TFile("out/"+dd4hep_fname+"_with_g4_info.root", "RECREATE")

        evt_tree = dd4hep_file.Get("evt")  # evt tree will be edited, meta tree will just be copied
        output_tree = evt_tree.CloneTree(0)

        num_entries = evt_tree.GetEntries()

        # did this event have an unmatched muon, that I later propagated thru geant4 simulation ?
        was_unmatched_val = array.array('f', [0.0])
        was_unmatched = output_tree.Branch("gen_prop_was_unmatched", was_unmatched_val, "gen_prop_was_unmatched/F")

        # if the event was unmatched, did it hit the nozzle ?
        hit_nozzle_val = array.array('f', [0.0])
        hit_nozzle = output_tree.Branch("gen_prop_hit_nozzle", hit_nozzle_val, "gen_prop_hit_nozzle/F")

        # sum_e
        sum_e_val = array.array('f', [0.0])
        sum_e = output_tree.Branch("gen_prop_sum_e", sum_e_val, "gen_prop_sum_e/F")

        # sum_e_minus_pz
        sum_e_minus_pz_val = array.array('f', [0.0])
        sum_e_minus_pz = output_tree.Branch("gen_prop_sum_e_minus_pz", sum_e_minus_pz_val, "gen_prop_sum_e_minus_pz/F")

        # eta_before_nozzle
        eta_before_nozzle_val  = array.array('f', [0.0])
        eta_before_nozzle = output_tree.Branch("gen_prop_eta_before_nozzle", eta_before_nozzle_val, "gen_prop_eta_before_nozzle/F")

        # eta_after_nozzle
        eta_after_nozzle_val  = array.array('f', [0.0])
        eta_after_nozzle = output_tree.Branch("gen_prop_eta_after_nozzle", eta_after_nozzle_val, "gen_prop_eta_after_nozzle/F")

        # energy_after_nozzle
        energy_after_nozzle_val  = array.array('f', [0.0])
        energy_after_nozzle = output_tree.Branch("gen_prop_energy_after_nozzle", energy_after_nozzle_val, "gen_prop_energy_after_nozzle/F")

        # smeared_e
        smeared_e_val  = array.array('f', [0.0])
        smeared_e = output_tree.Branch("gen_prop_smeared_e", smeared_e_val, "gen_prop_smeared_e/F")

        # smeared_pt
        smeared_pt_val = array.array('f', [0.0])
        smeared_pt = output_tree.Branch("gen_prop_smeared_pt", smeared_pt_val, "gen_prop_smeared_pt/F")

        for i in range(evt_tree.GetEntries()):  # this "i" is the event_idx
            evt_tree.GetEntry(i)

            if curr_csvline < len(prop_csv):
                curr_csventry = prop_csv[curr_csvline].strip().split(',')
            else: 
                curr_csventry = [0, curr_csvline, 0] # force option 1 (event did not go unmatched, so it was not propagated)

            if int(curr_csventry[1]) > i:  # basically if we missed it
                # if the event idx for this csv is greater than the one currently being processed
                # THEN:    this event did not go unmatched
                was_unmatched_val[0] = 0
                hit_nozzle_val[0] = 0
                sum_e_val[0] = 0
                sum_e_minus_pz_val[0] = 0
                eta_before_nozzle_val[0] = 0
                eta_after_nozzle_val[0] = 0
                energy_after_nozzle_val[0] = 0
                smeared_e_val[0] = 0
                smeared_pt_val[0] = 0

            else:
                was_unmatched_val[0] = 1
                hit_nozzle_val[0] = float(curr_csventry[19])
                sum_e_val[0] = float(curr_csventry[9])
                sum_e_minus_pz_val[0] = float(curr_csventry[10])
                eta_before_nozzle_val[0] = float(curr_csventry[14])
                eta_after_nozzle_val[0] = float(curr_csventry[15])
                energy_after_nozzle_val[0] = float(curr_csventry[16])
                smeared_e_val[0] = float(curr_csventry[17])
                smeared_pt_val[0] = float(curr_csventry[18])
                
                curr_csvline += 1
                num_unmatched += 1
                num_hits += int(curr_csventry[19])

            output_tree.Fill()

        print(f"Writing file: out/{dd4hep_fname}_with_g4_info.root with {num_entries} entries, {num_unmatched} unmatched muons, and {num_hits} hits on nozzle")

        output_tree.Write("",ROOT.TObject.kOverwrite)

        meta_tree = dd4hep_file.Get("meta")
        meta_tree_copy = meta_tree.CloneTree(-1)  # -1 copies all entries
        meta_tree_copy.SetName("meta")

        ROOT.SetOwnership(meta_tree, False)
        ROOT.SetOwnership(meta_tree_copy, False)

        # Write clone to output file
        output_file.cd()
        output_file.WriteTObject(meta_tree_copy, "meta")

        output_file.Close()
        dd4hep_file.Close()