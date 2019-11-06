'''
Main module to run the regDriver method

@author: Husen M. Umer
'''
from ParseCellInfo import parse_cellinfodict_to_populate_data, populate_cellinfo_dirs
from GetMotifMutScores import score_motifs_according_to_their_affect, file_len, \
        calculate_p_value_motifregions, get_number_of_mutations_per_sample_list_and_write_to_file
from AnnotateMuts import annotate_snvs_with_elements
from multiprocessing import Pool
import sys, os
from pybedtools import BedTool
import pybedtools
from TFExpression import parse_GTEx_metadafile,\
    get_expression_level_per_originType_per_TF,\
    get_expression_level_per_cell_per_TF
from ScoreAnnotatedMuts import getWeightsforAnnotations
from scipy import stats
import numpy as np
import shutil

def usage():
    print("***Usage***\n python regDriver.py [param=value] regDriver_params.conf_file (priority is given to the command line arguments then the input params file)")
    sys.exit(0)

def map_cellNames_to_originTypes(originType_cell_dict_input_file):
    
    tumor_cell_dict = {}
    with open(originType_cell_dict_input_file, 'r') as tumor_cell_dict_infile:
        lines = tumor_cell_dict_infile.readlines()
        for l in lines:
            if l.startswith('***'):
                break
            elif l.startswith('//'):
                continue
            elif '=' in l:
                sl = l.strip().split('=')
                if len(sl)==2:
                    if sl[0] not in tumor_cell_dict:
                        tumor_cell_dict[sl[0]] = sl[1].strip().split(',')
                    else:
                        tumor_cell_dict[sl[0]].extend(sl[1].strip().split(','))
    return tumor_cell_dict

def retreive_cell_elment_datasets(target_director_of_cells_info, list_of_file_tracks_to_add_cell_tracks = [], index_track_name_in_file_tracks=3):
    list_of_cell_tracks = []
    assay_type_cell_tracks_dict = {}
    tf_names = []
    cell_tracks_names_output_file = target_director_of_cells_info.split('/')[-1]+"_cellTracks.txt"
    tf_names_output_file = target_director_of_cells_info.split('/')[-1]+"_tfNames.txt"
    assay_type_cell_tracks_output_file = target_director_of_cells_info.split('/')[-1]+"_assayType_cells.txt"
    if os.path.exists(cell_tracks_names_output_file) and os.path.exists(tf_names_output_file) and os.path.exists(assay_type_cell_tracks_output_file):
        print("Using cell tracks listed in: " + os.path.abspath(cell_tracks_names_output_file))
        with open(cell_tracks_names_output_file, 'r') as cell_tracks_names_outfile:
            cell_tracks_line = cell_tracks_names_outfile.readline()
            list_of_cell_tracks = cell_tracks_line.strip().split(',')
        with open(tf_names_output_file, 'r') as tf_names_outfile:
            tf_names_line = tf_names_outfile.readline()
            tf_names = tf_names_line.strip().split(',')
        with open(assay_type_cell_tracks_output_file, 'r') as assay_type_cell_tracks_outfile:
            assay_type_lines = assay_type_cell_tracks_outfile.readlines()
            for a in assay_type_lines:
                if a.strip().split('\t')[0] not in assay_type_cell_tracks_dict.keys():
                    assay_type_cell_tracks_dict[a.strip().split('\t')[0]] = []
                assay_type_cell_tracks_dict[a.strip().split('\t')[0]].extend(a.strip().split('\t')[1].split(','))
    else:
        cell_dirs = os.listdir(target_director_of_cells_info)
        for cell_dir in cell_dirs:
            if os.path.isdir(target_director_of_cells_info+'/'+cell_dir):
                assay_types = os.listdir(target_director_of_cells_info+'/'+cell_dir)
                for assay_type in assay_types:
                    if assay_type == "DNase-seq":
                        list_of_cell_tracks.append(cell_dir+"#"+assay_type)
                        if assay_type not in assay_type_cell_tracks_dict.keys():
                            assay_type_cell_tracks_dict[assay_type] = []
                        assay_type_cell_tracks_dict[assay_type].append(cell_dir)
                    
                    elif assay_type == "ChromatinStates":
                        assay_type = "ChromHMM"
                        list_of_cell_tracks.append(cell_dir+"#"+assay_type)
                        if assay_type not in assay_type_cell_tracks_dict.keys():
                            assay_type_cell_tracks_dict[assay_type] = []
                        assay_type_cell_tracks_dict[assay_type].append(cell_dir)
                    
                    elif assay_type == "ChIP-seq":
                        factors = os.listdir(target_director_of_cells_info+'/'+cell_dir+'/'+assay_type)
                        if assay_type+'#'+cell_dir not in assay_type_cell_tracks_dict.keys():
                            assay_type_cell_tracks_dict[assay_type+'#'+cell_dir] = []
                        for factor in factors:
                            if ".bed4" in factor:
                                list_of_cell_tracks.append(cell_dir+"#"+factor.rstrip('.bed4'))
                                if "ChIP-seq" not in factor:
                                    tf_names.append(factor.strip('.bed4'))
                                    assay_type_cell_tracks_dict[assay_type+'#'+cell_dir].append(factor.strip('.bed4'))
        for file_track in list_of_file_tracks_to_add_cell_tracks:
            with open(file_track, 'r') as file_track_read:
                line = file_track_read.readline()
                print("reading: " + file_track)
                while line!="":
                    track_names = line.strip().split('\t')[index_track_name_in_file_tracks].split(',')
                    for track_name in track_names:
                        if '#'.join(track_name.strip().split('#')[0:2]) not in list_of_cell_tracks:
                            list_of_cell_tracks.append('#'.join(track_name.strip().split('#')[0:2]))#for tracks that contain cellName#assayType#Value just take cellName#assayType
                        if track_name.strip().split('#')[1] not in assay_type_cell_tracks_dict.keys():
                            assay_type_cell_tracks_dict[track_name.strip().split('#')[1]] = []
                        if track_name.strip().split('#')[0] not in assay_type_cell_tracks_dict[track_name.strip().split('#')[1]]:
                            assay_type_cell_tracks_dict[track_name.strip().split('#')[1]].append(track_name.strip().split('#')[0])
                    line = file_track_read.readline()
                    
        with open(cell_tracks_names_output_file, 'w') as cell_tracks_names_outfile:
            cell_tracks_names_outfile.write(','.join(list_of_cell_tracks))
        with open(tf_names_output_file, 'w') as tf_names_outfile: 
            tf_names_outfile.write(','.join(tf_names))
        with open(assay_type_cell_tracks_output_file, 'w') as assay_type_cell_tracks_outfile: 
            for k in assay_type_cell_tracks_dict.keys():
                assay_type_cell_tracks_outfile.write(k+ '\t' + ','.join(set(assay_type_cell_tracks_dict[k]))+'\n')
        
    return list_of_cell_tracks, tf_names, assay_type_cell_tracks_dict

def retreive_TFFamilyName_for_motifNames(TFFamily_matches_input_file):#TFFamilyName TF_name
    "Retrieves the TF family name for each TF name"
    
    motifTFName_TFNames_matches_dict = {}
    with open(TFFamily_matches_input_file, 'r') as TFFamily_matches_infile:
        lines = TFFamily_matches_infile.readlines()
        for line in lines:
            sl = line.strip().split('\t')
            if len(sl)>1:
                motifTF_name = sl[0].strip().upper()
                if motifTF_name not in motifTFName_TFNames_matches_dict.keys():
                    motifTFName_TFNames_matches_dict[motifTF_name] = []
                for s in sl:
                    if s.strip().upper()!="" and s.strip().upper() not in motifTFName_TFNames_matches_dict[motifTF_name]:
                        motifTFName_TFNames_matches_dict[motifTF_name].append(s.strip().upper())
            
    return motifTFName_TFNames_matches_dict

def annotate_mutations_file(mutations_file, elements_file, motif_sites_dir, motif_PFM_file, all_chromatin_makrs_all_cells_dir, tumor_cells_dict, 
                            tumor_cells_tracks_dict, list_of_cell_tracks, assay_type_cell_tracks_dict, motifTFName_TFNames_matches_dict, mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated_file, 
                            normal_expression_per_cell_per_TF, tumor_expression_per_cell_per_TF, index_sampleID_in_mutations_inputfile=8, remove_temp_files = True,
                            col_names=[], run_training= False, weights_per_param_dict = {}, take_abs_entropy_diff=True, log_base=10, header=True, only_unify_per_sample = False, window_overlap_to_merge = 20, chromatin_states_to_consider=[], use_estimates_from_simulation_set=False, mean_and_sd_from_the_simulations_indiv_sites_dict={}, retrieve_estimates_from_simulation_set = True, bed12_format_bool = True, params={}):
    print("a new process has started")
    mutations_motifs_intersected_output_file =  mutations_file.split('/')[-1] +"_overlappingmotifs"
    mutations_motifs_intersected_output_file_sorted =  mutations_file.split('/')[-1] +"_overlappingmotifssorted"
    no_motif_matching_file_was_found = True
    mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped = mutations_motifs_intersected_output_file_sorted+"_scoredmuts_allchromatin_makrs_all_cells_grouped"
    if not os.path.exists(mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped):
        mutations_file_obj = BedTool(mutations_file)
        mutations_elements = ""
        print("Operating on: " + mutations_file)
        #intersect the mutations with a preset of elements to search from (search domain,e.g a list of pre-difined regulatory regions or promoters) this is to minimize the searching domain and the number of mutations in the next processes
        if elements_file!="none":
            elements_file_obj = BedTool(elements_file)
            mutations_elements = mutations_file_obj.intersect(elements_file_obj, u =True, split=bed12_format_bool) 
        if mutations_elements != "":
            if os.stat(mutations_elements.fn).st_size>2:
                mutations_file_obj = mutations_elements
        #intersect mutations with motifs
        print("Finding overlap with motif sites")
        for motif_sites_file in os.listdir(motif_sites_dir):
            if motif_sites_file.split('.')[0] == mutations_file.split('/')[-1].split('.')[0].split('_')[0]: #to only intersect the mutations file with the motif file of the matching chromosome (chr1==chr1)
                no_motif_matching_file_was_found = False 
                motif_sites_file_obj = BedTool(motif_sites_dir+'/'+motif_sites_file)
                mutations_file_obj.intersect(motif_sites_file_obj, wo=True).saveas(mutations_motifs_intersected_output_file)
                #sort by chr,start,stop,sampleID
                os.system('sort -k1,1 -k2,2n -k3,3n -k' + str(index_sampleID_in_mutations_inputfile+1) + ',' + str(index_sampleID_in_mutations_inputfile+1) + ' '+ mutations_motifs_intersected_output_file + ' > ' + mutations_motifs_intersected_output_file_sorted)
                pybedtools.cleanup()
                #score each mut-motif
                print("Calculating entropy diff for the mutations found at motifs")
                mutations_motifs_scored_file = score_motifs_according_to_their_affect(motif_PFM_file, muts_mutated_motifs_input_file=mutations_motifs_intersected_output_file_sorted, index_mutations_info=0, index_mutated_motif_info=9, index_fromAllele = 3, index_toAllele = 4)
                if remove_temp_files:
                    if os.path.exists(mutations_motifs_intersected_output_file):
                        os.remove(mutations_motifs_intersected_output_file)
                        os.remove(mutations_motifs_intersected_output_file_sorted)
                    else:
                        print("No mutation was found at motif sites: " + mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated_file)
                #intersect then score mut-motifs with all chromatin signals
                mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped = mutations_motifs_scored_file+"_allchromatin_makrs_all_cells_grouped"
                all_chromatin_makrs_all_cells_files = os.listdir(all_chromatin_makrs_all_cells_dir)
                mutations_motifs_scored_file_obj = BedTool(mutations_motifs_scored_file)
                
                #check what chromosome is
                for all_chromatin_makrs_all_cells in all_chromatin_makrs_all_cells_files:
                    if all_chromatin_makrs_all_cells.split('.')[0] == mutations_file.split('/')[-1].split('.')[0].split('_')[0]: #to only intersect the mutations file with the chromain info of the matching chromosome (chr1==chr1)
                        all_chromatin_makrs_all_cells_obj = BedTool(all_chromatin_makrs_all_cells_dir + '/' + all_chromatin_makrs_all_cells)
                        print("Finding and Grouping overlaps between the mut-motifs and all chromatin signals")
                        mutations_motifs_scored_file_obj.intersect(all_chromatin_makrs_all_cells_obj, wo=True)\
                        .groupby(g=range(1, 18), c=21, o = ['distinct'])\
                        .saveas(mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped)
                        if not os.path.exists(mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped):
                            sys.exit("GroupBy error - sort before grouping")
                if remove_temp_files and os.path.exists(mutations_motifs_scored_file): 
                    os.remove(mutations_motifs_scored_file)
                pybedtools.cleanup()
    #annotate the mut-motifs with the chromatin signals
    print("Annotating the mutations")
    #use the expression level of each cell#TF from normal_expression_per_cell_per_TF, tumor_expression_per_cell_per_TF
    mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated = mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped+"_annotated"
    if run_training:
        mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated = mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated_file
    if (not os.path.exists(mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated)) and (not no_motif_matching_file_was_found):
        mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated = annotate_snvs_with_elements(mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped, tumor_cells_dict, tumor_cells_tracks_dict, list_of_cell_tracks, assay_type_cell_tracks_dict, motifTFName_TFNames_matches_dict, mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated, normal_expression_per_cell_per_TF, tumor_expression_per_cell_per_TF, col_names = col_names, run_training = run_training, weights_per_param_dict = weights_per_param_dict, take_abs_entropy_diff=take_abs_entropy_diff, log_base=log_base, header=header)
    if os.path.exists(mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped):
        if remove_temp_files:
            os.remove(mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped)
    if no_motif_matching_file_was_found:
        print("No matching-chr motif file was for: " + mutations_file  + " chromosome...")         
    
    if not run_training and os.path.exists(mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated):
    #unify the same mutations affecting different motifs in the same sample, merge mutations across samples and cancer projects
        if os.stat(mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated).st_size > 2:
            process_annotated_scored_mutations(mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated, mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated_file, only_unify_per_sample = False, window_overlap_to_merge = 20, header=header, elements_file=elements_file, chromatin_states_to_consider = chromatin_states_to_consider, use_estimates_from_simulation_set=use_estimates_from_simulation_set, mean_and_sd_from_the_simulations_indiv_sites_dict=mean_and_sd_from_the_simulations_indiv_sites_dict, retrieve_estimates_from_simulation_set = retrieve_estimates_from_simulation_set, bed12_format_bool = bed12_format_bool, params=params)
            if retrieve_estimates_from_simulation_set:
                shutil.move(mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated, mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated_file)
        if remove_temp_files and os.path.exists(mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated):
            os.remove(mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated)
    return mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated_file, mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated_file+"_withoutfilterunifiedpersample"

def get_col_names():
    col_names =  ['MutChr', 'MutStart', 'MutStop', 'Ref', 'Alt', 'Cancer-project', 'Mut_score', 'vcf_info', 'sample_file_id', 'MotifChr', 'MotifStart', 'MotifEnd', 'Motif_name', 'Motif_score', 'Motif_strand', 'Entropy_diff', 'MutMotif_position', 'normal_expression_level_of_this_motif', 'tumor_expression_level_of_this_motif',]
    #sep.join(cols_to_use_cobminedcells), #same_factor_overlaps +  sep + all_factors_same_cell_overlaps + sep + dnase_overlap + sep + contactingdomain_overlap + sep + loopdomain_overlap + sep + chromhmm_overlap_to_report + sep + replicationtiming_label_overlap_to_report + sep + cage_overlap
    ''''same_factor_overlaps', 'all_factors_same_cell_overlaps', 'dnase_overlap', 'contactingdomain_overlap', 'loopdomain_overlap', 'cage_overlap', 'chromhmm_overlap_to_report', 'replicationtiming_label_overlap_to_report',
    'same_cell_same_factor', 'same_cell_dnase', 'same_cell_chromatin_states', 'same_cell_replicationdomain', 'same_cell_contactingdomain', 'same_cell_loopdomain', 
    #Same cell Other factors
    'len_same_cell_other_factors_names', 'len_same_cell_other_factors_per', 'same_cell_factors_names', 'len_all_same_cell_other_factors',
    #the same factor in other cells
    'same_factor_overlaps_othercells_per', 'len_other_cells_same_factor', 'other_cells_same_factor', 'len_all_other_cells_same_factor', 
    #DNase1
    'dnase_overlap_othercells_per', 'len_other_cells_dnase', 'other_cells_dnase', 'len_alls_dnase_other_cells',
    #Chromatin states
    'len_other_cells_chromatinstates', 'other_cells_chromatinstates', 'len_all_chromatinstates_other_cells',
    #Replication domains
    'len_other_cells_replicationdomains', 'other_cells_replicationdomains', 'len_all_replicationdomains_other_cells',
    #Contacting domains
    'contactingdomain_overlap_othercells_per', 'len_other_cells_contactingdomains', 'other_cells_contactingdomains', 'len_all_contactingdomains_other_cells', 
    #Loop domains
    'loopdomain_overlap_othercells_per', 'len_other_cells_loopdomains', 'other_cells_loopdomains', 'len_all_loopdomains_other_cells',
    'CAGE_peak_expr_same_cell', 'CAGE_peak_expr_other_cells', 'WAScore']'''

    col_names.extend(["same_factor", "other_factors", "dnase", "contactingdomain", "loopdomain", "CAGE_expr", "chromatin_states", "replicationdomain"])
    col_names.extend(["same_cell_same_factor","same_cell_dnase","same_cell_chromatin_states","same_cell_replicationdomain","same_cell_contactingdomain","same_cell_loopdomain","same_cell_CAGE_expr","same_cell_other_factors_names","same_cell_other_factors", "other_cells_same_factor","other_cells_dnase","other_cell_chromatin_states","other_cells_replicationdomain","other_cells_contactingdomain","other_cells_loopdomain","other_cells_CAGE_expr"])
    col_names.append('WAScore')
    
    return col_names

def rank_scored_muts(tumor_type_file, tumor_type_file_ranked, index_score=61):
    #if the input was a file, split by tumor type, otherwise start processing
    tumor_type_file_readfile = open(tumor_type_file, 'r')
    line = tumor_type_file_readfile.readline()
    scores = []
    print("Reading scores for ranking")
    while line!="":
        scores.append(float(line.strip().split('\t')[index_score]))
        line = tumor_type_file_readfile.readline()
    tumor_type_file_readfile.close()
    
    #rank the scores
    ranks = (stats.rankdata(scores, "average")/len(scores)).tolist()
    
    print("Writing the ranks")
    tumor_type_file_readfile = open(tumor_type_file, 'r')
    line = tumor_type_file_readfile.readline().strip()
    tumor_type_file_ranked_write_file = open(tumor_type_file_ranked, 'w')
    line_number = 0
    while line!="":
        tumor_type_file_ranked_write_file.write(line + "\t" + str(ranks[line_number]) + '\n')
        line = tumor_type_file_readfile.readline().strip()
        line_number +=1
    tumor_type_file_readfile.close()
    tumor_type_file_ranked_write_file.close()

    return tumor_type_file_ranked

def process_annotated_scored_mutations(annotated_mutations_final_output_file_ranked, annotated_mutations_final_output_file_scored_merged, only_unify_per_sample = False, window_overlap_to_merge = 20, header=False, elements_file='none', chromatin_states_to_consider = [], use_estimates_from_simulation_set=True, mean_and_sd_from_the_simulations_indiv_sites_dict={}, retrieve_estimates_from_simulation_set = True, bed12_format_bool = True, params={}):
    
    skip_line = ""
    awk_cond_st=""
    if len(chromatin_states_to_consider)>0 and chromatin_states_to_consider[0]!='' and chromatin_states_to_consider[0]!='all':
        states_to_consider_cond = []
        for state in chromatin_states_to_consider:
            states_to_consider_cond.append("$26~ \""+ state + "\"")
        awk_cond_st = "if(" + ' || '.join(states_to_consider_cond) + ")"#{$44=$44+1.244;} "
    if header:
        skip_line = "if(NR>1) "#skip the header line from the annotated file if there exists one 
        
    #get overall score values from a simulation test without any filtering
    if retrieve_estimates_from_simulation_set:
        if not os.path.exists(annotated_mutations_final_output_file_scored_merged+"_withoutfilterunifiedpersample"):
            awk_stm = """awk 'BEGIN{FS=OFS="\t"}{""" + skip_line + awk_cond_st + """print $1,$2,$3,$4">"$5,$6,$9,$10":"$11"-"$12"#"$13"#"$15"#"$16"#"$17"#"$44,$26,$27,$44}' """ + annotated_mutations_final_output_file_ranked
            awk_stm += " | sort -k1,1 -k2,2n -k3,3n -k6,6 -k5,5 " + """ | groupBy -g 1,2,3,6 -c 4,5,7,8,9,10 -o distinct,distinct,distinct,distinct,distinct,mean | awk 'BEGIN{FS=OFS="\t"}{gsub(",","|"); print}' | sort -k1,1 -k2,2n -k3,3n """ 
            os.system(awk_stm + " > " + annotated_mutations_final_output_file_scored_merged+"_withoutfilterunifiedpersample")
            
        #mean_and_sd_from_the_simulations_indiv_sites_dict = get_mean_and_sd_from_simulation_sets(annotated_mutations_final_output_file_scored_merged+"_withoutfilterunifiedpersample", params['mean_and_sd_output_file']+"_withoutfilterunifiedpersample", score_index_simulation_file=int(params['score_index_simulation_file']))
    #get the mean of the score dist from the given simulation if there exists one
    else:
        #use mutations_elements to combine muts when an elements file is given instead of the d window_overlap_to_merge
        #unify muts per sample (combine a single mut afftecting several motifs in the same sample; groupBy mut_pos and sample_id)
        #min_indiv_score_to_consider = mean_and_sd_from_the_simulations_indiv_sites_dict['mean']#+(1.96*mean_and_sd_from_the_simulations_indiv_sites_dict['std']) 
        #if min_indiv_score_to_consider ==0:
        min_indiv_score_to_consider = 0#1.47 #:avg of single motifs genomewide (sanger sim) (FIX: make this dynamic)
        if 'mean' in mean_and_sd_from_the_simulations_indiv_sites_dict.keys():
            if mean_and_sd_from_the_simulations_indiv_sites_dict['mean'] >= 0 and mean_and_sd_from_the_simulations_indiv_sites_dict['mean']!='':
                min_indiv_score_to_consider = mean_and_sd_from_the_simulations_indiv_sites_dict['mean']
        
        #use the avg score dist from the given simulation to filterout muts with a lower score
        if not use_estimates_from_simulation_set:
            min_indiv_score_to_consider = 0
            
        awk_stm = """awk 'BEGIN{FS=OFS="\t"}{""" + skip_line + awk_cond_st + """print $1,$2,$3,$4">"$5,$6,$9,$10":"$11"-"$12"#"$13"#"$15"#"$16"#"$17"#"$44,$26,$27,$44}' """ + annotated_mutations_final_output_file_ranked
        awk_stm += " | sort -k1,1 -k2,2n -k3,3n -k6,6 -k5,5 " + """ | groupBy -g 1,2,3,6 -c 4,5,7,8,9,10 -o distinct,distinct,distinct,distinct,distinct,mean | awk 'BEGIN{FS=OFS="\t"}{gsub(",","|"); if($10>""" + str(min_indiv_score_to_consider) + """) print}' | sort -k1,1 -k2,2n -k3,3n """ 
        os.system(awk_stm + " > " + annotated_mutations_final_output_file_scored_merged+"_unifiedpersample")
        awk_stm=""
        elements_muts_temp = '.'.join(elements_file.split('/')[-1].split('.')[0:-1])+'_'+annotated_mutations_final_output_file_scored_merged+"_unifiedpersample_temp"
        if not only_unify_per_sample:
            if elements_file!="none" and os.path.exists(elements_file):
                if os.path.exists(annotated_mutations_final_output_file_scored_merged+"_unifiedpersample"):
                    if os.stat(annotated_mutations_final_output_file_scored_merged+"_unifiedpersample").st_size > 2:
                        number_cols_in_elements_file = 4
                        with open(elements_file, 'r') as elements_file_read:#get number of lines from the first line of the elements file; only the first four cols (pos, name) are used here 
                            number_cols_in_elements_file = len(elements_file_read.readline().strip().split('\t'))
                        
                        muts_unified_per_sample_obj = BedTool(annotated_mutations_final_output_file_scored_merged+"_unifiedpersample")
                        elements_file_obj = BedTool(elements_file)
                        elements_file_obj.intersect(muts_unified_per_sample_obj, wo =True, split=bed12_format_bool).saveas(elements_muts_temp)
                        #combine muts cross samples and blocks of the given element (sum of scores)
                        awk_stm = ("cut -f1-4," + str(number_cols_in_elements_file+1) + "- " + elements_muts_temp + """ | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$8,$9,$10,$5":"$6"-"$7"@"$8"@"$9"@"$10"@"$11"@"$12"@"$13"@"$14,$12,$13,$14}' | sort -k1,1 -k2,2n -k3,3n -k4,4 | groupBy -g 1,2,3,4 -c 5,5,6,7,8,9,10,11 -o collapse,count,collapse,distinct,collapse,distinct,distinct,sum""") 
                        #Input: element_chr, element_start,element_end, elementID, mut_chr, start, end, sampleID, refalt, tumorType, info, chromatinDomain, repliDomain, score
                    else:
                        print('no muts have a score larger than: ' + str(min_indiv_score_to_consider))
                    #intermediate (prepared for groupBy): element_chr, element_start,element_end, elementID, sampleID, tumorType, chromatin_label, repliLabel, score, mut_chr:start-end@sampleID@refalt@tumorType@info@chromatinDomain@repliDomain@score
                #final statement: element_chr, element_start,element_end, elementID, sampleID, sample_count, refToalt, tumorTypes, chromatin_label, repliLabel, sum_score, mut_chr:start-end@sampleID@refalt@tumorType@info@chromatinDomain@repliDomain@score
            else:
                #final statement: m_chr, m_start, m_end, mID, sampleID, sample_count, refToalt, tumorTypes, motifs_info, chromatin_label, repliLabel, scores, sum_score
                awk_stm = """mergeBed -i """ + annotated_mutations_final_output_file_scored_merged+"_unifiedpersample" + """ -d """ + str(window_overlap_to_merge) + """ -c 4,4,5,6,7,8,9,10,10 -o collapse,count,collapse,collapse,collapse,distinct,distinct,collapse,sum | awk 'BEGIN{FS=OFS="\t"}{split($6,nucl_split,">"); if(length(nucl_split[1])==length(nucl_split[2]) && ($3-$2+1 > length(nucl_split[1])) && nucl_split[1]!="-" && nucl_split[2]!="-"){$2=$2+1; $3=$3-1;} print $1,$2,$3,$1":"$2"-"$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}'""" #mergeBed adds one base to the start pos and one to the end, so they have to be extracted again 
        
        #print awk_stm + " > " + annotated_mutations_final_output_file_scored_merged
        os.system(awk_stm + " > " + annotated_mutations_final_output_file_scored_merged)
        if os.path.exists(elements_muts_temp):
            os.remove(elements_muts_temp)
        if os.path.exists(annotated_mutations_final_output_file_scored_merged+"_unifiedpersample"):
            os.remove(annotated_mutations_final_output_file_scored_merged+"_unifiedpersample")
    return annotated_mutations_final_output_file_scored_merged

def process_annotated_scored_mutations_rankedMutsbyPercentile(annotated_mutations_final_output_file_ranked, annotated_mutations_final_output_file_scored_merged, window_overlap_to_merge = 20):
    col_names_to_report = ['MutChr', 'MutStart', 'MutStop', 'Ref', 'Alt', 'Cancer-project', 'Mut_score', 'vcf_info', 'sample_file_id', 
                             'MotifChr', 'MotifStart', 'MotifEnd', 'Motif_name', 'Motif_score', 'Motif_strand', 'Entropy_diff', 'MutMotif_position',                        
                             'normal_expression_level_of_this_motif', 'tumor_expression_level_of_this_motif',
                             'same_factor_overlaps', 'all_factors_same_cell_overlaps', 'dnase_overlap', 'contactingdomain_overlap', 'loopdomain_overlap', 'cage_overlap', 
                             'chromhmm_overlap_to_report', 'replicationtiming_label_overlap_to_report','WAScore', 'rank']
    col_numbers_to_report = []
    col_number_for_filtering = 63
    filter_value = 0.95
    for c in col_names_to_report:
        if c in col_names:
            #get col number for awk
            col_numbers_to_report.append('$' + str(col_names.index(c)+1))
        else:
            del col_names_to_report[col_names_to_report.index(c)]
    
    col_names_to_sortBy = ['MutChr', 'MutStart', 'MutStop', 'Cancer-project', 'Mut_score', 'vcf_info', 'sample_file_id'] 
    col_numbers_to_sortBy = []
    
    numeric_columns_sortBy = ['MutStart', 'MutStop']
    for c in col_names_to_sortBy:
        if c in col_names:
            #get col number for awk
            if c in numeric_columns_sortBy:
                col_numbers_to_sortBy.append('-k' + str(col_names.index(c)+1) + ',' + str(col_names.index(c)+1)+'n')
            else:
                col_numbers_to_sortBy.append('-k' + str(col_names.index(c)+1) + ',' + str(col_names.index(c)+1))
        else:
            del col_names_to_sortBy[col_names_to_sortBy.index(c)]
    
    awk_stm = ("""awk 'BEGIN{FS=OFS="\t"}{if($""" + str(col_number_for_filtering) + ' >= ' + str(filter_value) + """){gsub(",","+"); print """ + ','.join(col_numbers_to_report[0:3]) + "," + col_numbers_to_report[8] + ',' + '"@"'.join(col_numbers_to_report[0:9]) + "," + col_numbers_to_report[12] + ',' + col_numbers_to_report[-2] + ',' + col_numbers_to_report[-1] + ',' + '"@"'.join(col_numbers_to_report[9::]) + """}}' """ + annotated_mutations_final_output_file_ranked +
               " | sort -k1,1 -k2,2n -k3,3n -k4,4 " +  
              """| groupBy -g 1,2,3,4 -c 5,6,7,8,9 -o distinct,collapse,max,max,collapse """ + #report each mutaton only once 
              " | sort -k1,1 -k2,2n -k3,3n | " + #merge muts withing a window 
              " mergeBed -d " + str(window_overlap_to_merge) + " -c 4,4,5,6,7,8,9 -o count,distinct,collapse,collapse,sum,sum,collapse | sort -k4,4nr | " +
              """ awk 'BEGIN{FS=OFS="\t"}{if($2!=$3){$2=$2+1; $3=$3-1;} print}' > """ + #mergeBed adds one base to the start pos and one to the end, so they have to be extracted again 
              annotated_mutations_final_output_file_scored_merged) 
    os.system(awk_stm)        

def get_mean_and_sd_from_simulation_sets(scored_simulation_input_file, mean_and_sd_output_file, score_index_simulation_file=-1):
    
    mean_std_median_dict = {}
    if os.path.exists(mean_and_sd_output_file):
        print("reading mean and sd from an existsing fiel: " + mean_and_sd_output_file)
        with open(mean_and_sd_output_file, 'r') as mean_and_sd_read:
            lines = mean_and_sd_read.readlines()
            for l in lines:
                if len(l.split('\t'))==2:
                    mean_std_median_dict[l.strip().split('\t')[0]] = float(l.strip().split('\t')[1])
    else:
        scores_from_simulations = []
        with open(scored_simulation_input_file, 'r') as scored_simulation_read:
            l = scored_simulation_read.readline()
            while l!="":
                scores_from_simulations.append(float(l.split('\t')[score_index_simulation_file].strip()))
                l = scored_simulation_read.readline()
        mean_std_median_dict['mean'] = np.mean(scores_from_simulations)
        mean_std_median_dict['std'] = np.std(scores_from_simulations)
        mean_std_median_dict['median'] = np.median(scores_from_simulations)
        if mean_and_sd_output_file!='none':
            with open(mean_and_sd_output_file, 'w') as mean_and_sd_write:
                for k in mean_std_median_dict.keys():
                    mean_and_sd_write.write(k + '\t' + str(mean_std_median_dict[k]) + '\n')
    return mean_std_median_dict

def get_pvalues_per_element_score(annotated_mutations_statcalc_output_file, annotated_mutations_statcalc_output_file_scores_significance, mean_and_sd_from_the_simulations, element_score_index_obs_file=-3):
    
    if 'mean' in mean_and_sd_from_the_simulations.keys() and 'std' in mean_and_sd_from_the_simulations.keys():
        with open(annotated_mutations_statcalc_output_file, 'r') as annotated_mutations_statcalc_read:
            with open(annotated_mutations_statcalc_output_file_scores_significance, 'w') as annotated_mutations_statcalc_output_file_scores_significance_write:
                read_line = annotated_mutations_statcalc_read.readline()
                while read_line != "":
                    score_current_line  = float(read_line.split('\t')[element_score_index_obs_file].strip())
                    z_score = (score_current_line - mean_and_sd_from_the_simulations['mean'])/mean_and_sd_from_the_simulations['std']
                    p_value = stats.norm.sf(z_score)
                    annotated_mutations_statcalc_output_file_scores_significance_write.write(read_line.strip() + '\t' + str(p_value) + '\n')
                    read_line = annotated_mutations_statcalc_read.readline()
    else:
        print("mean or std is not provided in the given file")
    return annotated_mutations_statcalc_output_file_scores_significance
            
def get_matching_cells_tracks_per_tumor(list_of_cell_tracks, tumor_cells_dict):
    tumor_cells_tracks_dict = {}
    for tumor in tumor_cells_dict:
        for cell_track in list_of_cell_tracks:
            if cell_track.split('#')[0] in tumor_cells_dict[tumor]:
                if tumor not in tumor_cells_tracks_dict:
                    tumor_cells_tracks_dict[tumor] = []
                tumor_cells_tracks_dict[tumor].append(cell_track)
    
    return tumor_cells_tracks_dict

def get_params(params_list):
    
    params = {}
    for arg in params_list:#priority is for the command line
        if '=' in arg: 
            if len(arg.strip().split('='))==2:
                if arg.split('=')[0] not in params.keys():
                    params[arg.strip().split('=')[0]] = arg.strip().split('=')[1]
    if 'param_file' in params:     
        with open(params['param_file'], 'r') as params_infile:
            params_from_file = params_infile.readlines()
            for line in params_from_file:
                if not line.startswith('//') and not line.startswith('#') and '=' in line:
                    if len(line.strip().split('='))==2:
                        if line.strip().split('=')[0] not in params.keys():
                            params[line.strip().split('=')[0]] = line.strip().split('=')[1].split('#')[0].split('//')[0]
    return params

def get_value(str):
    if 'true' in str.lower() or 'yes' in str.lower():
        return True
    else:
        return False

if __name__ == '__main__':
    
    params = get_params(sys.argv)
    if len(params.keys())==0:
        usage()
    print("Given params", params)
        
    run_training = get_value(params['run_training_arg'])    
    take_abs_entropy_diff = get_value(params['take_abs_entropy_diff_arg'])
    use_gene_expression_for_scoring = get_value(params['use_gene_expression_for_scoring_arg'])
    header = get_value(params['header_param'])
    rank_scores = get_value(params['rank_scores_param'])
    compute_significance = get_value(params['compute_significance_param'])
    compute_score_sig = get_value(params['compute_score_sig_param'])
    run_in_parallele = get_value(params['run_in_parallele_param'])
    remove_temp_files = get_value(params['remove_temp_files'])
    use_estimates_from_simulation_set = get_value(params['use_estimates_from_simulation_set'])
    only_unify_per_sample = False
    retrieve_estimates_from_simulation_set = get_value(params['retrieve_estimates_from_simulation_set'])
    bed12_format_bool = get_value(params['bed12_format_bool'])
    
    
    col_names = get_col_names()
    #if not run_training:
    #    col_names.append("final_mut_motif_tracks_score")
    weights_per_param_dict = {}
    if not run_training:
        if os.path.exists(params['weights_per_param_dict_arg_file']):
            with open(params['weights_per_param_dict_arg_file'], 'r') as weights_per_param_dict_arg_infile:
                param_weights_from_input_file  = weights_per_param_dict_arg_infile.readline().strip().split(',')
                for element in param_weights_from_input_file:
                    if len(element.split('='))==2:
                        if element.split('=')[0] not in weights_per_param_dict.keys():
                            weights_per_param_dict[element.strip().split('=')[0]] = float(element.strip().split('=')[1]) 
    
    list_of_file_tracks_to_add_cell_tracks = []
    list_of_file_tracks_to_add_cell_tracks.append(params['ContactingDomains_file_path'])
    list_of_file_tracks_to_add_cell_tracks.append(params['LoopDomains_file_path'])
    list_of_file_tracks_to_add_cell_tracks.append(params['ReplicationDomains_file_path'])
    list_of_file_tracks_to_add_cell_tracks.append(params['CAGE_expr_file_path'])
    
    Mutations_dir_list = []#if a directory is given then list its content otherwise put the file in the list
    if os.path.isdir(params['mutations_dir']):
        for f in os.listdir(params['mutations_dir']):
            Mutations_dir_list.append(params['mutations_dir']+"/"+f)
    elif os.path.isfile(params['mutations_dir']):#in such case mutations_dir is not a dir but a file
        Mutations_dir_list = [params['mutations_dir']]
    
    mean_and_sd_from_the_simulations_indiv_sites_dict = {'mean': 0, 'std': 0, 'median': 0}
    
    if not os.path.exists(params['annotated_mutations_final_output_file']):
        if use_estimates_from_simulation_set and not retrieve_estimates_from_simulation_set:
            if os.path.exists(params['scored_simulation_input_file']+"_withoutfilterunifiedpersample") and params['mean_and_sd_output_file']!='none':
                print('Getting mean and Std from: ' + params['scored_simulation_input_file']+"_withoutfilterunifiedpersample")
                mean_and_sd_from_the_simulations_indiv_sites_dict = get_mean_and_sd_from_simulation_sets(params['scored_simulation_input_file']+"_withoutfilterunifiedpersample", params['mean_and_sd_output_file']+"_withoutfilterunifiedpersample", score_index_simulation_file=int(params['score_index_simulation_file']))
            else: 
                print("using the default values for mean and sd (0,0) for single sites beacuse the score_per_site file doesn't exist")
        list_of_cell_tracks, list_tf_names_from_tracks, assay_type_cell_tracks_dict  = retreive_cell_elment_datasets(params['CellInfo_target_dir'], list_of_file_tracks_to_add_cell_tracks)
        tumor_cells_dict = map_cellNames_to_originTypes(params['TumorCellInfo_matches_dict'])
        tissue_cells_dict = map_cellNames_to_originTypes(params['TissueCellInfo_matches_dict'])
        motifTFName_TFNames_matches_dict = retreive_TFFamilyName_for_motifNames(params['TF_family_matches_file'])
        tumor_cells_tracks_dict = get_matching_cells_tracks_per_tumor(list_of_cell_tracks, tumor_cells_dict)
        #use list_of_cell_tracks to gather expression level of each TF in each cell
        #gene expression values from the normal samples input file 
        
        normal_expression_per_cell_per_TF = {}
        tf_names_to_extract_gene_expression_for = []#list_tf_names_from_tracks#get names of TFs from the TFFamily file and the dirs contaning ChIP-seq datasets
        for k in motifTFName_TFNames_matches_dict.keys():
            tf_names_to_extract_gene_expression_for.append(k)
        tf_names_to_extract_gene_expression_for = list(set(tf_names_to_extract_gene_expression_for))
        
        if os.path.exists(params['normal_gene_expression_inputfile']) and os.path.exists(params['metadata_samples_normal_gene_expression_file']):
            dict_tissue_type_sampleIDs = parse_GTEx_metadafile(params['metadata_samples_normal_gene_expression_file'])
            origin_gene_expression_values = get_expression_level_per_originType_per_TF(tf_names_to_extract_gene_expression_for, dict_tissue_type_sampleIDs, params['normal_gene_expression_inputfile'])
            print("Getting expression per cell#TF in normal samples")
            origin_gene_expression_values_outputfile = params['normal_gene_expression_inputfile'] + "_perCell_perTF"
            normal_expression_per_cell_per_TF = get_expression_level_per_cell_per_TF(tf_names_to_extract_gene_expression_for, origin_gene_expression_values, tissue_cells_dict, origin_gene_expression_values_outputfile)
        
        #gene expression values from the tumor samples input file 
        tumor_expression_per_cell_per_TF = {}
        if os.path.exists(params['tumor_gene_expression_inputfile']) and os.path.exists(params['metadata_samples_tumor_gene_expression_file']):
            dict_tissue_type_sampleIDs = parse_GTEx_metadafile(params['metadata_samples_tumor_gene_expression_file'])
            origin_gene_expression_values = get_expression_level_per_originType_per_TF(tf_names_to_extract_gene_expression_for, dict_tissue_type_sampleIDs, params['tumor_gene_expression_inputfile'])
            print("Getting expression per cell#TF in tumor samples")
            origin_gene_expression_values_outputfile = params['normal_gene_expression_inputfile'] + "_perCell_perTF"
            tumor_expression_per_cell_per_TF = get_expression_level_per_cell_per_TF(tf_names_to_extract_gene_expression_for, origin_gene_expression_values, tumor_cells_dict, origin_gene_expression_values_outputfile)
        
        mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated_files = []
        
        if not os.path.exists(params['all_chromatin_makrs_all_cells_combined_dir_path']):
            dict_cell_lines_info = parse_cellinfodict_to_populate_data(params['CellInfo_dict_file'], cell_names_start_with="#")
            populate_cellinfo_dirs(dict_cell_lines_info, params['CellInfo_target_dir'])
        
            os.mkdir(params['all_chromatin_makrs_all_cells_combined_dir_path'])
            print("Generating chromatin data for all the cells")
            os.system("cat " + params['CellInfo_target_dir'] + "/*/ChIP-seq/*ChIP-seq.bed4 " + params['CellInfo_target_dir'] + "/*/ChromatinStates/*_ChromatinStates.bed4 " + params['CellInfo_target_dir'] + "/*/DNase-seq/*_DNase-seq.bed4 | "  + """awk '{print $0 >> \"""" + params['all_chromatin_makrs_all_cells_combined_dir_path'] + """/"$1".bed4"}'""")
            if os.path.exists(params['CAGE_expr_file_path']):
                print("Appending CAGE peak expr data to the cell info data")
                os.system("cat " + params['CAGE_expr_file_path']+ """ | awk '{print $0 >> \"""" + params['all_chromatin_makrs_all_cells_combined_dir_path'] + """/"$1".bed4"}'""")
            if os.path.exists(params['ContactingDomains_file_path']):
                print("Appending Contacting domains to the cell info data")
                os.system("cat " + params['ContactingDomains_file_path'] + """ | awk '{print $0 >> \"""" + params['all_chromatin_makrs_all_cells_combined_dir_path'] + """/"$1".bed4"}'""")
            if os.path.exists(params['LoopDomains_file_path']):
                print("Appending Loop domains to the cell info data")
                os.system("cat " + params['LoopDomains_file_path'] + """ | awk '{print $0 >> \"""" + params['all_chromatin_makrs_all_cells_combined_dir_path'] + """/"$1".bed4"}'""")
            
            if os.path.exists(params['ReplicationDomains_file_path']):
                print("Appending Replication domains to the cell info data")
                os.system("cat " + params['ReplicationDomains_file_path'] + """ | awk '{print $0 >> \"""" + params['all_chromatin_makrs_all_cells_combined_dir_path'] + """/"$1".bed4"}'""")
            print("collected all the chromatin states.")
        
        if run_in_parallele and len(Mutations_dir_list)>1:
            p = Pool(int(params['number_processes_to_run_in_parallel']))
            for mutations_file in Mutations_dir_list:
                if "chr" not in mutations_file.split('/')[-1]:
                    continue 
                mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated_file = mutations_file.split('/')[-1] + "_annotatedscoredmerged"
                mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated_files.append(mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated_file)
                if not os.path.exists(mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated_file): 
                    #run in parallel for different mutations
                    print("Added to the queue: " + mutations_file)
                    p.apply_async(annotate_mutations_file, args=(mutations_file, params['elements_file'], params['motif_sites_dir'], params['motif_PFM_file'], params['all_chromatin_makrs_all_cells_combined_dir_path'], tumor_cells_dict, tumor_cells_tracks_dict, list_of_cell_tracks, assay_type_cell_tracks_dict, motifTFName_TFNames_matches_dict, mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated_file, normal_expression_per_cell_per_TF, tumor_expression_per_cell_per_TF, 8, remove_temp_files, col_names, run_training, weights_per_param_dict, take_abs_entropy_diff, int(params['log_base']), header, only_unify_per_sample, int(params['window_overlap_to_merge']), params['chromatin_states_to_consider'].split(','), use_estimates_from_simulation_set, mean_and_sd_from_the_simulations_indiv_sites_dict, retrieve_estimates_from_simulation_set, bed12_format_bool, params))
                else:
                    print("Using existing: " + mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated_file)
            p.close()
            p.join()
            print("finished all process")
        else:
            for mutations_file in Mutations_dir_list:
                if "chr" not in mutations_file.split('/')[-1]:
                    continue
                mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated_file = mutations_file.split('/')[-1] + "_annotatedscoredmerged"
                mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated_files.append(mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated_file)
                if not os.path.exists(mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated_file): 
                    print("start seq run")
                    annotate_mutations_file(mutations_file, params['elements_file'], params['motif_sites_dir'], params['motif_PFM_file'], params['all_chromatin_makrs_all_cells_combined_dir_path'], tumor_cells_dict, tumor_cells_tracks_dict, list_of_cell_tracks, assay_type_cell_tracks_dict, motifTFName_TFNames_matches_dict, mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated_file, normal_expression_per_cell_per_TF, tumor_expression_per_cell_per_TF, 8, remove_temp_files, col_names, run_training, weights_per_param_dict, take_abs_entropy_diff, int(params['log_base']), header, only_unify_per_sample, int(params['window_overlap_to_merge']), params['chromatin_states_to_consider'].split(','), use_estimates_from_simulation_set, mean_and_sd_from_the_simulations_indiv_sites_dict, retrieve_estimates_from_simulation_set, bed12_format_bool, params)
                else:
                    print("Using existing: " + mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated_file)
        if len(mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated_files) < 1:
            print('no mutation file is processed (check name (it should contain chr in the name) and content of the given directory or the given file)')
        #collect the results and merge them
        if params['annotated_mutations_final_output_file']!="none":
            if retrieve_estimates_from_simulation_set:
                with open(params['annotated_mutations_final_output_file']+"_withoutfilterunifiedpersample", 'w') as annotated_mutations_final_output_file_ranked_outfile_unifiedpersample:
                    for annotated_mutations_file in mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated_files:
                        if os.path.exists(annotated_mutations_file+"_withoutfilterunifiedpersample"):
                            if os.stat(annotated_mutations_file+"_withoutfilterunifiedpersample").st_size > 2:
                                with open(annotated_mutations_file+"_withoutfilterunifiedpersample", 'r') as ranked_file_infile_unifiedpersample:
                                    annotated_mutations_final_output_file_ranked_outfile_unifiedpersample.write(ranked_file_infile_unifiedpersample.read())
                            if remove_temp_files:
                                os.remove(annotated_mutations_file+"_withoutfilterunifiedpersample")
                                
                if os.stat(params['annotated_mutations_final_output_file']+"_withoutfilterunifiedpersample")>2:
                    #get mean and sd from the unified muts
                    overall_scores_input_file = params['annotated_mutations_final_output_file']+"_withoutfilterunifiedpersample"
                    mean_and_sd_from_the_simulations_indiv_sites_dict = get_mean_and_sd_from_simulation_sets(params['annotated_mutations_final_output_file']+"_withoutfilterunifiedpersample", 'none', score_index_simulation_file=int(params['score_index_simulation_file']))
                    #get element scores after filtering those muts that have a score less than mean from the previous
                    processed_files = []
                    for annotated_mutations_file in mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated_files:
                        if os.path.exists(annotated_mutations_file):
                            if os.stat(annotated_mutations_file).st_size > 2:
                                processed_file = process_annotated_scored_mutations(annotated_mutations_file, annotated_mutations_file+'_merged', only_unify_per_sample = only_unify_per_sample, window_overlap_to_merge = params['window_overlap_to_merge'], header=header, elements_file=params['elements_file'], chromatin_states_to_consider =params['chromatin_states_to_consider'].split(','), use_estimates_from_simulation_set=True, mean_and_sd_from_the_simulations_indiv_sites_dict=mean_and_sd_from_the_simulations_indiv_sites_dict, retrieve_estimates_from_simulation_set = False, bed12_format_bool = bed12_format_bool, params=params)
                                if os.path.exists(annotated_mutations_file):
                                    os.remove(annotated_mutations_file)
                                processed_files.append(processed_file)
                    mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated_files = processed_files
                    
            with open(params['annotated_mutations_final_output_file'], 'w') as annotated_mutations_final_outfile:
                #if header:
                    #annotated_mutations_final_outfile.write('\t'.join(col_names[0:]) + '\n')#write header only to the final combined file
                print("Generating a combined set of annotations for all mutation")
                for annotated_mutations_file in mutations_motifs_scored_all_chromatin_makrs_all_cells_grouped_annotated_files:
                    if os.path.exists(annotated_mutations_file):
                        if os.stat(annotated_mutations_file).st_size > 2:
                            print("Concatenating (size check): " + annotated_mutations_file)
                            with open(annotated_mutations_file, 'r') as annotated_mutations_infile:
                                annotated_mutations_final_outfile.write(annotated_mutations_infile.read())
                        if remove_temp_files:
                            os.remove(annotated_mutations_file) 
    else:
        print("Using the existing final annotated muts file: " + params['annotated_mutations_final_output_file'])
        #check intersection of annotated mutations with required sets and filter
    #filter the annotated mutations according to given conditions, merge the motifs to calculate p-values per window
    #needs improvement from here and optimization from above
    if not run_training:
        #ranking process is deprecated
        if rank_scores:
            annotated_mutations_final_output_file_ranked = params['annotated_mutations_final_output_file']+"_ranked"
            col_names.append('rank')
            if not os.path.exists(annotated_mutations_final_output_file_ranked):
                index_tumor_type = 5
                index_score=61
                ranked_files = []
                tumor_types_dir = params['annotated_mutations_final_output_file'] + "_splitForRanking/"
                #normalize the scores per cancer type
                if os.path.isfile(params['annotated_mutations_final_output_file']):
                    if not os.path.exists(tumor_types_dir):
                        os.system("mkdir " + tumor_types_dir)
                    os.system('awk {print $0 >> "' + tumor_types_dir + '"$'+ str(index_tumor_type+1)  + ' }') 
                tumor_type_files = os.listdir(tumor_types_dir)
                #start a pool of processes
                if run_in_parallele:
                    print("running in parallel - ranking")
                    k = Pool(int(params['number_processes_to_run_in_parallel']))
                    for tumor_type_file in tumor_type_files:
                        tumor_type_file_ranked = tumor_types_dir + tumor_type_file + "_ranked"
                        ranked_files.append(tumor_type_file_ranked)
                        if not os.path.exists(tumor_type_file_ranked): 
                            #run in parallel for different mutations
                            print("Added to the queue: " + tumor_types_dir + tumor_type_file)
                            k.apply_async(rank_scored_muts, args=(tumor_types_dir+tumor_type_file, tumor_type_file_ranked, index_score))
                    k.close()
                    k.join()
                else:
                    print("run in seq.")
                    for tumor_type_file in tumor_type_files:
                        tumor_type_file_ranked = tumor_type_file+"_ranked"
                        ranked_files.append(tumor_type_file_ranked)
                        if not os.path.exists(tumor_type_file_ranked): 
                            rank_scored_muts(tumor_types_dir+tumor_type_file, tumor_types_dir+tumor_type_file_ranked, index_score)
                
                #collect the results and merge them
                if params['annotated_mutations_final_output_file']!="none":
                    with open(annotated_mutations_final_output_file_ranked, 'w') as annotated_mutations_final_output_file_ranked_outfile:
                        print("Generating a combined file for ranked annotated-mutations")
                        for ranked_file in ranked_files:
                            if os.path.exists(ranked_file):
                                if os.stat(ranked_file).st_size > 2:
                                    print("Concatinating (size check): " + ranked_file)
                                    with open(ranked_file, 'r') as ranked_file_infile:
                                        annotated_mutations_final_output_file_ranked_outfile.write(ranked_file_infile.read())
                            
                print("Finished score ranking")
            else:
                print("Using the existing ranked file: " + annotated_mutations_final_output_file_ranked)
            #filter low scored muts
            #combine recurrent motifs (add scores of each mut-motif so that highly recurrent mutated motifs get higher scores than non-recurrent ones or ?)
            #print("Filtering, merging and grouping mut-motifs ")
            annotated_mutations_final_output_file_scored_merged = annotated_mutations_final_output_file_ranked + "rankedmerged"
            if not os.path.exists(annotated_mutations_final_output_file_scored_merged):
                process_annotated_scored_mutations_rankedMutsbyPercentile(annotated_mutations_final_output_file_ranked, annotated_mutations_final_output_file_scored_merged, only_unify_per_sample = False)
        #chr,start,stop,sampleID,mutinfo,motifName,score,rank,motifinfo
        #if not os.path.exists(annotated_mutations_final_output_file_filteredmerged):
        #    os.system('''awk '$1 ~ "chr"' ''' + annotated_mutations_final_output_file + " > " + params['annotated_mutations_final_output_file']_filteredmerged)#""" | awk 'BEGIN{FS=OFS="\t"}{if($16>=0.4 && $18 !~ "NA" && $18>0 && ($19=="NA" || $19>1.0) && ($20=="NA" || $20>1.0) && $21 !~ "NA" && $21>0 && $20 !~ "NA" && $22 !~ "Quies") print}' | awk 'BEGIN{FS=OFS="\t"}{print $10,$11,$12,$13,$14,$15,$1"%"$2"%"$3"%"$4"%"$5"%"$6"%"$7"%"$8"%"$9,$9}' | awk 'BEGIN{FS=OFS="\t"}{gsub(",","&", $7); print}' | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 -k6,6 | groupBy -g 1-6 -c 8,8,7 -o count,distinct,distinct > """ + annotated_mutations_final_output_file_filteredmerged)
        annotated_mutations_final_output_file_scored_merged = params['annotated_mutations_final_output_file']
        annotated_mutations_statcalc_output_file =  annotated_mutations_final_output_file_scored_merged + "_statcalc"
        if compute_significance:
        #calculate a p-value for each mutated motif, does the motif get such a score or a higher score at random? also take into account the mutation frequency of the samples and motif-length+others?
        #calculate FDR for each motif based on on all the reported motifs
            if not os.path.exists(annotated_mutations_statcalc_output_file):
                print("Calculating p-values")
                print("getting mut frequency per sample")
                sample_id_and_number_of_mutations_per_sample_dict = get_number_of_mutations_per_sample_list_and_write_to_file(Mutations_dir_list, "", index_sample_ids=8)
                number_of_elements_tested = file_len(annotated_mutations_final_output_file_scored_merged)
                if header:
                    number_of_elements_tested-=1
                calculate_p_value_motifregions(annotated_mutations_final_output_file_scored_merged, sample_id_and_number_of_mutations_per_sample_dict, mutated_regions_pval_outfile=annotated_mutations_statcalc_output_file, index_mutation_frequency=5, index_sample_ids=4, index_elment_start_coordinate=1, index_elment_stop_coordinate=2, genome_size=3100000000.0, total_number_tested_regions=number_of_elements_tested)
        else:
            annotated_mutations_statcalc_output_file =  annotated_mutations_final_output_file_scored_merged
        
        annotated_mutations_statcalc_output_file_scores_significance = annotated_mutations_statcalc_output_file + "_scoresSig"
        if compute_score_sig:
            if not os.path.exists(annotated_mutations_statcalc_output_file_scores_significance):
                print("Calculating p-values for the scores")
                mean_and_sd_from_the_simulations_dict = get_mean_and_sd_from_simulation_sets(params['scored_simulation_input_file'], params['mean_and_sd_output_file'], score_index_simulation_file=int(params['score_index_simulation_file']))
                print("mean and sd from: " + params['scored_simulation_input_file'])
                print(mean_and_sd_from_the_simulations_dict)
                get_pvalues_per_element_score(annotated_mutations_statcalc_output_file, annotated_mutations_statcalc_output_file_scores_significance, mean_and_sd_from_the_simulations_dict, element_score_index_obs_file=int(params['element_score_index_obs_file']))
        else:
            annotated_mutations_statcalc_output_file_scores_significance = annotated_mutations_statcalc_output_file
        #calculate significance of the given score for each mutated element based by comparing the score to the distribution of scores in the simulated sets.
    
    #run training set to get the odds ratio (weight) of each param
    else:
        #run logit regression on the combined dataset for a given list of column names, indexes, outcome_col name and index
        #write the ORs from the logit to a file in a single line
        #a subset of col_names
        #from the header: ["same_factor", "same_cell_other_factors", "dnase", "contactingdomain", "loopdomain", "CAGE_expr", "chromatin_states", "replicationdomain"]
        col_names_to_weight = ['Entropy_diff', 
                                'same_factor', 'other_factors', 'dnase', 'CAGE_expr']#, 'contactingdomain', 'loopdomain', 'CAGE_expr']#, 'chromatin_states', 'replicationdomain']
        if use_gene_expression_for_scoring:
            col_names_to_weight = ['Entropy_diff', 'normal_expression_level_of_this_motif',  
                                'same_factor', 'other_factors', 'dnase', 'CAGE_expr']#, 'contactingdomain', 'loopdomain', 'CAGE_expr']#, 'chromatin_states', 'replicationdomain']
        col_indexes_to_weight = []
        for c in col_names_to_weight:
            if c in col_names:
                col_indexes_to_weight.append(col_names.index(c))
            else:
                del col_names_to_weight[c]
        outcome_col_name = 'Mut_score'
        outcome_col_index = 6
        
        weights_per_param = getWeightsforAnnotations(params['annotated_mutations_final_output_file'],
                                 outcome_col_name = outcome_col_name,
                                 outcome_col_index = outcome_col_index, 
                                 list_of_cols_names = col_names_to_weight,
                                 list_of_cols_index = col_indexes_to_weight,
                                 take_abs_entropy_diff=take_abs_entropy_diff, log_base=int(params['log_base']))
        #write the logit output to a file to be used for non-training sets
        with open(params['weights_per_param_dict_arg_file'], 'w') as weights_per_param_dict_arg_outfile:
            out_line = []
            for param in range(0, len(weights_per_param)):
                out_line.append(weights_per_param.index[param] + "=" + str(weights_per_param[param]))
            if len(out_line)>0:
                weights_per_param_dict_arg_outfile.write(','.join(out_line)+'\n')
        
