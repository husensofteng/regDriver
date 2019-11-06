'''
Created on Mar 29, 2016

@author: Husen M. Umer
'''
import sys,os
from pybedtools import BedTool
import numpy as np
import math
from collections import Counter

'''The input file should have the following format:
9-col mutation-info
6-col motif-info
1-col entropy diff
1-col snv_pos at motif
4-col element-info
'''

def group_snvs_motifs_with_element_names(snv_motif_element_input_file, colnum_mut_info=1, colnum_motif_info=10, colnum_element_info=18):
    snv_motif_element_grouped = '.'.join(snv_motif_element_input_file.split('.')[0:-1])+"_grouped.bed" + str(colnum_motif_info+6+2)
    if not os.path.exists(snv_motif_element_grouped):
        bedTools_obj = BedTool(snv_motif_element_input_file)
        bedTools_obj.groupby(g=range(colnum_mut_info, colnum_element_info), c=colnum_element_info+3, o = ['distinct']).saveas(snv_motif_element_grouped) 
    
    return snv_motif_element_grouped

def annotate_snvs_with_elements(snv_motif_element_grouped_input_file, tumor_cells_dict, tumor_cells_tracks_dict, list_of_cell_tracks, assay_type_cell_tracks_dict, motifTFName_TFNames_matches_dict, 
                                snv_motif_element_grouped_output_file="", normal_expression_per_cell_per_TF = {}, tumor_expression_per_cell_per_TF = {}, 
                                index_mut_info=0, index_motif_info=9, index_element_info=15, col_names = [], run_training = True, weights_per_param_dict = {}, take_abs_entropy_diff=True, log_base=10, header=True):
    
    index_tumor_name = 5
    index_motif_name = 12
    overlaps_info_index = 17
    sep = '\t'
    if snv_motif_element_grouped_output_file=="":
        snv_motif_element_grouped_output_file = snv_motif_element_grouped_input_file + '_annotated'
    
    with open(snv_motif_element_grouped_input_file, 'r') as snv_motif_element_grouped_infile:
        with open(snv_motif_element_grouped_output_file, 'w') as snv_motif_element_grouped_outfile:
            if header:
                snv_motif_element_grouped_outfile.write(sep.join(col_names) + '\n')
            
            split_line = snv_motif_element_grouped_infile.readline().strip().split(sep)
            
            overlapping_motifs_same_mutation = []
            index_mutations_info = 0; index_mutated_motif_info = 9
            index_tf_normal_expr=17; index_tf_tumor_expr=18;
            index_same_factor=19;
            current_mutation = split_line[index_mutations_info:index_mutated_motif_info]
            write_last_mut_lines = False
            counter = 0
            sum = 0
            print("Mutation Annotation has Started")

            while len(split_line)>=overlaps_info_index:
                "a line every after processing each 100000 lines"
                counter+=1
                if counter==100000:
                    sum+=counter
                    print("#lines processed: " + str(sum) + " from " + snv_motif_element_grouped_input_file + " into " + snv_motif_element_grouped_output_file)
                    counter = 0
                #define and initialize the variables    
                same_factor = 'NA'
                same_cell_same_factor = ['NA']
                other_cells_same_factor = ['NA']#for other cell names
                
                other_factors = 'NA'
                same_cell_other_factors = ['NA']#names of other factors sames cell
                same_cell_other_factors_names = ['NA']
                
                dnase = 'NA'
                same_cell_dnase = ['NA']
                other_cells_dnase = ['NA']#for other cell names
                
                contactingdomain = 'NA'
                same_cell_contactingdomain = ['NA']
                other_cells_contactingdomain = ['NA'] 
                
                loopdomain = 'NA'
                same_cell_loopdomain = ['NA']
                other_cells_loopdomain = ['NA']
                
                chromatin_states = 'NA'
                same_cell_chromatin_states = ['NA']
                same_cell_chromatin_states_lables = ['NA']#just for counting
                other_cell_chromatin_states = ['NA']
                other_cell_chromatin_states_lables = ['NA']#just for counting
                
                replicationdomain = 'NA'
                same_cell_replicationdomain = ['NA']
                same_cells_replicationdomain_lables = ['NA']#just for counting
                other_cells_replicationdomain = ['NA']
                other_cells_replicationdomain_lables = ['NA']#just for counting
                
                CAGE_expr = 'NA'
                same_cell_CAGE_expr = ['NA']
                same_cell_CAGE_expr_values = 'NA'
                other_cells_CAGE_expr = ['NA']
                other_cells_CAGE_expr_values = 'NA'
                
                normal_expression_level_of_this_motif_from_itsFactors = 'NA'
                normal_expression_level_of_this_motif_from_overlappingFactors = ['NA']

                tumor_expression_level_of_this_motif_from_itsFactors = 'NA'
                tumor_expression_level_of_this_motif_from_overlappingFactors = ['NA']
                    
                overlap_track_names = split_line[overlaps_info_index].split(',')#get the track names from the overlapping column which is originally taken from the combined cell chromatin data bed4 file.
                tumor_name = split_line[index_tumor_name]#type of tumor which this mutation belong to
                factor_name_from_motif_name = split_line[index_motif_name].strip().split('_')[0].upper()
                
                if tumor_name in tumor_cells_dict.keys():#assign the values to 0 when the track was found otherwise let it be NA
                    cells_of_this_tumor = tumor_cells_dict[tumor_name]#from the TumorCell dict get the matching cells of this tumor 
                    for overlap_track in overlap_track_names:
                        overlap_track_split = overlap_track.split('#')
                        cell_name_from_overlap_track = overlap_track_split[0]
                        assay_type = overlap_track_split[1]
                        if overlap_track_split[0]+'#'+overlap_track_split[1] in tumor_cells_tracks_dict[tumor_name]:#{'PBCA-DE': ['SK-N-SH#ChromHMM#ChromatinStates', 'SK-N-SH#YY1',...], ...}
                            #a matching cell is found
                            #for tracks that have no value just showing in which of the matching cells it is found is sufficient
                            if assay_type == "DNase-seq":
                                if same_cell_dnase[0]=='NA':
                                    same_cell_dnase = []
                                same_cell_dnase.append(cell_name_from_overlap_track) 
                            
                            elif assay_type == "ContactingDomain":
                                if same_cell_contactingdomain[0]=="NA":
                                    same_cell_contactingdomain=[]
                                same_cell_contactingdomain.append(cell_name_from_overlap_track)
                            
                            elif assay_type == "LoopDomain":
                                if same_cell_loopdomain[0]=="NA":
                                    same_cell_loopdomain=[]
                                same_cell_loopdomain.append(cell_name_from_overlap_track)
                            #for those tracks that have a value we need to show their value in each of the matching cells
                            elif assay_type == "ChromHMM":
                                if same_cell_chromatin_states[0]=="NA":
                                    same_cell_chromatin_states=[]
                                    same_cell_chromatin_states_lables = []
                                same_cell_chromatin_states.append(overlap_track)
                                same_cell_chromatin_states_lables.append(overlap_track_split[2])
                                
                            elif assay_type == "RepliDomain":
                                if same_cell_replicationdomain[0]=="NA":
                                    same_cell_replicationdomain=[]
                                    same_cells_replicationdomain_lables = []
                                same_cell_replicationdomain.append(overlap_track)
                                same_cells_replicationdomain_lables.append(overlap_track_split[2])
                            
                            elif assay_type == "FANTOM":
                                if same_cell_CAGE_expr[0]=='NA':
                                    same_cell_CAGE_expr = []
                                same_cell_CAGE_expr.append(overlap_track)
                                if same_cell_CAGE_expr_values=='NA':
                                    same_cell_CAGE_expr_values = 0.0
                                same_cell_CAGE_expr_values += float(overlap_track_split[2])
                            else:#ChIP-seq #if assay_type != "ChromHMM" and assay_type != "DNase-seq" and assay_type != "ContactingDomain" and assay_type != "LoopDomain" and assay_type != "FANTOM" and assay_type != "RepliDomain":
                                if overlap_track.split('#')[1].upper() == factor_name_from_motif_name or overlap_track.split('#')[1].upper() in motifTFName_TFNames_matches_dict[factor_name_from_motif_name]:
                                    #the motif has a matching TF overlap
                                    if same_cell_same_factor[0] == 'NA':
                                        same_cell_same_factor = []
                                    same_cell_same_factor.append(cell_name_from_overlap_track)#save what cells had the signal
                                else:
                                    if same_cell_other_factors[0]=="NA":
                                        same_cell_other_factors = []
                                    same_cell_other_factors.append(overlap_track.split('#')[1])
                                    if same_cell_other_factors_names[0]=="NA":
                                        same_cell_other_factors_names = []
                                    same_cell_other_factors_names.append(overlap_track)#save what factors have binding at this place
                        else:#the tracks belong to a non-matching cell
                            if assay_type == "DNase-seq":
                                if other_cells_dnase[0]=='NA':
                                    other_cells_dnase = []
                                other_cells_dnase.append(cell_name_from_overlap_track) 
                            
                            elif assay_type == "ContactingDomain":
                                if other_cells_contactingdomain[0]=="NA":
                                    other_cells_contactingdomain=[]
                                other_cells_contactingdomain.append(cell_name_from_overlap_track)
                            
                            elif assay_type == "LoopDomain":
                                if other_cells_loopdomain[0]=="NA":
                                    other_cells_loopdomain=[]
                                other_cells_loopdomain.append(cell_name_from_overlap_track)
                            
                            elif assay_type == "ChromHMM":
                                if other_cell_chromatin_states[0]=="NA":
                                    other_cell_chromatin_states=[]
                                    other_cell_chromatin_states_lables = []
                                other_cell_chromatin_states.append(overlap_track)
                                other_cell_chromatin_states_lables.append(overlap_track_split[2])
                                
                            elif assay_type == "RepliDomain":
                                if other_cells_replicationdomain[0]=="NA":
                                    other_cells_replicationdomain=[]
                                    other_cells_replicationdomain_lables = []
                                other_cells_replicationdomain.append(overlap_track)
                                other_cells_replicationdomain_lables.append(overlap_track_split[2])
                                
                            elif assay_type == "FANTOM":
                                if other_cells_CAGE_expr[0]=='NA':
                                    other_cells_CAGE_expr = []
                                other_cells_CAGE_expr.append(overlap_track)
                                if other_cells_CAGE_expr_values=='NA':
                                    other_cells_CAGE_expr_values = 0.0
                                other_cells_CAGE_expr_values += float(overlap_track_split[2])
                                
                            else:#ChIP-seq
                                if overlap_track.split('#')[1].upper() == factor_name_from_motif_name or overlap_track.split('#')[1].upper() in motifTFName_TFNames_matches_dict[factor_name_from_motif_name]:
                                    #the motif has a matching TF overlap
                                    if other_cells_same_factor[0] == 'NA':
                                        other_cells_same_factor = []
                                    other_cells_same_factor.append(cell_name_from_overlap_track)#save what cells had the signal
                    
                    #Get expression level for the TF from the matching tumor cell lines
                    for cell in cells_of_this_tumor:
                        #Normal gene expression values
                        if cell+'#'+factor_name_from_motif_name in normal_expression_per_cell_per_TF.keys():
                            if normal_expression_per_cell_per_TF[cell+'#'+factor_name_from_motif_name]!='NA':
                                if normal_expression_level_of_this_motif_from_overlappingFactors[0] == 'NA':
                                    normal_expression_level_of_this_motif_from_overlappingFactors = []
                                normal_expression_level_of_this_motif_from_overlappingFactors.append(normal_expression_per_cell_per_TF[cell+'#'+factor_name_from_motif_name])
                        else:
                            for alt_factor in motifTFName_TFNames_matches_dict[factor_name_from_motif_name]:
                                if cell+'#'+alt_factor in normal_expression_per_cell_per_TF.keys():
                                    if normal_expression_per_cell_per_TF[cell+'#'+alt_factor]!='NA':
                                        if normal_expression_level_of_this_motif_from_overlappingFactors[0] == 'NA':
                                            normal_expression_level_of_this_motif_from_overlappingFactors = []
                                        normal_expression_level_of_this_motif_from_overlappingFactors.append(normal_expression_per_cell_per_TF[cell+'#'+alt_factor])
                        
                        #Tumor gene expression values        
                        if cell+'#'+factor_name_from_motif_name in tumor_expression_per_cell_per_TF.keys():
                            if tumor_expression_per_cell_per_TF[cell+'#'+factor_name_from_motif_name]!='NA':
                                if tumor_expression_level_of_this_motif_from_overlappingFactors[0] == 'NA':
                                    tumor_expression_level_of_this_motif_from_overlappingFactors = []
                                tumor_expression_level_of_this_motif_from_overlappingFactors.append(tumor_expression_per_cell_per_TF[cell+'#'+factor_name_from_motif_name])
                        else:
                            for alt_factor in motifTFName_TFNames_matches_dict[factor_name_from_motif_name]:
                                if cell+'#'+alt_factor in tumor_expression_per_cell_per_TF.keys():
                                    if tumor_expression_per_cell_per_TF[cell+'#'+alt_factor ]!='NA':
                                        if tumor_expression_level_of_this_motif_from_overlappingFactors[0] == 'NA':
                                            tumor_expression_level_of_this_motif_from_overlappingFactors = []
                                        tumor_expression_level_of_this_motif_from_overlappingFactors.append(tumor_expression_per_cell_per_TF[cell+'#'+alt_factor])
                    
                    #Get final values for each of the variables
                    #Normal gene expression
                    if normal_expression_level_of_this_motif_from_overlappingFactors[0]!='NA':
                        normal_expression_level_of_this_motif_from_itsFactors = np.mean(normal_expression_level_of_this_motif_from_overlappingFactors)    
                    #Tumor gene expression
                    if tumor_expression_level_of_this_motif_from_overlappingFactors[0]!='NA':
                        tumor_expression_level_of_this_motif_from_itsFactors = np.mean(tumor_expression_level_of_this_motif_from_overlappingFactors)    
                    #Same TF
                    same_tf_chipseq_in_cells_for_this_tumor = 0
                    for cell in cells_of_this_tumor:
                        if "ChIP-seq#"+cell in assay_type_cell_tracks_dict.keys():
                            if factor_name_from_motif_name in assay_type_cell_tracks_dict["ChIP-seq#"+cell]:
                                same_tf_chipseq_in_cells_for_this_tumor += 1
                            else:
                                for alt_factor in motifTFName_TFNames_matches_dict[factor_name_from_motif_name]:
                                    if alt_factor in assay_type_cell_tracks_dict["ChIP-seq#"+cell]:
                                        same_tf_chipseq_in_cells_for_this_tumor += 1
                                        break
                    same_tf_in_cells_for_this_tumor = 0
                    if same_cell_same_factor[0]!='NA':
                        same_factor=1.0
                        same_tf_in_cells_for_this_tumor = len(set(same_cell_same_factor))
                    elif other_cells_same_factor[0]!='NA':
                        if same_tf_chipseq_in_cells_for_this_tumor>0:#if chip-seq for the factor was available for the matching cell but no signal was found then report zero and don't consider the other cells
                            same_factor = 0
                            same_cell_same_factor[0]  = 'None'
                        else:#if no chip-seq was available for the factor in the matching cells then check the other cells
                            same_factor=len(set(other_cells_same_factor))*0.25
                            if same_factor>1.0:
                                same_factor=1.0
                    else:   # same_cell_same_factor[0]=='NA' and other_cells_same_factor[0]=='NA':
                        same_tf_other_cells = False 
                        for assay_name in assay_type_cell_tracks_dict.keys():
                            if "ChIP-seq#" in assay_name:
                                if factor_name_from_motif_name in assay_type_cell_tracks_dict[assay_name]:
                                    same_tf_other_cells += 1
                                else:
                                    for alt_factor in motifTFName_TFNames_matches_dict[factor_name_from_motif_name]:
                                        if alt_factor in assay_type_cell_tracks_dict[assay_name]:
                                            same_tf_other_cells += 1
                                            break
                        if same_tf_chipseq_in_cells_for_this_tumor>0:#if chip-seq for the factor was available for the matching cell but no signal was found then report zero and don't consider the other cells
                            same_factor = 0
                            same_cell_same_factor[0]  = 'None'
                        if same_tf_other_cells-same_tf_chipseq_in_cells_for_this_tumor>0:
                            other_cells_same_factor[0] = 'None'
                            same_factor = 0
                        else:
                            same_factor = 'NA'
                    
                    #Other TFs - same cell
                    if same_cell_other_factors[0]!='NA':
                        other_factors = len(set(same_cell_other_factors))
                        '''if len(set(same_cell_other_factors))>=1:
                            other_factors = 1#len(set(same_cell_other_factors))
                        else:
                            other_factors = 0'''
                    else:
                        num_tfs_for_this_tumor = 0
                        for cell in cells_of_this_tumor:
                            if "ChIP-seq#"+cell in assay_type_cell_tracks_dict.keys():
                                num_tfs_for_this_tumor+=len(assay_type_cell_tracks_dict["ChIP-seq#"+cell])
                        if num_tfs_for_this_tumor-same_tf_chipseq_in_cells_for_this_tumor>0:
                            same_cell_other_factors[0] = 'None'
                            same_cell_other_factors_names[0] = 'None'
                            other_factors = 0
                        else:
                            other_factors = 'NA'
                    #DNase1
                    if same_cell_dnase[0]!='NA':
                        dnase=1.0
                    elif other_cells_dnase[0]!='NA':
                        dnase=len(set(other_cells_dnase))/len(set(assay_type_cell_tracks_dict['DNase-seq']))
                        if dnase>1.0:
                            dnase=1.0
                    else:   #if same_cell_dnase[0]=='NA' and other_cells_dnase[0]=='NA':
                        num_dnase_in_cells_of_this_tumor = 0
                        num_dnase_in_other_cells = 0
                        if "DNase-seq" in assay_type_cell_tracks_dict.keys():
                            for cell in cells_of_this_tumor:
                                if cell in assay_type_cell_tracks_dict["DNase-seq"]:
                                    num_dnase_in_cells_of_this_tumor+=1
                            num_dnase_in_other_cells = len(assay_type_cell_tracks_dict["DNase-seq"])-num_dnase_in_cells_of_this_tumor
                        if num_dnase_in_cells_of_this_tumor>0 or num_dnase_in_other_cells>0: 
                            if num_dnase_in_cells_of_this_tumor>0:
                                same_cell_dnase[0] = 'None'
                                dnase = 0
                            if num_dnase_in_other_cells>0:
                                other_cells_dnase[0] = 'None'
                                dnase = 0
                        else:
                            dnase = 'NA'
                    #TAD contacting domains
                    if same_cell_contactingdomain[0]!='NA':
                        contactingdomain=1.0
                    elif other_cells_contactingdomain[0]!='NA':
                        contactingdomain=len(set(other_cells_contactingdomain))/len(set(assay_type_cell_tracks_dict['ContactingDomain']))
                        if contactingdomain>1.0:
                            contactingdomain=1.0
                    else:   #if same_cell_contactingdomain[0]=='NA' and other_cells_contactingdomain[0]=='NA':
                        num_contactingdomain_in_cells_of_this_tumor = 0
                        num_contactingdomain_in_other_cells = 0
                        if "ContactingDomain" in assay_type_cell_tracks_dict.keys():
                            for cell in cells_of_this_tumor:
                                if cell in assay_type_cell_tracks_dict["ContactingDomain"]:
                                    num_contactingdomain_in_cells_of_this_tumor+=1
                            num_contactingdomain_in_other_cells = len(assay_type_cell_tracks_dict["ContactingDomain"])-num_contactingdomain_in_cells_of_this_tumor
                        if num_contactingdomain_in_cells_of_this_tumor>0 or num_contactingdomain_in_other_cells>0:
                            if num_contactingdomain_in_cells_of_this_tumor>0:
                                same_cell_contactingdomain[0] = 'None'
                                contactingdomain = 0
                            if num_contactingdomain_in_other_cells>0:
                                other_cells_contactingdomain[0] = 'None'
                                contactingdomain = 0
                        else:
                            contactingdomain = 'NA'
                    #TAD Loops
                    if same_cell_loopdomain[0]!='NA':
                        loopdomain=1.0
                    elif other_cells_loopdomain[0]!='NA':
                        loopdomain=len(set(other_cells_loopdomain))/len(set(assay_type_cell_tracks_dict['LoopDomain']))
                        if loopdomain>1.0:
                            loopdomain=1.0
                    else:   #if same_cell_loopdomain[0]=='NA' and other_cells_loopdomain[0]=='NA':
                        num_loopdomain_in_cells_of_this_tumor = 0
                        num_loopdomain_in_other_cells = 0
                        if "LoopDomain" in assay_type_cell_tracks_dict.keys():
                            for cell in cells_of_this_tumor:
                                if cell in assay_type_cell_tracks_dict["LoopDomain"]:
                                    num_loopdomain_in_cells_of_this_tumor+=1
                            num_loopdomain_in_other_cells = len(assay_type_cell_tracks_dict["LoopDomain"])-num_loopdomain_in_cells_of_this_tumor
                        if num_loopdomain_in_cells_of_this_tumor>0 or num_loopdomain_in_other_cells>0:
                            if num_loopdomain_in_cells_of_this_tumor>0:
                                same_cell_loopdomain[0] = 'None'
                                loopdomain = 0
                            if num_loopdomain_in_other_cells>0:
                                other_cells_loopdomain[0] = 'None'
                                loopdomain = 0
                        else:
                            loopdomain = 'NA'
                    #CAGE
                    if same_cell_CAGE_expr_values!='NA':
                        CAGE_expr = same_cell_CAGE_expr_values/(len(same_cell_CAGE_expr)*1.0)
                    elif other_cells_CAGE_expr_values!='NA':
                        CAGE_expr= other_cells_CAGE_expr_values/len(set(assay_type_cell_tracks_dict['FANTOM']))
                    else:   #same_cell_CAGE_expr_values=='NA' and other_cells_CAGE_expr_values=='NA':
                        num_cage_expr_in_cells_of_this_tumor = 0
                        num_cage_expr_in_other_cells = 0
                        if "FANTOM" in assay_type_cell_tracks_dict.keys():
                            for cell in cells_of_this_tumor:
                                if cell in assay_type_cell_tracks_dict["FANTOM"]:
                                    num_cage_expr_in_cells_of_this_tumor+=1
                            num_cage_expr_in_other_cells = len(assay_type_cell_tracks_dict["FANTOM"])-num_cage_expr_in_cells_of_this_tumor
                        if num_cage_expr_in_cells_of_this_tumor>0 or num_cage_expr_in_other_cells>0:
                            if num_cage_expr_in_cells_of_this_tumor>0:
                                same_cell_CAGE_expr[0] = 'None'
                                CAGE_expr = 0
                            if num_cage_expr_in_other_cells>0:
                                other_cells_CAGE_expr[0] = 'None'
                                CAGE_expr = 0
                        else:
                            CAGE_expr = 'NA'
                    
                    #chromHMM
                    if same_cell_chromatin_states_lables[0]!='NA':
                        chromatin_states = Counter(same_cell_chromatin_states_lables).most_common(1)[0][0]
                    elif other_cell_chromatin_states_lables[0]!='NA':
                        chromatin_states = Counter(other_cell_chromatin_states_lables).most_common(1)[0][0]
                    else:
                        num_chromatin_states_in_cells_of_this_tumor = 0
                        num_chromatin_states_in_other_cells = 0
                        if "ChromatinStates" in assay_type_cell_tracks_dict.keys():
                            for cell in cells_of_this_tumor:
                                if cell in assay_type_cell_tracks_dict["ChromatinStates"]:
                                    num_chromatin_states_in_cells_of_this_tumor+=1
                            num_chromatin_states_in_other_cells = len(assay_type_cell_tracks_dict["ChromatinStates"])-num_chromatin_states_in_cells_of_this_tumor
                        if num_chromatin_states_in_cells_of_this_tumor>0 or num_chromatin_states_in_other_cells>0: 
                            if num_chromatin_states_in_cells_of_this_tumor>0:
                                same_cell_chromatin_states[0] = 'None'
                                same_cell_chromatin_states_lables[0] = 'None'
                                chromatin_states = 'None'
                            if num_chromatin_states_in_other_cells>0:
                                other_cell_chromatin_states[0] = 'None'
                                other_cell_chromatin_states_lables[0]='None'
                                chromatin_states = 'None'
                        else:
                            chromatin_states = 'NA'
                    #Replication domains
                    if same_cells_replicationdomain_lables[0]!='NA':
                        replicationdomain = Counter(same_cells_replicationdomain_lables).most_common(1)[0][0]
                    elif other_cells_replicationdomain_lables[0]!='NA':
                        replicationdomain = Counter(other_cells_replicationdomain_lables).most_common(1)[0][0]
                    else:
                        num_replicationdomain_in_cells_of_this_tumor = 0
                        num_replicationdomain_in_other_cells = 0
                        if "RepliDomain" in assay_type_cell_tracks_dict.keys():
                            for cell in cells_of_this_tumor:
                                if cell in assay_type_cell_tracks_dict["RepliDomain"]:
                                    num_replicationdomain_in_cells_of_this_tumor+=1
                            num_replicationdomain_in_other_cells = len(assay_type_cell_tracks_dict["RepliDomain"])-num_replicationdomain_in_cells_of_this_tumor
                        if num_replicationdomain_in_cells_of_this_tumor>0 or num_replicationdomain_in_other_cells>0: 
                            if num_replicationdomain_in_cells_of_this_tumor>0:
                                same_cell_replicationdomain[0] = 'None'
                                same_cells_replicationdomain_lables[0] = 'None'
                                replicationdomain = 'None'
                            if num_replicationdomain_in_other_cells>0:
                                other_cells_replicationdomain[0] = 'None'
                                other_cells_replicationdomain_lables[0]='None'
                                replicationdomain = 'None'
                        else:
                            replicationdomain = 'NA'                        
                        
                cols_to_use_cobminedcells = [str(same_factor), str(other_factors), str(dnase), str(contactingdomain), str(loopdomain), str(CAGE_expr), chromatin_states, replicationdomain]
                #report the calculated variables and entropy_diff,motif name info ,mut info    
                output_line = ('\t'.join(split_line[index_mut_info: index_mut_info+9]) + sep + sep.join(split_line[index_motif_info: index_motif_info+8])
                                                                                                     #TF, gene expr (norml|tumor)
                                                                                                     + sep + str(normal_expression_level_of_this_motif_from_itsFactors)
                                                                                                     + sep + str(tumor_expression_level_of_this_motif_from_itsFactors)
                                                                                                     #info for scoring
                                                                                                     + sep + sep.join(cols_to_use_cobminedcells)#same_factor_overlaps +  sep + all_factors_same_cell_overlaps + sep + dnase_overlap + sep + contactingdomain_overlap + sep + loopdomain_overlap + sep + chromhmm_overlap_to_report + sep + replicationtiming_label_overlap_to_report + sep + cage_overlap 
                                                                                                     #overlapping tracks from matching cells (TFBS, DNase1, chromatin state, replication label, contacting domain, loop domain)
                                                                                                     + sep + ','.join(set(same_cell_same_factor)) 
                                                                                                     + sep + ','.join(set(same_cell_dnase)) 
                                                                                                     + sep + ','.join(set(same_cell_chromatin_states))
                                                                                                     + sep + ','.join(set(same_cell_replicationdomain))
                                                                                                     + sep + ','.join(set(same_cell_contactingdomain))
                                                                                                     + sep + ','.join(set(same_cell_loopdomain))
                                                                                                     + sep + ','.join(set(same_cell_CAGE_expr))
                                                                                                     #Same cell Other factors
                                                                                                     + sep + ','.join(set(same_cell_other_factors_names))
                                                                                                     + sep + ','.join(set(same_cell_other_factors))
                                                                                                     #the same factor in other cells
                                                                                                     + sep + ','.join(set(other_cells_same_factor))
                                                                                                     #DNase1
                                                                                                     + sep + ','.join(set(other_cells_dnase))
                                                                                                     #Chromatin states
                                                                                                     + sep + ','.join(set(other_cell_chromatin_states))
                                                                                                     #Replication domains
                                                                                                     + sep + ','.join(set(other_cells_replicationdomain))
                                                                                                     #Contacting domains
                                                                                                     + sep + ','.join(set(other_cells_contactingdomain))
                                                                                                     #Loop domains
                                                                                                     + sep + ','.join(set(other_cells_loopdomain))
                                                                                                     #CAGE
                                                                                                     + sep + ','.join(set(other_cells_CAGE_expr))
                                                                                                     + '\n') 
                
                overlap_line = output_line.strip().split(sep)
                if overlap_line[index_mutations_info:index_mutated_motif_info] != current_mutation:#when a new mutation is processed write the content of the dict to the output file to save results of the previous mutation (from all the tf families, with one line for each family - only the one with the highest diff_freq)
                    #scan the dict to determine what to write from the overlapping motifs of this mut
                    annotated_mutations_output_line = []
                    if run_training:
                        annotated_mutations_output_lines = get_overlapping_motifs_to_write(overlapping_motifs_same_mutation, index_same_factor=19, index_other_factors=20, index_dnase=21, index_cage_expr=24)
                        if len(annotated_mutations_output_lines)>0:
                            annotated_mutations_output_line = get_one_motif(annotated_mutations_output_lines, index_entropy_diff=15, sep=sep)
                    else:
                        annotated_mutations_output_line = filter_and_score_overlapping_motifs_to_write(overlapping_motifs_same_mutation, sep, col_names, weights_per_param_dict, take_abs_entropy_diff=take_abs_entropy_diff, log_base=log_base)
                    if len(annotated_mutations_output_line)>0:
                        snv_motif_element_grouped_outfile.write('\n'.join(annotated_mutations_output_line)+'\n')
                    overlapping_motifs_same_mutation = []
                    current_mutation = overlap_line[index_mutations_info:index_mutated_motif_info]
                
                if overlap_line not in overlapping_motifs_same_mutation:
                    overlapping_motifs_same_mutation.append(overlap_line)
                
                line = snv_motif_element_grouped_infile.readline()
                split_line = line.strip().split(sep)
                if line=="":
                    write_last_mut_lines = True
                #snv_motif_element_grouped_outfile.write(output_line)
            if write_last_mut_lines:
                annotated_mutations_output_line = []
                if run_training:
                    annotated_mutations_output_lines = get_overlapping_motifs_to_write(overlapping_motifs_same_mutation, index_same_factor=19, index_other_factors=20, index_dnase=21, index_cage_expr=24)
                    if len(annotated_mutations_output_lines)>0:
                        annotated_mutations_output_line = get_one_motif(annotated_mutations_output_lines, index_entropy_diff=15, sep=sep)
                else:
                    annotated_mutations_output_line = filter_and_score_overlapping_motifs_to_write(overlapping_motifs_same_mutation, sep, col_names, weights_per_param_dict, take_abs_entropy_diff=take_abs_entropy_diff, log_base=log_base)
                if len(annotated_mutations_output_line)>0:
                    snv_motif_element_grouped_outfile.write('\n'.join(annotated_mutations_output_line) + '\n' )
                overlapping_motifs_same_mutation = []

    print("Finished the annotation process")
    return snv_motif_element_grouped_output_file

"""
Keep motifs that have a direct factor support, if non then keep motifs with indirect support (same factor other cells) and same cell DNase1, if none then keep motifs for which the factor chip-seq has not been performed in any cell but dnase from some cell and other factors from the same cell do exist. 
"""
def get_overlapping_motifs_to_write(overlapping_motifs_same_mutation, 
                               index_same_factor=19, index_other_factors=20, index_dnase=21, index_cage_expr=24):
    
    motifs_with_samefactor_support = []
    motifs_with_samefactor_othercells_support= []
    motifs_with_na_same_factor = []
    motifs_with_support_of_no_binding = []
    for mut_info in overlapping_motifs_same_mutation:    
        if mut_info[index_same_factor] == "NA":
            motifs_with_na_same_factor.append(mut_info)
        elif float(mut_info[index_same_factor])>=1.0:#direct support from the same factor or at least x other cells
            motifs_with_samefactor_support.append(mut_info)
        elif float(mut_info[index_same_factor]) > 0:#indirect support from the same factor other cells (value is between 0-1)
            motifs_with_samefactor_othercells_support.append(mut_info)
        elif float(mut_info[index_same_factor]) == 0.0 or mut_info[index_same_factor] == "0":
            motifs_with_support_of_no_binding.append(mut_info)

    if len(motifs_with_samefactor_support)>0:
        return motifs_with_samefactor_support
    elif len(motifs_with_samefactor_othercells_support)>0:
        return motifs_with_samefactor_othercells_support
    elif len(motifs_with_support_of_no_binding)>0:
        return motifs_with_support_of_no_binding
    elif len(motifs_with_na_same_factor)>0:
        return motifs_with_na_same_factor
        
def get_one_motif(annotated_mutations_output_lines, index_entropy_diff=15, sep='\t'):
    index_highest_entropy_diff = 0
    highest_entropy_diff = 0
    for i in range(0, len(annotated_mutations_output_lines)):
        if abs(float(annotated_mutations_output_lines[i][index_entropy_diff])) > highest_entropy_diff:
            highest_entropy_diff = abs(float(annotated_mutations_output_lines[i][index_entropy_diff]))
            index_highest_entropy_diff = i
    return [sep.join(annotated_mutations_output_lines[index_highest_entropy_diff])]
    
def filter_and_score_overlapping_motifs_to_write(overlapping_motifs_same_mutation, sep, col_names, weights_per_param_dict, index_normal_expression_level_of_this_motif=17, index_tumor_expression_level_of_this_motif=18, index_same_factor=19, index_all_factors_same_cell_overlaps=20, index_dnase=21, index_cage_overlap=24, take_abs_entropy_diff=True, log_base=10):
    
    overlapping_motifs_to_report = []
    for mut_motif_info in overlapping_motifs_same_mutation:
        Sum_params_with_logit_larger_than_zero = 0.0 #weighted accumulative score
        Sum_params_with_logit_smaller_than_zero = 0.0
        final_score_is_calculated = False
        filter = add_this_line_to_overlapping_motifs_check(mut_motif_info, index_normal_expression_level_of_this_motif, index_tumor_expression_level_of_this_motif, index_same_factor)
        if filter:
            for param in weights_per_param_dict.keys():
                if param in col_names:# and weights_per_param_dict[param]>0:
                    if mut_motif_info[col_names.index(param)]!="NA" and mut_motif_info[col_names.index(param)]!="None":
                        if float(mut_motif_info[col_names.index(param)]) != 0: 
                            final_score_is_calculated = True#if the value of this feature is not NA and not zero then use it to calculate its annotation score 
                            if param == 'normal_expression_level_of_this_motif' or param == 'tumor_expression_level_of_this_motif' or param == 'other_factors' or param == 'CAGE_expr':
                                if float(mut_motif_info[col_names.index(param)]) >=1: #log10 of x is in minus when x<1
                                    try:
                                        new_value = (math.log(float(mut_motif_info[col_names.index(param)]), log_base))
                                        if weights_per_param_dict[param]>0:                                            
                                            new_value*=weights_per_param_dict[param]
                                            Sum_params_with_logit_larger_than_zero+=new_value
                                        else:
                                            Sum_params_with_logit_smaller_than_zero+=new_value
                                        mut_motif_info[col_names.index(param)] = str(new_value)
                                    except ValueError:
                                        continue 
                            elif param == 'Entropy_diff':
                                new_value = abs(float(mut_motif_info[col_names.index(param)]))
                                if weights_per_param_dict[param]>0:
                                    new_value*=weights_per_param_dict[param]
                                    Sum_params_with_logit_larger_than_zero+=new_value
                                else:
                                    Sum_params_with_logit_smaller_than_zero+=new_value

                                if take_abs_entropy_diff:
                                    #keep the minus to differentiate between motif disruption and enhancement
                                    if float(mut_motif_info[col_names.index(param)])<0:
                                        new_value = new_value*-1
                                mut_motif_info[col_names.index(param)] = str(new_value)
                            else:
                                new_value = float(mut_motif_info[col_names.index(param)])
                                if weights_per_param_dict[param]>0:
                                    new_value*=weights_per_param_dict[param]
                                    Sum_params_with_logit_larger_than_zero+=new_value
                                else:
                                    Sum_params_with_logit_smaller_than_zero+=new_value
                                mut_motif_info[col_names.index(param)] = str(new_value)
                                
            if final_score_is_calculated:#in case the final score was computed then write the results out
                #use the log_coeffecients and report exp(final score)
                Sum_params_with_logit_smaller_than_zero = 0
                mut_motif_info.append("{0:.3f}".format(np.exp(Sum_params_with_logit_larger_than_zero)))#+Sum_params_with_logit_smaller_than_zero))#Ignore negative coefficients 
                overlapping_motifs_to_report.append(sep.join(mut_motif_info))
    return overlapping_motifs_to_report
   
def add_this_line_to_overlapping_motifs_check(mut_motif_info, index_normal_expression_level_of_this_motif, index_tumor_expression_level_of_this_motif, index_same_factor):
    tf_expression_check = False
    samefactor_check = False
    if mut_motif_info[index_normal_expression_level_of_this_motif]!="NA":
        if float(mut_motif_info[index_normal_expression_level_of_this_motif])>=1:
            tf_expression_check = True
    if mut_motif_info[index_tumor_expression_level_of_this_motif]!="NA": 
        if float(mut_motif_info[index_tumor_expression_level_of_this_motif])>=1:
            tf_expression_check = True
    if mut_motif_info[index_normal_expression_level_of_this_motif]=="NA" and mut_motif_info[index_tumor_expression_level_of_this_motif]=="NA":
        tf_expression_check = True
    if mut_motif_info[index_same_factor]!="NA":
        if float(mut_motif_info[index_same_factor])>0:    
            samefactor_check = True
    else:
        samefactor_check = True
    
    if tf_expression_check and samefactor_check:
        return True
    else:
        return False




