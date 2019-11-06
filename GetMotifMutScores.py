'''
Created on March 26, 2016

@author: Husen M. Umer
'''
import re
import string
import sys, os
from itertools import groupby
from scipy.stats import binom, hypergeom
from statsmodels.sandbox.stats.multicomp import multipletests
'''
Input file of containing mutations (9cols) followed by motifs (6 col)
Processes: identify the snv position and calcualte the entropy-change 
>unify the motifs 
>generate motif-position distribution matrix   
'''
def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False
'''
desc: reports motif-breaking info for mutations in each motif sites. It check the difference between nucleotide frequency of the ref allele and the mutated-to allele of each mutation in the PWM of the anchor motif.
in: motifs_PFM matrix (get it from ENCODE, for instance), mutated motifs infput file (contains mutations info), output file name
out: all mutated motifs with adding breaking info to the mutations, all mutations at the motifs with adding breaking info  
out: only mutated motifs have motif-breaking mutation (difference between ref and mut allele > given_Threshold) and the fourth returned file is the list of motif breaking mutations (diff>Threshold)
'''
def get_freq_per_motif(motif_PFM_input_file):
    "given a PEM file return a dict, a key for each tf and the freq as value"
    PFM_motifs_dict = {}
    PFM_motifs_infile = open(motif_PFM_input_file, 'r')
    PFM_motifs_lines = PFM_motifs_infile.readlines()
    PFM_motifs_infile.close()
    A = "A"
    C = "C"
    G = "G"
    T = "T"
    motif_info_sep = ' '
    motif_name = ""
    header_line = PFM_motifs_lines[0]
    PFM_fomrat = ''
    if 'MEME' in PFM_motifs_lines[0]: 
        PFM_fomrat = 'MEME'
        motif_info_sep = ' '
        freq_nucl_sep = '\t'
        for line in PFM_motifs_lines:
            if line.strip()!="" and not line.startswith('letter') and not line.startswith('URL'):
                if line.startswith('MOTIF'):
                    motif_name = line.strip().split(motif_info_sep)[2].upper()
                else:
                    if motif_name!="":#if it has been initialized
                        if motif_name not in PFM_motifs_dict:
                            PFM_motifs_dict[motif_name] = []
                        split_line = line.strip().split(freq_nucl_sep)
                        freq_per_allele = []
                        for s in split_line:
                            if s.strip()!="" and isfloat(s.strip()):
                                freq_per_allele.append(s.strip())
                        if len(freq_per_allele)==4: #freq of the 4 nucleotides
                            PFM_motifs_dict[motif_name].append({A: float(split_line[0]), C: float(split_line[1]), G: float(split_line[2]), T: float(split_line[3])})
        #END of reading motif info (ACGT frequencies at each position of the motif)
    else:
        PFM_fomrat = "ENCODE"
        for line in PFM_motifs_lines:
            if line.startswith('>'):
                if '_' in line.strip().split(motif_info_sep)[0].strip('>'):
                    motif_name = line.strip().split(motif_info_sep)[0].strip('>').split('_')[0].upper() + "_" + line.strip().split(motif_info_sep)[0].strip('>').split('_')[1]
                else:
                    motif_name = line.strip().split(motif_info_sep)[0].strip('>').upper()
            else:
                if motif_name not in PFM_motifs_dict:
                    PFM_motifs_dict[motif_name] = []
                split_line = line.strip().split(motif_info_sep)
                if len(split_line)==5: #starts with a letter and freq of the 4 nucleotides
                    PFM_motifs_dict[motif_name].append({A: float(split_line[1]), C: float(split_line[2]), G: float(split_line[3]), T: float(split_line[4])})
        #END of reading motif info (ACGT frequencies at each position of the motif)
    return PFM_motifs_dict, PFM_fomrat

def score_motifs_according_to_their_affect(motif_PFM_input_file, muts_mutated_motifs_input_file, mutated_motifs_scored_output_file="", index_mutations_info=0, index_mutated_motif_info=9, index_fromAllele = 3, index_toAllele = 4):
    "given a mutation mutated_motifs file (containing muts info (9cols) and followed by motifs info (6 cols)"
    if mutated_motifs_scored_output_file=="":
        mutated_motifs_scored_output_file = muts_mutated_motifs_input_file + "_scoredmuts"
    
    PFM_motifs_dict, PFM_fomrat = get_freq_per_motif(motif_PFM_input_file)
    
    motif_lengths_not_matched = {}
    motif_names_not_matched = {}
    with open(muts_mutated_motifs_input_file, 'r') as muts_mutated_motifs_inpfile:
        with open(mutated_motifs_scored_output_file, 'w') as scored_mut_motif_outfile:
            muts_mutated_motifs_line = muts_mutated_motifs_inpfile.readline()
            overlapping_motifs_same_mutation = {}
            split_line = muts_mutated_motifs_line.strip().split('\t')
            current_mutation = '\t'.join(split_line[index_mutations_info:index_mutated_motif_info])
            write_last_mut_lines = False
            while muts_mutated_motifs_line!="":
                if not ((int(split_line[index_mutations_info+2]) < int(split_line[index_mutated_motif_info+1].strip())) or (int(split_line[index_mutations_info+1]) > int(split_line[index_mutated_motif_info+2].strip()))): #if the mut start was greater than the motif end or the mut end was smaller than the motif start then it means that the mutation is outside of the motif and has no overlap with it, thus it needs to be ignored.
                    #Get motif_name, stand and length
                    motif_name = ""
                    if PFM_fomrat=="ENCODE":
                        motif_name = split_line[index_mutated_motif_info+3]
                    if PFM_fomrat=="MEME":
                        motif_name = split_line[index_mutated_motif_info+3].split('_')[0].upper()
                    motif_site_strand = split_line[index_mutated_motif_info+5]
                    motif_length = int(split_line[index_mutated_motif_info+2])-int(split_line[index_mutated_motif_info+1])+1 #because it is 1-based
                    #JUST FOR REPOTING
                    if motif_name not in PFM_motifs_dict.keys():
                        if motif_name not in motif_names_not_matched:
                            motif_names_not_matched[motif_name] = str(motif_length)
                        muts_mutated_motifs_line = muts_mutated_motifs_inpfile.readline()
                        split_line = muts_mutated_motifs_line.strip().split('\t')
                        if muts_mutated_motifs_line=="":
                            write_last_mut_lines = True
                        continue
                    if motif_length != len(PFM_motifs_dict[motif_name]):
                        if motif_name not in motif_lengths_not_matched:
                            motif_lengths_not_matched[motif_name] = "Motif_length:" + str(motif_length) + "-" + "PFM_motif_length:"+str(len(PFM_motifs_dict[motif_name]))
                        muts_mutated_motifs_line = muts_mutated_motifs_inpfile.readline()
                        split_line = muts_mutated_motifs_line.strip().split('\t')
                        if muts_mutated_motifs_line=="":
                            write_last_mut_lines = True
                        continue
                    #Get info about the current mut
                    snv_position_index = (int(split_line[index_mutations_info+1].strip())) - int(split_line[index_mutated_motif_info+1].strip()) #get index of the snv position in the motif (the position would be index+1) #SNV start position - Motif site start position 
                    from_allele = split_line[index_fromAllele].strip()
                    to_allele = split_line[index_toAllele].strip()
                    ref_allele_length = len(from_allele)
                    to_allele_length = len(to_allele)
                    if motif_site_strand == "-": 
                        #in order to read the motif vector from end to start for motif sites which are from the negative strand
                        snv_position_index = int(split_line[index_mutated_motif_info+2].strip()) - int(split_line[index_mutations_info+1].strip())
                        #PFM_motifs_dict[motif_name][snv_position_index]['T'] > PFM_motifs_dict[motif_name][snv_position_index]['A']:
                        from_allele = from_allele.translate(string.maketrans('ACGT','TGCA'))
                        to_allele = to_allele.translate(string.maketrans('ACGT','TGCA'))
                    
                    #beacuse the positions are stored as list of indeces so they start from 0, thus no need to +1
                    #if snv_position_index>=0: #since some mutations were overlabpping with motifs, but after unifying the motifs some of them were cut to unify them for this reason mutations that were located in those regions are not part of the reported motif so they lie outside of them
                    #    snv_position_index+=1 #because of BED indexing
                    #the idea is to convert the alleles to match the motif logo (if A has a higher freq than T in this position then show the A>x version of the allele otherwise show T>x version, this is true since the motif sequence can be found in either strand but since we only show the + strand of the motif we should conver all the alleles to match those of the motif 
                    #Calculate entropy diff
                    diff_freq = 0.0
                    if ref_allele_length==to_allele_length:#sign of SNV, MNV or single base indel 
                        if to_allele=='-' or from_allele=='-':
                            if to_allele =='-' and from_allele!='-':#single base deletion and no insertion
                                diff_freq = float(format(PFM_motifs_dict[motif_name][snv_position_index][from_allele], '.5f'))
                            if from_allele=='-':
                                diff_freq = 1.0
                        else:#no del and no ins
                            #if ref_allele_length==1 and to_allele_length==1:
                            #   ref_mut_freq = format(PFM_motifs_dict[motif_name][snv_position_index][from_allele], '.5f')
                            #  alt_mut_freq = format(PFM_motifs_dict[motif_name][snv_position_index][to_allele], '.5f')
                            # diff_freq  = float(ref_mut_freq) - float(alt_mut_freq)
                            #elif ref_allele_length>1 and to_allele_length>1:#multiple base mutation
                            for i in range(0,ref_allele_length):
                                if snv_position_index+i<motif_length:
                                    diff_freq += (PFM_motifs_dict[motif_name][snv_position_index+i][from_allele[i]] - PFM_motifs_dict[motif_name][snv_position_index+i][to_allele[i]])
                            if diff_freq>1.0:
                                diff_freq=1.0#max allowed score for any mut is set to 1.0
                            #diff_freq=1.0#in case of multiple base mutation set the entropy diff to 1
                    
                    elif ref_allele_length < to_allele_length: #sign of insertion
                        if to_allele_length==2:
                            if snv_position_index+1<motif_length:#if it is one base insertion just calculate the additional value of the insertion from the base after the indel position
                                ins_allele = to_allele[1]
                                diff_freq = 1.0 - float(format(PFM_motifs_dict[motif_name][snv_position_index+1][ins_allele], '.5f'))  
                        else:
                            diff_freq = 1.0#in case of multiple base insertion set the entropy diff to 1
                    elif ref_allele_length>to_allele_length:#sign of deletion
                        if snv_position_index+1<motif_length:#account for the base after the indel
                            if ref_allele_length==2:
                                del_allele = from_allele[1]#if it is one base del just calculate the deducted value of the deletion
                                diff_freq = float(format(PFM_motifs_dict[motif_name][snv_position_index+1][del_allele], '.5f'))  
                        else:
                            diff_freq = 1.0#in case of multiple base deletion set the entropy diff to 1
                    if '\t'.join(split_line[index_mutations_info:index_mutated_motif_info]) != current_mutation:#when a new mutation is processed write the content of the dict to the output file to save results of the previous mutation (from all the tf families, with one line for each family - only the one with the highest diff_freq)
                        for mut_motif in overlapping_motifs_same_mutation.keys():
                            scored_mut_motif_outfile.write(overlapping_motifs_same_mutation[mut_motif][1] + '\n')
                        overlapping_motifs_same_mutation = {}
                        current_mutation = '\t'.join(split_line[index_mutations_info:index_mutated_motif_info])
                        
                    mutation_motifFamily = '\t'.join(split_line[index_mutations_info:index_mutated_motif_info])+'\t'+motif_name.split('_')[0]
                    if mutation_motifFamily not in overlapping_motifs_same_mutation.keys():
                        overlapping_motifs_same_mutation[mutation_motifFamily] = [diff_freq, '\t'.join(split_line[index_mutations_info:index_mutated_motif_info]) + "\t" + '\t'.join(split_line[index_mutated_motif_info:index_mutated_motif_info+6]) + '\t' + str(diff_freq) + '\t' + str(snv_position_index+1)]
                    else:
                        if abs(diff_freq) > abs(overlapping_motifs_same_mutation[mutation_motifFamily][0]):#replace with the previous line of the same tf family if it had a higher diff_freq
                            overlapping_motifs_same_mutation[mutation_motifFamily] = [diff_freq, '\t'.join(split_line[index_mutations_info:index_mutated_motif_info]) + "\t" + '\t'.join(split_line[index_mutated_motif_info:index_mutated_motif_info+6]) + '\t' + str(diff_freq) + '\t' + str(snv_position_index+1)]
                        
                    #scored_mut_motif_outfile.write('\t'.join(split_line[index_mutations_info:index_mutated_motif_info]) + "\t" + '\t'.join(split_line[index_mutated_motif_info:index_mutated_motif_info+6]) + '\t' + str(diff_freq) + '\t' + str(snv_position_index) + '\n')    
                
                muts_mutated_motifs_line = muts_mutated_motifs_inpfile.readline()
                #if the end line (empty line) was reached then write the content of the previous mut because there will be no next iteration
                split_line = muts_mutated_motifs_line.strip().split('\t')
                if muts_mutated_motifs_line=="":
                    write_last_mut_lines = True
            if write_last_mut_lines:#if it was the last line then write it
                for mut_motif in overlapping_motifs_same_mutation.keys():
                    scored_mut_motif_outfile.write(overlapping_motifs_same_mutation[mut_motif][1] + '\n')
                    
                overlapping_motifs_same_mutation = {}
                
    if len(motif_lengths_not_matched)>0:
        print( "Motifs that their length didn't match those in the PFM file")
        print(motif_lengths_not_matched)
    if len(motif_names_not_matched)>0:
        print("Motifs that name didn't match any motif in the PFM file")
        print(motif_names_not_matched)
        
    return mutated_motifs_scored_output_file

#makes a vector for each motif by counting number of mutations at each of its positions

def calculate_p_value_motifregions(mutated_regions, sample_id_and_number_of_mutations_per_sample_dict, mutated_regions_pval_outfile="", index_mutation_frequency=7, index_sample_ids=8, index_elment_start_coordinate=10, index_elment_stop_coordinate=11, genome_size=3000000000.0, total_number_tested_regions=121667):
    
    mutated_regions_list = []
    if mutated_regions_pval_outfile=="":
        mutated_regions_pval_outfile = mutated_regions + "_calpvalues"
    with open(mutated_regions, 'r') as region_mutations_infile:
        mutated_regions_list = region_mutations_infile.readlines()
    mutated_regions_pval_outputfile = open(mutated_regions_pval_outfile, "w")
    reported_p_values = []#this holds p-values of all the regions, it will be used for p-value correction and after correction it is written to the output file as an additional column after the calculated p-value, the full list of p-values is need to make the correction test
    
    for line in mutated_regions_list:
        split_line = line.strip().split("\t")
        mutation_frequency=1
        if index_mutation_frequency!=-1:
            mutation_frequency = int(split_line[index_mutation_frequency])
        sample_ids = split_line[index_sample_ids].split(',')
        avg_proportion_of_mutations_in_the_samples_of_this_region = 0.0
        samples_counted = []
        
        for sample_id in sample_ids:
            if sample_id not in samples_counted:
                samples_counted.append(sample_id)
                if sample_id in sample_id_and_number_of_mutations_per_sample_dict.keys():
                    avg_proportion_of_mutations_in_the_samples_of_this_region += ((float(sample_id_and_number_of_mutations_per_sample_dict[sample_id]))/genome_size)
                #else:
                #    print("Key Error - sample ID not found in the initial mutations found: " + sample_id)
        p = avg_proportion_of_mutations_in_the_samples_of_this_region/(len(samples_counted)*1.0)
        n = (int(split_line[index_elment_stop_coordinate]) - int(split_line[index_elment_start_coordinate])) #* len(sample_id_and_number_of_mutations_per_sample_dict.keys()) #region length (end-start) multiplied by the total number of tested samples
        k = mutation_frequency
        p_val_of_this_region = 1 - (binom.cdf(k, n, p))
        reported_p_values.append(p_val_of_this_region)
    
    print(len(reported_p_values))
    while len(reported_p_values) < total_number_tested_regions:
        reported_p_values.append(1)
    print("correcting p-values for multiple testing")
    if len(reported_p_values)>0:
        significant_bool_report, corrected_p_values_array, alphacSidak, alphacBonf = multipletests(reported_p_values, alpha=0.05, method='fdr_bh', returnsorted=False) #returns 4 things: a boolean array contains True or False for each value meaning wether the value after correction compared to the given alpha is significant or not, an array of the values after correction, a single for corrected alpha for Sidak method, a single value for corrected alpha for Bonferroni method 
        corrected_p_values_list = corrected_p_values_array.tolist()
        for l in range(0, len(mutated_regions_list)):
            mutated_regions_pval_outputfile.write(mutated_regions_list[l].strip() + "\t" + str(reported_p_values[l]) + "\t" + str(corrected_p_values_list[l]) + "\n")
    return mutated_regions_pval_outfile

def get_number_of_mutations_per_sample_list_and_write_to_file(Mutations_dir_list, numberofmutationspersample_output_file, index_sample_ids=8):
    
    if numberofmutationspersample_output_file == "":
        numberofmutationspersample_output_file = "numebrofmutationspersample.txt"
    sample_id_and_number_of_mutations_per_sample_dict = {}  
    if not os.path.exists(numberofmutationspersample_output_file):
        print("Counting number of mutations per sample from the initial mutation file")
        for mutations_file in Mutations_dir_list:
            with open(mutations_file, "r") as mutations_infile:
                mutations_line = mutations_infile.readline()
                while mutations_line!="":
                    split_line = mutations_line.strip().split("\t")
                    sample_id_of_this_mutation = split_line[index_sample_ids].strip()
                    if sample_id_of_this_mutation not in sample_id_and_number_of_mutations_per_sample_dict.keys():
                        sample_id_and_number_of_mutations_per_sample_dict[sample_id_of_this_mutation] = 0
                    sample_id_and_number_of_mutations_per_sample_dict[sample_id_of_this_mutation] +=1
                    mutations_line = mutations_infile.readline()
        with open(numberofmutationspersample_output_file, 'w') as numberofmutationspersample_writefile:
            for sample_id in sample_id_and_number_of_mutations_per_sample_dict.keys():#write the sample and its number of mutations to the output file
                numberofmutationspersample_writefile.write(sample_id + "\t" + str(sample_id_and_number_of_mutations_per_sample_dict[sample_id]) +"\n")
    else:
        #in case the number of mutations per sample was already available then just read them and insert them to the returned list
        with open(numberofmutationspersample_output_file, 'r') as numberofmutationspersample_readfile:
            numberofmutationspersample_lines = numberofmutationspersample_readfile.readlines()
            for line in numberofmutationspersample_lines:
                sample_id_and_number_of_mutations_per_sample_dict[line.split('\t')[0].strip()] = int(line.split('\t')[1].strip()) # append sample id # append number of mutations in this sample id
    return sample_id_and_number_of_mutations_per_sample_dict


'''
in: mutated motifs 
out: unified mutated motifs
'''
def unify_motifs_sites_of_the_same_tfs(motif_sites_in_file, unified_motif_sites_out_file = "unified_knwon_motifs_in_peaks_thesameTF_in_regions_with_snp", snps_out_file = "unified_knwon_motifs_in_peaks_thesameTF_in_regions_with_snp.SNV", index_site_name=3, end_index_motif_info = 6, snp_info_exisits=True):
    "separates the given motif sites list according to the TF. A sublist for each TF, then it combines the overlapping motif sites of each TF to generate a unified set of motif sites for each TF"
    motif_sites_in = open(motif_sites_in_file, 'r') #input motif sites file, each line should start with 4 cols of (motif info): chr start stop motifName_knownN:TFName and any other followed column (info cols) will be written back. When two overlapping sites of the same TF are merged, their info cols will also be updated accordingly.
    
    unified_motif_sites_out_both_strands = open(unified_motif_sites_out_file, 'w')
    if snp_info_exisits:
        snps_in_sites = open(snps_out_file, 'w')
    sites = motif_sites_in.readlines()
    #handle the sites of every TF in the TFs_dict. Note:(TFs that belong to the same factor group defined in the mapping file, have the same motifs (the numbers in the brackets indicate how many bps should be added in the beginning and end of the mentioned motif.
    #when *NOALIGHNN is attached to a motif name it means that motif will not be aligned to the first motif thus it is reported as its without changing its length. 
    TFGroup_motif_names_dict = {'CTCF' : ['CTCF_known1', 0, 0, 'CTCF_known2', 1, 1], 
                      'FOXA' : ['FOXA_known4', 0, 0, 'FOXA_known1', 0, 3, 'FOXA_known2', 5, 2, 'FOXA_known3', 3, 2, 'FOXA_known5', 5, 2, 'FOXA_known6', 5, 1, 'FOXA_known7', 0, 1, 'FOX_1*NOALLIGN', 2, 4], #FOX_1 motif site sequences have the highest similarity to Jaspars FOXA2: MA0047.1 motif 
                      'HNF4' : ['HNF4_known12', 0, 0, 'HNF4_known13', 5, 4, 'HNF4_known1', 2, 2, 'HNF4_known2', 4, 5, 'HNF4_known3', 4, 4, 'HNF4_known4', 6, 4, 'HNF4_known5', 5, 5, 'HNF4_known6', 5, 5, 'HNF4_known7', 5, 5, 'HNF4_known8', 10, 4, 'HNF4_known9', 5, 4, 'HNF4_known10', 12, 5, 'HNF4_known11', 12, 5, 'HNF4_known14', 5, 5,  'HNF4_known15', 6, 0, 'HNF4_known16', 4, 3, 'HNF4_known17', 4, 3, 'HNF4_known18', 4, 3, 'HNF4_known19', 4, 5, 'HNF4_known20', 4, 3, 'HNF4_known21', 4, 4, 'HNF4_known22', 4, 3, 'HNF4_known23', 4, 4, 'HNF4_known24', 3, 4, 'HNF4_known25', 11, 4, 'HNF4_known26', 9, 1],
                      'MAF' : ['MAF_known7', 0, 0, 'MAF_known6', -1, 10, 'MAF_known5', -2, 8, 'MAF_known4', -1, 5, 'MAF_known3', 5, 8, 'MAF_known2', 4, 6, 'MAF_known8', -1, 10, 'MAF_known9', 3, 3, 'MAF_known10', 0, 9, 'MAF_known11', 2, 2, 'MAF_known12', 2, 1, 'MAF_known1*NOALLIGN', 0, 0, 'MAFF_1*NOALLIGN', 0, 0], #'MAF_known1', 0, 0: is not similar to the other motifs so it is considered as a seprate motif for its own. and ignored for the allignment
                      'AP1' : ['AP1_known6', 0, 0, 'AP1_known1', 1, 1, 'AP1_known2', 1, 1, 'AP1_known3', 1, 1, 'AP1_known4', 1, 1, 'AP1_known5', 2, 2, 'AP1_known7', 3, -2, 'AP1_known8', 2, 2, 'AP1_known9', 2, 3, 'AP1_known10', 3, 3, 'AP1_known11', 3, 3],
                      'MYC' : ['MYC_known9', 0, 0, 'MYC_known1', 3, 3, 'MYC_known2', 3, 3, 'MYC_known3', 3, 3, 'MYC_known4', 3, 3, 'MYC_known5', 4, 4, 'MYC_known6', 5, 5, 'MYC_known7', 6, 6, 'MYC_known8', 5, 5, 'MYC_known10', 7, 7, 'MYC_known11', 4, 4, 'MYC_known12', 6, 7, 'MYC_known13', 4, 4, 'MYC_known14', 4, 4, 'MYC_known15', 4, 6, 'MYC_known16', 5, 4, 'MYC_known17', 6, 7, 'MYC_known18', 5, 5, 'MYC_known19', 1, 3, 'MYC_known20', 7, -4, 'MYC_known21', 5, 5, 'MYC_known22', 5, 5],
                      'CEBPB' : ['CEBPB_known5', 0, 0, 'CEBPB_known1', 1, 3, 'CEBPB_known2', 1, 3, 'CEBPB_known3', 1, 4, 'CEBPB_known4', 1, 3, 'CEBPB_known6', 1, 5, 'CEBPB_known7', 1, 5, 'CEBPB_known8', 3, 5, 'CEBPB_known9', 3, 5, 'CEBPB_known10', 3, 5],
                      'SRF' : ['SRF_known6', 0, 0, 'SRF_known1', 1, 0, 'SRF_known2', 3, 2, 'SRF_known3', 1, 3, 'SRF_known4', 4, -3, 'SRF_known5', 5, -1, 'SRF_known7', 3, 4, 'SRF_known8', 3, 2, 'SRF_known9', 4, 3, 'SRF_known10', 2, 1],
                      'RXRA' : ['RXRA_known2', 0, 0, 'RXRA_known1', 0, 2, 'RXRA_known3', 2, 13, 'RXRA_known4', 2, 13, 'RXRA_known5', 2, 1, 'RXRA_known6', 2, 4, 'RXRA_known7', 2, 6, 'RXRA_known8', 4, 2, 'RXRA_known9', -2, 8, 'RXRA_known10', 3, 6, 'RXRA_known11', 3, 6, 'RXRA_known12', 3, 6, 'RXRA_known13', 3, 6, 'RXRA_known14', 3, 6, 'RXRA_known15', 4, 7],
                      'SP1' : ['SP1_known1', 0, 0,'SP1_known2', -2, -1, 'SP1_known4', 0, 0, 'SP1_known5', -2, -1, 'SP1_known6', -1, 1, 'SP1_known7', -1, 1, 'SP1_known8', -1, 0, 'SP1_known9', -1, 0, 'SP1_known3*NOALIGN', 0, 0 ],#ignored SP1_known3 since it was very different#'SP1' : ['SP1_known1', -1, -1,'SP1_known2', -3, -2, 'SP1_known3', -6, 0, 'SP1_known4', -1, -1, 'SP1_known5', -3, -2, 'SP1_known6', -2, 0, 'SP1_known7', -2, 0, 'SP1_known8', -2, -1, 'SP1_known9', -2, -1 ], #just get the core motif and trim the rest
                      'ETS' : ['ETS_known17', 0, 0, 'ETS_known1', -3, 5, 'ETS_known2', -3, 7, 'ETS_known3', 0, 8, 'ETS_known4', 0, 3, 'ETS_known5', -1, 5, 'ETS_known6', 0, 6, 'ETS_known7', -2, 10, 'ETS_known8', 2, 10, 'ETS_known9', 1, 6, 'ETS_known10', -4, 5, 'ETS_known11', 0, 8, 'ETS_known12', 0, 8, 'ETS_known13', 8, -7, 'ETS_known14', 0, 8, 'ETS_known15', 0, 0, 'ETS_known16', 0, 8, 'ETS_known18', 0, 8, 'ETS_1', 1, 5, 'ETS_2', 2, 8], #shifting ETS_known13 motif might be incorrect so it is not changed since it is opposite to known17 motif - if you still want to report one motif for ETS (instead of the current two known17 & known13) then change 'known13',0,0 to 'known13',-7,8. ETS_2 sequence contains 'AGGAAGT' from the fasta sequences
                      'REST' : ['REST_known1', 0, 0,  'REST_known2', 0, 0,  'REST_known3', 1, 1,  'REST_known4', 0, 0, ],
                      'BHLHE40' : ['BHLHE40_known2', 0, 0, 'BHLHE40_known1', 4, 4, 'BHLHE40_known3', 6, 6, 'BHLHE40_known2', 6, 6],
                      'NRF1' : ['NRF1_known2', 0, 0, 'NRF1_known1', 2, 0],
                      'RFX5' : ['RFX5_known2', 0, 0, 'RFX5_known1', 1, 0, 'RFX5_known3', 2, 2, 'RFX5_known4', 9, 0, 'RFX5_known5', 4, -1, 'RFX5_known6', 1, 1, 'RFX5_known7', 2, 1, 'RFX5_known8', 1, 1, 'RFX5_known9', 1, 1],
                      'YY1' : ['YY1_known2', 0, 0, 'YY1_known1', 1, 2, 'YY1_known3', 3, 8, 'YY1_known4', 5, 6, 'YY1_known5', 2, 7, 'YY1_known6', 5, 9, 'YY1_known7', 2, 7],
                      'ELF1' : ['ELF1_known1', 0, 0, 'ELF1_known2', 0, 0, 'ELF1_known3', 0, 0],
                      'TCF7L2' : ['TCF7L2_known5', 0, 0, 'TCF7L2_known1', 4, 5, 'TCF7L2_known2', 5, 6, 'TCF7L2_known3', 4, 2, 'TCF7L2_known4', 3, 4, 'TCF7L2_known6', 0, 0, 'TCF7L2_known7', 1, 1],
                      #'TCF7' : ['TCF7L2_known5', 0, 0,  'TCF7_1', 1, 0, 'TCF7_2', 5, 1], #length of TCF7_1 is 16 and its fasta sequences are like: AGTTCTTTGATCTCCA. length of TCF7_2 is 11 and its fasta sequences are like: CTTTGATCCTC
                      'ATF3' : ['ATF3_known16', 0, 0, 'ATF3_known1', 1, 1, 'ATF3_known2', 4, 4, 'ATF3_known3', 0, 4, 'ATF3_known4', -2, 3, 'ATF3_known5', 3, -2, 'ATF3_known6', 2, 2, 'ATF3_known7', 2, 2, 'ATF3_known8', 3, 1, 'ATF3_known9', 1, 1, 'ATF3_known10', 1, 4, 'ATF3_known11', 7, 3, 'ATF3_known12', 2, 0, 'ATF3_known13', 1, 4, 'ATF3_known14', 3, 4, 'ATF3_known15', 4, 4],
                      'NFIC' : ['NFIC_1', 0, 0, 'NFIC_2*NOALLIGN', 0, 0, 'NFIC_3*NOALLIGN', 0, 0, 'NR2E1..NFIC_1*NOALLIGN', 0, 0], #'NFIC_2', 'NFIC_3'], #_1 have sites with length 13 and 28. _2 & _3 have different sequence variations so it is hard to make a consensus motif for this TF
                      'MYB' : ['MYB_4', 0, 0, 'MYB_1*NOALLIGN', 0, 0, 'MYB_2*NOALLIGN', 0, 0, 'MYB_3*NOALLIGN', 0, 0, 'MYB_4*NOALLIGN', 0, 0, 'MYB_5*NOALLIGN', 0, 0, 'MYB_6*NOALLIGN', 0, 0, 'MYBL2_1*NOALLIGN', 0, 0, ],
                      'ESRRA': ['ESRRA_known4', 0, 0, 'ESRRA_known1', 3, -2, 'ESRRA_known2', -1, 7, 'ESRRA_known3', 4, 5, 'ESRRA_known5', -1, 4, 'ESRRA_known6', 4, -1, 'ESRRA_known7', 1, 8, 'ESRRA_known8', -6, 9, 'ESRRA_known9', -8, 9, 'ESRRA_known10', -6, 9, 'ESRRA_known11', 1, 8],
                      'TEAD4' : ['TEAD4_1', 0, 0], #change it to ['TEAD4_1', 1, 2] to match it to the Jaspar's TEAD4_1 motif
                      'MXI1' : ['MXI1_known1', 0, 0],
                      'TATA' : ['TATA_known6', 0, 0, 'TATA_known1', 2, 4, 'TATA_known2', 3, -2, 'TATA_known3', 4, 4, 'TATA_known4', 6, 3, 'TATA_known5', 3, -2],
                      'NR2C2' : ['NR2C2_known1', 0, 0],
                      'SREBP' : ['SREBP_known4', 0, 0, 'SREBP_known1', 2, -1, 'SREBP_known2', 2, -1, 'SREBP_known3', 5, 0, 'SREBP_known5', 3, -1, 'SREBP_known6', 3, -1],
                      'HSF' : ['HSF_known3', 0, 0, 'HSF_known1', 5, 2, 'HSF_known2', 1, 3, 'HSF_known4', 1, 3, 'HSF_known5', 1, 3],
                      'ZBTB7A' : ['ZBTB7A_known2', 0, 0, 'ZBTB7A_known1', 0, 6, 'ZBTB7A_known3', -1, 4, 'ZBTB7A_known4', 0, 3],
                      'CEBPD' : ['CEBPD_1', 0, 0, 'CEBPD_2', 1, 1], #binds to CCAAT but has no motif in the databases so I generated one from sequences of its sites
                      'ARID3A' : ['ARID3A_2', 0, 0],
                      'EP300' : ['EP300_known1', 0, 0],
                      'IRF' : ['IRF_known14', 0, 0, 'IRF_known1', 7, 1, 'IRF_known2', 7, 1, 'IRF_known3', 0, 6, 'IRF_known4', 5, -2, 'IRF_known5', 3, 6, 'IRF_known6', 4, 10, 'IRF_known7', 2, 4, 'IRF_known8', 2, 8, 'IRF_known9', 2, 7, 'IRF_known10', 1, 2, 'IRF_known11', 1, 6, 'IRF_known12', 7, -1, 'IRF_known13', 1, 5, 'IRF_known15', 6, 1, 'IRF_known16', 10, 0, 'IRF_known17', 6, 1, 'IRF_known18', 9, -5, 'IRF_known19', 6, 1, 'IRF_known20', 6, 1, 'IRF_known21', 5, 1],
                      'TCF12' : ['TCF12_known1', 0, 0],
                      }
    
    TFs_dict = {} #this dict holds all sites of each TF where key=TF_name and every key has a list attached to it containing all the motif sites of the key TF.
    for site in sites:
        split_site = site.strip().split('\t')
        TF_name = split_site[index_site_name].split(':')[1]
        TF_group_name = split_site[index_site_name].split(':')[0].split('_')[0]
        motif_name = split_site[index_site_name].split(':')[0]
        #sometimes the motifs of the same TF groups are written with different names such as FOX and FOXA but since only one of this names are accounted for in the matrix thus they need to renamed to match the name in the matrix.
        if TF_group_name == "FOX": TF_group_name = "FOXA"
        if TF_group_name == "MAFF": TF_group_name = "MAF"
        if TF_group_name == "MYBL2": TF_group_name = "MYB"
        if TF_group_name == "NR2E1..NFIC": TF_group_name = "NFIC"
        if TF_group_name in TFGroup_motif_names_dict.keys(): #only add sites of TFs that have their TF group name in the distance matrix and their motif name in the list of TF group motifs in the distance matrix
            if motif_name in TFGroup_motif_names_dict[TF_group_name] or motif_name+"*NOALLIGN" in TFGroup_motif_names_dict[TF_group_name]: #if the speicified motif was listen in the unified motifs distance matrix then include it for further analysis otherwise ignore it 
                if TF_name not in TFs_dict:
                    TFs_dict[TF_name] = [] #make a new list for each TF to hold all its sites
                TFs_dict[TF_name].append(site)
    #reset the variables
    split_site = []
    TF_name = ""
    TF_group_name = ""
    motif_name = ""
    all_snps = []
    number_of_unified_sites = 0
    number_unified_sites_per_TF = 0
    list_of_unified_TFs = []
    print("Unifying motif sites.....")
    for TF in TFs_dict.keys():
        TFj_motifs = [] #this will be filled with the motif names of each TF, motif-name, distance-from-start and end of-the-largest-motif
        site_snps = []
        site_snp = ""
        number_unified_sites_per_TF = 0
        if len(TFs_dict[TF])==1:
            split_site_i = TFs_dict[TF][0].strip().split('\t') #current site
            if snp_info_exisits:
                site_snp = ','.join(split_site_i[end_index_motif_info:-1])
                if split_site_i[24:-1] not in all_snps:
                    all_snps.append(split_site_i[end_index_motif_info:-1])
                if site_snp not in site_snps:
                    site_snps.append(site_snp)
            if snp_info_exisits:
                unified_motif_sites_out_both_strands.write('\t'.join(split_site_i[0:end_index_motif_info])+ "\t" + ';'.join(site_snps) + "\n")
            else: 
                unified_motif_sites_out_both_strands.write('\t'.join(split_site_i) + "\n")
            site_snps = []
            site_snp = ""    
            number_of_unified_sites+=1
            number_unified_sites_per_TF+=1
            continue
        for i in range(1, len(TFs_dict[TF])):
            split_site_i = TFs_dict[TF][i].strip().split('\t') #current site
            split_site_iminus1 = TFs_dict[TF][i-1].strip().split('\t') #previous site
            TF_name = split_site_i[index_site_name].split(':')[1]
            TF_group_name = split_site_i[index_site_name].split(':')[0].split('_')[0]
            if TF_group_name == "FOX": TF_group_name = "FOXA"
            if TF_group_name == "MAFF": TF_group_name = "MAF"
            if TF_group_name == "MYBL2": TF_group_name = "MYB"
            if TF_group_name == "NR2E1..NFIC": TF_group_name = "NFIC"
            motif_name_i = split_site_i[index_site_name].split(':')[0]
            motif_name_iminus_1 = split_site_iminus1[index_site_name].split(':')[0]
            if TF_group_name in TFGroup_motif_names_dict.keys():
                TFj_motifs = TFGroup_motif_names_dict[TF_group_name]
            else:
                TFj_motifs = []
            #if the TF was in the to-be unified TF set then do otherwise ignore it
            if len(TFj_motifs)>0:
                #for each site make a list of its mutations (uniq) and write that list to the final unified site
                #the snp which is located in each site is written at the end of the line (index -2:-6) chr1    3245103 3245103 1       MU2296257    1       C->T    SNV
                if snp_info_exisits:
                    if split_site_iminus1[end_index_motif_info:-1] not in all_snps:
                        all_snps.append(split_site_iminus1[end_index_motif_info:-1])
                    if split_site_i[end_index_motif_info:-1] not in all_snps:
                        all_snps.append(split_site_i[end_index_motif_info:-1])
                
                #in cases of full overlap: 
                #if the current site was wider (larger start & end) than the previous one ignore the previous site
                if split_site_i[0] == split_site_iminus1[0] and int(split_site_i[1]) <= int(split_site_iminus1[1]) and int(split_site_i[2]) >= int(split_site_iminus1[2]): #if the previous site was fully located in the current one then ignore the previous one
                    if snp_info_exisits:
                        site_snp = ','.join(split_site_iminus1[end_index_motif_info:-1])
                        if site_snp not in site_snps:
                            site_snps.append(site_snp)
                    else:
                        pass
                elif split_site_i[0] == split_site_iminus1[0] and int(split_site_iminus1[1]) <= int(split_site_i[1]) and int(split_site_iminus1[2]) >= int(split_site_i[2]):
                    #if the previous site was wider than the current one the update the current site to the previous line, because in the next round this site will be compared to the next one
                    if snp_info_exisits:
                        site_snp = ','.join(split_site_i[end_index_motif_info:-1])
                        if site_snp not in site_snps:
                            site_snps.append(site_snp)
                    TFs_dict[TF][i] = '\t'.join(split_site_iminus1)
                    
                #in cases of no overlap, if its the largest motif site report it, otherwise extend it to make it the same as the largest motif, change its name and then report it 
                #if the previous site was on a different chr or located before or after the current site then its a different site since there is no overlap between the two sites
                elif split_site_i[0] != split_site_iminus1[0] or (split_site_i[0] == split_site_iminus1[0] and \
                                                                 (int(split_site_iminus1[1]) < int(split_site_i[1]) and int(split_site_iminus1[2]) < int(split_site_i[1]))):
                    if snp_info_exisits:
                        site_snp = ','.join(split_site_iminus1[end_index_motif_info:-1])
                        if site_snp not in site_snps:
                            site_snps.append(site_snp)
                    
                    
                    #get motif name and it is distance from the largest motif (if the distance was zero then the site itself is belonging to the largest motif so need to change it 
                    if motif_name_iminus_1+"*NOALLIGN" not in TFj_motifs:
                        distance_of_start_motif_site_from_start_the_largest_motif = TFj_motifs[TFj_motifs.index(motif_name_iminus_1)+1]
                        distance_of_end_motif_site_from_end_the_largest_motif = TFj_motifs[TFj_motifs.index(motif_name_iminus_1)+2]
                    
                        #if distance_of_start_motif_site_from_start_the_largest_motif != 0 or distance_of_end_motif_site_from_end_the_largest_motif != 0: #check if its needs to be extended
                        #if split_site_iminus1[5].strip() == "+" or split_site_iminus1[5].strip() == ".":
                        split_site_iminus1[1] = str(int(split_site_iminus1[1])-distance_of_start_motif_site_from_start_the_largest_motif) #extend its start to match the widest motif site (the first element in the matrix)
                        split_site_iminus1[2] = str(int(split_site_iminus1[2])+distance_of_end_motif_site_from_end_the_largest_motif) #extend its end to match the widest motif site (the first element in the matrix)
                        split_site_iminus1[3] = TFj_motifs[0]+ ":" + TF_name #set name of the motif to the largest motif since its coordinates are maximized to match that one.
                    
                    if snp_info_exisits:
                        unified_motif_sites_out_both_strands.write('\t'.join(split_site_iminus1[0:end_index_motif_info])+ "\t" + ';'.join(site_snps) + "\n")
                    else: 
                        unified_motif_sites_out_both_strands.write('\t'.join(split_site_iminus1) + "\n")
                    site_snps = []
                    site_snp = ""
                    number_of_unified_sites+=1 #just to keep the count (how many unified motif sites for all TFs are reported)
                    number_unified_sites_per_TF+=1 #just to keep the count (how many unified motif sites for this TF are reported)
                #in cases of partial overlap, just keep the largest one until there is no more overlap (which will be reported by the previous elif condition) (simple version)
                #strategies: either report site of the largest one and ignore the partially overlapping short ones or 
                #merge all the overlapping sites (of the same and other motifs)
                #here the first strategy is applied to avoid false positive site since there are not many of such cases
                
                elif split_site_i[0] == split_site_iminus1[0] and \
                (  
                 (int(split_site_i[1]) > int(split_site_iminus1[1]) and int(split_site_i[2]) > int(split_site_iminus1[2]) and int(split_site_i[1]) < int(split_site_iminus1[2])) \
                 or (int(split_site_i[1]) < int(split_site_iminus1[1]) and int(split_site_i[2]) < int(split_site_iminus1[2])  and int(split_site_i[2]) > int(split_site_iminus1[1]))): #the current site's start and end are shifted from start and end of the previous site
                    if int(int(split_site_i[2])-int(split_site_i[1])) < int(int(split_site_iminus1[2])-int(split_site_iminus1[1])):
                        if snp_info_exisits:
                            site_snp = ','.join(split_site_i[end_index_motif_info:-1])
                            if site_snp not in site_snps:
                                site_snps.append(site_snp)
                        TFs_dict[TF][i] = '\t'.join(split_site_iminus1) #if the previous motif was bigger (end-start) then set it to the current site
                    else:
                        pass#the current one will be automatically set to the previous in the next round
                
                #to write the last site from the list, when it reached the last site no need to keep it just extend it and write it
                if i == len(TFs_dict[TF])-1:
                    if snp_info_exisits:
                        site_snp = ','.join(split_site_i[end_index_motif_info:-1])
                        if site_snp not in site_snps:
                            site_snps.append(site_snp)
                    
                    if motif_name_iminus_1+"*NOALLIGN" not in TFj_motifs:
                        #get motif name and it is distance from the largest motif (if the distance was zero then the site itself is belonging to the largest motif so need to change it 
                        distance_of_start_motif_site_from_start_the_largest_motif = TFj_motifs[TFj_motifs.index(motif_name_i)+1]
                        distance_of_end_motif_site_from_end_the_largest_motif = TFj_motifs[TFj_motifs.index(motif_name_i)+2]
                        
                        split_site_i[1] = str(int(split_site_i[1])-distance_of_start_motif_site_from_start_the_largest_motif) #extend its start to match the widest motif site (the first element in the matrix)
                        split_site_i[2] = str(int(split_site_i[2])+distance_of_end_motif_site_from_end_the_largest_motif) #extend its end to match the widest motif site (the first element in the matrix)
                        split_site_i[3] = TFj_motifs[0]+ ":" + TF_name #set name of the motif to the largest motif since its coordinates are maximized to match that one.
                    if snp_info_exisits:
                        unified_motif_sites_out_both_strands.write('\t'.join(split_site_i[0:end_index_motif_info])+ "\t" + ';'.join(site_snps) + "\n")
                    else:
                        unified_motif_sites_out_both_strands.write('\t'.join(split_site_i)+ "\n")
                    site_snps = []
                    site_snp = ""
                    number_of_unified_sites+=1
                    number_unified_sites_per_TF+=1
        list_of_unified_TFs.append([TF, number_unified_sites_per_TF])
    print("Total no. unified sites (mutated) of all TFs: " + str(number_of_unified_sites))
    #print("[Unified TFs, number of sites]:")
    #print(list_of_unified_TFs)
    print("********")
    
    unified_motif_sites_out_both_strands.close()
    
    if snp_info_exisits:    
        for snp in all_snps:
            snps_in_sites.write('\t'.join(snp)+"\n")
        snps_in_sites.close()
        return unified_motif_sites_out_file, snps_out_file
    else:
        return unified_motif_sites_out_file


def file_len(fname):
    i = -1
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def get_SNV_position_in_motif_sites(motif_sites_SNV_file, positions_output_file ="", report_on_negative_strand = False, index_motif_site_name=3, SNV_info_index=24, index_sampleIDs=9 ,min_number_snps_for_significant_motif_site = 10, report_motif_per_factor_group=False, extension="_PerTF"):
    #output_flie_name = '.'.join(motif_sites_SNV_file.split('.')[0:-1]) + extension +"_SNV_positions.csv"
    output_flie_name =  '.'.join(positions_output_file.split('.')[0:-1]) + extension + "_positions.csv"
    output_file = open(output_flie_name, 'w')
    #output_file_significant = open('.'.join(motif_sites_SNV_file.split('.')[0:-1]) + extension + "_SNV_significant" + str(min_number_snps_for_significant_motif_site) + "_positions.csv", 'w')
    sites_SNV = open(motif_sites_SNV_file, 'r')
    
    motifs_info_dict = {} #each key will belong to a unique motif id and stores all the info needed to make the matrix for that motif (a list to store frequency of mutations per position, another list to store all unique mutation ids per position,.... each list index+1 corresponding to a motif position
    motifs_lines = sites_SNV.readlines()
    
    for motif_line in motifs_lines:
        split_motif_line = motif_line.split('\t')
        motif_length = int(split_motif_line[index_motif_site_name-1].strip()) - (int(split_motif_line[index_motif_site_name-2].strip())) #get motif length (end- (start+1)) because the start is written as zero-based so the first bp should be ignored. Although the actual motif length is this calculated length +1 but since we are using and index based method to store the number of mutations per positions so we should use length-1.  
        motif_name = split_motif_line[index_motif_site_name].strip() + "(" + str(motif_length+1) +")" #count the motifs of each TF separately, ex count FOXA_x:FOX1 independently than FOXA_x:FOXA2
        motif_site_strand = split_motif_line[index_motif_site_name+2].strip()
        if motif_name not in motifs_info_dict:
            motifs_info_dict[motif_name] = {}
            #a motif name is set as key and a list containing 5 sub lists stores information related to this motif/key 
            #initialize the starting values for each sublist, in order to define the indeces for each position of the motif
            
            #add the first list to hold frequency of mutations per position, and total frequency of this motif at the last index
            motifs_info_dict[motif_name]["frequency"] = []
            motifs_info_dict[motif_name]["AlleleChanges"] = [] #add the second list to hold allele changes and their frequency per position, and all the observed allele changes at this motif
            motifs_info_dict[motif_name]["UnifiedAlleleChanges"] = []
            motifs_info_dict[motif_name]["SampleIDs"] = [] #add the third list to hold unique sample IDs per position, and all the sampleIDs found at this motif
            motifs_info_dict[motif_name]["No.UniqueSamples"] = [] #add the fourth list to hold the number of unique sample IDs that have a mutation at this position, and the total number of sampleIDs with mutation at this motif
            motifs_info_dict[motif_name]["No.RecurrentMuts"] = [] #add the fifth list to hold number of recurrent mutations per position and the total number of recurrent mutations at this motif
            motifs_info_dict[motif_name]["MutIds"] = [] #add the sixth list to hold unique mutationIDs found at each position, and all mutationIDs at this motif at the last index
            motifs_info_dict[motif_name]["MutPositions"] = [] #add the seventh list to hold unique genomic positoin of each mutation per each position, and the genomic positions of all the mutations at this motif 
            motifs_info_dict[motif_name]["MotifPositions"] = [] #add the eight list to hold genomic positions of the mutated sites (report no duplicates), and the last index to contain genomic positions of all the mutated sites of this motif 
            motifs_info_dict[motif_name]["No.UniqueMotifs"] = [] #add the ninth list to hold the number of unique mutated sites per position, and the last index to contain the number of unique mutated sites per motif 
            for l in range(0,motif_length+2):
                motifs_info_dict[motif_name]["frequency"].append(0)
                motifs_info_dict[motif_name]["AlleleChanges"].append([])
                motifs_info_dict[motif_name]["UnifiedAlleleChanges"].append([])
                motifs_info_dict[motif_name]["SampleIDs"].append([])
                motifs_info_dict[motif_name]["No.UniqueSamples"].append(0)
                motifs_info_dict[motif_name]["No.RecurrentMuts"].append(0)
                motifs_info_dict[motif_name]["MutIds"].append([])
                motifs_info_dict[motif_name]["MutPositions"].append([])
                motifs_info_dict[motif_name]["MotifPositions"].append([])
                motifs_info_dict[motif_name]["No.UniqueMotifs"].append(0)
          
        motif_genomic_position = "{}:{}-{}".format(split_motif_line[0], split_motif_line[1], split_motif_line[2])
        snps_at_this_motif = split_motif_line[SNV_info_index].split(';') #example of a snp info line: chr19,21325153,21325153,.,ssm_v2_29140,1,C->T,SNV,HX9T
        for snp in snps_at_this_motif:
            snp_inf = snp.split(',')
            variant_covers_the_whole_motif = False
            number_of_extra_bp_from_the_begining = 0
            if (int(snp_inf[1].strip()) < int(split_motif_line[index_motif_site_name-2].strip())):
                number_of_extra_bp_from_the_begining = int(split_motif_line[index_motif_site_name-2].strip()) - int(snp_inf[1].strip()) 
                if split_motif_line[5] == "-": 
                    number_of_extra_bp_from_the_begining = int(snp_inf[1].strip()) - int(split_motif_line[index_motif_site_name-2].strip())
                    #in cases where the motif is located on minus strand and the variant covers the whole motif (snp_start<motif_start and snp_end>motif_end) then consider the variants as located at position zero of the motif
                    if (int(snp_inf[2].strip()) > int(split_motif_line[index_motif_site_name-1].strip())):
                        #number_of_extra_bp_from_the_begining = int(split_motif_line[index_motif_site_name-1].strip()) - int(snp_inf[2].strip()) 
                        variant_covers_the_whole_motif = True
            snv_position = (int(snp_inf[1].strip())) - int(split_motif_line[index_motif_site_name-2].strip())  + number_of_extra_bp_from_the_begining#SNV start position (+1) - Motif site start position (+ 1) (because the motifs and mutations are zero-based so the first base should be ignored))
            if not report_on_negative_strand:
                if split_motif_line[5] == "-": #in order to read the motif vector from end to start for motif sites which are from the negative strand
                    snv_position = int(split_motif_line[index_motif_site_name-1].strip()) - int(snp_inf[1].strip())  + number_of_extra_bp_from_the_begining
                    #to put the variant at position zero of the motif not end (if this conditiona and the previous one for calculating number_of_extra_bp_from_the_begining were not given the variant position in the motif would become the last position instead of the begining as it is now with these conditions
                    
                    if variant_covers_the_whole_motif:
                        snv_position = 0#int(split_motif_line[index_motif_site_name-1].strip()) - int(snp_inf[2].strip())  + number_of_extra_bp_from_the_begining #motif-end - snv-end and add the extra bases, variants that fully cover a full motif start from position 0
                    
            #snv_position = (int(snp_inf[1].strip())) - int(split_motif_line[index_motif_site_name-2].strip())  #SNV start position (+1) - Motif site start position (+ 1) (because the motifs and mutations are zero-based so the first base should be ignored))
            mutation_genomic_position = "{}:{}-{}".format(snp_inf[0], snp_inf[1], snp_inf[2])
            mutation_id = snp_inf[4].strip()
            frequency_of_this_mutation = int(snp_inf[5].strip()) #chr, start, stop, strand, name, frquency
            sample_ids = snp_inf[index_sampleIDs].strip().replace(".somatic", "KAN").split('.')
            from_to_allele = snp_inf[6].strip()
            
            if snv_position <= motif_length and snv_position>=0: #this condition makes the multi-indels disappear from the position matrix and the final motif position figure. some of the mutations are located outside of the unified motifs (they were located in a motif but after unifying some motifs were chopped off so these mutations are no longer considered to be located within the motif thus should not be considered 
                motifs_info_dict[motif_name]["frequency"][snv_position]+=frequency_of_this_mutation
                
                split_from_to_allele = from_to_allele.split('.')
                for allele in split_from_to_allele: #to consider cases such as: A>C.C>T|A frequency=4 then A>C is written twice while C>T and C>A are writen once each
                    if allele.startswith("del"):
                        if len(allele)>5:
                            allele = "multi_del"
                    elif allele.startswith("ins"):
                        if len(allele)>5: #delAC is ok but for larger deletions just write multi_del
                            allele = "multi_ins"
                            
                    motifs_info_dict[motif_name]["AlleleChanges"][snv_position].append(allele)
                    motifs_info_dict[motif_name]["AlleleChanges"][-1].append(allele)
                
                for sample_id in sample_ids:
                    motifs_info_dict[motif_name]["SampleIDs"][snv_position].append(sample_id)
                motifs_info_dict[motif_name]["MutIds"][snv_position].append(mutation_id)
                motifs_info_dict[motif_name]["MutPositions"][snv_position].append(mutation_genomic_position)
                motifs_info_dict[motif_name]["MotifPositions"][snv_position].append(motif_genomic_position)
                if frequency_of_this_mutation>1:
                    motifs_info_dict[motif_name]["No.RecurrentMuts"][snv_position]+=1
                
                motifs_info_dict[motif_name]["frequency"][-1] += frequency_of_this_mutation
                for sample_id in sample_ids:
                    motifs_info_dict[motif_name]["SampleIDs"][-1].append(sample_id)
                motifs_info_dict[motif_name]["MutIds"][-1].append(mutation_id)
                motifs_info_dict[motif_name]["MutPositions"][-1].append(mutation_genomic_position)
                motifs_info_dict[motif_name]["MotifPositions"][-1].append(motif_genomic_position)
                if frequency_of_this_mutation>1:
                    motifs_info_dict[motif_name]["No.RecurrentMuts"][-1]+=1
                
    #Update the dicutionaries to keep the unique IDs and count number of occurences of MutIds and SampleIds, frequency of allele changes
    for motif in motifs_info_dict.keys():
        for position_index in range(0, len(motifs_info_dict[motif]["SampleIDs"])):
            motifs_info_dict[motif]["SampleIDs"][position_index] = list(set(motifs_info_dict[motif]["SampleIDs"][position_index]))
            motifs_info_dict[motif]["No.UniqueSamples"][position_index] = len(motifs_info_dict[motif]["SampleIDs"][position_index])
        
        for position_index in range(0, len(motifs_info_dict[motif]["MotifPositions"])):
            motifs_info_dict[motif]["MotifPositions"][position_index] = list(set(motifs_info_dict[motif]["MotifPositions"][position_index]))
            motifs_info_dict[motif]["No.UniqueMotifs"][position_index] = len(motifs_info_dict[motif]["MotifPositions"][position_index])
        
        for position_index in range(0, len(motifs_info_dict[motif]["MutIds"])):
            motifs_info_dict[motif]["MutIds"][position_index] = list(set(motifs_info_dict[motif]["MutIds"][position_index]))
            motifs_info_dict[motif]["MutPositions"][position_index] = list(set(motifs_info_dict[motif]["MutPositions"][position_index]))
        
        #count frequency of each allele, replace ----- (indels) to a single -
        for position_index in range(0, len(motifs_info_dict[motif]["AlleleChanges"])):
            motifs_info_dict[motif]["UnifiedAlleleChanges"][position_index] = [key + ":" + str(len(list(group))) for key, group in groupby(sorted(motifs_info_dict[motif]["AlleleChanges"][position_index]))]
            
    sep = "\t"
    for motif_name in motifs_info_dict:
        output_file.write(motif_name + sep + sep.join(str(x) for x in motifs_info_dict[motif_name]["frequency"])+"\n")
        output_file.write("UnifiedAlleleChanges" + sep + re.sub("--*", "-", sep.join(','.join(x) for x in motifs_info_dict[motif_name]["UnifiedAlleleChanges"]))+"\n")
        output_file.write("AlleleChanges" + sep + sep.join(','.join(x) for x in motifs_info_dict[motif_name]["AlleleChanges"])+"\n")
        output_file.write("SampleIDs" + sep + sep.join(','.join(x) for x in motifs_info_dict[motif_name]["SampleIDs"])+"\n")
        output_file.write("No.UniqueSamples" + sep + sep.join(str(x) for x in motifs_info_dict[motif_name]["No.UniqueSamples"])+"\n")
        output_file.write("No.RecurrentMuts" + sep + sep.join(str(x) for x in motifs_info_dict[motif_name]["No.RecurrentMuts"])+"\n")
        output_file.write("MutIds" + sep + sep.join(','.join(x) for x in motifs_info_dict[motif_name]["MutIds"])+"\n")
        output_file.write("MutPositions" + sep + sep.join(','.join(x) for x in motifs_info_dict[motif_name]["MutPositions"])+"\n")
        output_file.write("MotifPositions" + sep + sep.join(','.join(x) for x in motifs_info_dict[motif_name]["MotifPositions"])+"\n") 
        output_file.write("No.UniqueMotifs" + sep + sep.join(str(x) for x in motifs_info_dict[motif_name]["No.UniqueMotifs"])+"\n")
    output_file.close()


if __name__ == '__main__':
    
    mutated_elements_file = sys.argv[1]
    mutations_file = sys.argv[2]
    
    sample_id_and_number_of_mutations_per_sample_dict = get_number_of_mutations_per_sample_list_and_write_to_file(mutations_file)
    number_of_elements_tested = file_len(mutated_elements_file)
    calculate_p_value_motifregions(mutated_elements_file, sample_id_and_number_of_mutations_per_sample_dict, mutated_regions_pval_outfile="", index_mutation_frequency=6, index_sample_ids=7, index_elment_start_coordinate=1, index_elment_stop_coordinate=2, genome_size=3000000000.0, total_number_tested_regions=number_of_elements_tested)
    
    if len(sys.argv) != 3:
        print("Usage: motif_PFM_input_file muts_mutated_motifs_input_file (muts-9col and motifs-6col")
        exit()
        
    '''motif_PFM_input_file="motifs_converted_tab_to_space.PFM" #sys.argv[1]
    muts_mutated_motifs_input_file="temp.bed" #sys.argv[2] 
    score_motifs_according_to_their_affect(motif_PFM_input_file=motif_PFM_input_file, muts_mutated_motifs_input_file=muts_mutated_motifs_input_file, mutated_motifs_scored_output_file="", index_mutations_info=0, index_mutated_motif_info=9, index_fromAllele = 3, index_toAllele = 4)
    '''
    #score_motifs_according_to_their_affect(motif_PFM_input_file=sys.argv[1], mutated_motifs_input_file=sys.argv[2], mutated_motifs_scored_output_file="", mutations_output_file="", index_mutations_info=-1, mutation_info_sep=',', mutations_sep=';', motif_breakness_diff_freq_threshold=0.4, index_fromToAllele = 6, index_mutationFrequency=5, index_sampleID = 8)
    