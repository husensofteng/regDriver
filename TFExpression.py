'''
Created on May 3, 2016

@author: Husen M. Umer
'''
import os
import numpy as np

def parse_GTEx_metadafile(sample_info_input_file, sample_name_index=0, tissue_type_index=5, index_UseME_col_indicator=14):
    
    dict_tissue_type_sampleIDs = {}
    with open(sample_info_input_file, 'r') as sample_info_infile:
        sample_info_lines = sample_info_infile.readlines()
        for line in sample_info_lines:
            split_line = line.strip().split('\t')
            if len(split_line)>index_UseME_col_indicator and len(split_line)>tissue_type_index and len(split_line)>index_UseME_col_indicator>=sample_name_index:
                if split_line[tissue_type_index]!="" and split_line[index_UseME_col_indicator]=="USE ME":
                    if split_line[tissue_type_index] not in dict_tissue_type_sampleIDs:
                        dict_tissue_type_sampleIDs[split_line[tissue_type_index]] = []
                    dict_tissue_type_sampleIDs[split_line[tissue_type_index]].append(split_line[sample_name_index])    
    return dict_tissue_type_sampleIDs

'''
@desc For each origin#Gene in list_of_cell_tracks report a list of expression values (all samples) detected in the corresponding tissue/tumor-types (aka. origin) 
Input params: 
a list containing names TFs (TF) for which expression levels to be reported,
dict with origin types as key and sample ids as values from the gene expression file. Indicates which samples should be used for each originType
expression file path
normalization to be applied for the gene values from the expression file (eg. log2),
gene name index in in the expression file, 
sample columns start index in the expression file, 
sample columns end index in the expression file, 
the summary statistic to be applied (mean, median) for each sample set,
Output:
a dictionary with each Tissue#TF(Gene) as a key and a list of expression levels from the samples belonging to that tissue or tumor type as values   
a file with the contents of the dictionary, if such file exisits the dictionary is re-bulit based on it to avoid re-processing the whole expression file next time
'''
def get_expression_level_per_originType_per_TF(list_of_tf_names, dict_tissue_type_sampleIDs, expression_inputfile,
                                     normalization_method_exp='NA',
                                     gene_name_index_expr=1,
                                     sample_cols_start_index_expr=2, sample_cols_end_index_expr=8554, 
                                     num_header_lines_expr=2, line_number_containing_sampleIDs_expr=3):
    
    expression_output_file = expression_inputfile + "_perTissue_perTF.txt"
    origin_gene_expression_values = {} #origin#gene as keys (origin can be tissue type): gene expression levels for the samples (genes in rows and samples in cols), #limit this list to contain only for genes given in the cell#TF list  
    
    if os.path.exists(expression_output_file):
        with open(expression_output_file, 'r') as expression_outfile:
            lines = expression_outfile.readlines()
            for line in lines:
                if len(line.split('\t'))==2:
                    origin_gene_expression_values[line.strip().split('\t')[0]] = []
                    for x in line.strip().split('\t')[1].split(','):
                        origin_gene_expression_values[line.strip().split('\t')[0]].append(float(x))
                else:
                    print(expression_output_file + " should have a correct format (originType#TFName \t GeneValue1,GeneValue2,GeneValueN") 
    else:
        sample_IDs = []
        TFs_extract_expression = []
        for tf_name in list_of_tf_names:
            if tf_name not in TFs_extract_expression:
                TFs_extract_expression.append(tf_name.split('(')[0])
                #add TF names after removing dashes and splitting them by :: (in case they were in the name)
                if '-' in tf_name:#in case - was in the name remove it
                    TFs_extract_expression.append(tf_name.replace('-',''))
                if '::' in tf_name:
                    s_tf_name = tf_name.split('::')
                    for s in s_tf_name:
                        if s not in TFs_extract_expression:
                            TFs_extract_expression.append(s)
        
        with open(expression_inputfile, 'r') as expression_infile:
            expression_line = expression_infile.readline()
            line_number_counter = 1
            print("Processing gene expression values for TFs from:" + expression_inputfile)
            while expression_line!="":
                #ignore info header lines
                if line_number_counter <= num_header_lines_expr:
                    expression_line = expression_infile.readline()
                    line_number_counter+=1
                    continue
                
                split_expression_line = expression_line.strip().split('\t')
                #read the sampleIDs from the header
                if line_number_counter == line_number_containing_sampleIDs_expr:
                    sample_IDs = split_expression_line[sample_cols_start_index_expr:sample_cols_end_index_expr+1]
                    expression_line = expression_infile.readline()
                    line_number_counter+=1
                    continue
                
                #start reading the gene expression values
                if split_expression_line[gene_name_index_expr].upper() in TFs_extract_expression:
                    for tissue_type in dict_tissue_type_sampleIDs.keys():
                        if tissue_type+"#"+split_expression_line[gene_name_index_expr].upper() not in origin_gene_expression_values:
                            origin_gene_expression_values[tissue_type+"#"+split_expression_line[gene_name_index_expr].upper()] = []
                            
                            for sampleID in dict_tissue_type_sampleIDs[tissue_type]:
                                if sampleID in sample_IDs:
                                    if normalization_method_exp=='NA':
                                        origin_gene_expression_values[tissue_type+"#"+split_expression_line[gene_name_index_expr]].append(
                                            float(split_expression_line[sample_IDs.index(sampleID)+sample_cols_start_index_expr])) #+sample_cols_start_index_expr because the cols at the begining of the line are gene names and gene info
                                    
                                    elif normalization_method_exp=='log2':
                                        if normalization_method_exp=='NA':
                                            try:
                                                origin_gene_expression_values[tissue_type+"#"+split_expression_line[gene_name_index_expr].upper()].append(
                                                    np.log2(float(split_expression_line[sample_IDs.index(sampleID)+sample_cols_start_index_expr]))) #+sample_cols_start_index_expr because the cols at the begining of the line are gene names and gene info
                                            except ValueError:
                                                origin_gene_expression_values[tissue_type+"#"+split_expression_line[gene_name_index_expr].upper()].append(
                                                    float(split_expression_line[sample_IDs.index(sampleID)+sample_cols_start_index_expr])) #+sample_cols_start_index_expr because the cols at the begining of the line are gene names and gene info
                                            
                line_number_counter+=1            
                expression_line = expression_infile.readline()
        
        #write the results to avoid processing the whole file for next time
        with open(expression_output_file, 'w') as expression_outfile:
            for originType in origin_gene_expression_values.keys():
                if len(origin_gene_expression_values[originType])>0:#only write those tissue#tf that have values from at least one sample
                    expression_outfile.write(originType + "\t" + ','.join('%.3f' % s for s in origin_gene_expression_values[originType]) + "\n")
    
    return origin_gene_expression_values

'''
@desc For each TF in list_of_tf_names report a summary of expression values detected in the corresponding tissue(s)
Input params: 
a list containing names of TFs for which expression levels to be reported,
a dict mapping cell names to tissue names that are used in the expression file #to know from which tissues obtain the expression values,
normalization to be applied for the gene values from the expression file (eg. log2),
the summary statistic to be applied (mean, median) for each sample set
Output:
a dictionary with each cell#TF as a key and normal expression level and tumor expression level as values   
'''

def get_expression_level_per_cell_per_TF(list_of_tf_names, origin_gene_expression_values, originType_cells_mapping_dict,
                                         origin_gene_expression_values_outputfile, summary_statistics_to_apply='mean'):
    
    expression_per_cell_per_TF = {}
    if not os.path.exists(origin_gene_expression_values_outputfile):
        #return NA for the cell#TFs that have no expression levels in the input file
        for tf_name in list_of_tf_names:
            #collect the expression of this gene from all the tissues that this cell belongs to
            for originType in originType_cells_mapping_dict.keys():
                for cell in originType_cells_mapping_dict[originType]:
                    cell_tf = cell+'#'+tf_name
                    if cell_tf not in expression_per_cell_per_TF.keys():
                        if originType+'#'+tf_name in origin_gene_expression_values.keys():
                            expression_per_cell_per_TF[cell_tf] = origin_gene_expression_values[originType+'#'+tf_name]
                            if summary_statistics_to_apply == 'mean':
                                expression_per_cell_per_TF[cell_tf] = np.mean(expression_per_cell_per_TF[cell_tf])
                            elif summary_statistics_to_apply == 'median':
                                expression_per_cell_per_TF[cell_tf] = np.median(expression_per_cell_per_TF[cell_tf])
                        
                        elif originType+'#'+tf_name.replace('-', '') in origin_gene_expression_values.keys():
                            expression_per_cell_per_TF[cell_tf] = origin_gene_expression_values[originType+'#'+tf_name.replace('-', '')]
                            if summary_statistics_to_apply == 'mean':
                                expression_per_cell_per_TF[cell_tf] = np.mean(expression_per_cell_per_TF[cell_tf])
                            elif summary_statistics_to_apply == 'median':
                                expression_per_cell_per_TF[cell_tf] = np.median(expression_per_cell_per_TF[cell_tf])
                        
                        elif '::' in tf_name:
                            mean_complexes = 0.0
                            for s in tf_name.split('::'):
                                if s in origin_gene_expression_values.keys():
                                    value_s = origin_gene_expression_values[originType+'#'+s]
                                    if summary_statistics_to_apply == 'mean':
                                        value_s = np.mean(value_s)
                                    elif summary_statistics_to_apply == 'median':
                                        value_s = np.median(value_s = np.mean(value_s))
                                    mean_complexes+=value_s
                            expression_per_cell_per_TF[cell_tf] = mean_complexes
                
        with open(origin_gene_expression_values_outputfile, 'w') as origin_gene_expression_values_outfile:
            for cell_tf in expression_per_cell_per_TF.keys():
                    origin_gene_expression_values_outfile.write(cell_tf + '\t' + str(expression_per_cell_per_TF[cell_tf])+'\n')
    else:
        with open(origin_gene_expression_values_outputfile, 'r') as origin_gene_expression_values_infile:
            lines = origin_gene_expression_values_infile.readlines()
            for l in lines:
                sl = l.strip().split('\t')
                if sl[1]=="NA":
                    expression_per_cell_per_TF[sl[0].strip()] = sl[1]
                else:
                    expression_per_cell_per_TF[sl[0].strip()] = float(sl[1])
    return expression_per_cell_per_TF

