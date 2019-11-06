'''
Created on May 3, 2016

@author: Husen M. Umer
'''
import os,sys
import numpy as np
#import regDriver #for testing purpose in this module only -- can safely be removed after removing the main function of this module

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
                    print expression_output_file + " should have a correct format (originType#TFName \t GeneValue1,GeneValue2,GeneValueN" 
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
            print "Processing gene expression values for TFs from:" + expression_inputfile
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

if __name__ == '__main__':
    sample_info_input_file = sys.argv[1]
    expression_inputfile = sys.argv[2]
    TissueCellInfo_matches_inputfile = sys.argv[3]
    list_of_cell_tracks = ['A549#ATF3','A549#BCL3','A549#BHLHE40','A549#CEBPB','A549#CREB1','A549#CTCF','A549#E2F6','A549#ELF1','A549#EP300','A549#ETS1','A549#FOSL2','A549#FOXA1','A549#FOXA2','A549#GABPA','A549#GATA3','A549#JUND','A549#MAX','A549#MYC','A549#NR3C1','A549#PBX3','A549#POLR2A','A549#POLR2AphosphoS2','A549#RAD21','A549#REST','A549#SIN3A','A549#SIX5','A549#SP1','A549#SREBF2','A549#TAF1','A549#TCF12','A549#TEAD4','A549#USF1','A549#YY1','A549#ZBTB33','GM12878#ATF2','GM12878#ATF3','GM12878#BATF','GM12878#BCL11A','GM12878#BCL3','GM12878#BCLAF1','GM12878#BHLHE40','GM12878#BRCA1','GM12878#CEBPB','GM12878#CEBPZ','GM12878#CHD1','GM12878#CHD2','GM12878#CREB1','GM12878#CREM','GM12878#CTCF','GM12878#CUX1','GM12878#E2F4','GM12878#EBF1','GM12878#EGR1','GM12878#ELF1','GM12878#ELK1','GM12878#EP300','GM12878#ESRRA','GM12878#ETS1','GM12878#ETV6','GM12878#EZH2','GM12878#FOS','GM12878#FOXM1','GM12878#GABPA','GM12878#IKZF1','GM12878#IRF3','GM12878#IRF4','GM12878#JUND','GM12878#KAT2A','GM12878#MAFK','GM12878#MAX','GM12878#MAZ','GM12878#MEF2A','GM12878#MEF2C','GM12878#MTA3','GM12878#MXI1','GM12878#MYC','GM12878#NFATC1','GM12878#NFE2','GM12878#NFIC','GM12878#NFYA','GM12878#NFYB','GM12878#NR2C2','GM12878#NRF1','GM12878#PAX5','GM12878#PBX3','GM12878#PML','GM12878#POLR2A','GM12878#POLR2AphosphoS2','GM12878#POLR2AphosphoS5','GM12878#POLR3G','GM12878#POU2F2','GM12878#RAD21','GM12878#RCOR1','GM12878#RELA','GM12878#REST','GM12878#RFX5','GM12878#RUNX3','GM12878#RXRA','GM12878#SIN3A','GM12878#SIX5','GM12878#SMAD5','GM12878#SMC3','GM12878#SP1','GM12878#SPI1','GM12878#SREBF1','GM12878#SREBF2','GM12878#SRF','GM12878#STAT1','GM12878#STAT3','GM12878#STAT5A','GM12878#SUPT20H','GM12878#TAF1','GM12878#TBL1XR1','GM12878#TBP','GM12878#TCF12','GM12878#TCF3','GM12878#TCF7','GM12878#USF1','GM12878#USF2','GM12878#WRNIP1','GM12878#YY1','GM12878#ZBED1','GM12878#ZBTB33','GM12878#ZEB1','GM12878#ZNF143','GM12878#ZNF274','GM12878#ZNF384','GM12878#ZZZ3','HCT116#ATF3','HCT116#CBX3','HCT116#CEBPB','HCT116#CTCF','HCT116#EGR1','HCT116#ELF1','HCT116#FOSL1','HCT116#JUND','HCT116#MAX','HCT116#POLR2A','HCT116#POLR2AphosphoS5','HCT116#RAD21','HCT116#REST','HCT116#SIN3A','HCT116#SP1','HCT116#SRF','HCT116#TCF7L2','HCT116#TEAD4','HCT116#USF1','HCT116#YY1','HCT116#ZBTB33','HEK293#CTCF','HEK293#ELK4','HEK293#POLR2A','HEK293#TCF7L2','HEK293#TRIM28','HEK293#ZNF263','HeLa-S3#BDP1','HeLa-S3#BRCA1','HeLa-S3#BRF1','HeLa-S3#BRF2','HeLa-S3#CEBPB','HeLa-S3#CHD1','HeLa-S3#CHD2','HeLa-S3#CTCF','HeLa-S3#E2F1','HeLa-S3#E2F4','HeLa-S3#E2F6','HeLa-S3#ELK1','HeLa-S3#ELK4','HeLa-S3#EP300','HeLa-S3#EZH2','HeLa-S3#FOS','HeLa-S3#GABPA','HeLa-S3#GTF2F1','HeLa-S3#GTF3C2','HeLa-S3#HCFC1','HeLa-S3#IRF3','HeLa-S3#JUN','HeLa-S3#JUND','HeLa-S3#KAT2A','HeLa-S3#MAFK','HeLa-S3#MAX','HeLa-S3#MAZ','HeLa-S3#MXI1','HeLa-S3#MYC','HeLa-S3#NFYA','HeLa-S3#NFYB','HeLa-S3#NR2C2','HeLa-S3#NRF1','HeLa-S3#POLR2A','HeLa-S3#POLR2AphosphoS2','HeLa-S3#POLR3A','HeLa-S3#PRDM1','HeLa-S3#RAD21','HeLa-S3#RCOR1','HeLa-S3#REST','HeLa-S3#RFX5','HeLa-S3#SMARCA4','HeLa-S3#SMARCB1','HeLa-S3#SMARCC1','HeLa-S3#SMARCC2','HeLa-S3#SMC3','HeLa-S3#SREBF2','HeLa-S3#STAT1','HeLa-S3#STAT3','HeLa-S3#SUPT20H','HeLa-S3#TAF1','HeLa-S3#TBP','HeLa-S3#TCF7L2','HeLa-S3#USF2','HeLa-S3#ZKSCAN1','HeLa-S3#ZNF143','HeLa-S3#ZNF274','HeLa-S3#ZZZ3','HepG2#ARID3A','HepG2#ATF3','HepG2#BHLHE40','HepG2#BRCA1','HepG2#CBX1','HepG2#CEBPB','HepG2#CEBPD','HepG2#CEBPZ','HepG2#CHD2','HepG2#CREB1','HepG2#CTCF','HepG2#ELF1','HepG2#EP300','HepG2#ESRRA','HepG2#EZH2','HepG2#FOSL2','HepG2#FOXA1','HepG2#FOXA2','HepG2#GABPA','HepG2#HCFC1','HepG2#HDAC2','HepG2#HNF4A','HepG2#HNF4G','HepG2#HSF1','HepG2#IRF3','HepG2#JUN','HepG2#JUND','HepG2#MAFF','HepG2#MAFK','HepG2#MAX','HepG2#MAZ','HepG2#MBD4','HepG2#MXI1','HepG2#MYBL2','HepG2#MYC','HepG2#NFIC','HepG2#NR2C2','HepG2#NR2F2','HepG2#NR3C1','HepG2#NRF1','HepG2#POLR2A','HepG2#POLR2AphosphoS2','HepG2#POLR2AphosphoS5','HepG2#PPARGC1A','HepG2#RAD21','HepG2#RCOR1','HepG2#REST','HepG2#RFX5','HepG2#RXRA','HepG2#SIN3A','HepG2#SMC3','HepG2#SP1','HepG2#SP2','HepG2#SREBF1','HepG2#SREBF2','HepG2#SRF','HepG2#TAF1','HepG2#TBP','HepG2#TCF12','HepG2#TCF7L2','HepG2#TEAD4','HepG2#TFAP4','HepG2#USF1','HepG2#USF2','HepG2#YY1','HepG2#ZBTB33','HepG2#ZBTB7A','HepG2#ZEB1','HepG2#ZHX2','HepG2#ZKSCAN1','HepG2#ZNF274','IMR-90#CEBPB','IMR-90#CHD1','IMR-90#CTCF','IMR-90#MAFK','IMR-90#MAZ','IMR-90#MXI1','IMR-90#POLR2A','IMR-90#RAD21','IMR-90#RCOR1','IMR-90#RFX5','Ishikawa#CEBPB','Ishikawa#CREB1','Ishikawa#CTCF','Ishikawa#EGR1','Ishikawa#EP300','Ishikawa#ESR1','Ishikawa#FOXA1','Ishikawa#FOXM1','Ishikawa#MAX','Ishikawa#NFIC','Ishikawa#NR3C1','Ishikawa#POLR2A','Ishikawa#RAD21','Ishikawa#REST','Ishikawa#SRF','Ishikawa#TAF1','Ishikawa#TCF12','Ishikawa#TEAD4','Ishikawa#USF1','Ishikawa#YY1','Ishikawa#ZBTB7A','K562#ARID3A','K562#ATF1','K562#ATF3','K562#BACH1','K562#BCLAF1','K562#BDP1','K562#BHLHE40','K562#BRF1','K562#BRF2','K562#CBX1','K562#CBX2','K562#CBX3','K562#CBX5','K562#CBX8','K562#CCNT2','K562#CEBPB','K562#CEBPD','K562#CEBPZ','K562#CHD1','K562#CHD2','K562#CHD4','K562#CHD7','K562#CREB1','K562#CREB3','K562#CREBBP','K562#CREM','K562#CTCF','K562#CTCFL','K562#CUX1','K562#DDX20','K562#DIDO1','K562#E2F4','K562#E2F6','K562#EGR1','K562#ELF1','K562#ELK1','K562#EP300','K562#ETS1','K562#ETV1','K562#ETV6','K562#EZH2','K562#FOS','K562#FOSL1','K562#GABPA','K562#GATA1','K562#GATA2','K562#GTF2B','K562#GTF2F1','K562#GTF3C2','K562#HCFC1','K562#HDAC1','K562#HDAC2','K562#HDAC6','K562#HDAC8','K562#HINFP','K562#HMGN3','K562#ID3','K562#ILK','K562#IRF1','K562#IRF9','K562#JUN','K562#JUNB','K562#JUND','K562#KAT2B','K562#KDM1A','K562#KDM5B','K562#KLF1','K562#KLF13','K562#MAFF','K562#MAFG','K562#MAFK','K562#MAX','K562#MAZ','K562#MEF2A','K562#MITF','K562#MXI1','K562#MYC','K562#NCOR1','K562#NFE2','K562#NFE2L1','K562#NFYA','K562#NFYB','K562#NR2C2','K562#NR2F2','K562#NR4A1','K562#NRF1','K562#PHF8','K562#PML','K562#POLR2A','K562#POLR2AphosphoS2','K562#POLR2AphosphoS5','K562#POLR3A','K562#POLR3G','K562#PTRF','K562#PTTG1','K562#PYGO2','K562#RAD21','K562#RBBP5','K562#RCOR1','K562#RELA','K562#REST','K562#RFX5','K562#RNF2','K562#SAP30','K562#SETDB1','K562#SIN3A','K562#SIRT6','K562#SIX5','K562#SMAD5','K562#SMARCA4','K562#SMARCB1','K562#SMC3','K562#SP1','K562#SP2','K562#SPI1','K562#SRF','K562#STAT1','K562#STAT2','K562#STAT5A','K562#SUZ12','K562#TAF1','K562#TAF7','K562#TAL1','K562#TBL1XR1','K562#TBP','K562#TCF7','K562#TEAD2','K562#TEAD4','K562#TFDP1','K562#THAP1','K562#TRIM28','K562#TSC22D4','K562#UBTF','K562#USF1','K562#USF2','K562#WHSC1','K562#XRCC4','K562#YY1','K562#ZBED1','K562#ZBTB33','K562#ZBTB7A','K562#ZC3H11A','K562#ZKSCAN1','K562#ZMIZ1','K562#ZNF143','K562#ZNF24','K562#ZNF263','K562#ZNF274','K562#ZNF384','LNCaP clone FGC#CTCF','MCF 10A#E2F4','MCF 10A#FOS','MCF 10A#MYC','MCF 10A#POLR2A','MCF 10A#STAT3','MCF-7#CEBPB','MCF-7#CTCF','MCF-7#E2F1','MCF-7#EGR1','MCF-7#ELF1','MCF-7#ELK1','MCF-7#EP300','MCF-7#FOSL2','MCF-7#FOXM1','MCF-7#GABPA','MCF-7#GATA3','MCF-7#HCFC1','MCF-7#HDAC2','MCF-7#JUND','MCF-7#MAX','MCF-7#MYC','MCF-7#NR2F2','MCF-7#PML','MCF-7#POLR2A','MCF-7#RAD21','MCF-7#REST','MCF-7#SIN3A','MCF-7#SRF','MCF-7#TAF1','MCF-7#TCF12','MCF-7#TCF7L2','MCF-7#TEAD4','MCF-7#ZNF217','Panc1#POLR2AphosphoS5','Panc1#REST','Panc1#SIN3A','Panc1#TCF7L2','SK-N-SH#CTCF','SK-N-SH#ELF1','SK-N-SH#EP300','SK-N-SH#FOSL2','SK-N-SH#FOXM1','SK-N-SH#GABPA','SK-N-SH#GATA3','SK-N-SH#JUND','SK-N-SH#MAX','SK-N-SH#MEF2A','SK-N-SH#MXI1','SK-N-SH#NFIC','SK-N-SH#NRF1','SK-N-SH#PBX3','SK-N-SH#POLR2AphosphoS5','SK-N-SH#RAD21','SK-N-SH#REST','SK-N-SH#RFX5','SK-N-SH#RXRA','SK-N-SH#SIN3A','SK-N-SH#SMC3','SK-N-SH#TAF1','SK-N-SH#TCF12','SK-N-SH#TEAD4','SK-N-SH#USF1','SK-N-SH#YY1','SK-N-SH#ZBTB33','T47D#CTCF','T47D#EP300','T47D#ESR1','T47D#FOXA1','T47D#GATA3','T47D#JUND','U2OS#SETDB1','U2OS#TRIM28','astrocyte#CTCF','astrocyte#EZH2','epithelial cell of esophagus#CTCF','keratinocyte#CTCF','keratinocyte#EZH2','osteoblast#CTCF','osteoblast#EP300']
    
    print "Reading metadafile"
    dict_tissue_type_sampleIDs = parse_GTEx_metadafile(sample_info_input_file)
    print "Getting expression level per originType#TF"
    origin_gene_expression_values = get_expression_level_per_originType_per_TF(list_of_cell_tracks, dict_tissue_type_sampleIDs, expression_inputfile)
    print "Mapping cells and originTypes" 
    originType_cells_mapping_dict = regDriver.map_cellNames_to_originTypes(TissueCellInfo_matches_inputfile)
    print "Getting expression per cell#TF"
    expression_per_cell_per_TF = get_expression_level_per_cell_per_TF(list_of_cell_tracks, origin_gene_expression_values, originType_cells_mapping_dict)
    print expression_per_cell_per_TF
    