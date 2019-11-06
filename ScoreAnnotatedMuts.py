'''
Created on May 18, 2016

@author: Husen M. Umer
'''
import pandas as pd
import statsmodels.api as sm
import numpy as np
#import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
#from patsy.contrasts import Treatment

'''
Given a list of annotated mutations whose scores are known find the weight for each annotation.
input: 
- a tab separated file containing the scores for each mutation and a score attribute. A header line should exist including names of the annotations
- a list containing names of annotaton to use in the regression model
process: run a linear regression model to find the coefficient for each of the annotations
output: return a dict containing the weights for each annotation
'''

def getWeightsforAnnotations(annoated_scored_muts_input_file, outcome_col_name="", outcome_col_index=0, list_of_cols_names=[], list_of_cols_index=[], take_abs_entropy_diff=True, log_base=10):
    
    all_names_cols_list = [outcome_col_name]
    all_names_cols_list.extend(list_of_cols_names)
    all_index_cols_list = [outcome_col_index]
    all_index_cols_list.extend(list_of_cols_index)
    pcsv = pd.read_csv(annoated_scored_muts_input_file, header=0, sep='\t', na_filter=True, na_values=["NA", "None"], usecols=all_index_cols_list, names=all_names_cols_list)#list the columns to be read from hear and the decision
    
    df = pd.DataFrame(pcsv).dropna()
    if take_abs_entropy_diff:
        if 'Entropy_diff' in list_of_cols_names:
            df['Entropy_diff'] = abs(df['Entropy_diff'])
    df[outcome_col_name] = np.where(df[outcome_col_name]==True, 1, 0)
    if 'normal_expression_level_of_this_motif' in list_of_cols_names:
        df['normal_expression_level_of_this_motif']= np.where(df['normal_expression_level_of_this_motif']==0, 1, df['normal_expression_level_of_this_motif'])
        df['normal_expression_level_of_this_motif'] = np.log(np.array(df['normal_expression_level_of_this_motif']))/np.log(log_base)

    if 'other_factors' in list_of_cols_names:
        df['other_factors']= np.where(df['other_factors']==0, 1, df['other_factors'])
        df['other_factors'] = np.log(np.array(df['other_factors']))/np.log(log_base)
    
    #plot_corr(df)
    if 'CAGE_expr' in list_of_cols_names:
        df['CAGE_expr']= np.where(df['CAGE_expr']==0, 1, df['CAGE_expr'])
        df['CAGE_expr'] = np.log(np.array(df['CAGE_expr']))/np.log(log_base)
    
    if 'chromatin_states' in list_of_cols_names:
        df = pd.get_dummies(df, columns=['chromatin_states'], sparse=True)
        list_of_cols_names = list(df.columns)
        del list_of_cols_names[list_of_cols_names.index(outcome_col_name)]

    if 'replicationdomain' in list_of_cols_names:
        df = pd.get_dummies(df, columns=['replicationdomain'], sparse=True)
        list_of_cols_names = list(df.columns)
        del list_of_cols_names[list_of_cols_names.index(outcome_col_name)]
    df['intercept']  = 1
    list_of_cols_names.append('intercept')
    X = df[list_of_cols_names]
    y  = df[outcome_col_name]
    
    #error here
    logit = sm.Logit(y,X).fit(method='bfgs', full_output=True)#, maxiter=1000000000
    print logit.summary()
    
    for k in X.columns:
        kk = X[k]
        print y.corr(kk)

    print "log odds ratio is reported: logit.params:"
    return logit.params#think of negative values (excluded) ((weightXi[without exp, just use the result from here logit.params]*valueXi)+(weightXi+1*valueXi+1))
    #return np.exp(logit.params)# think of values less than one (should be excluded, may be). (weightXi^valueXi)*weightXi+1*exp(valueXi+1)))

def plot_corr(df,size=10):
    '''Function plots a graphical correlation matrix for each pair of columns in the dataframe.

    Input:
        df: pandas DataFrame
        size: vertical and horizontal size of the plot'''

    corr = df.corr()
    fig, ax = plt.subplots(figsize=(size, size))
    sm = ax.matshow(corr)
    cbar = plt.colorbar(sm, alpha=0.05, aspect=16, shrink=0.4)
    cbar.solids.set_edgecolor("face")

    plt.xticks(range(len(corr.columns)), corr.columns);
    plt.yticks(range(len(corr.columns)), corr.columns);
    
def get_col_names():
    col_names =  ['MutChr', 'MutStart', 'MutStop', 'Ref', 'Alt', 'Cancer-project', 'Mut_score', 'vcf_info', 'sample_file_id', 
                  'MotifChr', 'MotifStart', 'MotifEnd', 'Motif_name', 'Motif_score', 'Motif_strand', 'Entropy_diff', 'MutMotif_position',                        
                  'normal_expression_level_of_this_motif', 'tumor_expression_level_of_this_motif',]
    col_names.extend(["same_factor", "other_factors", "dnase", "contactingdomain", "loopdomain", "CAGE_expr", "chromatin_states", "replicationdomain"])
    col_names.extend(["same_cell_same_factor","same_cell_dnase","same_cell_chromatin_states","same_cell_replicationdomain","same_cell_contactingdomain","same_cell_loopdomain","same_cell_CAGE_expr","same_cell_other_factors_names","same_cell_other_factors", "other_cells_same_factor","other_cells_dnase","other_cell_chromatin_states","other_cells_replicationdomain","other_cells_contactingdomain","other_cells_loopdomain","other_cells_CAGE_expr"])
    col_names.append('WAScore')
    
    return col_names


if __name__ == '__main__':
    outcome_col_name = 'Mut_score'
    outcome_col_index = 6
    col_names = get_col_names()
    col_names_to_weight = ['Entropy_diff', 'normal_expression_level_of_this_motif',  
                                'same_factor_overlaps', 'all_factors_same_cell_overlaps', 'dnase_overlap', 'contactingdomain_overlap', 'loopdomain_overlap', 'cage_overlap']#, 'chromhmm_overlap_to_report', 'replicationtiming_label_overlap_to_report']
    col_indexes_to_weight = []
    for c in col_names_to_weight:
        if c in col_names:
            col_indexes_to_weight.append(col_names.index(c))
        else:
            del col_names_to_weight[c]
    outcome_col_name = 'Mut_score'
    outcome_col_index = 6
    take_abs_entropy_diff = False
    log_base =10
    k = getWeightsforAnnotations('~/Documents/Group-Projects/PCAWG/regDriverTest/MauranoTrainingAnnotatedNoDiscMotifs', 
                             outcome_col_name = outcome_col_name,
                             outcome_col_index = outcome_col_index, 
                             list_of_cols_names = col_names_to_weight,
                             list_of_cols_index = col_indexes_to_weight,
                             take_abs_entropy_diff=take_abs_entropy_diff, log_base=log_base)
    k.name = "ORs"
    for p in range(0, len(k)):
        print k.index[p] + "=" + str(k[p])
    print k.index[0]
