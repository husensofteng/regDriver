'''
Created on Jun 11, 2016

@author: Husen M. Umer
'''
import os,sys
import pandas as pd
import statsmodels.api as sm
import numpy as np
from scipy import stats

def get_project_codes_metadata_file(metadata_input_file, index_tumor_type=1, index_matching_code_from_snv_file=25, index_donor_id=3, sep='\t'):
    #index_donor_id is not used 
    dict_code_tumor_type = {}
    with open(metadata_input_file, 'r') as metadata_infile:
        inline = metadata_infile.readline().strip().split(sep)
        while len(inline)>1:
            codes = inline[index_matching_code_from_snv_file].split(',')
            for c in codes:
                if c not in dict_code_tumor_type:
                    dict_code_tumor_type[c] = inline[index_tumor_type]
            inline = metadata_infile.readline().strip().split(sep)
            
    return dict_code_tumor_type

def mege_files(input_dir, merged_output_flie_name, dict_code_tumor_type, index_filter=6, index_info=7, sep='\t'):
    
    list_dir = os.listdir(input_dir)
    with open(merged_output_flie_name, 'w') as merged_outflie:
        for infile in list_dir:
            if infile.split('.')[-1]=='vcf':
                print "merging: " + input_dir+'/'+infile
                tumor_type = dict_code_tumor_type[infile.split('.')[0]]
                with open(input_dir+'/'+infile, 'r') as input_datafile:
                    inline = input_datafile.readline()
                    while inline!="":
                        split_inline = inline.strip().split(sep)
                        write_line  = True
                        if not inline.startswith('#'):
                            split_inline = inline.strip().split(sep)
                            if split_inline[index_filter]=="LOWSUPPORT":
                                write_line = False
                            elif "1000genomes_AF" in split_inline[index_info]:
                                info_cols = split_inline[index_info].split(';')
                                for col in info_cols:
                                    if col.startswith('1000genomes_AF'):
                                        scol = col.split('=')[1].split(',')
                                        for c in scol:
                                            if float(c)<0.01:
                                                write_line = False
                            if write_line:
                                out_line = "chr" + split_inline[0] + sep + split_inline[1] + sep + split_inline[1] + sep + split_inline[3] + sep + split_inline[4] + sep + tumor_type  + sep + split_inline[5] + sep + split_inline[7] + sep + infile.split('.')[0] + '\n'
                                merged_outflie.write(out_line)
                        inline = input_datafile.readline()
    return merged_output_flie_name

def rank_scores(input_file, scores_index):
    
    pcsv = pd.read_csv(input_file, header=0, sep='\t', names = ['proj', 'scores'])
    df = pd.DataFrame(pcsv)
    with open('tout', 'w') as tout_write:
        for k in range(0, len(df['scores'])):
            ranks = stats.rankdata(df['scores'][k].split(','), "average")/len(df['scores'][k].split(','))
            tout_write.write(','.join(str(x) for x in ranks.tolist())+'\n')
    return ranks
    
def testf(infile):
    print infile
    f = open(infile, 'r+')
    l = f.readline()
    while l!="":
        k = "ddd" + l
        f.write(k)
        l = f.readline()
    
def get_primary_sample_per_donor(mutations_file_MAF_input_file, mutations_file_MAF_output_file, selected_primary_tumors_inputfile, sample_per_donor_input_file, 
                                 index_tumor_selected_primary_tumors=1, index_donorID_sample_per_donor_input_file = 1, index_tumorID_sample_per_donor_input_file=2,
                                 chr=1,start=2,end=3,ref=7,alt=9,proj=44, other_cols=[16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28], tumor=12, normal=13,donor=45,#index of cols to extract from the MAF file
                                 mut_sep = '\t', other_cols_sep='|'):
    
    sample_per_donor_read_file = open(sample_per_donor_input_file, 'r')
    sample_per_donor_lines = sample_per_donor_read_file.readlines()
    
    selected_primary_tumors_infile = open(selected_primary_tumors_inputfile, 'r')
    selected_primary_tumors_lines = selected_primary_tumors_infile.readlines()
    selected_primary_tumors = []
    for l in selected_primary_tumors_lines:
        if l.strip().split('\t')[index_tumor_selected_primary_tumors] not in selected_primary_tumors: 
            selected_primary_tumors.append(l.strip().split('\t')[index_tumor_selected_primary_tumors])
    donor_tumor_id_dict = {}    
    for l in sample_per_donor_lines:
        if l.strip().split('\t')[index_donorID_sample_per_donor_input_file] not in donor_tumor_id_dict.keys(): 
            donor_tumor_id_dict[l.strip().split('\t')[index_donorID_sample_per_donor_input_file]] = []
        if l.strip().split('\t')[index_tumorID_sample_per_donor_input_file] not in donor_tumor_id_dict[l.strip().split('\t')[index_donorID_sample_per_donor_input_file]]:
            donor_tumor_id_dict[l.strip().split('\t')[index_donorID_sample_per_donor_input_file]].append(l.strip().split('\t')[index_tumorID_sample_per_donor_input_file]) 
    
    for d in donor_tumor_id_dict.keys():
        for t in donor_tumor_id_dict[d]:
            if t in selected_primary_tumors:
                donor_tumor_id_dict[d]=[t]
            break
    print "extracting lines..."
    with open(mutations_file_MAF_input_file, 'r') as mutations_file_MAF_read_file:
        with open(mutations_file_MAF_output_file, 'w') as mutations_file_MAF_outfile:
            mut_line = mutations_file_MAF_read_file.readline()
            while mut_line!="":
                split_mut_line = mut_line.strip().split(mut_sep)
                if split_mut_line[donor] in donor_tumor_id_dict.keys():
                    if split_mut_line[tumor] in donor_tumor_id_dict[split_mut_line[donor]]:
                        mutations_file_MAF_outfile.write("chr"+split_mut_line[chr] + mut_sep + split_mut_line[start] + mut_sep + split_mut_line[end] + mut_sep + split_mut_line[ref] + mut_sep + split_mut_line[alt] + mut_sep + split_mut_line[proj] + mut_sep + 
                                                         other_cols_sep.join(split_mut_line[other_cols[0]::]) + mut_sep + split_mut_line[tumor]+"<>"+split_mut_line[normal] + mut_sep + split_mut_line[donor]  + '\n')
                mut_line = mutations_file_MAF_read_file.readline()
    print "done, check: " + mutations_file_MAF_output_file 
    
    
def find_motifTFname_expTFName_matches(motif_tfnames_input_file, exp_tfnames_input_file):
    motif_tfnames_lines = open(motif_tfnames_input_file, 'r').readlines() 
    exp_tfnames_lines = open(exp_tfnames_input_file, 'r').readlines()
    motif_tfnames = []
    exp_tfnames = []
    exp_tfnames_motif_tfnames_matches = {}
    for l in motif_tfnames_lines:
        if l.strip() not in motif_tfnames:
            motif_tfnames.append(l.strip())
    for l in exp_tfnames_lines:
        if l.strip() not in exp_tfnames:
            sl = l.strip().split('\t')
            exp_tfnames.append(sl)
    for motif_tfname in motif_tfnames:
        for exp_tfname in exp_tfnames:
            if len(motif_tfname)==1:
                if motif_tfname in exp_tfname:
                    if motif_tfname not in exp_tfnames_motif_tfnames_matches.keys():
                        exp_tfnames_motif_tfnames_matches[motif_tfname] = []
                    if '\t'.join(exp_tfname) not in exp_tfnames_motif_tfnames_matches[motif_tfname]:
                            exp_tfnames_motif_tfnames_matches[motif_tfname].append('\t'.join(exp_tfname))
            else:
                for exp_tfname_and_group in exp_tfname:
                    if exp_tfname_and_group.upper().find(motif_tfname.upper())>=0 or motif_tfname.upper().find(exp_tfname_and_group.upper())>=0:
                        if motif_tfname not in exp_tfnames_motif_tfnames_matches.keys():
                            exp_tfnames_motif_tfnames_matches[motif_tfname] = []
                        n=0
                        for exp_tfname_and_group2 in exp_tfname:
                            if n==1:
                                exp_tfname_and_group2+="#G"#to indicate its TF family name and not just TF name in the TF names file the first col is for TF name and second col is for TF family (if exists)
                            if exp_tfname_and_group2 not in exp_tfnames_motif_tfnames_matches[motif_tfname]:#to avoid duplicate writing
                                exp_tfnames_motif_tfnames_matches[motif_tfname].append(exp_tfname_and_group2)
                            n+=1
    
    for exp_tfname in exp_tfnames:
        tf_exists = False
        for k in exp_tfnames_motif_tfnames_matches.keys():
            if exp_tfname[0] in exp_tfnames_motif_tfnames_matches[k]:
                tf_exists = True
                continue
        if not tf_exists:
            exp_tfnames_motif_tfnames_matches[exp_tfname[0]] = ["#E"]
    
    for tfname in motif_tfnames:
        if tfname not in exp_tfnames_motif_tfnames_matches.keys():
            exp_tfnames_motif_tfnames_matches[tfname] = ["#M"]
      
    for k in exp_tfnames_motif_tfnames_matches.keys():
        print k + "\t" + '\t'.join(exp_tfnames_motif_tfnames_matches[k])
              
def get_dict_TFs_per_motif(mapping_file=sys.argv[1]):
    
    mapping_infile = open(mapping_file, 'r')
    mapping_outfile = open(mapping_file+'_TFNameMappingDict', 'w')
    mapping_lines = mapping_infile.readlines()
    mapping_dict = {}
    for l in mapping_lines:
        sl = l.strip().split('\t')
        if len(sl)>=2:# and ('#ME' in sl or '#MG' in sl or '#M' in sl or '#E' in sl):
            if sl[0].upper() not in mapping_dict.keys():
                mapping_dict[sl[0].upper()]=[]
            for s in sl[1::]:
                if s.strip().upper() not in mapping_dict[sl[0].upper()] and s.strip().upper()!="" and s.strip().upper()!="#G":
                    mapping_dict[sl[0].upper()].append(s.strip().upper())
    for k in mapping_dict:
        mapping_outfile.write(k+ '\t' + '\t'.join(mapping_dict[k])+'\n')
        
    return mapping_dict
#takes an input file where TF names from JASPAR and HOCO. are listed 
#returns a list of JASPAr TF names from the given dict (just checks the dict keys) (jaspar_motifs)
#it can also generate a list (tf_names_not_present) for motifs that are present in  HOCO. but not in JASPAR (uncomment the first two conditions) in the last main loop
def compare_resources(mapping_dict, tf_names_input_file):
    tf_names_infile = open(tf_names_input_file, 'r')
    tf_names_lines = tf_names_infile.readlines()
    tf_names = {"JASPAR": [], "HOCOMOCO":[]}
    tf_names_not_present = []
    
    for l in tf_names_lines:
        if l.strip().split('\t')[1]=="JASPAR":
            tf_names["JASPAR"].append(l.strip().split('\t')[0].upper())
        elif l.strip().split('\t')[1]=="HOCOMOCO":
            tf_names["HOCOMOCO"].append(l.strip().split('\t')[0].upper())
        else:
            print "Noneee"
    jaspar_motifs=[]
    for k in mapping_dict.keys():
        present = False
        if True:#("#ME" in mapping_dict[k] or "#MG" in mapping_dict[k]):# or "#M" in mapping_dict[k] ):#it is in the list of motifs but not in the list of tf motifs
            if k in tf_names["JASPAR"]:# or k in tf_names["HOCOMOCO"]:
                jaspar_motifs.append(mapping_dict[k])
                
                for l in mapping_dict[k]:
                    if "#" not in l:
                        if l in tf_names["JASPAR"]:
                            present = True
                if not present:
                    tf_names_not_present.append(mapping_dict[k])
            
    for k in jaspar_motifs:
        print '\t'.join(k)
    #print tf_names_not_present

def filter_hyper_mutated_samples(mutations_input_file, mutations_output_file="", index_sample_id=-1, max_number_of_muts_to_remove_sample=100000):
    
    number_muts_sample = {}
    if mutations_output_file=="":
        mutations_output_file = '.'.join(mutations_input_file.split('.')[0:-1]) + '_filteredsamples' + str(max_number_of_muts_to_remove_sample)+ '.' + mutations_input_file.split('.')[-1]
    if os.path.exists(mutations_output_file):
        return mutations_output_file
    print mutations_output_file
    n=0
    with open(mutations_input_file, 'r') as mutations_infile:
        line = mutations_infile.readline()
        while line!='':
            if line.strip().split('\t')[index_sample_id] not in number_muts_sample.keys():
                number_muts_sample[line.strip().split('\t')[index_sample_id]]=0
            number_muts_sample[line.strip().split('\t')[index_sample_id]]+=1
            line = mutations_infile.readline()
            n+=1
            if n%1000000==0:
                print "processed: ", n, " lines"
    with open(mutations_output_file, 'w') as mutations_outfile:
        with open(mutations_input_file, 'r') as mutations_infile:
            line = mutations_infile.readline()
            while line!='':
                if number_muts_sample[line.strip().split('\t')[index_sample_id]] <= max_number_of_muts_to_remove_sample:
                    mutations_outfile.write(line)
                line = mutations_infile.readline()


def get_muts_from_borad_simulations(muts_input_file, samples_hist_matchings_input_file):
    
    muts_output_file = '.'.join(muts_input_file.split('.')[0:-1])+'_extracted'
    print muts_output_file
    samples_hist_dict = {}
    samples_count_dict = {}
    with open(samples_hist_matchings_input_file, 'r') as samples_hist_matchings_infile:
        samples_hist = samples_hist_matchings_infile.readlines()
        for s in samples_hist:
            if s.strip().split('\t')[0] not in samples_hist_dict.keys():
                samples_hist_dict[s.strip().split('\t')[0]] = s.strip().split('\t')[1]
    #with open(muts_input_file, 'w') as muts_outfile:
    print len(samples_hist_dict.keys())
    with open(muts_input_file, 'r') as muts_infile:
        mut_line = muts_infile.readline()#skip the header
        n=0
        while mut_line != '':
            mut_line = muts_infile.readline()#skip the header
            split_mut_line = mut_line.strip().split('\t')
            if len(split_mut_line)>2:
                if split_mut_line[2] in samples_hist_dict.keys():
                    if split_mut_line[2] not in samples_count_dict.keys() :
                        samples_count_dict[split_mut_line[2]] = 0
                    else:
                        samples_count_dict[split_mut_line[2]]+=1
            n+=1
            if n%1000000==0:
                print "processed ", n
            
    print len(samples_count_dict.keys())
    with open(muts_output_file, 'w') as muts_outfile:
        with open(muts_input_file, 'r') as muts_infile:
            n=0
            mut_line = muts_infile.readline()#skip the header
            i = 1
            print mut_line
            while mut_line != '':
                mut_line = muts_infile.readline()
                split_mut_line = mut_line.strip().split('\t')
                if len(split_mut_line[2])>2:
                    if split_mut_line[2] in samples_count_dict.keys(): 
                        if samples_count_dict[split_mut_line[2]] < 100000:
                            muts_outfile.write('chr'+split_mut_line[0] + '\t' + split_mut_line[1] + '\t' + split_mut_line[1] + '\t' + split_mut_line[3] + '\t' + split_mut_line[4] + '\t' + samples_hist_dict[split_mut_line[2]] + '\t' + 'B' + str(i) + '\t' + 'SNP' + '\t' + split_mut_line[2] + '\n')
                    i += 1
                n+=1
                if n%1000000==0:
                    print "processed ", n
            
            print 'number of lines written' + str(i)
    
    return 
    
def normalize_scores(input_dir, output_dir):
    for f in os.listdir("results"):
        if f.split('.')[-1]=='txt_statcalc_scoresSig':
            os.system("""awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$15}' """ + "results" + '/' + f + ' > ' + output_dir + '/' + f)
            #os.system("cat t " + "regDriver_coreProm/"+f + """ | awk 'BEGIN{FS=OFS="\t"}{if($5==1000) $5="NA"; print}' | sort | uniq | sort -k4 | groupBy -g 1,2,3,4 -c 5 -o collapse | awk 'BEGIN{FS=OFS="\t"}{if($5~","); gsub("NA,", ""); gsub(",NA", ""); print}' > """ + "regDriver_coreProm/"+f.replace("_statcalc_scoresSig", ""))
            

def generate_5col_bed(input_dir, output_dir, element_file, file_keyword="_statcalc_scoresSig"):
    
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    fs = os.listdir(input_dir)
    for f in fs:
        if file_keyword in f:
            os.system("cat " + element_file + ' ' + input_dir+'/'+f + """ | awk 'BEGIN{FS=OFS="\t"}{if(NF==12) $5="NA"; else if(NF==15) $5=$15; print $1,$2,$3,$4,$5}' | sort | uniq | sort -k4 | groupBy -g 1,2,3,4 -c 5 -o collapse | awk 'BEGIN{FS=OFS="\t"}{if($5~","); gsub("NA,", ""); gsub(",NA", ""); print}' > """ + output_dir + '/' + f.replace("_statcalc_scoresSig", ""))
    print "tar -czf " + output_dir + '.tar.gz' + ' ' + output_dir
    os.system('tar -czf ' + output_dir+'.tar.gz' + ' ' + output_dir)
    

if __name__ == '__main__':
    '''metadata_input_file = sys.argv[1]
    input_dir = sys.argv[2]
    merged_output_flie_name = sys.argv[3]
    
    dict_code_tumor_type = get_project_codes_metadata_file(metadata_input_file, index_tumor_type = 1, index_matching_code_from_snv_file = 25, index_donor_id=3, sep='\t')
    mege_files(input_dir, merged_output_flie_name, dict_code_tumor_type, index_filter=6, index_info=7, sep='\t')
    '''
    #rank_scores(input_file="~/Documents/Group-Projects/PCAWG/regDriverTest/t11", scores_index=1)
    #testf(infile="/Users/husensofteng/Documents/Group-Projects/PCAWG/regDriverTest/tu")
    
    #get_primary_sample_per_donor(mutations_file_MAF_input_file=sys.argv[1], mutations_file_MAF_output_file=sys.argv[2], selected_primary_tumors_inputfile=sys.argv[3], sample_per_donor_input_file=sys.argv[4])
    
    #find_motifTFname_expTFName_matches(motif_tfnames_input_file=sys.argv[1], exp_tfnames_input_file=sys.argv[2])
    #find_motifTFname_expTFName_matches(motif_tfnames_input_file="/Users/husensofteng/Downloads/ENCODE_motifs_tf_names.txt2", exp_tfnames_input_file="/Users/husensofteng/Downloads/ENCODE_TF_factor_group_name_metadataENCODETFChIPSeq21Sep2016_TFnames_map.tsv")
    
    #mapping_dict = get_dict_TFs_per_motif(mapping_file=sys.argv[1])#MotifTF_ChIPSeqTF_name_mapping.tsv (manually generated dict file)
    #compare_resources(mapping_dict, tf_names_input_file=sys.argv[2])#temp_tf_names.txt
    
    #filter_hyper_mutated_samples(mutations_input_file=sys.argv[1])
    #get_muts_from_borad_simulations(muts_input_file=sys.argv[1], samples_hist_matchings_input_file=sys.argv[2])
    
    generate_5col_bed(input_dir=sys.argv[1], output_dir=sys.argv[2], element_file=sys.argv[3], file_keyword="_statcalc_scoresSig")
    