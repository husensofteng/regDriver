'''
Created on Mar 23, 2016

@author: Husen M. Umer
'''

import os, sys
from urllib.request import urlopen
import gzip
import pybedtools
from pybedtools import BedTool
import shutil

def generate_list_of_accession_info(encode_metadata_inputfile, biosample_name_to_extract, assay_type, accepted_file_formats=[], default_target_name="Unknown", assembly="hg19",
                                    index_file_accession=0, index_file_format=1, index_output_type=2, index_assay_type=4, 
                                    index_biosample_name=6, index_target=15, index_assembly=40):#in the new version of metadata from ENCODE assembly is given on col41 while in the older version it was on col40
    
    metadata_lines = []
    targets_dict = {}
    
    with open(encode_metadata_inputfile, 'r') as encode_metadata_infile:
        metadata_lines = encode_metadata_infile.readlines()
    for line in metadata_lines:
        split_line = line.strip().split('\t')
        if biosample_name_to_extract==split_line[index_biosample_name] and assay_type==split_line[index_assay_type] and split_line[index_file_format] in accepted_file_formats and split_line[index_assembly]==assembly:
            target_name = split_line[index_target].replace("-human", "").replace("eGFP-", "").replace("HA-", "").replace("FLAG-", "") #all the factor names end with -human so there is no need to keep them; some factor names are written as eGFP-factorName or HA-factorName to indidicate the antibody type so they have to be removed to combine the same factors regardless of their antibody type 
            if target_name=="":
                target_name=default_target_name
            if target_name not in targets_dict.keys():
                targets_dict[target_name] = {}
            file_format_output_type = split_line[index_file_format]+"_"+split_line[index_output_type]
            if file_format_output_type not in targets_dict[target_name].keys():
                targets_dict[target_name][file_format_output_type]=[]
            targets_dict[target_name][file_format_output_type].append(split_line[index_file_accession])
            
    return targets_dict
    
def parse_cellinfodict_to_populate_data(cellinfodict_inputfile, cell_names_start_with="#"):
    cellinfo_lines = []
    with open(cellinfodict_inputfile, 'r') as cellinfodict_infile:
        cellinfo_lines = cellinfodict_infile.readlines()
    
    dict_cell_lines_info = {}
    cell_name = ""
    for line in cellinfo_lines:
        if line.startswith('//'):#skip lines start with // researved for comments in the config file
            continue
        elif line.startswith('***'):
            break
        elif line.startswith(cell_names_start_with):
            cell_name=line.strip().strip(cell_names_start_with)
            if cell_name not in dict_cell_lines_info.keys():
                dict_cell_lines_info[cell_name]={}
        elif (not line.startswith(cell_names_start_with)) and "=" in line:
            split_line = line.strip().split('=')            
            if len(split_line)==2:
                assay_type = split_line[0]
                assay_info = split_line[1].strip().split(',')
                if assay_info[0] != "":
                    dict_cell_lines_info[cell_name][assay_type]=assay_info
    return dict_cell_lines_info
    
def select_ENCODE_datasets_per_factor(ENCODE_accession_codes_dict, highest_priority_output_types=['optimal idr thresholded peaks', 'conservative idr thresholded peaks'], flag_indicating_high_qualtity_peaks="high", flag_indicating_low_qualtity_peaks="low"):
 
    selected_datasets_per_factor = {}
    for factor in ENCODE_accession_codes_dict.keys():
        if factor not in selected_datasets_per_factor.keys():
            selected_datasets_per_factor[factor] = []
        else:
            print("Warning: duplicate Factors:" + factor)
        
        for output_type in ENCODE_accession_codes_dict[factor].keys():
            quality_flag= flag_indicating_low_qualtity_peaks
            for highest_priority_output_type in highest_priority_output_types:
                if output_type==highest_priority_output_type:
                    quality_flag = flag_indicating_high_qualtity_peaks
                    break
            if quality_flag in selected_datasets_per_factor[factor]:
                selected_datasets_per_factor[factor][selected_datasets_per_factor[factor].index(quality_flag)+1].extend(ENCODE_accession_codes_dict[factor][output_type])
            else:
                selected_datasets_per_factor[factor].append(quality_flag)
                selected_datasets_per_factor[factor].append(ENCODE_accession_codes_dict[factor][output_type]) 
    return selected_datasets_per_factor
    
def download_and_unify_datasets(cell_name, assay_type, assay_info_dict, target_cellinfo_dirs_path, 
                                number_of_votes_from_highquality_datasets=1, 
                                number_of_votes_from_lowquality_datasets=2, 
                                number_of_files_to_consider_from_highquality_datasets='all', 
                                number_of_files_to_consider_from_lowquality_datasets='all', 
                                dont_consider_low_quality_datasets_when_highquality_datasets_available=True, 
                                consider_peak_score_from_peak_file = True, peak_score_index=6):
    
    current_dir = os.getcwd()
    final_dataset_of_this_assay_cell = cell_name+"_"+assay_type+".bed4"
    final_datasets_of_this_assay_cell = []
    if not os.path.exists(target_cellinfo_dirs_path+'/'+cell_name):
        os.mkdir(target_cellinfo_dirs_path+'/'+cell_name)
    if not os.path.exists(target_cellinfo_dirs_path+'/'+cell_name+'/'+assay_type):
        os.mkdir(target_cellinfo_dirs_path+'/'+cell_name+'/'+assay_type)
    os.chdir(target_cellinfo_dirs_path+'/'+cell_name+'/'+assay_type)
    print(target_cellinfo_dirs_path+'/'+cell_name+'/'+assay_type)
    if not os.path.exists(final_dataset_of_this_assay_cell):
        for factor in assay_info_dict.keys():
            peak_score_from_peak_file_exists = True
            final_dataset = factor+".bed4"
            if os.path.exists(final_dataset): #if the final merged file of this factor was already available then no need to do any more operations
                final_datasets_of_this_assay_cell.append(final_dataset)
                continue
            list_of_high_quality_datasets_from_this_factor = []
            list_of_low_quality_datasets_from_this_factor = []
            print(assay_info_dict[factor])
            if 'high' in assay_info_dict[factor]:
                if number_of_files_to_consider_from_highquality_datasets=='all':
                    list_of_high_quality_datasets_from_this_factor=assay_info_dict[factor][assay_info_dict[factor].index("high")+1]
                else:
                    for i in range(0, number_of_files_to_consider_from_highquality_datasets):
                        if i < len(assay_info_dict[factor][assay_info_dict[factor].index("high")+1]):
                            list_of_high_quality_datasets_from_this_factor.append(assay_info_dict[factor][assay_info_dict[factor].index("high")+1][i])
                        else:
                            break
            if 'low' in assay_info_dict[factor]:
                if 'high' in assay_info_dict[factor] and dont_consider_low_quality_datasets_when_highquality_datasets_available:
                    pass
                #in case no high quality dataset was available or it was specifically asked to include even with the availablility of high quality datasets then use the low quality datasets as well
                else:
                    if number_of_files_to_consider_from_lowquality_datasets=='all':
                        list_of_low_quality_datasets_from_this_factor=assay_info_dict[factor][assay_info_dict[factor].index("low")+1]
                    else:
                        for i in range(0, number_of_files_to_consider_from_lowquality_datasets):
                            if i < len(assay_info_dict[factor][assay_info_dict[factor].index("low")+1]):
                                list_of_low_quality_datasets_from_this_factor.append(assay_info_dict[factor][assay_info_dict[factor].index("low")+1][i])
                            else:
                                break
            
            #process the datasets from the high quality list
            list_of_high_quality_peakfiles_from_this_factor = []
            final_dataset_high_quality = factor+"_high" + ".bed"
            final_dataset_high_quality_name = open(final_dataset_high_quality, 'w')
            for dataset in list_of_high_quality_datasets_from_this_factor:
                dataset_name = ""
                if "ENCFF" in dataset:
                    dataset_path = "https://www.encodeproject.org/files/"+dataset+"/@@download/"+dataset+".bed.gz"
                    dataset_name = factor+"_"+dataset+".bed"
                    if not os.path.exists(dataset_name):
                        if not os.path.exists(dataset_name+".gz"):
                            downloaded_obj = urlopen(dataset_path)
                            print("downloading.... " + dataset_path)
                            with open(os.path.basename(dataset_name+".gz"), 'wb') as local_file:
                                local_file.write(downloaded_obj.read())
                        with gzip.open(dataset_name+".gz", 'rb') as dataset_name_zip, open(dataset_name, 'w') as dataset_name_unzipped:
                            dataset_name_unzipped.write(dataset_name_zip.read())
                        #os.system("gunzip " + dataset_name+".gz")
                elif dataset.startswith("http://") or dataset.startswith("ftp://"):
                    dataset_path = dataset
                    dataset_name = factor+"_"+dataset.strip().split('/')[-1]
                    dataset_name_unzipped = dataset_name
                    if "." in dataset_name: 
                        if dataset_name.split('.')[-1]=="gz":
                            dataset_name_unzipped = '.'.join(dataset_name.split('.')[0:-1])
                            
                    if os.path.exists(dataset_name_unzipped):#this could be the gzip or the unzipped file
                        dataset_name = '.'.join(dataset_name.split('.')[0:-1])
                    else:
                        if not os.path.exists(dataset_name):
                            downloaded_obj = urlopen(dataset_path)
                            print("downloading.... " + dataset_path)
                            with open(dataset_name, 'wb') as local_file:
                                local_file.write(downloaded_obj.read())
                        if "." in dataset_name: 
                            if dataset_name.split('.')[-1]=="gz":
                                with gzip.open(dataset_name, 'rb') as dataset_name_unzip_read, open(dataset_name_unzipped, 'wb') as dataset_name_unzipped_write:
                                    dataset_name_unzipped_write.write(dataset_name_unzip_read.read())
                                #os.system("gunzip " + dataset_name)
                                dataset_name = '.'.join(dataset_name.split('.')[0:-1])
                else:#path to a local file
                    dataset_path = dataset
                    dataset_name = factor+"_"+dataset.strip().split('/')[-1]
                    if not os.path.exists(dataset_name):
                        shutil.copy(dataset_name, "./")
                        if "." in dataset_name:
                            if dataset_name.split('.')[-1]=="gz":
                                with gzip.open(dataset_name, 'rb') as dataset_name_unzip_read, open('.'.join(dataset_name.split('.')[0:-1]), 'wb') as dataset_name_unzip_write:
                                    dataset_name_unzip_write.write(dataset_name_unzip_read.read())
                                #os.system("gunzip " + dataset_name)
                                dataset_name = '.'.join(dataset_name.split('.')[0:-1])
                                
                if dataset_name!="":
                    dataset_sort_bedtools = pybedtools.BedTool(dataset_name)
                    sorting_result = dataset_sort_bedtools.sort() 
                    list_of_high_quality_peakfiles_from_this_factor.append(sorting_result.fn)
            #Combine all high quality peak files into one
            peak_score_from_peak_file_exists = True
            if len(list_of_high_quality_peakfiles_from_this_factor)!=0:
                print(cell_name + ": high: " + assay_type + ":" + factor + ": "  + ','.join(list_of_high_quality_peakfiles_from_this_factor))
                #merge the high quality datasets
                if len(list_of_high_quality_peakfiles_from_this_factor)==1:
                    if assay_type == "ChromatinStates":
                        final_dataset_high_quality = list_of_high_quality_peakfiles_from_this_factor[0]
                    else:
                        merged_output = open(list_of_high_quality_peakfiles_from_this_factor[0], 'r').readlines()
                        try:#in case the line didn't have col index or the value of col index was not convertable to float then it's an indication of no score availabbility
                            peak_score = float(merged_output[0][peak_score_index])
                        except (IndexError, ValueError) as e:
                            repr( e )
                            peak_score_from_peak_file_exists = False
                        if peak_score_from_peak_file_exists and consider_peak_score_from_peak_file:
                            for line in merged_output:
                                final_dataset_high_quality_name.write('\t'.join(line.strip().split('\t')[0:3]) + '\t' + line.strip().split('\t')[peak_score_index] +"\n")
                        else:
                            for line in merged_output:
                                final_dataset_high_quality_name.write('\t'.join(line.strip().split('\t')[0:3]) +"\n")
                        final_dataset_high_quality_name.close()
                elif len(list_of_high_quality_peakfiles_from_this_factor)>1:
                    if assay_type == "ChromatinStates":
                        #write all the files into one
                        with open(final_dataset_high_quality, 'w') as concatenated_file_write: 
                            for file_name in list_of_high_quality_peakfiles_from_this_factor:
                                with open(file_name, 'r') as infile: 
                                    concatenated_file_write.write(infile.read())
                    else:
                        bedTools_obj = BedTool()
                        merging_all = bedTools_obj.multi_intersect(i=list_of_high_quality_peakfiles_from_this_factor).filter(lambda x: int(x[3]) >= number_of_votes_from_highquality_datasets).sort().merge()
                        
                        #in case the line didn't have col index or the value of col index was not convertable to float then it's an indication of no score availabbility
                        for file_i in  list_of_high_quality_peakfiles_from_this_factor:#check if all the files have peak scores
                            with open(file_i) as read_file_i:
                                try:
                                    h = read_file_i.readline()
                                    peak_score = float(h.strip().split('\t')[peak_score_index])
                                except (IndexError, ValueError) as e:
                                    repr( e )
                                    peak_score_from_peak_file_exists = False
                        if peak_score_from_peak_file_exists and consider_peak_score_from_peak_file:
                            #tmp_dir = './tmp_dir_to_remove_{}'.format(list_of_high_quality_peakfiles_from_this_factor[0].split('/')[-1])
                            #os.makedirs(tmp_dir)
                            list_of_high_quality_peakfiles_from_this_factor_updated = []
                            for i_file in list_of_high_quality_peakfiles_from_this_factor:
                                with open(i_file, 'r') as ifile, open(i_file + "_tmp", 'w') as ofile:
                                    for line in ifile.readlines():
                                        ofile.write('\t'.join(line.strip().split('\t')[0:3]) + '\t{}\n'.format(line.strip().split('\t')[peak_score_index]))
                                list_of_high_quality_peakfiles_from_this_factor_updated.append(i_file + "_tmp")
                            peak_score_index_updated = 3
                            #os.system('cp ' + merging_all.fn + ' . ' )
                            merging_all = merging_all.intersect(list_of_high_quality_peakfiles_from_this_factor_updated, wo=True).sort().groupby(g=[1,2,3], c=peak_score_index_updated+1+4, o=['mean'])#4 cols from the mergeBed and one extra from the intersection then it follows the cols from each file
                            #os.system('cp ' + merging_all.fn + ' . ' )
                            for l in list_of_high_quality_peakfiles_from_this_factor_updated:
                                os.remove(l)
                        merged_output = open(merging_all.fn, 'r').readlines()
                        for line in merged_output:
                            final_dataset_high_quality_name.write('\t'.join(line.strip().split('\t')[0::]) +"\n")
                        final_dataset_high_quality_name.close()
            #handling peak files from low quality datasets    
            list_of_low_quality_peakfiles_from_this_factor = []
            final_dataset_low_quality = factor+"_low" + ".bed"
            final_dataset_low_quality_name = open(final_dataset_low_quality, 'w')
            for dataset in list_of_low_quality_datasets_from_this_factor:
                dataset_name = ""
                if "ENCFF" in dataset:
                    dataset_path = "https://www.encodeproject.org/files/"+dataset+"/@@download/"+dataset+".bed.gz"
                    dataset_name = factor+"_"+dataset+".bed"
                    if not os.path.exists(dataset_name):
                        if not os.path.exists(dataset_name+".gz"):
                            downloaded_obj = urlopen(dataset_path)
                            print("downloading.... " + dataset_path)
                            with open(os.path.basename(dataset_name+".gz"), 'wb') as local_file:
                                local_file.write(downloaded_obj.read())
                        with gzip.open(dataset_name+".gz", 'rb') as dataset_name_zip, open(dataset_name, 'w') as dataset_name_unzipped:
                            dataset_name_unzipped.write(dataset_name_zip.read())
                        #os.system("gunzip " + dataset_name+".gz")
                elif dataset.startswith("http://") or dataset.startswith("ftp://"):
                    dataset_path = dataset
                    dataset_name = factor+"_"+dataset.strip().split('/')[-1]
                    dataset_name_unzipped = dataset_name
                    if "." in dataset_name: 
                        if dataset_name.split('.')[-1]=="gz":
                            dataset_name_unzipped = '.'.join(dataset_name.split('.')[0:-1])
                            
                    if os.path.exists(dataset_name_unzipped):#this could be the gzip or the unzipped file
                        dataset_name = '.'.join(dataset_name.split('.')[0:-1])
                    else:
                        if not os.path.exists(dataset_name):
                            downloaded_obj = urlopen(dataset_path)
                            print("downloading.... " + dataset_path)
                            with open(dataset_name, 'wb') as local_file:
                                local_file.write(downloaded_obj.read())
                        if "." in dataset_name: 
                            if dataset_name.split('.')[-1]=="gz":
                                with gzip.open(dataset_name, 'rb') as dataset_name_unzip_read, open(dataset_name_unzipped, 'wb') as dataset_name_unzipped_write:
                                    dataset_name_unzipped_write.write(dataset_name_unzip_read.read())
                                #os.system("gunzip " + dataset_name)
                                dataset_name = '.'.join(dataset_name.split('.')[0:-1])
                else:#path to a local file
                    dataset_path = dataset
                    dataset_name = factor+"_"+dataset.strip().split('/')[-1]
                    if not os.path.exists(dataset_name):
                        shutil.copy(dataset, "./"+dataset_name)
                        if "." in dataset_name:
                            if dataset_name.split('.')[-1]=="gz":
                                with gzip.open(dataset_name, 'rb') as dataset_name_unzip_read, open('.'.join(dataset_name.split('.')[0:-1]), 'wb') as dataset_name_unzip_write:
                                    dataset_name_unzip_write.write(dataset_name_unzip_read.read())
                                #os.system("gunzip " + dataset_name)
                                dataset_name = '.'.join(dataset_name.split('.')[0:-1])
                                
                if dataset_name!="":
                    dataset_sort_bedtools = pybedtools.BedTool(dataset_name)
                    sorting_result = dataset_sort_bedtools.sort() 
                    list_of_low_quality_peakfiles_from_this_factor.append(sorting_result.fn)
            #Combine all low quality peak files into one
            peak_score_from_peak_file_exists = True
            if len(list_of_low_quality_peakfiles_from_this_factor)!=0:
                print(cell_name + ": low: " + assay_type + ":" + factor + ": "  + ','.join(list_of_low_quality_peakfiles_from_this_factor))
                #merge the low quality datasets
                if len(list_of_low_quality_peakfiles_from_this_factor)==1:
                    if assay_type == "ChromatinStates":
                        final_dataset_low_quality = list_of_low_quality_peakfiles_from_this_factor[0]
                    else:
                        merged_output = open(list_of_low_quality_peakfiles_from_this_factor[0], 'r').readlines()
                        try:#in case the line didn't have col index or the value of col index was not convertable to float then it's an indication of no score availabbility
                            peak_score = float(merged_output[0][peak_score_index])
                        except (IndexError, ValueError) as e:
                            repr( e )
                            peak_score_from_peak_file_exists = False
                        if peak_score_from_peak_file_exists and consider_peak_score_from_peak_file:
                            for line in merged_output:
                                final_dataset_low_quality_name.write('\t'.join(line.strip().split('\t')[0:3]) + '\t' + line.strip().split('\t')[peak_score_index] + "\n")
                        else:
                            for line in merged_output:
                                final_dataset_low_quality_name.write('\t'.join(line.strip().split('\t')[0:3]) + "\n")
                        final_dataset_low_quality_name.close()
                
                elif len(list_of_low_quality_peakfiles_from_this_factor)>1:
                    if assay_type == "ChromatinStates":
                        #write all the files into one
                        with open(final_dataset_low_quality, 'w') as concatenated_file_write: 
                            for file_name in list_of_low_quality_peakfiles_from_this_factor:
                                with open(file_name, 'r') as infile: 
                                    concatenated_file_write.write(infile.read())
                    else:
                        bedTools_obj = BedTool()
                        merging_all = bedTools_obj.multi_intersect(i=list_of_low_quality_peakfiles_from_this_factor).filter(lambda x: int(x[3]) >= number_of_votes_from_lowquality_datasets).sort().merge()
                        for file_i in  list_of_low_quality_peakfiles_from_this_factor:#check if all the files have peak scores
                            with open(file_i) as read_file_i:
                                try:
                                    h = read_file_i.readline()
                                    peak_score = float(h.strip().split('\t')[peak_score_index]) 
                                except (IndexError, ValueError) as e:
                                    repr( e )
                                    peak_score_from_peak_file_exists = False
                        if peak_score_from_peak_file_exists and consider_peak_score_from_peak_file:
                            list_of_low_quality_datasets_from_this_factor_updated = []
                            for i_file in list_of_low_quality_peakfiles_from_this_factor:
                                with open(i_file, 'r') as ifile, open(i_file + "_tmp", 'w') as ofile:
                                    for line in ifile.readlines():
                                        ofile.write('\t'.join(line.strip().split('\t')[0:3]) + '\t{}\n'.format(line.strip().split('\t')[peak_score_index]))
                                list_of_low_quality_datasets_from_this_factor_updated.append(i_file + "_tmp")
                            peak_score_index_updated = 3
                            merging_all = merging_all.intersect(list_of_low_quality_datasets_from_this_factor_updated, wo=True).sort().groupby(g=[1,2,3], c=peak_score_index_updated+1+4, o=['mean'])#4 cols from the mergeBed and one extra from the intersection then it follows the cols from each file
                            for l in list_of_low_quality_datasets_from_this_factor_updated:
                                os.remove(l)
                        merged_output = open(merging_all.fn, 'r').readlines()
                        for line in merged_output:
                            final_dataset_low_quality_name.write('\t'.join(line.strip().split('\t')[0::]) +"\n")
                        final_dataset_low_quality_name.close()
            
            #Combine results of low and high quality peak files and merge them with adding the factor name
            peak_score_from_peak_file_exists = True
            final_file = ""
            if os.stat(final_dataset_high_quality).st_size==0 and os.stat(final_dataset_low_quality).st_size==0:
                continue
            else:
                merge_final_lines = []
                highlow_combined  = "highlow_combined"
                os.system("cat " + final_dataset_high_quality + " " + final_dataset_low_quality + " > " + highlow_combined)
                if assay_type == "ChromatinStates":#because the chromatinstates are defined for all genome bins merging them would cause create 25 regions only since all the bins are starting consequentively 
                    final_file = highlow_combined
                else:
                    with open(highlow_combined) as read_file_i:
                        try:
                            h = read_file_i.readline()
                            peak_score = float(h.strip().split('\t')[3])
                        except (IndexError, ValueError) as e:
                            repr( e )
                            peak_score_from_peak_file_exists = False
                    highlow_combined_obj = BedTool(highlow_combined)
                    merge_final = ""
                    if peak_score_from_peak_file_exists and consider_peak_score_from_peak_file:
                        merge_final = highlow_combined_obj.sort().merge(c=4, o='mean')
                    else:
                        merge_final = highlow_combined_obj.sort().merge()
                    final_file = merge_final.fn
                with open(final_file, 'r') as merge_final_read:
                    merge_final_lines = merge_final_read.readlines()
                    with open(final_dataset, 'w') as final_dataset_writer:  
                        if assay_type == "ChromatinStates":
                            for line in merge_final_lines:
                                final_dataset_writer.write('\t'.join(line.strip().split('\t')[0:3]) + '\t' + cell_name+"#ChromHMM#"+line.strip().split('\t')[3].replace(" ", "-") + '\n')
                        elif assay_type == "ChIP-seq":
                            if peak_score_from_peak_file_exists and consider_peak_score_from_peak_file:
                                for line in merge_final_lines:
                                    peak_score = "#"+str(line.strip().split('\t')[3])
                                    final_dataset_writer.write('\t'.join(line.strip().split('\t')[0:3]) + '\t' + cell_name+"#TFBinding#"+factor.replace(" ", "-")+peak_score + '\n')
                            else:
                                for line in merge_final_lines:
                                    peak_score = ""
                                    final_dataset_writer.write('\t'.join(line.strip().split('\t')[0:3]) + '\t' + cell_name+"#TFBinding#"+factor.replace(" ", "-")+peak_score + '\n')
                        else:
                            if peak_score_from_peak_file_exists and consider_peak_score_from_peak_file:
                                for line in merge_final_lines:
                                    peak_score = "#"+str(line.strip().split('\t')[3])
                                    final_dataset_writer.write('\t'.join(line.strip().split('\t')[0:3]) + '\t' + cell_name+"#"+factor.replace(" ", "-")+peak_score + '\n')
                            else:
                                peak_score = ""
                                for line in merge_final_lines:
                                    final_dataset_writer.write('\t'.join(line.strip().split('\t')[0:3]) + '\t' + cell_name+"#"+factor.replace(" ", "-")+peak_score + '\n')
                final_datasets_of_this_assay_cell.append(final_dataset)
                os.remove(highlow_combined)
            os.remove(factor+"_high" + ".bed")
            os.remove(factor+"_low" + ".bed")
        #combine peak files of all the factors into one
        with open(final_dataset_of_this_assay_cell, 'w') as final_dataset_of_this_cell_out:
            for peak_file in final_datasets_of_this_assay_cell:
                with open(peak_file, 'r') as infile:
                    final_dataset_of_this_cell_out.write(infile.read())
    
    os.chdir(current_dir)
    return final_dataset_of_this_assay_cell, final_datasets_of_this_assay_cell
    
def populate_cellinfo_dirs(dict_cell_lines_info, target_cellinfo_dirs_path):
    #for a given cell name and assay_type make subdirs and retreive the accession codes to download the files
    summary_file_name = "cell_info_summary"
    summary_file = open("cell_info_summary", 'w')
    final_dataset_of_this_assay_cell = ""
    final_datasets_of_this_assay_cell = []
    for cell_name in dict_cell_lines_info:
        for assay_type in dict_cell_lines_info[cell_name]:
            final_dataset_of_this_assay_cell = cell_name+"_"+assay_type+".bed4"
            if os.path.exists(target_cellinfo_dirs_path + "/" + cell_name + "/" + assay_type + "/" + final_dataset_of_this_assay_cell):#if the the combined file of this assay in this cell has been created there is no need to go any further
                continue
            ENCODE_accession_codes_dict = {}
            if assay_type=="ChIP-seq" or assay_type=="DNase-seq" or assay_type=="ChromatinStates" or assay_type=="Repli-seq":
                selected_datasets_per_factor_dict = {}
                selected_datasets_per_factor_dict_from_metadatafile = {}
                #generate the list of dataset IDs per factor in each assay type
                for datasource in dict_cell_lines_info[cell_name][assay_type]:
                    if "metadataENCODE" in datasource:
                        if assay_type=="ChIP-seq":
                            ENCODE_accession_codes_dict = generate_list_of_accession_info(datasource, cell_name, assay_type, accepted_file_formats=['bed narrowPeak'])#'bed broadPeak', 
                            selected_datasets_per_factor_dict_from_metadatafile = select_ENCODE_datasets_per_factor(ENCODE_accession_codes_dict, highest_priority_output_types= ['bed narrowPeak_optimal idr thresholded peaks', 'bed narrowPeak_conservative idr thresholded peaks'], flag_indicating_high_qualtity_peaks="high", flag_indicating_low_qualtity_peaks="low")
                        elif assay_type=="DNase-seq":
                            ENCODE_accession_codes_dict = generate_list_of_accession_info(datasource, cell_name, assay_type, accepted_file_formats=['bed narrowPeak'], default_target_name=assay_type)
                            selected_datasets_per_factor_dict_from_metadatafile = select_ENCODE_datasets_per_factor(ENCODE_accession_codes_dict, highest_priority_output_types= ['bed narrowPeak_peaks'], flag_indicating_high_qualtity_peaks="high", flag_indicating_low_qualtity_peaks="low")
                        elif assay_type=="Repli-seq":
                            ENCODE_accession_codes_dict = generate_list_of_accession_info(datasource, cell_name, assay_type, accepted_file_formats=['bigWig'], default_target_name=assay_type)
                            selected_datasets_per_factor_dict_from_metadatafile = select_ENCODE_datasets_per_factor(ENCODE_accession_codes_dict, highest_priority_output_types= ['bigWig_percentage normalized signal'], flag_indicating_high_qualtity_peaks="high", flag_indicating_low_qualtity_peaks="low")
                        if len(selected_datasets_per_factor_dict)==0:
                            selected_datasets_per_factor_dict=selected_datasets_per_factor_dict_from_metadatafile
                        else:
                            for factor in selected_datasets_per_factor_dict_from_metadatafile.keys():
                                if factor not in selected_datasets_per_factor_dict.keys():
                                    selected_datasets_per_factor_dict[factor]=selected_datasets_per_factor_dict_from_metadatafile[factor]
                                else:
                                    for items in range(0,len(selected_datasets_per_factor_dict_from_metadatafile[factor]),2):
                                        quality_level_from_the_metadatENCODE_dict=selected_datasets_per_factor_dict_from_metadatafile[factor][items]
                                        if quality_level_from_the_metadatENCODE_dict in selected_datasets_per_factor_dict[factor]:
                                            selected_datasets_per_factor_dict[factor][selected_datasets_per_factor_dict[factor].index(quality_level_from_the_metadatENCODE_dict)+1].extend(selected_datasets_per_factor_dict_from_metadatafile[factor][items+1])
                    elif "#" in datasource: #FactorName#high|low#list_of_files|URLs-separated by semicolon]]
                        datasource_split = datasource.strip().split('#')
                        factor_name_from_otherdatasource = datasource_split[0]
                        quality_level_from_otherdatasource = datasource_split[1]
                        source_from_otherdatasource = datasource_split[2].split(';')
                        if factor_name_from_otherdatasource not in selected_datasets_per_factor_dict.keys():
                            selected_datasets_per_factor_dict[factor_name_from_otherdatasource] = []
                        if quality_level_from_otherdatasource in selected_datasets_per_factor_dict[factor_name_from_otherdatasource]:
                            selected_datasets_per_factor_dict[factor_name_from_otherdatasource][selected_datasets_per_factor_dict[factor_name_from_otherdatasource].index(quality_level_from_otherdatasource)+1].extend(source_from_otherdatasource)
                        else:
                            selected_datasets_per_factor_dict[factor_name_from_otherdatasource].append(quality_level_from_otherdatasource)
                            selected_datasets_per_factor_dict[factor_name_from_otherdatasource].append(source_from_otherdatasource)
                
                if assay_type=="ChIP-seq":
                    final_dataset_of_this_assay_cell, final_datasets_of_this_assay_cell = download_and_unify_datasets(cell_name, assay_type, selected_datasets_per_factor_dict, target_cellinfo_dirs_path, number_of_votes_from_highquality_datasets=1, number_of_votes_from_lowquality_datasets=2, number_of_files_to_consider_from_highquality_datasets='all', number_of_files_to_consider_from_lowquality_datasets='all', dont_consider_low_quality_datasets_when_highquality_datasets_available=True, consider_peak_score_from_peak_file = True, peak_score_index=6)
                elif assay_type=="DNase-seq":
                    final_dataset_of_this_assay_cell, final_datasets_of_this_assay_cell = download_and_unify_datasets(cell_name, assay_type, selected_datasets_per_factor_dict, target_cellinfo_dirs_path, number_of_votes_from_highquality_datasets=1, number_of_votes_from_lowquality_datasets=2, number_of_files_to_consider_from_highquality_datasets='all', number_of_files_to_consider_from_lowquality_datasets='all', dont_consider_low_quality_datasets_when_highquality_datasets_available=True, consider_peak_score_from_peak_file = True, peak_score_index=6)
                elif assay_type=="Repli-seq":
                    download_and_unify_datasets(cell_name, assay_type, selected_datasets_per_factor_dict, target_cellinfo_dirs_path, number_of_votes_from_highquality_datasets=1, number_of_votes_from_lowquality_datasets=1, number_of_files_to_consider_from_highquality_datasets=1, number_of_files_to_consider_from_lowquality_datasets=1, dont_consider_low_quality_datasets_when_highquality_datasets_available=True, consider_peak_score_from_peak_file = True, peak_score_index=6)
                elif assay_type == "ChromatinStates":
                    final_dataset_of_this_assay_cell, final_datasets_of_this_assay_cell = download_and_unify_datasets(cell_name, assay_type, selected_datasets_per_factor_dict, target_cellinfo_dirs_path, number_of_votes_from_highquality_datasets=1, number_of_votes_from_lowquality_datasets=1, number_of_files_to_consider_from_highquality_datasets=1, number_of_files_to_consider_from_lowquality_datasets='all', dont_consider_low_quality_datasets_when_highquality_datasets_available=True, consider_peak_score_from_peak_file = False)
            if final_dataset_of_this_assay_cell!="":
                summary_file.write(final_dataset_of_this_assay_cell+"\t"+str(len(final_datasets_of_this_assay_cell)) + '\t' +','.join(final_datasets_of_this_assay_cell)+'\n')
    summary_file.close()
    return summary_file_name


if __name__ == '__main__':
    
    if len(sys.argv)==3:
        cellinfodict_inputfile=sys.argv[1] # path to CellInfoDict
        target_cellinfo_dirs_path = sys.argv[2]#"path to PanCanRegulators
    else:
        print("Usage: python ParseCellInfo.py CellInfoDict_input_file target_cellinfo_dirs_path")  
    dict_cell_lines_info =  parse_cellinfodict_to_populate_data(cellinfodict_inputfile, cell_names_start_with="#")
    populate_cellinfo_dirs(dict_cell_lines_info, target_cellinfo_dirs_path)
    
    
    
    
    