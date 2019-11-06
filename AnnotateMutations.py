'''
Created on Mar 31, 2017

@author: husensofteng
'''
'''
Annotates a list of genomic coordinates using overlapping_tracks of cell lines corresponding to the tumor type given in each line
'''
import os, sys, shutil
from pybedtools import BedTool, set_tempdir
from collections import Counter

temp_dir = 'tmp_pybedtoos'
if not os.path.exists(temp_dir):
    os.mkdir(temp_dir)  
set_tempdir(temp_dir)

def get_muts_tracks_info(muts_input_file, tracks_dir, muts_dir_out, split_muts_file_by_chr=False):
    
    muts_tracks_files = []
    
    tracks_files = [x for x in os.listdir(tracks_dir) if x.endswith('.bed')]
    if split_muts_file_by_chr:
        muts_files = []
        chr_ext = "." + tracks_files[0].split('.')[-1]
        if os.path.exists(muts_dir_out):
            muts_tracks_files = [muts_dir_out+'/'+x for x in os.listdir(muts_dir_out) if x.endswith('_overlapping_tracks.bed10')]
            if len(muts_tracks_files)>0:
                return muts_tracks_files
        else:
            os.mkdir(muts_dir_out)
            muts_files = [muts_dir_out+'/'+x for x in os.listdir(muts_dir_out) if x.endswith(chr_ext)]
            if len(muts_files)<=0:
                os.system("""awk '{{print $0 >> "{muts_dir}/"$1"{chr_ext}"}}' {muts_file} """.format(muts_dir=muts_dir_out, chr_ext=chr_ext, muts_file=muts_input_file))
                muts_files = [muts_dir_out+'/'+x for x in os.listdir(muts_dir_out) if x.endswith(chr_ext)]
        
        print('muts_files:  ', muts_files)
        print('tracks_files: ', tracks_files)
        for muts_file in muts_files:
            if muts_file.split('/')[-1] in tracks_files:
                muts_tracks_file = muts_file+"_overlapping_tracks.bed10"
                if not os.path.exists(muts_tracks_file):
                    print("Intersecting and Grouping: ", muts_tracks_file)
                    BedTool(muts_file).intersect(BedTool(tracks_dir+'/'+tracks_files[tracks_files.index(muts_file.split('/')[-1])]), wo=True, loj=True).groupby(g=[1,2,3,4,5,6,7,8,9], c=13, o=['collapse']).saveas(muts_tracks_file)
                muts_tracks_files.append(muts_tracks_file)
    else:
        for tracks_file in tracks_files:
            if not os.path.exists(muts_dir_out):
                os.mkdir(muts_dir_out)
                
            muts_tracks_file = muts_dir_out+'/'+tracks_file+"_overlapping_tracks.bed10"
            if not os.path.exists(muts_tracks_file):
                print("Intersecting and Grouping: ", muts_tracks_file)
                BedTool(muts_input_file).intersect(BedTool(tracks_dir+'/'+tracks_file), wo=True, loj=True).groupby(g=[1,2,3,4,5,6,7,8,9], c=13, o=['collapse']).saveas(muts_tracks_file)
            muts_tracks_files.append(muts_tracks_file)
    print('muts_tracks_files: ', muts_tracks_files)
    return muts_tracks_files

def annotate_muts(muts_tracks_file, muts_tracks_ouput_file, tissue_cell_assays, matching_cell_name_representative_dict, cancer_type_index = 9, tracks_info_index = 3):
    with open(muts_tracks_file, 'r') as mut_tracks_ifile, open(muts_tracks_ouput_file, 'w') as mut_tracks_ofile:
        l = mut_tracks_ifile.readline().strip().split('\t')
        while len(l)>cancer_type_index:
            annotations = {}
            tracks = l[tracks_info_index].split(',')
            for track in tracks:
                #print(track, l[cancer_type_index])
                cellname = track.split('#')[0]
                try:
                    cellname = matching_cell_name_representative_dict[cellname][0]
                except KeyError:
                    continue
                
                try:
                    if cellname in tissue_cell_assays[l[cancer_type_index]]:
                        try:
                            annotations[track.split('#')[1]].append(float(track.split('#')[2]))#annotations[track.split('#')[1]].append(cellname+"#"+'#'.join(track.split('#')[1:]))
                        except ValueError:
                            annotations[track.split('#')[1]].append(track.split('#')[2])
                        except KeyError:
                            try:
                                annotations[track.split('#')[1]] = [float(track.split('#')[2])]
                            except ValueError:
                                annotations[track.split('#')[1]] = [track.split('#')[2]]
                            except IndexError:
                                annotations[track.split('#')[1]] = [1.0]
                        except IndexError:
                                annotations[track.split('#')[1]] = [1.0]
                except KeyError:
                    print("Cancer type not found: ", l[cancer_type_index])
                    print(track)
                    print(cellname)
                    print(l[cancer_type_index])
                    print(tissue_cell_assays[l[cancer_type_index]])
                    break
                
            for anno in sorted(annotations.keys()):
                if anno=='TFBinding':
                    annotations[anno] = ';'.join(set(annotations[anno]))
                    continue
                try:
                    annotations[anno] = str(sum(annotations[anno])/(len(annotations[anno])*1.0))
                except (ValueError, TypeError):
                    annotations[anno] = Counter(annotations[anno]).most_common(1)[0][0]
            
            annotations_combined = [x+":"+annotations[x] for x in sorted(annotations.keys())]
            if len(annotations_combined)==0:
                annotations_combined = ['NaN']
            mut_tracks_ofile.write('\t'.join(l[0:9]) + '\t' + '|'.join(annotations_combined) + '\n')
            l = mut_tracks_ifile.readline().strip().split('\t')
            
    return muts_tracks_ouput_file
    

def retreive_key_values_from_dict_file(dict_input_file, key_value_sep='=', values_sep=','):#TFFamilyName TF_name
    "Retrieves the key and its values"
    key_values_dict = {}
    value_key_dict = {}
    with open(dict_input_file, 'r') as dict_input_file_infile:
        lines = dict_input_file_infile.readlines()
        for line in lines:
            if line.startswith('//') or line.startswith('#'):# or '=' not in line:
                continue
            sl = line.strip().split(key_value_sep)
            key_value = sl[0].strip()
            if key_value not in key_values_dict.keys():
                key_values_dict[key_value] = []
            if key_value not in value_key_dict:
                value_key_dict[key_value]=[key_value]
            if len(sl)>1:
                for s in sl[1].split(values_sep):
                    if s.strip()!="" and s.strip() not in key_values_dict[key_value]:
                        key_values_dict[key_value].append(s.strip())
                    if s.strip()!="":
                        if s.strip() not in value_key_dict:
                            value_key_dict[s.strip()]=[]
                        value_key_dict[s.strip()].append(key_value)
    return key_values_dict, value_key_dict
    
def get_tissue_cell_mappings(tissue_cell_mappings_file, matching_cell_name_representative_dict, key_value_sep='=', values_sep=',', cell_assay_sepe=':'):
    
    tissue_cell_assays = {}    
    with open(tissue_cell_mappings_file, 'r') as ifile:
        lines = ifile.readlines()
        for line in lines:
            if line.startswith('//') or line.startswith('#') or '=' not in line:
                continue
            sl = line.strip().split(key_value_sep)
            key_value = sl[0]
            if key_value not in tissue_cell_assays.keys():
                tissue_cell_assays[key_value] = []
            
            for s in sl[1].split(values_sep):
                try:
                    cell = matching_cell_name_representative_dict[s.split(':')[0]][0]#get the rep name of the cell
                    if cell not in tissue_cell_assays[key_value]:
                        tissue_cell_assays[key_value].append(cell)
                except KeyError:
                    continue
    return tissue_cell_assays


def get_annotated_muts(muts_input_file, tracks_dir, muts_out):
    cell_names_to_use = 'CellNamesDict'
    tissue_cell_mappings_file='TissueCellMatches'
    representative_cell_name_matchings_dict, matching_cell_name_representative_dict = retreive_key_values_from_dict_file(cell_names_to_use)
    tissue_cell_assays = get_tissue_cell_mappings(tissue_cell_mappings_file=tissue_cell_mappings_file,
                                                  matching_cell_name_representative_dict=matching_cell_name_representative_dict, 
                                                  key_value_sep='=', values_sep=',', cell_assay_sepe=':') 
    
    if not os.path.exists(muts_out):
        muts_tracks_files = get_muts_tracks_info(muts_input_file=muts_input_file, tracks_dir=tracks_dir, muts_dir_out=muts_out+'_tmp')
        with open(muts_out, 'w') as muts_ofile:
            for muts_tracks_file in muts_tracks_files:
                muts_tracks_ouput_file = muts_tracks_file+"_annotated" 
                print("Annotating: ", muts_tracks_ouput_file)
                annotate_muts(muts_tracks_file, muts_tracks_ouput_file, tissue_cell_assays, matching_cell_name_representative_dict, cancer_type_index = 5, tracks_info_index = 9)
                with open(muts_tracks_ouput_file, 'r') as muts_tracks_ouput_ifile:
                    muts_ofile.write(muts_tracks_ouput_ifile.read())
        #if os.path.exists(muts_out+'_tmp'):
        #    shutil.rmtree(muts_out+'_tmp')
    return muts_out
    
if __name__ == '__main__':
    
    muts_input_file = sys.argv[1]
    tracks_dir = sys.argv[2]
    muts_out = sys.argv[3]#muts_input_file+'_annotations'
    get_annotated_muts(muts_input_file, tracks_dir, muts_out)