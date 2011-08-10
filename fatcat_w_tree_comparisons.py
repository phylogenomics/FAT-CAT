import os
import glob
import subprocess
from Bio import SeqIO
from validation_scripts import compare_two_trees, PID

#Global variables for original destinations
ORIGINAL_DATASET_FOLDER = os.path.join('FAT-CAT-1_dataset_validation', 'original_dataset')
NEWICK_FILES_DIR = os.path.join(ORIGINAL_DATASET_FOLDER, 'gpcr_trees')
SEQUENCE_FILES_DIR = os.path.join(ORIGINAL_DATASET_FOLDER, 'insert_seqs')
ALN_FILES_DIR = os.path.join(ORIGINAL_DATASET_FOLDER, 'gpcr_mafft_alignments')
ORIGINAL_MSA = os.path.join(ORIGINAL_DATASET_FOLDER, 'gpcr_nr95gc50_mafft.fasta')

FINAL_DESTINATION = 'FAT-CAT_Output'

def make_output_folder(folder_name):
    try:
        os.mkdir(folder_name)
    except:
        pass

def get_leaf_name_from_seq_file(seq_file):
    handle = open(seq_file, "rU")
    record = SeqIO.read(handle, "fasta")
    recordid = record.id
    return recordid

def run_comparisons(orig_tree, final_tree, inserted_leaf_name):
    rf_metric = 'NA'; rf_distance = 'NA'; sopd_distance = 'NA',
    try:
        rf_metric = (compare_two_trees.calculate_rf_metric(
            (orig_tree, final_tree), "/"))/2
    except:
        pass
    try:
        rf_distance = compare_two_trees.calculate_rf_distance(
            (orig_tree, final_tree), "/")
    except:
        pass
    try:
        sopd_distance = compare_two_trees.calculate_distance_of_one_leaf_between_trees(
            (orig_tree,final_tree),inserted_leaf_name)
    except:
            pass
    return rf_metric, rf_distance, sopd_distance

def main():
    final_string = []
    hpid = PID.Highest_PID_Match(ORIGINAL_MSA)
    final_string.append('Deleted Sequence\tRF metric\tRF distance\tSOPD\tRF metric\tRF distance\tSOPD\tRF metric\tRF distance\tSOPD\tRF metric\tRF distance\tSOPD\tClosest Sequence (pwid)\tPWID score')
    for newickFile in glob.glob(NEWICK_FILES_DIR+ '/*'):
        file_name = os.path.basename(newickFile)
        file_name, ext = os.path.splitext(file_name)
        if not ext == '.newick': continue
        output_folder = file_name
        #make_output_folder(output_folder)
        tree_file = '%s/%s.newick' % (NEWICK_FILES_DIR, file_name) # Repeated for flow
        msa_file = '%s/%s.aln' % (ALN_FILES_DIR, file_name)
        seq_file = '%s/%s.seq' % (SEQUENCE_FILES_DIR, file_name)
        subprocess.check_call(['/opt/local/bin/python', 'fatcat-1.py', '-t', tree_file,
                               '-a', msa_file, '-s', seq_file,
                               '-o', output_folder,'-v'])
        orig_tree = 'gpcr_original.newick'
        final_tree = '%s_fat_cat/final/tree_plus_o.newick' % file_name
        tree_plus = '%s_fat_cat/tree/tree_plus.newick' % file_name
        tree_n =  '%s_fat_cat/tree/tree_n.newick' % file_name
        tree_s ='%s_fat_cat/tree/tree_s.newick' % file_name
        inserted_leaf_name = get_leaf_name_from_seq_file(seq_file)
        closest_seq_details = hpid.get_hpid_by_id(inserted_leaf_name)
        this_file_string = file_name;
        for tree in [final_tree,tree_plus,tree_n,tree_s]:
            rf_metric, rf_distance, sopd_distance = run_comparisons(orig_tree, tree, inserted_leaf_name)
            for var in [rf_metric, rf_distance, sopd_distance]:
                this_file_string = this_file_string+'\t'+str(var)
        this_file_string += ("\t%s\t%s" %
                             (closest_seq_details['Highest scored id'],
                              closest_seq_details['Highest score']))
        final_string.append(this_file_string)
        print 'done'
    final_string = '\n'.join(final_string)
    f = open("FATCAT_STATS", 'w'); f.write(final_string); f.close()
    print 'Wrapper done'
    
if __name__ == "__main__":
    main()
    
