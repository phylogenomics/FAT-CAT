'''
A utility function to return percent id based closest match to every sequence in the msa.

Call this function by
>> from validation_scripts import get_percent_identity as gpi

To run the whole alignment file
>> output = gpi.main("gpcr_nr95gc50_mafft.fasta")
The output should look like this
[(id of seq1 in msa, id of closest id match to seq1 in msa, percent identity),
...
]

To run a single id against the whole file
>> output = gpi.get_top_scoring_percent_identity_by_id(r"gi|110832828|sp|Q95136.2|DRD1_BOVIN/1-446", "gpcr_nr95gc50_mafft.fasta")

'''

import os
import os
import glob
import subprocess
from Bio import AlignIO

def get_top_scoring_percent_identity_by_sequence(seq, others_alignio_obj):
    """
    Given a sequence and other records (from alignio) as arguments,
    return the top_scoring_id and the top_score
    """
    top_scoring_id = ""
    top_score = 0.0
    for other in others_alignio_obj:
        other_pct_identity = get_percent_identity_two_sequences(seq, other.seq)
        if other_pct_identity > top_score:
            top_score = other_pct_identity
            top_scoring_id = other.id
    return top_scoring_id, top_score

def get_percent_identity_two_sequences(seq1, seq2):
    """
    Given two sequences, ignore characters where there is an _ in each sequence,
    Sum all remaining to total (total),
    Sum perfect matches (matches)
    return % match
    """
    matches = 0.0; total = 0.0
    for ind, seq1char in enumerate(seq1):
        seq2char = seq2[ind]
        if seq1char == "_" and seq2char == "_":
            continue
        total += 1.0
        if seq1char == seq2char:
            matches += 1.0
    return (matches/total * 100)

def get_top_scoring_percent_identity_by_id(seq_id, aln_file):
    '''
    Given an input sequence_id list containing [top_scoring_id, top_score]
    '''
    alignment = AlignIO.read(open(aln_file), "fasta")
    this_record = seq_id; top_scoring_id = ""; top_score = 0.0
    for record in alignment:
        if record.id == this_record:
            this_seq = record.seq
            other_records = [other for other in alignment if not other.id == this_record]
            top_scoring_id, top_score = get_top_scoring_percent_identity_by_sequence(this_seq, other_records)
    return [top_scoring_id, top_score]

def main(aln_file):
    '''
    Given an input file, return a list of lists. Inner lists have [test_id, top_scoring_id, top_score]
    '''
    outer_list = []
    alignment = AlignIO.read(open(aln_file), "fasta")
    for record in alignment:
        this_record = record.id
        this_seq = record.seq
        other_records = [other for other in alignment if not other.id == this_record]
        top_scoring_id, top_score = get_top_scoring_percent_identity_by_sequence(this_seq, other_records)
        outer_list += [[this_record, top_scoring_id, top_score]]
    return outer_list

if __name__ == "__main__":
    alignment_file = '/home/ajith/Documents/Write/Programming/Python/SjolanderLab/FAT_CAT/gpcr_nr95gc50_mafft.fasta'
    #main(alignment_file)
    print get_top_scoring_percent_identity_by_id(r"gi|1168246|sp|P35348.2|ADA1A_HUMAN/1-466", alignment_file)
