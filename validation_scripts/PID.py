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
>> output = gpi.get_hpid_by_id(r"gi|110832828|sp|Q95136.2|DRD1_BOVIN/1-446", "gpcr_nr95gc50_mafft.fasta")

'''

import os
import os
import glob
import subprocess
from Bio import AlignIO

class Highest_PID_Match():
    '''
    Base class. Call this with
    >> from validation_scripts import PID
    '''
    def __init__(self, aln_file):
        """
        Initial parameters are the overall alignment file. Call this with
        >> hpid = PID.Highest_PID_Match('gpcr_nr95gc50_mafft.fasta')
        """
        self.aln_file = aln_file
        self.highest_pid_dict = self.get_pid_dict()

    def get_pid_dict(self):
        '''
        Given an input file, return a dictionary with the ids as keys
        '''
        return_dict = {}
        alignment = AlignIO.read(open(self.aln_file), "fasta")
        for record in alignment:
            this_record = record.id
            this_seq = record.seq
            other_records = [other for other in alignment if not other.id == this_record]
            top_scoring_id, top_score = self.get_hpid_by_sequence(this_seq, other_records)
            return_dict[this_record] = {'Highest scored id' : top_scoring_id,
                                        'Highest score' : top_score}
        return return_dict
    
    def get_hpid_by_sequence(self, seq, others_alignio_obj):
        """
        Given a sequence and other records (from alignio) as arguments,
        return the top_scoring_id and the top_score
        """
        top_scoring_id = ""
        top_score = 0.0
        for other in others_alignio_obj:
            other_pct_identity = self.get_pid_two_seqs(seq, other.seq)
            if other_pct_identity > top_score:
                top_score = other_pct_identity
                top_scoring_id = other.id
        return top_scoring_id, top_score

    def get_pid_two_seqs(self, seq1, seq2):
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

    def get_hpid_by_id(self, seq_id):
        '''
        Given an input sequence_id return a dict containing the Highest scored id and the highest score
        '''
        return self.highest_pid_dict[seq_id]



if __name__ == "__main__":
    alignment_file = '/home/ajith/Documents/Write/Programming/Python/SjolanderLab/FAT_CAT/gpcr_nr95gc50_mafft.fasta'
    #main(alignment_file)
    hpid = Highest_PID_Match(alignment_file)
    print hpid.get_hpid_by_id('gi|110832828|sp|Q95136.2|DRD1_BOVIN/1-446')
