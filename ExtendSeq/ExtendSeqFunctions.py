# -*- coding: utf8 -*-
"""
ExtendSeq - Kazuki Nakamae 2020
Software pipeline for classifing mutation type predicted by InDelPhi
https://github.com/KazukiNakamae/ExtendSeq
(c) 2020 Hiroshima University. All Rights Reserved.
"""

from ExtendSeq.ExtendSeqShared import *

### MODULES ############################
from Bio.Seq import Seq

### FUNCTIONS ############################

def create_multi_fasta( seq_pddf ):
    """Create a multiple fasta from the pandas dataframe.

    Make fasta format per row and merge them.

    Args:
        seq_pddf: An panda dataframe that contains the label and sequence.

    Returns:
        NULL

    """
    seq_list = []
    for row in seq_pddf.itertuples( name=None ):
        record = SeqRecord( Seq(row[2]), id=str(row[1]), description="" )
        seq_list.append(record)
        temp1 = Seq(row[2])
        temp2 = str(row[1])
    SeqIO.write(seq_list, ".ExtendSeq_temp/output.fasta", "fasta") # make output in the special hidden directory

def execute_blat( ref, seq_len ):
    """Execute Blat.

    Execute Blat using output.fasta. The output file is output.psl

    Args:
        ref: A reference sequence file (.fa).

    Returns:
        the value of subprocess.check_call()

    """
    return( subprocess.check_call(["blat", "-minMatch=0", "-minScore=" + str(seq_len), str(ref), ".ExtendSeq_temp/output.fasta", ".ExtendSeq_temp/output.psl"]) )

def extend_seqs( bed_pddf, extended_left_length, extended_right_length, chr_len_dict ):
    """Extend sequence in the input table.

    Modified values into the coordinates. The value change acoording to the strand that the sequence is on

    Args:
        bed_pddf: A pandas dataframe that contains bed format data.
        extended_left_length: The length that the software is to extend base toward the left.
        extended_right_length: The length that the software is to extend base toward the right.
        chr_len_dict: The dictionary: The key is a chromosome name. The value is a length of the chromosome.

    Returns:
        Modified pandas dataframe

    """
    def f_extend_seq( row, extended_left_length, extended_right_length, chr_len_dict ):
        if row[ "strand" ] == "-":
            extended_left_length, extended_right_length = extended_right_length, extended_left_length
        row[ "chromStart" ] = row[ "chromStart" ] - extended_left_length
        row[ "chromEnd" ] = row[ "chromEnd" ] + extended_right_length

        if row["chromStart"] < 0:
            warn("A extended coordinate is over the lenght of the chromosome. The start position is forced to be 0")
            row["chromStart"] = 0
        if row[ "chromEnd" ] > chr_len_dict[ row[ "chrom" ] ]:
            warn("A extended coordinate is over the lenght of the chromosome. The end position is forced to be " + chr_len_dict[ row[ "strand" ] ])
            row[ "chromEnd" ] = chr_len_dict[ row[ "chrom" ] ]
        return(row)
        
    bed_modified_pddf = bed_pddf.apply(f_extend_seq, axis='columns', extended_left_length=extended_left_length, extended_right_length=extended_right_length, chr_len_dict=chr_len_dict)
    return( bed_modified_pddf )
    
def clean_files():
    """clean immediate files.

    Execute rm -r to remove the immediate files

    Args:

    Returns:
        True: process properly ends.

    """
    try:
        # Remove the hidden directory.
        if os.path.exists(".ExtendSeq_temp") and ( subprocess.check_call([ "rm", "-r", ".ExtendSeq_temp" ]) != 0 ):
            raise ExecutionException( "Cannot remove the .ExtendSeq_temp directory." )
        return( True )

    except ExecutionException as e:
        #traceback.print_exc(file = sys.stdout)
        error("Execution error. Please check your input.\n\nERROR: %s" % e)
        sys.exit(1)
        

