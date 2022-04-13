# -*- coding: utf8 -*-
"""
ExtendSeq - Kazuki Nakamae 2020
Software pipeline for extending sequence based on a genome reference
https://github.com/KazukiNakamae/ExtendSeq
(c) 2020 Hiroshima University. All Rights Reserved.
"""

__version__ = "1.0.0"

from ExtendSeq.ExtendSeqShared import *
# from MaChIAto.MaChIAtoClassifier import MaChIAtoClassifier, MaChIAtoPEClassifier
from ExtendSeq.ExtendSeqFunctions import create_multi_fasta, execute_blat, extend_seqs, clean_files

def main():
    try:
        print("  \n===\ExtendSeq/===")
        print("""
        -Software pipeline for extending sequence based on a genome reference-
        """)
        print("Author: Kazuki Nakamae at Sakuma Tetsushi and Takashi Yamamoto lab, Hiroshima Univ, Japan")
        print("For support contact kazukinakamae@gmail.com")
        info("ExtendSeq Version :" + __version__)

        parser = argparse.ArgumentParser(description="ExtendSeq Parameters",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        # Arguments
        parser.add_argument("-i", "--input", type=str,  help="This parameter allows for the specification of the input file that contains the original sequence.", required=True)
        parser.add_argument("-o", "--output", type=str,  help="This parameter allows for the specification of the output file that contains the extended sequence (default: output.csv).", default="output.csv", required=False)
        parser.add_argument("-s", "--sequence_column", type=str,  help="This parameter allows for the specification of the column name that contains sequences.", required=True)
        parser.add_argument("-u", "--unique_column", type=str,  help="This parameter allows for the specification of the column name that contains reference labels. Each label should be unique to merge processed values with the inout table.", required=True)
        parser.add_argument("-l", "--extended_left_length", type=int,  help="This parameter allows for the specification of the length that the software is to extend base toward the left  (default: 0).", default=0, required=False)
        parser.add_argument("-r", "--extended_right_length", type=int,  help="This parameter allows for the specification of the length that the software is to extend base toward the right direction (default: 0).", default=0, required=False)
        args = parser.parse_args()

        # Read the input table
        info("Read table...")
        input_pddf = pd.read_csv( args.input )

        # Make a hidden directory that is to contain the immediate files
        info("Prepare the system...")
        if not os.path.exists(".ExtendSeq_temp"):
            if subprocess.check_call([ "mkdir", ".ExtendSeq_temp" ]) != 0:
                raise ExecutionException( "Cannot make the .ExtendSeq_temp directory." )

        # Convert the table to FASTA file
        info("Extract sequences...")
        if not os.path.exists(".ExtendSeq_temp/output.fasta"):
            create_multi_fasta( input_pddf.loc[:, [ args.unique_column, args.sequence_column ] ] )

        # Perform Blat
        info("Search the coordinates of each sequence...")
        seq_len = int( input_pddf.loc[ 1, [ args.sequence_column ] ].str.len() )
        if not os.path.exists(".ExtendSeq_temp/output.psl"):
            if execute_blat( "./hg38/hg38.fa", seq_len ) != 0:
                raise ExecutionException('Cannot execute blat.')
        
        # Convert psl into bed
        if not os.path.exists(".ExtendSeq_temp/output.bed"):
            if subprocess.check_call([ "psltobed", ".ExtendSeq_temp/output.psl", ".ExtendSeq_temp/output.bed" ]) != 0:
                raise ExecutionException( "Cannot execute psltobed." )

        # Read bed files
        info("Read the coordinates...")
        bed_name_tpl = ("chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")
        bed_pddf = pd.read_table( ".ExtendSeq_temp/output.bed", names=bed_name_tpl )

        # Get Maximum length of each chromosome
        info("Estimate the length of reference sequence...")
        chr_len_dict = {}
        with open("./hg38/hg38.fa", "rU") as fh:
            for seq_record in SeqIO.parse(fh, "fasta"):
                chr_len_dict[seq_record.id] = len(seq_record.seq)

        # Modify bed format data
        info("Extend the coordinates...")
        bed_modified_pddf = extend_seqs( bed_pddf, args.extended_left_length, args.extended_right_length, chr_len_dict )

        # Write the modified bed file
        info("Extract the extended coordinates...")
        bed_modified_pddf.to_csv( ".ExtendSeq_temp/modified.bed", header=False, index=False, sep="\t" )

        # Make sequence file based on the modified bed file
        info("Extract the extended sequences...")
        if subprocess.check_call([ "bedtools", "getfasta", "-fi", "./hg38/hg38.fa", "-bed", ".ExtendSeq_temp/modified.bed", "-fo", ".ExtendSeq_temp/modified.tsv", "-name", "-tab", "-s" ]) != 0:
            raise ExecutionException( "Cannot execute bedtools getfasta." )
        
        # Read the modified bed data as pandas dataframe
        info("Read the extended sequences...")
        extended_pddf = pd.read_table( ".ExtendSeq_temp/modified.tsv", names=("name", "extended_sequence") )

        # Extract number for the identification
        extended_pddf[ "name" ] = pd.DataFrame( extended_pddf[ "name" ].str.extractall( "(^[0-9]+)\:\:.*" ) ).reset_index()[ 0 ]
        
        # Change type
        # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
        extended_pddf[ "name" ] = extended_pddf[ "name" ].astype("int64")

        # Remove duplicates (multi hit sequence)
        info("Remove the non-specific sequences in this process...")
        extended_pddf = extended_pddf.drop_duplicates(subset=["name"])

        # Merge and remove entry that does not exist in the extended_pddf
        info("Assemble the extended sequences...")
        merge_pddf = input_pddf.merge( extended_pddf, left_on=args.unique_column, right_on="name" )

        # Write output table
        info("Write the result...")
        merge_pddf.to_csv( args.output, index=False )

        # Remove the immediate files
        info("The process is ending...")
        if clean_files():
            info("Done.")

        sys.exit(0)
    except Exception as e:
        # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
        error("Error. Please check your input.\n\nERROR: %s" % e)
        # clean_files()
        sys.exit(-1)
    except InputException as e:
        #traceback.print_exc(file = sys.stdout)
        error("Input error. Please check your input.\n\nERROR: %s" % e)
        print("""
        usage:
        python ExtendSeq.py -i <INPUT> -o <OUTPUT> -n <NAME COLUMN> -s <SEQUENCE COLUMN> -r <REFERENCE COLUMN> -l <EXTEND LEFT LENGTH> -r <EXTENDED RIGHT LENGTH>
        """)
        clean_files()
        sys.exit(1)
    except ExecutionException as e:
        #traceback.print_exc(file = sys.stdout)
        error("Execution error. Please check your input.\n\nERROR: %s" % e)
        clean_files()
        sys.exit(1)
    except UnexpectedException as e:
        #traceback.print_exc(file = sys.stdout)
        error("Unexpected error. Please get touch in the author.[Kazuki Nakamae : kazukinakamae@gmail.com]\n\nERROR: %s" % e)
        # clean_files()
        sys.exit(1)
