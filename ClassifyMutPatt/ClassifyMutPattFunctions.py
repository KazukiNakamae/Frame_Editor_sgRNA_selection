# -*- coding: utf8 -*-
"""
ClassifyMutPatt - Kazuki Nakamae 2020
Software pipeline for classifing mutation type predicted by InDelPhi
https://github.com/KazukiNakamae/ClassifyMutPatt
(c) 2020 Hiroshima University. All Rights Reserved.
"""

from ClassifyMutPatt.ClassifyMutPattShared import *

### MODULES ############################
from Bio.Seq import Seq

### FUNCTIONS ############################
    
def clean_files():
    """clean immediate files.

    Execute rm -r to remove the immediate files

    Args:

    Returns:
        True: process properly ends.

    """
    try:
        # Remove the hidden directory.
        if os.path.exists(".ClassifyMutPatt_temp") and ( subprocess.check_call([ "rm", "-r", ".ClassifyMutPatt_temp" ]) != 0 ):
            raise ExecutionException( "Cannot remove the .ClassifyMutPatt_temp directory." )
        return( True )

    except ExecutionException as e:
        #traceback.print_exc(file = sys.stdout)
        error("Execution error. Please check your input.\n\nERROR: %s" % e)
        sys.exit(1)
        

