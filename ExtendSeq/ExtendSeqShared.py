# -*- coding: utf8 -*-
"""
ExtendSeq - Kazuki Nakamae 2020
Software pipeline for classifing mutation type predicted by InDelPhi
https://github.com/KazukiNakamae/ExtendSeq
(c) 2020 Hiroshima University. All Rights Reserved.
"""

### MODULES ############################
import sys
import os
import re
import argparse
import pandas as pd
import numpy as np
import logging
import warnings
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

### MODULE SETTING ############################

# logging
logging.basicConfig(level=logging.INFO,
	format='%(levelname)-5s @ %(asctime)s:\n\t %(message)s \n',
	datefmt='%a, %d %b %Y %H:%M:%S',
	stream=sys.stderr,
	filemode="w"
	)
error   = logging.critical
warn    = logging.warning
debug   = logging.debug
info    = logging.info
"""
from Bio import pairwise2
import csv
import datetime
from GPyOpt.methods import BayesianOptimization
import matplotlib.pyplot as plt

import math
import os
# import re
import regex as re
from tqdm import tqdm



"""
### EXCEPTIONS ############################
class ReferenceRegexException(Exception):
    pass
class UnexpectedException(Exception):
  	pass
class InputException(Exception):
  	pass
class ExecutionException(Exception):
  	pass
