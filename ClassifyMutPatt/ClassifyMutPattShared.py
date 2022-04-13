# -*- coding: utf8 -*-
"""
ClassifyMutPatt - Kazuki Nakamae 2020
Software pipeline for classifing mutation type predicted by InDelPhi
https://github.com/KazukiNakamae/ClassifyMutPatt
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
sys.path.append('./inDelphi-HEK293-model/')
import hek293_inDelphi
hek293_inDelphi.init_model(celltype = 'HEK293')
sys.path.append('./inDelphi-U2OS-model/')
import u2os_inDelphi
u2os_inDelphi.init_model(celltype = 'U2OS')
sys.path.append('./inDelphi-HCT116-model/')
import hct116_inDelphi
hct116_inDelphi.init_model(celltype = 'HCT116')
sys.path.append('./inDelphi-K562-model/')
import k562_inDelphi
k562_inDelphi.init_model(celltype = 'K562')
import torch
import torch.nn as nn
kldiv = nn.KLDivLoss( reduction="sum" ) # function calculating KL divergence

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
