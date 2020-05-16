# merge background (reanalysis) and observation (regression estimates)

import numpy as np
import auxiliary as au
from matplotlib import pyplot as plt
from scipy import io
import os
import sys
import h5py
import time
import random
import datetime
from optimal_interpolation import OImerge

