import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


twiss = pd.read_table('twiss_NICA_Gold.tbl', comment = '@', sep=" ")

print(twiss)
