# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 12:21:28 2023

@author: repha
"""

import pandas as pd
import numpy as np
import seaborn as sns
import time
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog
root = tk.Tk()
root.withdraw()
root.attributes("-topmost", True)
directory = filedialog.askdirectory()
root.destroy()

kansas = pd.read_csv(directory + '/Deliverables_Ti02_New.csv')

sns.scatterplot(data=kansas, x="SF", y="BlN",hue = "Slip Plane", palette="deep")
plt.ylabel('BlN Intensity (px)')
plt.xlabel('Schmid Factor')
plt.title('Ti_02')
plt.ylim(0,3.1)
plt.xlim(0,0.53)
plt.show()






