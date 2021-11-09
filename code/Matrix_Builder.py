
#   Y Matrix Builder
#   Michael Stickels
#   10/28/2021
#


# Imports
#import csv
import numpy as np


with open("Line_Data.csv") as file_name:
    array = np.loadtxt(file_name, delimiter=",")

print(array)