
#   Power Flow Solver
#   Michael Stickels
#   10/28/2021
#
#   All values are in Per Unit unless another unit is stated


# Imports
import numpy as np
import pandas as pd


# Constants & Parameters
S_BASE = 100 #MVA


# Read in excel input file
input = pd.read_excel("data/system_basecase.xlsx", sheet_name=None)
#print(input)
busData = input['BusData']
lineData = input['LineData']

# number of buses
num_Busses = busData.shape[0]
#print(num_Busses)

# Initialize admittance matrix
matrix_Y_real = np.zeros((num_Busses,num_Busses))
matrix_Y_imaginary = np.zeros((num_Busses,num_Busses))
#print(matrix_Y)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Build Y Matrix
# Y = G + jB

# Input line admittances (off-diagonal values) from input data
for ind in lineData.index:

    i = lineData['From'][ind] - 1
    j = lineData['To'][ind] - 1

    # Real part -> G(i,j) = (-R)/(R^2+X^2)
    matrix_Y_real[i,j] = (-lineData['Rtotal, p.u.'][ind])/(lineData['Rtotal, p.u.'][ind]**2 + lineData['Xtotal, p.u.'][ind]**2)

    # Imaginary part -> B(i,j) = (X)/(R^2+X^2)
    matrix_Y_imaginary[i,j] = (lineData['Xtotal, p.u.'][ind])/(lineData['Rtotal, p.u.'][ind]**2 + lineData['Xtotal, p.u.'][ind]**2)

# Sum rows/cols for imtermidiate bus admittance
for i in np.arange(0,matrix_Y_real.shape[0]):

    # Real part -> ?
    matrix_Y_real[i,i] = matrix_Y_real.sum(axis=1)[i]

    # Imaginary part -> ?
    matrix_Y_imaginary[i,i] = matrix_Y_imaginary.sum(axis=1)[i]

# Add shunt admittances to diagonals (bus admittances)
for ind in lineData.index:

    i = lineData['From'][ind] - 1
    j = lineData['To'][ind] - 1

    # Real part -> G(i,i) = G(i,i) + B_i/2
    matrix_Y_real[i,i] = matrix_Y_real[i,i] + lineData['Btotal, p.u.'][ind]/2
    matrix_Y_real[j,j] = matrix_Y_real[j,j] + lineData['Btotal, p.u.'][ind]/2

    # Imaginary part -> B(i,j) = (X)/(R^2+X^2)
    matrix_Y_imaginary[i,i] = matrix_Y_imaginary[i,i] + lineData['Btotal, p.u.'][ind]/2
    matrix_Y_imaginary[j,j] = matrix_Y_imaginary[j,j] + lineData['Btotal, p.u.'][ind]/2

# Export matrix Y to csv files to check
pd.DataFrame(matrix_Y_real).to_csv("output/matrix_Y_real.csv")
pd.DataFrame(matrix_Y_imaginary).to_csv("output/matrix_Y_imaginary.csv")

# print(matrix_Y_real)
# print(matrix_Y_imaginary)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Determine P and Q Equations at Buses



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Build Jacobian Matrix



# Invert Jacobian Matrix



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Calculate P and Q Mismatches



# Log Convergence For All Values to CSV



# Exit Function:    If P < 0.1 && Q < 0.1


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Determine P and Q At Generators



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Calculate Line Currents


# Calculate Line P, Q and S


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Build Ratings Matrix [Use Common Structure For All Matrices]


# Check Node Voltages Against Ratings


# Check Line Power Flows Against Ratings


# Check Generator Power Against Ratings     ???????????


# Flag Violations of Operating Limits / Settings


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
