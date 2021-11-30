
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
    matrix_Y_real[j,i] = (-lineData['Rtotal, p.u.'][ind])/(lineData['Rtotal, p.u.'][ind]**2 + lineData['Xtotal, p.u.'][ind]**2)

    # Imaginary part -> B(i,j) = (X)/(R^2+X^2)
    matrix_Y_imaginary[i,j] = (lineData['Xtotal, p.u.'][ind])/(lineData['Rtotal, p.u.'][ind]**2 + lineData['Xtotal, p.u.'][ind]**2)
    matrix_Y_imaginary[j,i] = (lineData['Xtotal, p.u.'][ind])/(lineData['Rtotal, p.u.'][ind]**2 + lineData['Xtotal, p.u.'][ind]**2)

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

# Need to combine Y-matrices together????

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Determine P and Q Equations at Buses

# Initialize matrix_PQ_equations (needs to be single dimension with P_n in top half and Q_n in bottom half)
matrix_PQ_equations = np.zeros((1,num_Busses*2))

# For rows i in Y matricies:
    # matrix_PQ_equations[1,i] = v_k * v_i * (matrix_Y_real[i,j] * cos (theta) + matrix_Y_imaginary[i,j] * sin (theta))                  ## P EQUATIONS  (v_k = this bus; v_i = opposite bus)
    # matrix_PQ_equations[1,i] += v_k * v_i * (matrix_Y_real[i,j] * cos (theta) + matrix_Y_imaginary[i,j] * sin (theta))                 ## P SUMATION
    # matrix_PQ_equations[1,i + num_Busses] = v_k * v_i * (matrix_Y_real[i,j] * sin (theta) - matrix_Y_imaginary[i,j] * cos (theta))     ## Q EQUATIONS
    # matrix_PQ_equations[1,i + num_Busses"] += v_k * v_i * (matrix_Y_real[i,j] * sin (theta) - matrix_Y_imaginary[i,j] * cos (theta))   ## Q SUMATION

# print(matrix_PQ_equations)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Initialize Jacobian matrix
matrix_Jacobian = np.zeros((num_Busses*2,num_Busses*2))

# For H quadrant:
# Off diagonals: 
# matrix_Jacobian[i,j] = v_k * v_i * (matrix_Y_real[i,j] * sin (theta) - matrix_Y_imaginary[i,j] * cos (theta))                 ## (v_k = this bus; v_i = opposite bus)

# Collector: collector += v_k * v_i * (matrix_Y_imaginary[i,j] * cos (theta) - matrix_Y_real[i,j] * sin (theta))                 ## (v_k = this bus; v_i = opposite bus)

# Diagonals: 
# matrix_Jacobian[i,i] = collector



# For M quadrant:
# Off diagonals:
# matrix_Jacobian[i,j + num_Busses] = v_k * (matrix_Y_real[i,j] * cos (theta) - matrix_Y_imaginary[i,j] * sin (theta))       ## (v_k = this bus; v_i = opposite bus)

# Collector: collector2 += v_i * (matrix_Y_real[i,j] * cos (theta) - matrix_Y_imaginary[i,j] * sin (theta))                  ## (v_k = this bus; v_i = opposite bus)

# Diagonals:
# matrix_Jacobian[i,i + num_Busses] = collector2 + 2 * matrix_Y_real[i,i] * v_k                  ## (v_k = this bus; v_i = opposite bus)



# For N quadrant:
# Off diagonals:
# matrix_Jacobian[i + num_Busses, j] = v_k * v_i * ((-1) * matrix_Y_real[i,j] * cos (theta) - matrix_Y_imaginary[i,j] * sin (theta))

# Collector: collector3 += v_k * v_i * (* matrix_Y_real[i,j] * cos (theta) + matrix_Y_imaginary[i,j] * sin (theta))

# Diagonals:
# matrix_Jacobian[i + num_Busses, i] = collector3

# For L quadrant:
# Off diagonals:
# matrix_Jacobian[i + num_Busses, j + num_Busses] = v_k * (matrix_Y_real[i,j] * sin (theta) - matrix_Y_imaginary[i,j])

# Collector: collector4 += v_i * (matrix_Y_real[i,j] * sin (theta) - matrix_Y_imaginary[i,j])

# Diagonals:
# matrix_Jacobian[ i+ num_Busses, i + num_Busses] = collector4 - 2 * matrix_Y_imaginary * v_k

# print(matrix_Jacobian)

# Invert Jacobian Matrix

matrix_Inverse_Jacobian = np.linalg.inv(matrix_jacobian)

# print(matrix_Inverse_Jacobian)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Calculate P and Q Mismatches

# Initialize matrix_PQ_given (needs to be single dimension with P_n in top half and Q_n in bottom half)
# Read P and Q data from csv and input to matrix_PQ_given
matrix_PQ_given = np.zeros((1,num_busses*2))

# Combine PQ matrices via loop
matrix_PQ_equations[i,j] = matrix_PQ_equations[i,j] - matrix_PQ_given[i,j]   ## matrix_PQ_equations comes from two sections above.

# Initialize matrix_deltaTheta_deltaV (needs to be same shape as matrix_PQ_equations)
matrix_deltaTheta_deltaV = np.zeros((1,num_Busses*2))

# Flat Start (set all values in matrix_deltaTheta_deltaV to 1.0 and 0.0, respectively)

# Update matrix_PQ_quantities after flat start

# Initialize matrix_PQ_mismatch and matrix_convergence_log (needs 4 columns)
matrix_PQ_mismatch = np.zeros((1,num_busses*2))
matrix_convergence_log = np.zeros((4,num_Busses*2))

# Determine largest mismatch in both P and Q
# Add largest mismatch of both P and Q and ij values to matrix_convergence_log

# While any values of matrix_PQ_mismatch > 0.1:
    # matrix_deltaTheta_deltaV =  matrix_Inverse_jacobian * matrix_PQ_quantities * (-1)
    # matrix_deltaTheta_deltaV += matrix_deltaTheta_deltaV
    # matrix_PQ_quantities = matrix_PQ_equations with matrix_deltaTheta_deltaV plugged in
    # determine largest mismatch in both P and Q
    # add largest mismatch of both P and Q and ij values to matrix_convergence_log

# Print matrix_convergence_log
# Print matrix_PQ_quantities
# Print 


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
