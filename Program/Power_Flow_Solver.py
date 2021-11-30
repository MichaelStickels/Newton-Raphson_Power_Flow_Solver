
#   Power Flow Solver
#   Michael Stickels
#   10/28/2021
#
#   All values are in Per Unit unless another unit is stated


# Imports
import numpy as np
from numpy.core.numeric import empty_like
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

    # Imaginary part -> B(i,i) = G(i,i) + B_i/2
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

def P_k_Equation(k, v_t_mat):
    p_temp = 0

    for i in np.arange(num_Busses):
        
        p_temp += v_t_mat[k + num_Busses] * v_t_mat[i + num_Busses] * (matrix_Y_real[k,i] * np.cos(v_t_mat[k] - v_t_mat[i]) + matrix_Y_imaginary[k,i] * np.sin(v_t_mat[k] - v_t_mat[i]))
    
    return p_temp

def Q_k_Equation(k, v_t_mat):
    q_temp = 0

    for i in np.arange(num_Busses):
        
        q_temp += v_t_mat[k + num_Busses] * v_t_mat[i + num_Busses] * (matrix_Y_real[k,i] * np.sin(v_t_mat[k] - v_t_mat[i]) - matrix_Y_imaginary[k,i] * np.cos(v_t_mat[k] - v_t_mat[i]))
    
    
    return q_temp



# For rows i in Y matricies:
    # matrix_PQ_equations[i,1] = v_k * v_i * (matrix_Y_real[i,j] * cos (theta) + matrix_Y_imaginary[i,j] * sin (theta))                  ## P EQUATIONS  (v_k = this bus; v_i = opposite bus)
    # matrix_PQ_equations[i,1] += v_k * v_i * (matrix_Y_real[i,j] * cos (theta) + matrix_Y_imaginary[i,j] * sin (theta))                 ## P SUMATION
    # matrix_PQ_equations[i + num_Busses,1] = v_k * v_i * (matrix_Y_real[i,j] * sin (theta) - matrix_Y_imaginary[i,j] * cos (theta))     ## Q EQUATIONS
    # matrix_PQ_equations[i + num_Busses,1] += v_k * v_i * (matrix_Y_real[i,j] * sin (theta) - matrix_Y_imaginary[i,j] * cos (theta))   ## Q SUMATION

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Initialize Jacobian matrix
# matrix_J = np.zeros((num_Busses*2,num_Busses*2))

# For H quadrant:
# Off diagonals: 
# matrix_J[i,j] = v_k * v_i * (matrix_Y_real[i,j] * sin (theta) - matrix_Y_imaginary[i,j] * cos (theta))                 ## (v_k = this bus; v_i = opposite bus)

# Collector: collector += v_k * v_i * (matrix_Y_imaginary[i,j] * cos (theta) - matrix_Y_real[i,j] * sin (theta))                 ## (v_k = this bus; v_i = opposite bus)

# Diagonals: 
# matrix_J[i,i] = collector

def H_quadrant_equation(k, i, v_t_mat):
    
    H_temp = 0

    if(k == i):

        for x in np.arange(num_Busses):
            if(x != k):
                
                H_temp += v_t_mat[k + num_Busses] * v_t_mat[x + num_Busses] * (-matrix_Y_real[k,x] * np.sin(v_t_mat[k] - v_t_mat[x]) + matrix_Y_imaginary[k,x] * np.cos(v_t_mat[k] - v_t_mat[x]))

    else:

        H_temp = v_t_mat[k + num_Busses] * v_t_mat[i + num_Busses] * (matrix_Y_real[k,i] * np.sin(v_t_mat[k] - v_t_mat[i]) - matrix_Y_imaginary[k,i] * np.cos(v_t_mat[k] - v_t_mat[i]))
    
    return H_temp




# For M quadrant:
# Off diagonals:
# matrix_J[i,j + num_Busses] = v_k * (matrix_Y_real[i,j] * cos (theta) - matrix_Y_imaginary[i,j] * sin (theta))       ## (v_k = this bus; v_i = opposite bus)

# Collector: collector2 += v_i * (matrix_Y_real[i,j] * cos (theta) - matrix_Y_imaginary[i,j] * sin (theta))                  ## (v_k = this bus; v_i = opposite bus)

# Diagonals:
# matrix_J[i,i + num_Busses] = collector2 + 2 * matrix_Y_real[i,i] * v_k                  ## (v_k = this bus; v_i = opposite bus)

def M_quadrant_equation(k, i, v_t_mat):
    
    M_temp = 0

    if(k == i):

        for x in np.arange(num_Busses):
            if(x != k):
                
                M_temp += v_t_mat[x + num_Busses] * (matrix_Y_real[k,x] * np.cos(v_t_mat[k] - v_t_mat[x]) + matrix_Y_imaginary[k,x] * np.sin(v_t_mat[k] - v_t_mat[x]))

        M_temp +=  2 * matrix_Y_real[k,k] * v_t_mat[k + num_Busses]

    else:

        M_temp = v_t_mat[k + num_Busses] * (matrix_Y_real[k,i] * np.cos(v_t_mat[k] - v_t_mat[i]) + matrix_Y_imaginary[k,i] * np.sin(v_t_mat[k] - v_t_mat[i]))
    
    return M_temp



# For N quadrant:
# Off diagonals:
# matrix_J[i + num_Busses, j] = v_k * v_i * ((-1) * matrix_Y_real[i,j] * cos (theta) - matrix_Y_imaginary[i,j] * sin (theta))

# Collector: collector3 += v_k * v_i * (* matrix_Y_real[i,j] * cos (theta) + matrix_Y_imaginary[i,j] * sin (theta))

# Diagonals:
# matrix_J[i + num_Busses, i] = collector3

def N_quadrant_equation(k, i, v_t_mat):
    
    N_temp = 0

    if(k == i):

        for x in np.arange(num_Busses):
            if(x != k):
                
                N_temp += v_t_mat[k + num_Busses] * v_t_mat[x + num_Busses] * (-matrix_Y_real[k,x] * np.cos(v_t_mat[k] - v_t_mat[x]) - matrix_Y_imaginary[k,x] * np.sin(v_t_mat[k] - v_t_mat[x]))

    else:

        N_temp = v_t_mat[k + num_Busses] * v_t_mat[i + num_Busses] * (-matrix_Y_real[k,i] * np.cos(v_t_mat[k] - v_t_mat[i]) - matrix_Y_imaginary[k,i] * np.sin(v_t_mat[k] - v_t_mat[i]))
    
    return N_temp

# For L quadrant:
# Off diagonals:
# matrix_J[i + num_Busses, j + num_Busses] = v_k * (matrix_Y_real[i,j] * sin (theta) - matrix_Y_imaginary[i,j])

# Collector: collector4 += v_i * (matrix_Y_real[i,j] * sin (theta) - matrix_Y_imaginary[i,j])

# Diagonals:
# matrix_J[ i+ num_Busses, i + num_Busses] = collector4 - 2 * matrix_Y_imaginary * v_k


def L_quadrant_equation(k, i, v_t_mat):
    
    L_temp = 0

    if(k == i):

        for x in np.arange(num_Busses):
            if(x != k):
                
                L_temp += v_t_mat[x + num_Busses] * (matrix_Y_real[k,x] * np.sin(v_t_mat[k] - v_t_mat[x]) - matrix_Y_imaginary[k,x] * np.cos(v_t_mat[k] - v_t_mat[x]))

        L_temp += -2 * matrix_Y_imaginary[k,k] * v_t_mat[k + num_Busses]

    else:

        L_temp = v_t_mat[k + num_Busses] * (matrix_Y_real[k,i] * np.sin(v_t_mat[k] - v_t_mat[i]) - matrix_Y_imaginary[k,i] * np.cos(v_t_mat[k] - v_t_mat[i]))
    
    return L_temp


# print(matrix_J)

# Invert Jacobian Matrix

# matrix_Inverse_Jacobian = np.linalg.inv(matrix_J)

# print(matrix_Inverse_Jacobian)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Function to calculate PQ matrix
def PQ_Calculate(v_t_mat):
    pq_mat = np.zeros((num_Busses*2, 1))

    for k in np.arange(num_Busses):

        pq_mat[k] = P_k_Equation(k, v_t_mat)

        pq_mat[k + num_Busses] = Q_k_Equation(k, v_t_mat)

    return pq_mat

# Function to calculate J matrix
#
#   J = (H M)
#       (N L)
#
def J_Calculate(v_t_mat):
    J_mat = np.zeros((num_Busses*2, num_Busses*2))

    for a in np.arange(num_Busses):

        for b in np.arange(num_Busses):

            J_mat[a,b] = H_quadrant_equation(a, b, v_t_mat)
            J_mat[a,b+num_Busses] = M_quadrant_equation(a, b, v_t_mat)
            J_mat[a+num_Busses,b] = N_quadrant_equation(a, b, v_t_mat)
            J_mat[a+num_Busses,b+num_Busses] = L_quadrant_equation(a, b, v_t_mat)

    return J_mat



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# Initialize starting point for VT Matrix
V_T_matrix = np.zeros((num_Busses*2, 1), dtype=float)

for y in np.arange(num_Busses):
    V_T_matrix[y + num_Busses][0] = busData['V Set'][y]

# print(V_T_matrix)

# Set mismatch target
acceptable_mismatch = 0.1


# Initialize mismatch to arbitrarily high value
max_mismatch = float("inf")


# Calculate initial PQ Matrix
PQ_matrix = PQ_Calculate(V_T_matrix)

max_iterations = 20
iteration = 1


while(max_mismatch >= acceptable_mismatch and iteration <= max_iterations - 1):

    # Build Jacobian
    J_matrix = J_Calculate(V_T_matrix)

    # Invert Jacobian
    J_inverse = np.linalg.inv(J_matrix)

    # Calculate corrections
    delta_VT_matrix = np.matmul(-J_inverse, PQ_matrix)

    # Update V and T
    V_T_new = V_T_matrix + delta_VT_matrix

    # Calculate Mismatch
    PQ_new = PQ_Calculate(V_T_new)

    PQ_mismatch = np.abs(PQ_new - PQ_matrix)

    max_mismatch = np.amax(PQ_mismatch)
    print('Max Mismatch:', max_mismatch)

    PQ_matrix = PQ_new
    V_T_matrix = V_T_new

    iteration += 1

print(iteration)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Calculate P and Q Mismatches

# Initialize matrix_PQ_given (needs to be single dimension with P_n in top half and Q_n in bottom half)
# Read P and Q data from csv and input to matrix_PQ_given
# matrix_PQ_given = np.zeros((1,num_Busses*2))

# Combine PQ matrices via loop
# matrix_PQ_equations[i,j] = matrix_PQ_equations[i,j] - matrix_PQ_given[i,j]   ## matrix_PQ_equations comes from two sections above.

# Initialize matrix_deltaTheta_deltaV (needs to be same shape as matrix_PQ_equations)
# matrix_deltaTheta_deltaV = np.zeros((1,num_Busses*2))

# Flat Start (set all values in matrix_deltaTheta_deltaV to 1.0 and 0.0, respectively)

# Update matrix_PQ_quantities after flat start

# Initialize matrix_PQ_mismatch and matrix_convergence_log (needs 4 columns)
# matrix_PQ_mismatch = np.zeros((1,num_Busses*2))
# matrix_convergence_log = np.zeros((4,num_Busses*2))

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
