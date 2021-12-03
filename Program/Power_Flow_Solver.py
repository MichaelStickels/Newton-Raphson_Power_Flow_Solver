
#   Power Flow Solver
#   Michael Stickels, Mike Mayhew
#   10/28/2021
#
#   All values are in Per Unit unless another unit is stated




# @@@@@@@@@@@@@@@
#
# Assumes that bus data is in order of bus number
# Assumes that bus 1 is the slack bus and has no load attached
# Assumes that all other buses have a load attached and a given P_load and Q_load
#
# @@@@@@@@@@@@@@@



# Imports
import numpy as np
import pandas as pd


# Constants & Parameters
S_BASE = 100 # MVA
acceptable_mismatch = 0.1 # Highest allowable mismatch for power flow calculation



# Read in excel input file
input = pd.read_excel("data/system_basecase.xlsx", sheet_name=None)
busData = input['BusData']
lineData = input['LineData']



# number of buses
total_Busses = busData.shape[0]

# number of busses where a P load is given ('P MW' > 0)
P_Busses = busData[busData['P MW'] > 0].count()["P MW"]

# number of busses with a generator ('P Gen' > 0)
gen_Busses = busData[busData['P Gen'] > 0].count()["P Gen"]


# Print bus info
print('Total Busses:', total_Busses)
print('Active Load Busses:', P_Busses)
print('Generator Busses (not including slack bus):', gen_Busses)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Build Y Matrix
# Y = G + jB

# Initialize empty admittance matrix
matrix_Y_real = np.zeros((total_Busses,total_Busses))
matrix_Y_imaginary = np.zeros((total_Busses,total_Busses))

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
    matrix_Y_real[i,i] = matrix_Y_real.sum(axis=1)[i] * -1

    # Imaginary part -> ?
    matrix_Y_imaginary[i,i] = matrix_Y_imaginary.sum(axis=1)[i] * -1

# Add shunt admittances to diagonals (bus admittances)
for ind in lineData.index:

    i = lineData['From'][ind] - 1
    j = lineData['To'][ind] - 1

    # Imaginary part -> B(i,i) = G(i,i) + B_i/2
    matrix_Y_imaginary[i,i] = matrix_Y_imaginary[i,i] + lineData['Btotal, p.u.'][ind]/2
    matrix_Y_imaginary[j,j] = matrix_Y_imaginary[j,j] + lineData['Btotal, p.u.'][ind]/2


# Export matrix Y to csv files for troubleshooting
pd.DataFrame(matrix_Y_real).to_csv("output/matrix_Y_real.csv")
pd.DataFrame(matrix_Y_imaginary).to_csv("output/matrix_Y_imaginary.csv")



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Build matrix of given active and reactive power for mismatch calculation
given_PQ = np.zeros((P_Busses * 2 - gen_Busses, 1))

for x in np.arange(P_Busses):
    
    # given active power (P)
    given_PQ[x] = ((busData[busData['Type'] != 'S'])['P MW'].to_numpy())[x] - ((busData.loc[busData['Type'] != 'S'])['P Gen'].to_numpy())[x]


for x in np.arange(P_Busses - gen_Busses):

    # given reactive power (Q)
    given_PQ[x + P_Busses] = ((busData[busData['Type'] == 'D'])['Q MVAr'].to_numpy())[x]


# Initialize starting point for VT Matrix
#   All Thetas set to 0
#   Bus voltage set to given value from input
V_T_matrix = np.zeros((total_Busses * 2, 1), dtype=float)

for y in np.arange(total_Busses):

    V_T_matrix[y + P_Busses + 1][0] = busData['V Set'][y]




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Determine P and Q Equations at Buses

# Calculate P_k
def P_k_Equation(k, v_t_mat):
    p_temp = 0
    k = int(k)

    for i in np.arange(P_Busses):
        
        p_temp += v_t_mat[int(k + total_Busses)] * v_t_mat[int(i + total_Busses)] * (matrix_Y_real[k][i] * np.cos(v_t_mat[k] - v_t_mat[i]) + matrix_Y_imaginary[k][i] * np.sin(v_t_mat[k] - v_t_mat[i]))
    
    p_temp += busData['P MW'][k] - busData['P Gen'][k]

    return p_temp

# Calculate Q_k
def Q_k_Equation(k, v_t_mat):
    q_temp = 0
    k = int(k)

    for i in np.arange(total_Busses):
   
        q_temp += v_t_mat[int(k + total_Busses)] * v_t_mat[int(i + total_Busses)] * (matrix_Y_real[k][i] * np.sin(v_t_mat[k] - v_t_mat[i]) - matrix_Y_imaginary[k][i] * np.cos(v_t_mat[k] - v_t_mat[i]))
    
    q_temp += busData['Q MVAr'][k]
    
    return q_temp



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Define functions for calculating Jacobian values
#
#   J = (H M)
#       (N L)
#

def H_quadrant_equation(k, i, v_t_mat):
    
    H_temp = 0

    if(k == i):

        for x in np.arange(P_Busses):
            if(x != k):
                
                H_temp += v_t_mat[k + P_Busses] * v_t_mat[x + P_Busses] * (-matrix_Y_real[k,x] * np.sin(v_t_mat[k] - v_t_mat[x]) + matrix_Y_imaginary[k,x] * np.cos(v_t_mat[k] - v_t_mat[x]))

    else:

        H_temp = v_t_mat[k + P_Busses] * v_t_mat[i + P_Busses] * (matrix_Y_real[k,i] * np.sin(v_t_mat[k] - v_t_mat[i]) - matrix_Y_imaginary[k,i] * np.cos(v_t_mat[k] - v_t_mat[i]))
    
    print(v_t_mat[k + P_Busses])
    return H_temp


def M_quadrant_equation(k, i, v_t_mat):
    
    M_temp = 0

    if(k == i):

        for x in np.arange(P_Busses):
            if(x != k):
                
                M_temp += v_t_mat[x + P_Busses] * (matrix_Y_real[k,x] * np.cos(v_t_mat[k] - v_t_mat[x]) + matrix_Y_imaginary[k,x] * np.sin(v_t_mat[k] - v_t_mat[x]))

        M_temp +=  2 * matrix_Y_real[k,k] * v_t_mat[k + P_Busses]

    else:

        M_temp = v_t_mat[k + P_Busses] * (matrix_Y_real[k,i] * np.cos(v_t_mat[k] - v_t_mat[i]) + matrix_Y_imaginary[k,i] * np.sin(v_t_mat[k] - v_t_mat[i]))
    
    return M_temp


def N_quadrant_equation(k, i, v_t_mat):
    
    N_temp = 0

    if(k == i):

        for x in np.arange(P_Busses):
            if(x != k):
                
                N_temp += v_t_mat[k + P_Busses] * v_t_mat[x + P_Busses] * (-matrix_Y_real[k,x] * np.cos(v_t_mat[k] - v_t_mat[x]) - matrix_Y_imaginary[k,x] * np.sin(v_t_mat[k] - v_t_mat[x]))

    else:

        N_temp = v_t_mat[k + P_Busses] * v_t_mat[i + P_Busses] * (-matrix_Y_real[k,i] * np.cos(v_t_mat[k] - v_t_mat[i]) - matrix_Y_imaginary[k,i] * np.sin(v_t_mat[k] - v_t_mat[i]))
    
    return N_temp


def L_quadrant_equation(k, i, v_t_mat):
    
    L_temp = 0

    if(k == i):

        for x in np.arange(P_Busses):
            if(x != k):
                
                L_temp += v_t_mat[x + P_Busses] * (matrix_Y_real[k,x] * np.sin(v_t_mat[k] - v_t_mat[x]) - matrix_Y_imaginary[k,x] * np.cos(v_t_mat[k] - v_t_mat[x]))

        L_temp += -2 * matrix_Y_imaginary[k,k] * v_t_mat[k + P_Busses]

    else:

        L_temp = v_t_mat[k + P_Busses] * (matrix_Y_real[k,i] * np.sin(v_t_mat[k] - v_t_mat[i]) - matrix_Y_imaginary[k,i] * np.cos(v_t_mat[k] - v_t_mat[i]))
    
    return L_temp



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# Function to add delta V and T to VT Matrix
def update_VT(v_t_calc_mat, v_t):

    # print(v_t_calc_mat)
    # print(v_t)

    # Put changes into full VT matrix
    for x in np.arange(P_Busses):
        
        v_t[int(v_t_calc_mat[x][0])] = v_t_calc_mat[x][1]


    # print(v_t)

    for x in np.arange(P_Busses - gen_Busses):

        # print(x)
        # print(int(v_t_calc_mat[x + P_Busses][0]))
        # print(v_t_calc_mat[x + P_Busses][1])
        # print(x + P_Busses)
        # print(v_t[int(v_t_calc_mat[x + P_Busses][0])])

        v_t[int(v_t_calc_mat[x + P_Busses][0] + P_Busses + 1)] = v_t_calc_mat[x + P_Busses][1]

    # print(v_t)

    return v_t



# Function to calculate PQ matrix
def PQ_Calculate(v_t_calc_mat, v_t_mat):

    
    pq_mat = np.zeros((P_Busses*2 - gen_Busses, 1))

    for y in np.arange(P_Busses):

        pq_mat[y] = P_k_Equation(y + 1, v_t_mat)

    for z in np.arange(P_Busses - gen_Busses):

        pq_mat[z + P_Busses] = Q_k_Equation(v_t_calc_mat[z + P_Busses][0], v_t_mat)

    return pq_mat



# Function to calculate J matrix
#
#   J = (H M)
#       (N L)
#
def J_Calculate(v_t_mat):

    J_mat = np.zeros((P_Busses * 2 - gen_Busses, P_Busses * 2 - gen_Busses))

    for a in np.arange(P_Busses):

        for b1 in np.arange(P_Busses):

            J_mat[a,b1] = H_quadrant_equation(a + 1, b1 + 1, v_t_mat)

        for b2 in np.arange(P_Busses - gen_Busses):

            J_mat[a,b2 + P_Busses] = M_quadrant_equation(a + 1, b2 + 1, v_t_mat)
            
    
    for a in np.arange(P_Busses - gen_Busses):

        for b1 in np.arange(P_Busses):

            J_mat[a + P_Busses,b1] = N_quadrant_equation(a + 1, b1 + 1, v_t_mat)

        for b2 in np.arange(P_Busses - gen_Busses):

            J_mat[a + P_Busses,b2 + P_Busses] = L_quadrant_equation(a + 1, b2 + 1, v_t_mat)
            

    return J_mat



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Set start values for power flow calculation


# gen_Busses = busData[busData['P Gen'] > 0].count()["P Gen"]

# Initialize mismatch to arbitrarily high value
max_mismatch = acceptable_mismatch + 10


# some fun parameters for a rainy day
max_iterations = 20
iteration = 1


# Define start points
PQ_matrix = np.copy(given_PQ)

V_T_calc_matrix = np.zeros((P_Busses * 2 - gen_Busses, 2))

for a in np.arange(P_Busses):
    
    V_T_calc_matrix[a][0] = a + 1

for a in np.arange(P_Busses - gen_Busses):
    
    V_T_calc_matrix[a + P_Busses][0] = ((busData[busData['Type'] == 'D'])['Bus #'].to_numpy())[a] - 1

    V_T_calc_matrix[a + P_Busses][1] = 1

# print(V_T_matrix)
# print(V_T_calc_matrix)


# Dataframe to save convergence history
convergence_history = np.array(['Iterations:', 'Max Mismatch:'])



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Newton-Raphson algorithm

while(max_mismatch >= acceptable_mismatch and iteration < max_iterations):

    # print(PQ_matrix)
    # print(V_T_matrix)

    # Build Jacobian
    J_matrix = J_Calculate(V_T_matrix)
    
    # Save Jacobian to CSV and halt
    # pd.DataFrame(J_matrix).to_csv("output/Jacobian.csv")
    # exit()

    # Invert Jacobian
    J_inverse = np.linalg.inv(J_matrix)

    # Calculate corrections
    delta_VT_matrix = np.matmul(-J_inverse, PQ_matrix)

    # print(delta_VT_matrix)


    # Update V and T
    V_T_calc_new = np.copy(V_T_calc_matrix)
    for a in np.arange(P_Busses * 2 - gen_Busses):
        V_T_calc_new[a][1] += delta_VT_matrix[a]

    # print(V_T_calc_new)

    V_T_new = update_VT(V_T_calc_new, V_T_matrix)

    # print(V_T_new)
    # exit()

    # Calculate Mismatch
    PQ_new = PQ_Calculate(V_T_calc_new, V_T_new)

    PQ_mismatch = np.abs(PQ_new)

    max_mismatch = np.amax(PQ_mismatch)

    print('Iteration:', iteration)
    print('Max Mismatch:', max_mismatch)

    convergence_history = np.vstack([convergence_history, [iteration, max_mismatch]])
    
    # print(V_T_new)

    PQ_matrix = PQ_new
    V_T_matrix = V_T_new

    iteration += 1


# export covnergence history CSV
pd.DataFrame(convergence_history).to_csv("output/Convergence History.csv")

# export final bus results
pd.DataFrame(PQ_matrix).to_csv("output/Final PQ Matrix.csv")
pd.DataFrame(V_T_matrix).to_csv("output/Final VT Matrix.csv")


print(PQ_matrix)
print(V_T_matrix)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Calculate Line Flow P, Q and S

# Funtion to calculate active power in line k to i
def active_line_power(k, i, v_t_mat):

    p = v_t_mat[k + total_Busses - 1] * v_t_mat[k + total_Busses - 1] * (matrix_Y_real[k - 1][i - 1] * np.cos(v_t_mat[k - 1] - v_t_mat[i -1]) + matrix_Y_imaginary[k - 1][i - 1] * np.sin(v_t_mat[k - 1] - v_t_mat[i -1]))

    return p

# Function to calculate reactive power in line k to i
def reactive_line_power(k, i, v_t_mat):

    q = v_t_mat[k + total_Busses - 1] * v_t_mat[k + total_Busses - 1] * (matrix_Y_real[k - 1][i - 1] * np.sin(v_t_mat[k - 1] - v_t_mat[i -1]) - matrix_Y_imaginary[k - 1][i - 1] * np.cos(v_t_mat[k - 1] - v_t_mat[i -1]))

    return q


# line_flows = [P, Q, S], export to CSV
line_flows = np.zeros((len(lineData), 3))

for a in np.arange(len(lineData)):

    line_flows[a][0] = active_line_power(lineData['From'][a], lineData['To'][a], V_T_matrix)

    line_flows[a][1] = reactive_line_power(lineData['From'][a], lineData['To'][a], V_T_matrix)

    line_flows[a][2] = np.sqrt(line_flows[a][0]**2 + line_flows[a][1]**2)

pd.DataFrame(line_flows).to_csv("output/Line Power Flows.csv")



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Build Ratings Matrix [Use Common Structure For All Matrices]



# Check Node Voltages Against Ratings



# Check Line Power Flows Against Ratings



# Flag Violations of Operating Limits / Settings



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
