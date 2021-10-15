#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 19:21:10 2021

@author: mikewang
"""


import time
import matplotlib.pyplot as plt
import numpy as np
#%%
#############################################################################################################
#                                    original choleski decompoistion                                        #
#############################################################################################################
def choleski_direst_method(A, b):
    # implement determine
    L = choleski_decomp_o(A)
    L_T = two_d_transpose(L)
    y = find_y(L,b)
    x = find_x(L_T,y)
    return x


def choleski_decomp_o(A): 
    dim_A = len(A)
    L_array = zeros(dim_A, dim_A)
    for j in range(dim_A):
        L_array[j][j] = np.sqrt(A[j][j] - L_square_sum(L_array,j))
        for i in range(j+1, dim_A):
            L_array[i][j] = (A[i][j] - L_multi_sum(L_array,j,i))/L_array[j][j]
    return L_array

def find_y(L, b):
    dim_L = L.shape
    y = zeros(dim_L[0],1)
    for i in range(dim_L[0]):
        y[i][0] = (b[i][0]- L_y_sum(L, y, i))/L[i][i]
    return y

def find_x(L_T, y):
    dim_L_T = L_T.shape
    x = zeros(dim_L_T[0],1)
    for i in range(dim_L_T[0]-1,-1,-1):
        x[i][0] = (y[i][0] -  L_x_sum(L_T, x, i, dim_L_T))/L_T[i][i]
    return x





#-------------------------------- Helper method for choleski_decomp_o --------------------------------
def L_square_sum(L_array, j):
    square_sum = 0
    for h in range(j):
        if j-1 < 0:
            square_sum = 0
        else:
            square_sum = square_sum + L_array[j][h] * L_array[j][h]
    return square_sum

def L_multi_sum(L_array, j , i):
    square_sum = 0
    for h in range(j):
        if j-1 < 0:
            square_sum = 0
        else:
            square_sum = square_sum + L_array[i][h] * L_array[j][h]
    return square_sum




#-------------------------------- Helper method for choleski_direst_method --------------------------------
def L_y_sum(L, y, i):
    sumation = 0
    count = 0
    for j in range(i):
        Ly = L[i][j]*y[j]
        if count == 0:
            sumation = Ly
        if count > 0:
            sumation = sumation + Ly
        count += 1
    return sumation

def L_x_sum(L_T, x, i, dim_L_T):
    sumation = 0
    count = 0
    for j in range(i+1, dim_L_T[0]):
        Lx = L_T[i][j]*x[j][0]
        if count == 0:
            sumation = Lx
        if count > 0:
            sumation = sumation + Lx
        count += 1
    return sumation

#############################################################################################################
#                                           general helper method                                           #
#############################################################################################################
# construct a all zero entry array with shape (dim1, dim2)
def zeros(dim1, dim2):
    overall_array = []
    for m in range(dim1):
        row_list = []
        for n in range(dim2):
            row_list.append(0.0)
        if len(overall_array) == 0:
            overall_array = np.array(row_list)[None]
        elif len(overall_array) > 0:
            overall_array = np.concatenate((overall_array, np.array(row_list)[None]), axis =0)
    return np.array(overall_array)

# sum of the all entries inside an array
def sum_array(array):
    sum_a = 0.0
    count = 0
    array_linear = array.flatten()
    for i in array_linear:
        if count == 0:
            sum_a = i
        elif count > 0:
            sum_a = sum_a + i
        count =+ 1
    return sum_a

# computes the summation of two 2-d array
def two_d_matrix_summation(A, B):
    dim_A = A.shape
    dim_B = A.shape
    sumed_array = zeros(A.shape[0], A.shape[1])
    if dim_A != dim_B:
        print("************** the dimension of the two entry array does not match! ***************")
    else: 
        for i in range(dim_A[0]):
            for j in range(dim_B[1]):
                sumed_array[i][j] = A[i][j]+B[i][j]
    return sumed_array

# computes the dot product of two 2-d array
def two_d_dot_product(A, B):
    dim_A = A.shape
    dim_B = B.shape
    resulted_matrix = zeros(dim_A[0], dim_B[1])
    if dim_A[1] != dim_B[0]:
        print("************** the inner size of the two array don't match! ***************")
    for i in range(dim_A[0]):
        for j in range(dim_B[1]):
            mult_inid_array = []
            A_row = A[i,:]
            B_colomn = B[:,j]
            for m in range(len(A_row)):
                mult_inid_array.append(A_row[m]*B_colomn[m])
            resulted_matrix[i][j]=sum_array(np.array(mult_inid_array))
    return resulted_matrix

# compute the transposation of an 2-d array
def two_d_transpose(A):
    dim_A = A.shape
    transposed_matrix = zeros(dim_A[1],dim_A[0])
    for i in range(dim_A[0]):
        for j in range(dim_A[1]):
            transposed_matrix[j][i]=A[i][j]
    return transposed_matrix


'''
============================================ Question 2 ========================================
Question setting: 
    Take a regular N by 2N finite-difference mesh and replace each horizontal and 
    vertical line by a 1 k resistor. This forms a linear, resistive network.
================================================================================================
'''
'''
Q2 a) Using the program you developed in question 1, find the resistance, R, between the node 
    at the bottom left corner of the mesh and the node at the top right corner of the mesh, 
    for N = 2, 3, …, 10. (You will probably want to write a small program that generates the 
    input file needed by the network analysis program. Constructing by hand the incidence matrix
    for a 200-node network is rather tedious). 
'''
print("------------------------ Q2 part a -----------------------------")
#############################################################################################################
#                                mesh resistence generator function                                         #
#############################################################################################################
'''
Here we are considering N as the number of nodes
'''
def mesh_resistence_generator(N,resistence,test_voltage):
    start_time = time.time()
    AyA, b = AyA_b_generator(N,resistence,test_voltage)
    v_n = choleski_direst_method(AyA,b)
    test_resistence_voltage = v_n[0]
    overall_mesh_resistence = resistence*(test_resistence_voltage/(test_voltage-test_resistence_voltage))
    time_used = time.time()-start_time
    return overall_mesh_resistence, time_used




def AyA_b_generator(N,resistence,test_voltage):
    resistence = float(resistence)
    voltage = float(test_voltage)
    A = A_generator(N)
    y = y_generator(N,resistence)
    J = J_generator(N)
    E = E_generator(N,voltage)
    AyA = two_d_dot_product(A,two_d_dot_product(y,two_d_transpose(A)))
    b = two_d_dot_product(A,(J-two_d_dot_product(y,E)))
    return AyA, b

#-------------------------------- Helper method to generate the mesh --------------------------------
def A_generator(N):
    num_row_A = N*2*N -1 #since we are taking one node grounded
    num_coloumn_A = ((N-1)*2*N)+((2*N-1)*N)+1
    num_vertical_branch_row = N-1
    num_vertical_branch_coloumn = 2*N
    num_horziontal_branch_row = N
    num_horziontal_branch_coloumn = 2*N-1
    
    A_temp = zeros(num_row_A+1,num_coloumn_A)
    A = zeros(num_row_A,num_coloumn_A)
    A_temp[0][0] = -1 # the testing branch

    '''
     iterate the branch in a seqence:
         0, the testing branch
         1, from bottom up, start from the left bottom conner vertical branch
         2, interates up-ward until reaching the top
         3, move to it's right neighbouring horziontal branches also starts from the buttom until reach the top
         4, move to it's right neighbouring vertical branches
         5, repeat step 1 to 4 for 2*N-1 times
         6, interates the right most vertical branches from the bottum to the top
    '''
    #print(A_temp)
    for m in range(num_horziontal_branch_coloumn):
        starting_colomn_coord_in_A = m*(num_vertical_branch_row + num_horziontal_branch_row)+1
        starting_row_coord_in_A = m*N
        for i in range(num_vertical_branch_row):
            A_temp[starting_row_coord_in_A+i][starting_colomn_coord_in_A+i] = 1
            A_temp[starting_row_coord_in_A+i+1][starting_colomn_coord_in_A+i] = -1 
        for j in range(num_horziontal_branch_row):
            A_temp[starting_row_coord_in_A+j][starting_colomn_coord_in_A+num_vertical_branch_row+j] = 1
            A_temp[starting_row_coord_in_A+j+N][starting_colomn_coord_in_A+num_vertical_branch_row+j] = -1
    # for the right most colomn branch 
    m = m+1
    starting_colomn_coord_in_A = m*(num_vertical_branch_row + num_horziontal_branch_row)+1
    starting_row_coord_in_A = m*N
    for i in range(num_vertical_branch_row):
              A_temp[starting_row_coord_in_A+i][starting_colomn_coord_in_A+i] = 1
              A_temp[starting_row_coord_in_A+i+1][starting_colomn_coord_in_A+i] = -1
    A = A_temp[:num_row_A,:]
    return A

# generate y matrix
def y_generator(N, resistence):
    num_resistor = ((N-1)*2*N)+((2*N-1)*N)+1
    y = zeros(num_resistor,num_resistor)
    for i in range(num_resistor):
        y[i][i] = 1/resistence
    return y

def J_generator(N):
    num_resistor = ((N-1)*2*N)+((2*N-1)*N)+1
    J = zeros(num_resistor,1)
    return J

def E_generator(N,voltage):
    num_resistor = ((N-1)*2*N)+((2*N-1)*N)+1
    E= zeros(num_resistor,1)
    E[0,0] = voltage
    return E

resistence = 1000
test_voltage = 10
# N = 2 
r_N_2, time_N_2 = mesh_resistence_generator(2,resistence,test_voltage)
print("overall resistence for a N*2N mesh resistor for N = 2 is " + str(r_N_2))
print("time it takes " + str(time_N_2))
# N = 3 
r_N_3, time_N_3 = mesh_resistence_generator(3,resistence,test_voltage)
print("overall resistence for a N*2N mesh resistor for N = 3 is " + str(r_N_3))
print("time it takes " + str(time_N_3))
# N = 4 
r_N_4, time_N_4 = mesh_resistence_generator(4,resistence,test_voltage)
print("overall resistence for a N*2N mesh resistor for N = 4 is " + str(r_N_4))
print("time it takes " + str(time_N_4))
# N = 5 
r_N_5, time_N_5 = mesh_resistence_generator(5,resistence,test_voltage)
print("overall resistence for a N*2N mesh resistor for N = 5 is " + str(r_N_5))
print("time it takes " + str(time_N_5))
# N = 6 
r_N_6, time_N_6 = mesh_resistence_generator(6,resistence,test_voltage)
print("overall resistence for a N*2N mesh resistor for N = 6 is " + str(r_N_6))
print("time it takes " + str(time_N_6))
# N = 7 
r_N_7, time_N_7 = mesh_resistence_generator(7,resistence,test_voltage)
print("overall resistence for a N*2N mesh resistor for N = 7 is " + str(r_N_7))
print("time it takes " + str(time_N_7))
# N = 8 
r_N_8, time_N_8 = mesh_resistence_generator(8,resistence,test_voltage)
print("overall resistence for a N*2N mesh resistor for N = 8 is " + str(r_N_8))
print("time it takes " + str(time_N_8))
# N = 9 
r_N_9, time_N_9 = mesh_resistence_generator(9,resistence,test_voltage)
print("overall resistence for a N*2N mesh resistor for N = 9 is " + str(r_N_9))
print("time it takes " + str(time_N_9))
# N = 10 
r_N_10, time_N_10 = mesh_resistence_generator(10,resistence,test_voltage)
print("overall resistence for a N*2N mesh resistor for N = 10 is " + str(r_N_10))
print("time it takes " + str(time_N_10))








#%%
'''
    b) In theory, how does the computer time taken to solve this problem increase with N, 
    for large N? Are the timings you observe for your practical implementation consistent 
    with this? Explain your observations.
'''
print("------------------------ Q2 part b -----------------------------")
# collecting the experiment data for time usage
time_list = []
resistence_list = []
for i in range(2,11):
    r, t = mesh_resistence_generator(i,resistence,test_voltage)
    time_list.append(t)
    resistence_list.append(r)
experiment_time_result_array = np.array(time_list)
print(experiment_time_result_array)


plot_horziontal_value = np.array(range(2,11))
best_fit_coef = np.polyfit(plot_horziontal_value,experiment_time_result_array,6)
best_fit_function = np.poly1d(best_fit_coef)
theory_fit_coef = best_fit_coef.copy()
theory_fit_coef[1:] = 0
theory_fit_coef[0] = np.abs(theory_fit_coef[0])
theory_function = np.poly1d(theory_fit_coef)
#plt.plot(plot_horziontal_value, experiment_time_result_array,label="experiment")
xp = np.linspace(2,11,1000)
_ = plt.plot(plot_horziontal_value,experiment_time_result_array, '*', label="experiment")
_ = plt.plot( xp, best_fit_function(xp), '--',label="best-fit")
#_ = plt.plot( xp, theory_function(xp), '-',label="theory")
#plt.plot( plot_horziontal_value, theory_time_result_array, label="theory")
plt.title('Best-fit and experiment time usage with respect to N (Q2 part b)')
plt.xlabel('N')
plt.ylabel('time used (s)')
plt.legend()
plt.show()

print("the equation of the best fitting line is: ")
print(best_fit_function)


#%%
'''
    c) Modify your program to exploit the sparse nature of the matrices to save computation time.
    What is the half-bandwidth b of your matrices? In theory, how does the computer time taken to 
    solve this problem increase now with N, for large N? Are the timings you for your practical 
    sparse implementation consistent with this? Explain your observations.
'''
#############################################################################################################
#                                    look ahead choleski decompoistion                                      #
#############################################################################################################
def choleski_look_ahead_method(A, b):
    L,y = choleski_decomp_look_ahead(A,b)
    L_T = two_d_transpose(L)
    #y = forward_elimination(L,b)
    x = find_x(L_T,y)
    return x
    
    
    
def choleski_decomp_look_ahead(A,b):
    dim_A = len(A)
    L = zeros(dim_A, dim_A)
    #print(A)
    for j in range(dim_A):
        #print("L[{}][{}] = np.sqrt(A[{}][{}])".format(j,j,j,j))
        L[j][j] = np.sqrt(A[j][j])
        b[j,0] = b[j,0]/L[j][j]
        #print("for loop 1: ")
        #print(L)
        for i in range(j+1, dim_A):
            #print("for loop 2: ")
            #print("L[{}][{}]=A[{}][{}]/L[{}][{}]".format(i,j,i,j,j,j))
            L[i][j]=A[i][j]/L[j][j]
            b[i,0]=b[i,0]-L[i][j]*b[j]
            #print(L)
            for k in range(j+1, i+1):
                #print("for loop 3: ")
                #print("A[{}][{}] = A[{}][{}] - L[{}][{}]*L[{}][{}]".format(i, k, i, k, i, j ,k, j))
                A[i][k] = A[i][k] - L[i][j]*L[k][j]
                #print(A)
    return L, b


def forward_elimination(L,b):
    dim_L = len(L)
    for j in range(dim_L):
        b[j,0] = b[j,0]/L[j][j]
        for i in range(j+1,dim_L):
            b[i,0]=b[i,0]-L[i][j]*b[j]
    return b



#############################################################################################################
#                    look ahead choleski decompoistion (half_bandwidth)                                     #
#############################################################################################################
def choleski_look_ahead_half_bandwidth_method(A, b):
    bandwidth = find_bandwidth(A)
    L,y = choleski_decomp_look_ahead_half_bandwidth(A, b, bandwidth)
    L_T = two_d_transpose(L)
    #y = forward_elimination(L,b)
    x = find_x(L_T,y)
    return x, bandwidth
    
def choleski_decomp_look_ahead_half_bandwidth(A, b, bandwidth):
    dim_A = len(A)
    L = zeros(dim_A, dim_A)
    #print(A)
    for j in range(dim_A):
        #print("L[{}][{}] = np.sqrt(A[{}][{}])".format(j,j,j,j))
        L[j][j] = np.sqrt(A[j][j])
        b[j,0] = b[j,0]/L[j][j]
        #print("for loop 1: ")
        #print(L)
        for i in range(j+1, dim_A):
            #print("for loop 2: ")
            #print("L[{}][{}]=A[{}][{}]/L[{}][{}]".format(i,j,i,j,j,j))
            if i > j+bandwidth:
                break
            else:
                L[i][j]=A[i][j]/L[j][j]
                b[i,0]=b[i,0]-L[i][j]*b[j]
                #print(L)
                for k in range(j+1, i+1):
                    #print("for loop 3: ")
                    #print("A[{}][{}] = A[{}][{}] - L[{}][{}]*L[{}][{}]".format(i, k, i, k, i, j ,k, j))
                    A[i][k] = A[i][k] - L[i][j]*L[k][j]
                    #print(A)
    return L, b

def mesh_resistence_generator_half_bandwidth(N,resistence,test_voltage):
    start_time = time.time()
    AyA, b = AyA_b_generator(N,resistence,test_voltage)
    v_n, bandwidth = choleski_look_ahead_half_bandwidth_method(AyA,b)
    test_resistence_voltage = v_n[0]
    overall_mesh_resistence = resistence*(test_resistence_voltage/(test_voltage-test_resistence_voltage))
    time_used = time.time()-start_time
    return overall_mesh_resistence, time_used, bandwidth


def find_bandwidth(A):
    dim_A = A.shape[0]
    bandwidth = 0
    for i in range(dim_A):
        if A[i][0] != 0 and (A[i+1:,0] == 0).all():
            bandwidth = i+1
    return bandwidth 



resistence = 1000
test_voltage = 10
# data collections
# N = 2 
r_N_2, time_N_2, bandwidth_N_2 = mesh_resistence_generator_half_bandwidth(2,resistence,test_voltage)
print("overall resistence for a N*2N mesh resistor for N = 2 is " + str(r_N_2))
print("time it takes " + str(time_N_2))
print("the bandwidth calculated for AyA = " + str(bandwidth_N_2))
# N = 3 
r_N_3, time_N_3, bandwidth_N_3 = mesh_resistence_generator_half_bandwidth(3,resistence,test_voltage)
print("overall resistence for a N*2N mesh resistor for N = 3 is " + str(r_N_3))
print("time it takes " + str(time_N_3))
print("the bandwidth calculated for AyA = " + str(bandwidth_N_3))
# N = 4 
r_N_4, time_N_4, bandwidth_N_4 = mesh_resistence_generator_half_bandwidth(4,resistence,test_voltage)
print("overall resistence for a N*2N mesh resistor for N = 4 is " + str(r_N_4))
print("time it takes " + str(time_N_4))
print("the bandwidth calculated for AyA = " + str(bandwidth_N_4))
# N = 5 
r_N_5, time_N_5, bandwidth_N_5 = mesh_resistence_generator_half_bandwidth(5,resistence,test_voltage)
print("overall resistence for a N*2N mesh resistor for N = 5 is " + str(r_N_5))
print("time it takes " + str(time_N_5))
print("the bandwidth calculated for AyA = " + str(bandwidth_N_5))
# N = 6 
r_N_6, time_N_6, bandwidth_N_6 = mesh_resistence_generator_half_bandwidth(6,resistence,test_voltage)
print("overall resistence for a N*2N mesh resistor for N = 6 is " + str(r_N_6))
print("time it takes " + str(time_N_6))
print("the bandwidth calculated for AyA = " + str(bandwidth_N_6))
# N = 7 
r_N_7, time_N_7, bandwidth_N_7 = mesh_resistence_generator_half_bandwidth(7,resistence,test_voltage)
print("overall resistence for a N*2N mesh resistor for N = 7 is " + str(r_N_7))
print("time it takes " + str(time_N_7))
print("the bandwidth calculated for AyA = " + str(bandwidth_N_7))
# N = 8 
r_N_8, time_N_8, bandwidth_N_8 = mesh_resistence_generator_half_bandwidth(8,resistence,test_voltage)
print("overall resistence for a N*2N mesh resistor for N = 8 is " + str(r_N_8))
print("time it takes " + str(time_N_8))
print("the bandwidth calculated for AyA = " + str(bandwidth_N_8))
# N = 9 
r_N_9, time_N_9, bandwidth_N_9 = mesh_resistence_generator_half_bandwidth(9,resistence,test_voltage)
print("overall resistence for a N*2N mesh resistor for N = 9 is " + str(r_N_9))
print("time it takes " + str(time_N_9))
print("the bandwidth calculated for AyA = " + str(bandwidth_N_9))
# N = 10 
r_N_10, time_N_10, bandwidth_N_10 = mesh_resistence_generator_half_bandwidth(10,resistence,test_voltage)
print("overall resistence for a N*2N mesh resistor for N = 10 is " + str(r_N_10))
print("time it takes " + str(time_N_10))
print("the bandwidth calculated for AyA = " + str(bandwidth_N_10))

time_list = []
resistence_list = []
bandwidths = []
for i in range(2,11):
    r, t, bandwidth= mesh_resistence_generator_half_bandwidth(i,resistence,test_voltage)
    time_list.append(t)
    resistence_list.append(r)
    bandwidths.append(bandwidth)
experiment_time_result_array_half_bandwidth = np.array(time_list)
experiment_resistence_result_array_half_bandwidth = np.array(resistence_list)
experiment_bandwidth_result_array_half_bandwidth = np.array(bandwidths)

print(experiment_time_result_array_half_bandwidth)


plot_horziontal_value = np.array(range(2,11))
best_fit_coef = np.polyfit(plot_horziontal_value,experiment_time_result_array,6)
best_fit_function = np.poly1d(best_fit_coef)
best_fit_coef_half_bandwidth = np.polyfit(plot_horziontal_value,experiment_time_result_array_half_bandwidth,2)
best_fit_half_bandwidth_function = np.poly1d(best_fit_coef_half_bandwidth)
#plt.plot(plot_horziontal_value, experiment_time_result_array,label="experiment")
xp = np.linspace(2,11,1000)
_ = plt.plot(plot_horziontal_value,experiment_time_result_array, '*', label="experiment")
_ = plt.plot(plot_horziontal_value,experiment_time_result_array_half_bandwidth, '.', label="experiment (half_bandwidth)")
_ = plt.plot( xp, best_fit_function(xp), '--',label="best-fit")
_ = plt.plot( xp, best_fit_half_bandwidth_function(xp), '-',label="best-fit (half_bandwidth)")
#_ = plt.plot( xp, theory_function(xp), '-',label="theory")
#plt.plot( plot_horziontal_value, theory_time_result_array, label="theory")
plt.title('Best-fit and experiment time usage with respect to N (Q2 part c)')
plt.xlabel('N')
plt.ylabel('time used (s)')
plt.legend()
plt.show()

print("the equation of the best fitting line is: ")
print(best_fit_function)
print("the equation of the best fitting line for half_bandwidth is: ")
print(best_fit_half_bandwidth_function)




#%%
'''
    d) Plot a graph of R versus N. Find a function R(N) that fits the curve reasonably well and 
        is asymptotically correct as N tends to infinity, as far as you can tell.
'''
print("------------------------ Q2 part d -----------------------------")
resistence_result_array = np.array(resistence_list)
best_fit_coef_r = np.polyfit(plot_horziontal_value,resistence_result_array[:,0],4)
best_fit_function_r = np.poly1d(best_fit_coef_r)
xp = np.linspace(2,11,1000)
_ = plt.plot(plot_horziontal_value,resistence_result_array, '*', label="result")
_ = plt.plot( xp, best_fit_function_r(xp), '--',label="best_fit")
plt.title('R(N) best fit curve')
plt.xlabel('N')
plt.ylabel('resistence (Ohm)')
plt.legend()
plt.show()

print("the equation of R(N) = ")
print(best_fit_function_r)


#%% testing field
A = A_generator(2)
y = y_generator(4, 1000)
J = J_generator(4)
E = E_generator(4, 10)



yA = two_d_dot_product(y,two_d_transpose(A))
AyA = two_d_dot_product(A, yA)
AyA_bandwidth = find_bandwidth(AyA)
b = two_d_dot_product(A,(J-two_d_dot_product(y,E)))
AyA2, b2 = AyA_b_generator(3,1000,10)

v = choleski_direst_method(AyA,b)
v_test_branch = v[0]
overall_resistence = 1000*((10-v_test_branch)/v_test_branch)
overall_resistence = mesh_resistence_generator(2,1000,10)



