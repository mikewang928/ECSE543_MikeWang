#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 15:18:29 2021

@author: mikewang
"""

import numpy as np

'''
============================================ Question 1 ========================================
'''

#%% define Choleski Decompsition
'''
Q1 a) Write a program to solve the matrix equation Ax=b by Choleski decomposition. 
    A is a real, symmetric, positive-definite matrix of order n.
'''
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
            #print("i = " + str(i))
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
    #print("-------------- L_square_sum: j-1 = " +str(j-1))
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
# 



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
    #print(sumed_array.shape)
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

print(print("=================== Q1 a) ======================="))
# n=3
M_L_n_3 = np.array([[3, 0, 0],[2, 1, 0], [1, 2, 3]])
M_L_n_3_T = two_d_transpose(M_L_n_3)
A_n_3 = two_d_dot_product(M_L_n_3, two_d_transpose(M_L_n_3))
result = choleski_decomp_o(A_n_3)
print(result)
#%%
'''
Q1 b) Construct some small matrices (n = 2, 3, 4, or 5) to test the program. 
    Remember that the matrices must be real, symmetric and positive-definite. 
    Explain how you chose the matrices.
'''
print("=================== Q1 b) =======================")
# n=2
M_L_n_2 = np.array([[3, 0],[2, 1]])
M_L_n_2_T = two_d_transpose(M_L_n_2)
A_n_2 = two_d_dot_product(M_L_n_2, two_d_transpose(M_L_n_2))

# n=3
M_L_n_3 = np.array([[3, 0, 0],[2, 1, 0], [1, 2, 3]])
M_L_n_3_T = two_d_transpose(M_L_n_3)
A_n_3 = two_d_dot_product(M_L_n_3, two_d_transpose(M_L_n_3))

# n=4
M_L_n_4 = np.array([[3, 0, 0, 0],[2, 1, 0, 0], [1, 2, 3, 0], [2, 3, 4, 5]])
M_L_n_4_T = two_d_transpose(M_L_n_4)
A_n_4 = two_d_dot_product(M_L_n_4, two_d_transpose(M_L_n_4))

# n=5
M_L_n_5 = np.array([[3, 0, 0, 0, 0],[2, 1, 0, 0, 0], [1, 2, 3, 0, 0], [2, 3, 4, 5, 0], [4, 5, 6, 7, 8]])
M_L_n_5_T = two_d_transpose(M_L_n_5)
A_n_5 = two_d_dot_product(M_L_n_5, two_d_transpose(M_L_n_5))



#%% Test Q1 a) program with matrices generated in Q1 b)
'''
Q1 c): Test the program you wrote in (a) with each small matrix you built in (b) 
    in the following way: invent an x, multiply it by A to get b, then give A and b 
    to your program and check that it returns x correctly.
'''
print("=================== Q1 c) =======================")
#------------------ n=2 ------------------
print("------------------ n=2 ------------------")
x_n_2 = np.random.rand(2)[None].reshape(2,1)
b_n_2 = two_d_dot_product(A_n_2, x_n_2)
x_test_n_2 = choleski_direst_method(A_n_2, b_n_2)
if abs(x_test_n_2[0][0]- x_n_2[0][0]) < 1e-5 and abs(x_test_n_2[1][0]- x_n_2[1][0]) < 1e-5 :
    print("n=2 case success!")
else: 
    print ("n=2 case fail!")

#------------------ n=3 ------------------
print("------------------ n=3 ------------------")
x_n_3 = np.random.rand(3)[None].reshape(3,1)
b_n_3 = two_d_dot_product(A_n_3, x_n_3)
x_test_n_3 = choleski_direst_method(A_n_3, b_n_3)
if abs(x_test_n_3[0][0]- x_n_3[0][0]) < 1e-5 and abs(x_test_n_3[1][0]- x_n_3[1][0]) < 1e-5 and abs(x_test_n_3[2][0]- x_n_3[2][0]) < 1e-5 :
    print("n=3 case success!")
else: 
    print ("n=3 case fail!")


#------------------ n=4 ------------------
print("------------------ n=4 ------------------")
x_n_4 = np.random.rand(4,1)[None].reshape(4,1)
x_n_4 = np.array([[1.], [2.], [3.], [4.]])
b_n_4 = two_d_dot_product(A_n_4, x_n_4)
x_test_n_4 = choleski_direst_method(A_n_4, b_n_4)
if abs(x_test_n_4[0][0]- x_n_4[0][0]) < 1e-5 and abs(x_test_n_4[1][0]- x_n_4[1][0]) < 1e-5 and abs(x_test_n_4[2][0]- x_n_4[2][0]) < 1e-5 and abs(x_test_n_4[3][0]- x_n_4[3][0]) < 1e-5:
    print("n=4 case success!")
else: 
    print ("n=4 case fail!")



#------------------ n=5 ------------------
print("------------------ n=5 ------------------")
x_n_5 = np.random.rand(5,1)[None].reshape(5,1)
b_n_5 = two_d_dot_product(A_n_5, x_n_5)
x_test_n_5 = choleski_direst_method(A_n_5, b_n_5)
if abs(x_test_n_5[0][0]- x_n_5[0][0]) < 1e-5 and abs(x_test_n_5[1][0]- x_n_5[1][0]) < 1e-5 and abs(x_test_n_5[2][0]- x_n_5[2][0]) < 1e-5 and abs(x_test_n_5[3][0]- x_n_5[3][0]) and abs(x_test_n_5[4][0]- x_n_5[4][0])< 1e-5:
    print("n=5 case success!")
else: 
    print ("n=5 case fail!")

#%%
'''
Q1 d) Write a program that reads from a file a list of network branches (Jk, Rk, Ek) 
        and a reduced incidence matrix, and finds the voltages at the nodes of the network. 
        Use the code from part (a) to solve the matrix problem. Explain how the data is organized 
        and read from the file. Test the program with a few small networks that you can check by hand.
        Compare the results for your test circuits with the analytical results you obtained by hand. 
        Cleary specify each of the test circuits used with a labeled schematic diagram.
'''
print("=================== Q1 d) =======================")
#------------------------- circuit 1 -------------------------
print("------------------- circuit 1 -------------------")
A_Q1_d_1 = np.array([[-1. , 1.]])
A_Q1_d_1_t = two_d_transpose(A_Q1_d_1)
y_Q1_d_1 = np.array([[0.1, 0.],[0., 0.1]])
J_Q1_d_1 = np.array([[0.],[0.]])
E_Q1_d_1 = np.array([[10.],[0.]])
# since A*y*transpose(A)*v_n = A*(J-y*E)
# A*y*transpose(A)
AyAt_Q1_d_1 = two_d_dot_product(two_d_dot_product(A_Q1_d_1, y_Q1_d_1), A_Q1_d_1_t)
print("A*y*transpose(A) = ")
print(AyAt_Q1_d_1)
# A*(J-y*E)
rh_Q1_d_1 = two_d_dot_product(A_Q1_d_1,(J_Q1_d_1-two_d_dot_product(y_Q1_d_1,E_Q1_d_1)))
print("A*(J-y*E) = ")
print(rh_Q1_d_1)
v_n_Q1_d_1 = choleski_direst_method(AyAt_Q1_d_1, rh_Q1_d_1)
print("node voltage = " )
print(v_n_Q1_d_1)
if v_n_Q1_d_1[0][0] - 5. < 1e-5:
    print("program from a) circuit case 1 calculation correct!")
else: 
    print("program from a) circuit case 1 calculation NOT correct!")

print("------------------- circuit 2 -------------------")
A_Q1_d_2 = np.array([[-1. , -1.]])
A_Q1_d_2_t = two_d_transpose(A_Q1_d_2)
y_Q1_d_2 = np.array([[0.1, 0.],[0., 0.1]])
J_Q1_d_2 = np.array([[-10.],[0.]])
E_Q1_d_2 = np.array([[0.],[0.]])
# since A*y*transpose(A)*v_n = A*(J-y*E)
# A*y*transpose(A)
AyAt_Q1_d_2 = two_d_dot_product(two_d_dot_product(A_Q1_d_2, y_Q1_d_2), A_Q1_d_2_t)
print("A*y*transpose(A) = ")
print(AyAt_Q1_d_2)
# A*(J-y*E)
rh_Q1_d_2 = two_d_dot_product(A_Q1_d_2,(J_Q1_d_2-two_d_dot_product(y_Q1_d_2,E_Q1_d_2)))
print("A*(J-y*E) = ")
print(rh_Q1_d_2)
v_n_Q1_d_2 = choleski_direst_method(AyAt_Q1_d_2, rh_Q1_d_2)
print("node voltage = " )
print(v_n_Q1_d_2)
if v_n_Q1_d_2[0][0] - 50. < 1e-5:
    print("program from a) circuit case 2 calculation correct!")
else: 
    print("program from a) circuit case 2 calculation NOT correct!")

print("------------------- circuit 3 -------------------")
A_Q1_d_3 = np.array([[-1. , -1.]])
A_Q1_d_3_t = two_d_transpose(A_Q1_d_3)
y_Q1_d_3 = np.array([[0.1, 0.],[0., 0.1]])
J_Q1_d_3 = np.array([[0.],[-10.]])
E_Q1_d_3 = np.array([[10.],[0.]])
# since A*y*transpose(A)*v_n = A*(J-y*E)
# A*y*transpose(A)
AyAt_Q1_d_3 = two_d_dot_product(two_d_dot_product(A_Q1_d_3, y_Q1_d_3), A_Q1_d_3_t)
print("A*y*transpose(A) = ")
print(AyAt_Q1_d_3)
# A*(J-y*E)
rh_Q1_d_3 = two_d_dot_product(A_Q1_d_3,(J_Q1_d_3-two_d_dot_product(y_Q1_d_3,E_Q1_d_3)))
print("A*(J-y*E) = ")
print(rh_Q1_d_3)
v_n_Q1_d_3 = choleski_direst_method(AyAt_Q1_d_3, rh_Q1_d_3)
print("node voltage = " )
print(v_n_Q1_d_3)
if v_n_Q1_d_3[0][0] - 55. < 1e-5:
    print("************** program from a) circuit case 3 calculation correct!")
else: 
    print("************** program from a) circuit case 3 calculation NOT correct!")


print("------------------- circuit 4 -------------------")
A_Q1_d_4 = np.array([[-1. , -1., 0., 1.],[0., 0., -1., -1.]])
A_Q1_d_4_t = two_d_transpose(A_Q1_d_4)
y_Q1_d_4 = np.array([[0.1, 0., 0., 0.],[0., 0.1, 0., 0.], [0., 0., 0.2, 0.],[0., 0., 0., 0.2]])
J_Q1_d_4 = np.array([[0.],[0.],[-10.],[0.]])
E_Q1_d_4 = np.array([[10.],[0.],[0.],[0.]])
# since A*y*transpose(A)*v_n = A*(J-y*E)
# A*y*transpose(A)
AyAt_Q1_d_4 = two_d_dot_product(two_d_dot_product(A_Q1_d_4, y_Q1_d_4), A_Q1_d_4_t)
print("A*y*transpose(A) = ")
print(AyAt_Q1_d_4)
# A*(J-y*E)
rh_Q1_d_4 = two_d_dot_product(A_Q1_d_4,(J_Q1_d_4-two_d_dot_product(y_Q1_d_4,E_Q1_d_4)))
print("A*(J-y*E) = ")
print(rh_Q1_d_4)
v_n_Q1_d_4 = choleski_direst_method(AyAt_Q1_d_4, rh_Q1_d_4)
print("node voltage = " )
print(v_n_Q1_d_4)
if v_n_Q1_d_4[0][0] - 20. < 1e-5 and v_n_Q1_d_4[1][0] - 35. < 1e-5:
    print("************** program from a) circuit case 4 calculation correct!")
else: 
    print("************** program from a) circuit case 4 calculation NOT correct!")


print("------------------- circuit 5 -------------------")
A_Q1_d_5 = np.array([[-1. , 1.,  1., 0., 0., 0.],[0., -1., 0., 1., 1., 0.],[0., 0., -1., -1., 0., 1.]])
A_Q1_d_5_t = two_d_transpose(A_Q1_d_5)
y_Q1_d_5 = np.array([[0.05, 0., 0., 0., 0., 0.],[0., 0.1, 0., 0., 0., 0.], [0., 0., 0.1, 0., 0., 0.],[0., 0., 0., 1/30, 0., 0.],[0., 0., 0., 0., 1/30, 0.],[0., 0., 0., 0., 0., 1/30]])
J_Q1_d_5 = np.array([[0.],[0.],[0.],[0.],[0.],[0.]])
E_Q1_d_5 = np.array([[10.],[0.],[0.],[0.],[0.],[0.]])
# since A*y*transpose(A)*v_n = A*(J-y*E)
# A*y*transpose(A)
AyAt_Q1_d_5 = two_d_dot_product(two_d_dot_product(A_Q1_d_5, y_Q1_d_5), A_Q1_d_5_t)
print("A*y*transpose(A) = ")
print(AyAt_Q1_d_5)
# A*(J-y*E)
rh_Q1_d_5 = two_d_dot_product(A_Q1_d_5,(J_Q1_d_5-two_d_dot_product(y_Q1_d_5,E_Q1_d_5)))
print("A*(J-y*E) = ")
print(rh_Q1_d_5)
v_n_Q1_d_5 = choleski_direst_method(AyAt_Q1_d_5, rh_Q1_d_5)
print("node voltage = " )
print(v_n_Q1_d_5)
if v_n_Q1_d_5[0][0] - 5. < 1e-5 and v_n_Q1_d_5[1][0] - 3.75 < 1e-5 and v_n_Q1_d_5[2][0] - 3.75 < 1e-5:
    print("************** program from a) circuit case 5 calculation correct!")
else: 
    print("************** program from a) circuit case 5 calculation NOT correct!")











#%%
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













#%% test field 
'''
temp = zeros(2,4)
test_list = [[[1, 1 , 1, 1, 1],[1, 1 , 1, 1, 1]],[[1, 1 , 1, 1, 1],[1, 1 , 1, 1, 1]]]
test_list_1 = [[1, 1 , 1, 1, 1],[1, 1 , 1, 1, 1]]
test_list_2 = [[2, 2 , 2, 2, 2],[2, 2 , 2, 2, 2]]
test_array = np.array(test_list)
test_array_1 = np.array(test_list_1)
test_array_2 = np.array(test_list_2)
test_array_3 = two_d_transpose(test_array_2)
temp_sum = sum_array(np.array(test_list))
temp_two_d_matrix_sum = two_d_matrix_summation(test_array_1,test_array_2)
temp_two_d_doc_product = two_d_dot_product(test_array_1,test_array_3)


# test choleski_decomp_o
temp_n_2 = choleski_decomp_o(A_n_2)
temp_n_3 = choleski_decomp_o(A_n_3)
temp_n_4 = choleski_decomp_o(A_n_4)
temp_n_5 = choleski_decomp_o(A_n_5)
'''
