# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 18:20:16 2021

@author: wsycx
"""

import time 
import matplotlib.pyplot as plt
import numpy as np
import math
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
        #print('i: ' + str(i))
        if count == 0:
            sum_a = i
        elif count > 0:
            sum_a = sum_a + i
        count =+ 1
        #print("sum_a: " + str(sum_a))
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
            #print("A_row: ")
            #print(A_row)
            #print("B_column")
            #print(B_colomn)
            for m in range(len(A_row)):
                mult_inid_array.append(A_row[m]*B_colomn[m])
                #print("A_row[" + str(m) + "]*B_colomn["+str(m)+"]: " + str(A_row[m]*B_colomn[m]))
                #print(np.array(mult_inid_array))
            resulted_matrix[i][j]=sum_array(np.array(mult_inid_array))
            #print(sum_array(np.array(mult_inid_array)))
            #print(resulted_matrix[i][j])
    return resulted_matrix

# compute the transposation of an 2-d array
def two_d_transpose(A):
    dim_A = A.shape
    transposed_matrix = zeros(dim_A[1],dim_A[0])
    for i in range(dim_A[0]):
        for j in range(dim_A[1]):
            transposed_matrix[j][i]=A[i][j]
    return transposed_matrix

#%%
'''
============================================ Question 3 ========================================
Question setting: 
    Figure 1 shows the cross-section of an electrostatic problem with translational symmetry: 
        a coaxial cable with a square outer conductor and a rectangular inner conductor. The
        inner conductor is held at 15 volts and the outer conductor is grounded.
================================================================================================
'''
'''
Q3 a) Write a computer program to find the potential at the nodes of a regular mesh in the 
    air between the conductors by the method of finite differences. Use a five-point difference
    formula. Exploit at least one of the planes of mirror symmetry that this problem has. Use an
    equal node-spacing, h, in the x and y directions. Solve the matrix equation by successive
    over-relaxation (SOR), with SOR parameter omega. Terminate the iteration when the magnitude
    of the residual at each free node is less than 10−5 
'''


def mesh_generator(mesh_length,mesh_inner_length,mesh_inner_height,inner_voltage,outter_voltage,h):
    mesh = zeros(mesh_length,mesh_length)
    # the top right conner of the mesh as a voltage of 15 volts
    for i in range(mesh_inner_height):
        for j in range(mesh_length-mesh_inner_length,mesh_length):
            mesh[i][j] = inner_voltage
    # set the Neuman conditions
    rateofChangeX = inner_voltage*h/(length_outter/2 - length_inner/2)
    #print(rateofChangeX)
    rateofChangeY = inner_voltage*h/(length_outter/2 - height_inner/2)
    #print(rateofChangeY)
    for x in  range(mesh_length-mesh_inner_length-1,0,-1):
        mesh[0][x] = mesh[0][x+1] - rateofChangeX
    for y in range(mesh_inner_height, mesh_length-1):
        mesh[y][mesh_length-1] = mesh[y-1][mesh_length-1] - rateofChangeY
    return mesh


def SOR_solver(mesh,w,x,y,mesh_inner_length,mesh_inner_height,h):
    i = 0 
    #print("--------------- iteration "  + str(i) + "---------------")
    #print(mesh)
    mesh_computed = mesh
    ll_x, ll_y = find_coord_low_left(x, y, h)
    while computeMaxResidual(mesh,mesh_inner_length,mesh_inner_height) > MIN_RESIDUAL:
        #print("in!")
        i += 1
        #print("--------------- iteration "  + str(i) + "---------------")
        mesh_computed = step_SOR(mesh,w,mesh_inner_length,mesh_inner_height)
        #print(mesh_copmuted)
        #print(computeMaxResidual(mesh))
    x_y_value = mesh_computed[ll_x][ll_y]
    return i, x_y_value


def computeMaxResidual(mesh,mesh_inner_length,mesh_inner_height):
    MaxResidual = 0.0
    mesh_length = len(mesh)
    for y in range (1,mesh_length - 1):
        for x in range (1,mesh_length - 1):
            if x < mesh_length-mesh_inner_length or y > mesh_inner_height-1:
                Residual =  mesh[y][x-1] + mesh[y][x+1] + mesh[y-1][x] + mesh[y+1][x] - 4 * mesh[y][x]
                Residual = math.fabs(Residual)
                if Residual > MaxResidual:
                    MaxResidual = Residual
                
    return MaxResidual


def step_SOR(mesh,w,mesh_inner_length,mesh_inner_height):
    #print("in step_SOR: ")
    mesh_length = len(mesh)
    for y in range (1,mesh_length - 1):
        for x in range (1,mesh_length - 1):
            if x < mesh_length-mesh_inner_length or y > mesh_inner_height-1:
                 #print("("+str(y)+","+str(x)+")")
                 mesh[y][x] = (1 - w) * mesh[y][x] + (w/4) * (mesh[y][x-1] + mesh[y][x+1] + mesh[y-1][x] + mesh[y+1][x])
                 #print(mesh)
    return mesh
'''
def find_coord_low_left(x, y, h):
    length_outter = 0.2
    
    if x > length_outter/2:
        x_lower_left = length_outter-x
    else: 
        x_lower_left = x
    if y < length_outter/2:
        y_lower_left = length_outter-y
    else:
        y_lower_left = y
    x_lower_left_transform = x_lower_left
    y_lower_left_transform = y_lower_left - length_outter/2 
    x_mesh = int(x_lower_left_transform/h)
    y_mesh = int(y_lower_left_transform/h)
    return x_mesh, y_mesh
'''
def find_coord_low_left(x,y,h):
    length_outter = 0.2
    if x > length_outter/2:
        x_lower_left = length_outter-x
    else: 
        x_lower_left = x    
    if y > length_outter/2:
        y_lower_left = length_outter-y
    else:
        y_lower_left = y
    x_mesh = int(x_lower_left/h)
    y_mesh = int(y_lower_left/h)
    return x_mesh, y_mesh    
        
print("++++++++++++++++++++++++++++++++++++ part a) +++++++++++++++++++++++++++++++++++++")
# test 
# question constances
h = 0.02
length_inner = 0.08
height_inner = 0.04 
length_outter = 0.2

inner_voltage = 15.0
outter_voltage = 0.0 

MIN_RESIDUAL = 1e-5

# due to symmetry we will only consider the lower left quarter of the overall function. 
mesh_length =  int(length_outter/(2*h))+1
mesh_inner_length = int(length_inner/(2*h))+1
mesh_inner_height = int(height_inner/(2*h))+1
mesh_temp = mesh_generator(mesh_length,mesh_inner_length,mesh_inner_height,inner_voltage,outter_voltage,h)
SOR_solver(mesh_temp,1.3, 0.06, 0.04,mesh_inner_length,mesh_inner_height,h)
find_coord_low_left(0.06, 0.04, h)

#%%
'''
Q3 b) With h = 0.02, explore the effect of varying omiga. For 10 values of omiga between 1.0 and
    2.0, tabulate the number of iterations taken to achieve convergence, and the corresponding 
    value of potential at the point (x ,y) = (0.06, 0.04). Plot a graph of number of iterations 
    versus omiga. 
'''
print("++++++++++++++++++++++++++++++++++++ part b) +++++++++++++++++++++++++++++++++++++")
# question constances
h = 0.02
length_inner = 0.08
height_inner = 0.04 
length_outter = 0.2

inner_voltage = 15.0
outter_voltage = 0.0 

MIN_RESIDUAL = 1e-5

# due to symmetry we will only consider the lower left quarter of the overall function. 
mesh_length =  int(length_outter/(2*h))+1
mesh_inner_length = int(length_inner/(2*h))+1
mesh_inner_height = int(height_inner/(2*h))+1


x = 0.06 
y = 0.04
#mesh_temp = mesh_generator(mesh_length,mesh_inner_length,mesh_inner_height,inner_voltage,outter_voltage)
iterations = []
x_y_values = []
omega = []
for i in range(10,20):
    mesh_temp = mesh_generator(mesh_length,mesh_inner_length,mesh_inner_height,inner_voltage,outter_voltage,h)
    iteration, x_y_value =  SOR_solver(mesh_temp,0.1*i,x, y, mesh_inner_length,mesh_inner_height, h)
    iterations.append(iteration)
    x_y_values.append(x_y_value)
    omega.append(0.1*i)


plt.plot(omega, iterations)
plt.plot(omega, iterations,'*')
plt.title("number of iterations versus omega")
plt.xlabel("omega")
plt.ylabel("number of iterations")
plt.legend()
plt.show()
print("see plot: 'number of iterations versus omega'")




#%%
'''
Q3 c) With an appropriate value of omiga, chosen from the above experiment, explore the effect 
    of decreasing h on the potential. Use values of h = 0.02, 0.01, 0.005, etc, and both tabulate 
    and plot the corresponding values of potential at (x, y) = (0.06, 0.04) versus 1/h. What do you
    think is the potential at (0.06, 0.04), to three significant figures? Also, tabulate and plot 
    the number of iterations versus 1/h. Comment on the properties of both plots.
'''
print("++++++++++++++++++++++++++++++++++++ part c) +++++++++++++++++++++++++++++++++++++")
omega = 1.3 
x, y = 0.06, 0.04
h_list = []
iterations_h = []
x_y_values_h = []
for h_mul in range(0,6):
    h = 0.02/(2**h_mul)
    print(1/h)
    length_inner = 0.08
    height_inner = 0.04 
    length_outter = 0.2
    
    inner_voltage = 15.0
    outter_voltage = 0.0 
    
    MIN_RESIDUAL = 1e-5
    
    # due to symmetry we will only consider the lower left quarter of the overall function. 
    mesh_length =  int(length_outter/(2*h))+1
    mesh_inner_length = int(length_inner/(2*h))+1
    mesh_inner_height = int(height_inner/(2*h))+1
    mesh_temp = mesh_generator(mesh_length,mesh_inner_length,mesh_inner_height,inner_voltage,outter_voltage,h)
    iteration, x_y_value =  SOR_solver(mesh_temp,omega,x, y, mesh_inner_length,mesh_inner_height, h)
    iterations_h.append(iteration)
    x_y_values_h.append(x_y_value)
    h_list.append(1/h)
    
plt.plot(h_list,x_y_values_h)    
plt.plot(h_list,x_y_values_h,"*")
plt.title("values of potential at (x, y) = (0.06, 0.04) versus 1/h")
plt.xlabel("1/h")
plt.ylabel("values of potential at (x, y) = (0.06, 0.04)")
plt.legend()
plt.show()
print("see plot: 'values of potential at (x, y) = (0.06, 0.04) versus 1/h'")


plt.plot(h_list,iterations_h)
plt.plot(h_list,iterations_h,'*')
plt.title("number of iterations versus 1/h")
plt.xlabel("1/h")
plt.ylabel("number of iterations")
plt.legend()
plt.show()
print("see plot: 'number of iterations versus 1/h'")
#%%
'''
Q3 d) Use the Jacobi method to solve this problem for the same values of h used in part (c). 
    Tabulate and plot the values of the potential at (x, y) = (0.06, 0.04) versus 1/h and the 
    number of iterations versus 1/h. Comment on the properties of both plots and compare to those 
    of SOR.
'''
print("++++++++++++++++++++++++++++++++++++ part d) +++++++++++++++++++++++++++++++++++++")
def Jacobi_solver(mesh,w,x,y,mesh_inner_length,mesh_inner_height,h):
    i = 0 
    #print("--------------- iteration "  + str(i) + "---------------")
    #print(mesh)
    mesh_computed = mesh
    ll_x, ll_y = find_coord_low_left(x, y, h)
    while computeMaxResidual(mesh,mesh_inner_length,mesh_inner_height) > MIN_RESIDUAL:
        #print("in!")
        i += 1
        #print("--------------- iteration "  + str(i) + "---------------")
        mesh_computed = step_Jacobi(mesh,w,mesh_inner_length,mesh_inner_height)
        #print(mesh_copmuted)
        #print(computeMaxResidual(mesh))
    x_y_value = mesh_computed[ll_x][ll_y]
    return i, x_y_value

def step_Jacobi(mesh,w,mesh_inner_length,mesh_inner_height):
    #print("in step_SOR: ")
    mesh_length = len(mesh)
    mesh_old = mesh
    for y in range (1,mesh_length - 1):
        for x in range (1,mesh_length - 1):
            if x < mesh_length-mesh_inner_length or y > mesh_inner_height-1:
                 #print("("+str(y)+","+str(x)+")")
                 mesh[y][x] = (1/4)*(mesh_old[y][x-1] + mesh_old[y][x+1] + mesh_old[y-1][x] + mesh_old[y+1][x])
                 #print(mesh)
    return mesh



omega = 1.3 
x, y = 0.06, 0.04
h_list_J = []
iterations_h_J = []
x_y_values_h_J = []
for h_mul in range(0,6):
    h_J = 0.02/(2**h_mul)
    print(1/h_J)
    length_inner = 0.08
    height_inner = 0.04 
    length_outter = 0.2
    
    inner_voltage = 15.0
    outter_voltage = 0.0 
    
    MIN_RESIDUAL = 1e-5
    
    # due to symmetry we will only consider the lower left quarter of the overall function. 
    mesh_length =  int(length_outter/(2*h_J))+1
    mesh_inner_length = int(length_inner/(2*h_J))+1
    mesh_inner_height = int(height_inner/(2*h_J))+1
    mesh_temp = mesh_generator(mesh_length,mesh_inner_length,mesh_inner_height,inner_voltage,outter_voltage,h_J)
    iteration_J, x_y_value_J =  Jacobi_solver(mesh_temp,omega,x, y, mesh_inner_length,mesh_inner_height, h_J)
    iterations_h_J.append(iteration_J)
    x_y_values_h_J.append(x_y_value_J)
    h_list_J.append(1/h_J)
    

plt.plot(h_list_J,x_y_values_h_J)    
plt.plot(h_list_J,x_y_values_h_J,"*")
plt.title("Jacobi values of potential at (x, y) = (0.06, 0.04) versus 1/h")
plt.xlabel("1/h")
plt.ylabel("Jacobi values of potential at (x, y) = (0.06, 0.04)")
plt.legend()
plt.show()
print("see plot: 'Jacobi values of potential at (x, y) = (0.06, 0.04) versus 1/h'")


plt.plot(h_list_J,iterations_h_J)
plt.plot(h_list_J,iterations_h_J,'*')
plt.title("Jacobi number of iterations versus 1/h")
plt.xlabel("1/h")
plt.ylabel("Jacobi number of iterations")
plt.legend()
plt.show()
print("see plot: 'Jacobi number of iterations versus 1/h'")

#%%

'''
Q3 e) Modify the program you wrote in part (a) to use the five-point difference formula derived 
    in class for non-uniform node spacing. An alternative to using equal node spacing, h, is to 
    use smaller node spacing in more “difficult” parts of the problem domain. Experiment with a 
    scheme of this kind and see how accurately you can compute the value of the potential at 
    (x, y) = (0.06, 0.04) using only as many nodes as for the uniform case h = 0.01 in part (c).
'''
print("++++++++++++++++++++++++++++++++++++ part e) +++++++++++++++++++++++++++++++++++++")

import math
#######################################################################################################
#Generates the initial mesh, taking into considering the boundary conditions
def nonuniform_mesh_generator(vertical_lines,horizontal_lines):
    cableHeight = 0.1
    cableWidth = 0.1
    coreHeight = 0.02
    coreWidth = 0.04
    corePot = 15.0    
    #Create the mesh, with Dirchlet conditions 
    vertical_lines_reversed = vertical_lines[::-1]
    mesh = [[corePot if x >= cableWidth-coreWidth-1e-5 and y >= cableWidth-coreHeight-1e-5 else 0.0 for x in vertical_lines] for y in horizontal_lines[::-1]]  
    #update the mesh to take into account the Neuman conditions
    rateofChangeX = corePot/(cableWidth - coreWidth)
    rateofChangeY = corePot/(cableHeight - coreHeight)    
    print(np.array(mesh))
    for x in range (len(vertical_lines)):
        if (vertical_lines[x] < cableWidth-coreWidth-1e-5):        
            mesh[0][x] = corePot - rateofChangeX * (cableWidth-coreWidth- vertical_lines[x])            
    for y in range (len(horizontal_lines)):
        if (horizontal_lines[y] > coreHeight):        
            mesh[y][len(vertical_lines)-1] = corePot - rateofChangeY * (horizontal_lines[y] - coreHeight)
    return np.array(mesh) 


def nonuniform_SOR_solver(mesh,w,x,y,horizontal_lines,vertical_lines):
    i = 0 
    #print("--------------- iteration "  + str(i) + "---------------")
    #print(mesh)
    mesh_computed = mesh
    ll_x, ll_y = find_coord_low_left_nonuniform(x, y, horizontal_lines,vertical_lines)
    while nonuniform_computeMaxResidual(mesh,horizontal_lines,vertical_lines) > MIN_RESIDUAL:
        #print("in!")
        i += 1
        #print("--------------- iteration "  + str(i) + "---------------")
        mesh_computed = nonuniform_step_SOR(mesh,w,vertical_lines,horizontal_lines)
        #print(mesh_copmuted)
        #print(computeMaxResidual(mesh))
    x_y_value = mesh_computed[ll_x][ll_y]
    return i, x_y_value

def nonuniform_step_SOR(mesh,w,vertical_lines,horizontal_lines):
    #print("in step_SOR: ")
    cableWidth = 0.1
    coreHeight = 0.02
    coreWidth = 0.04
    for y in range (1,len(horizontal_lines) - 1):
        for x in range (1,len(vertical_lines) - 1):
            if vertical_lines[x] < cableWidth-coreWidth-1e-5 or horizontal_lines[y] > coreHeight:
                 #print("("+str(y)+","+str(x)+")")
                 alpha1 = vertical_lines[x] - vertical_lines[x-1]
                 alpha2 = vertical_lines[x+1] - vertical_lines[x]
                 beta1 = horizontal_lines[y+1] - horizontal_lines[y]
                 beta2 = horizontal_lines[y] - horizontal_lines[y-1]
                 mesh[y][x] = (mesh[y][x-1]/(alpha1 * (alpha1 + alpha2)) + mesh[y][x+1]/(alpha2 * (alpha1 + alpha2)) + \
                             mesh[y-1][x]/(beta1 * (beta1 + beta2)) + mesh[y+1][x]/(beta2 * (beta1 + beta2))) / \
                             (1/(alpha1 * alpha2) + 1/(beta1 * beta2)) 
                 #print(mesh)
    return mesh

def nonuniform_computeMaxResidual(mesh,horizontal_lines,vertical_lines):
    cableWidth = 0.1
    coreHeight = 0.02
    coreWidth = 0.04
    maxRes = 0 
    for y in range (1,len(horizontal_lines) - 1):
        for x in range (1,len(vertical_lines) - 1):
            if vertical_lines[x] < cableWidth-coreWidth-1e-5 or horizontal_lines[y] > coreHeight:
                alpha1 = vertical_lines[x] - vertical_lines[x-1]
                alpha2 = vertical_lines[x+1] - vertical_lines[x]
                beta1 = horizontal_lines[y+1] - horizontal_lines[y]
                beta2 = horizontal_lines[y] - horizontal_lines[y-1]
                res = (mesh[y][x-1]/(alpha1 * (alpha1 + alpha2)) + mesh[y][x+1]/(alpha2 * (alpha1 + alpha2)) + mesh[y-1][x]/(beta1 * (beta1 + beta2)) + mesh[y+1][x]/(beta2 * (beta1 + beta2))) - (1/(alpha1 * alpha2) + 1/(beta1 * beta2))*mesh[y][x]
                res = math.fabs(res)
                if (res > maxRes):
                    #Updates variable with the biggest residue amongst the free point
                    maxRes = res
    return maxRes          

def find_coord_low_left_nonuniform(x, y, horizontal_lines,vertical_lines):
    length_outter = 0.2
    if x > length_outter/2:
        x_lower_left = length_outter-x
    else: 
        x_lower_left = x    
    if y > length_outter/2:
        y_lower_left = length_outter-y
    else:
        y_lower_left = y
    #print(x_lower_left)
    #print(y_lower_left)
    x_coord = np.where(np.array(vertical_lines)==x_lower_left)
    y_coord = np.where(np.array(horizontal_lines)==y_lower_left)
    return x_coord[0][0], y_coord[0][0]



horizontal_lines = [0.00, 0.020, 0.032, 0.04, 0.055, 0.065, 0.074, 0.082, 0.089, 0.096, 0.1]
vertical_lines = [0.00, 0.020, 0.032, 0.044, 0.055, 0.06, 0.074, 0.082, 0.089, 0.096, 0.1]

#horizontal_lines = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
#vertical_lines = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]

initialMesh = nonuniform_mesh_generator(horizontal_lines,vertical_lines)
mesh = nonuniform_step_SOR(initialMesh,1.1,vertical_lines,horizontal_lines)
x_index, y_index = find_coord_low_left_nonuniform(0.06, 0.04,horizontal_lines,vertical_lines )
iterations_non_uniform , x_y_value_non_uniform = nonuniform_SOR_solver(initialMesh,1.1,x,y,horizontal_lines,vertical_lines)
print("iterations_non_uniform = " + str(iterations_non_uniform))
print("potential at (" + str(x_index)+ ","+ str(y_index) + ") is " + str(x_y_value_non_uniform))






#%%
# test 
# question constances
h = 0.01
length_inner = 0.08
height_inner = 0.04 
length_outter = 0.2

inner_voltage = 15.0
outter_voltage = 0.0 

MIN_RESIDUAL = 1e-5

# due to symmetry we will only consider the lower left quarter of the overall function. 
mesh_length =  int(length_outter/(2*h))+1
mesh_inner_length = int(length_inner/(2*h))+1
mesh_inner_height = int(height_inner/(2*h))+1
mesh_temp = mesh_generator(mesh_length,mesh_inner_length,mesh_inner_height,inner_voltage,outter_voltage,h)
iteration, x_y_value = SOR_solver(mesh_temp,1.1, 0.06, 0.04,mesh_inner_length,mesh_inner_height,h)
print(iteration)
print(x_y_value)
print(find_coord_low_left(0.06, 0.04, h))
#print(find_coord_low_left_2(0.06, 0.04, h))


mesh_temp_J = mesh_generator(mesh_length,mesh_inner_length,mesh_inner_height,inner_voltage,outter_voltage,h)
iteration_J, x_y_value_J = Jacobi_solver(mesh_temp_J,1.7, 0.06, 0.04,mesh_inner_length,mesh_inner_height,h)
print(iteration_J)
print(x_y_value_J)

# Operating System List
systems = ['Windows', 'macOS', 'Linux']
print('Original List:', systems)

# List Reverse
systems.reverse()


# updated list
print('Updated List:', systems)




