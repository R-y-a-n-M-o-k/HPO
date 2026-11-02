import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import cm

##########
# Much of the code/functions are taken from the HPO_PES/make_HPO_pes.py file,
# with modifications to fit the DMS data
##########

def get_coeffs_orders(filename):
    """Obtain the linear coefficients of PES
       obtains associated orders, torsion argument for given coefficients
       filename is the name of the txt file containing fit data
       
       Assumes up to 10th order functions in fit"""
    
    params = []
    orders = []

    with open(filename, 'r') as file:
        #skip unnecessary lines 
        line = file.readline()
        while line[0] != 'f':
            line = file.readline()
    
        #start getting orders, parameters
        line_split = line.split()
        while line_split[0][0] == 'f':
            orders_string = line_split[0][1:] #orders as a string (eg. 010)
        
            #case 1: 4-digit order string (eg. 0100 is 0 10 0)
            # 10th order only
            if len(orders_string) == 4:

                #initialise orders values
                orders_line = [0.0]*3
                index = orders_string.index('1') #where the 1 is
                orders_line[index] = 10.0

                #append the orders and corresponding parameter
                param = float(line_split[-1])
                params.append(param)
                orders.append(orders_line)

            elif len(orders_string) == 3:
                orders_line = list(orders_string)
                #convert order strings to float values
                for j in range(len(orders_line)):
                    orders_line[j] = float(orders_line[j])

                #append the orders and corresponding parameter
                param = float(line_split[-1])
                params.append(param)
                orders.append(orders_line)

            ############
            line = file.readline()
            line_split = line.split()
    
    return (np.array(params), np.array(orders))

#########
def get_non_linear_p(filename):
    """Obtain the non-linear coefficients of PES fit from 
       given txt file, filename."""
    
    non_linear_params = []
    with open(filename, 'r') as file:

        #nonlinear param lines do not start with f
        line = file.readline()
        while line[0] != 'f':
            line_split = line.split()

            #get nonlinear param values (2nd column)
            non_lin_param = float(line_split[-1])
            non_linear_params.append(non_lin_param)
            line = file.readline()
        
    return np.array(non_linear_params)

def get_dms(linear_params, nonlin_eq, orders):
    """computes the dms function for a given component direction,
       assumes fitting is taylor expansion centered at equilibrium geometry
    """
    def dms(geometries):
        """Fitted dms function (R^3 -> R) (one component)
           - Assumes geometries is a n by 6 ndarray of input values, 
           - Assumes angles are given in degrees, (function converts to radians)
           Outputs value of component of dms for those given geometries.
        """
        #initialise DMS values to constant term of fit
        dms_values = np.full(np.shape(geometries)[0], 0.0)

        #set up variable arg to be modified to calculate pes values
        #start with input geometries
        arg = geometries.copy()
        nonlin_params = nonlin_eq.copy()

        #convert deg to rad, change angle poh
        arg[:, 2] = arg[:, 2] * np.pi/180
        nonlin_params[2] = nonlin_params[2] * np.pi/180

        #compute one sum-product term per loop
        for j in range(len(orders[:,0])):
            #initialise product parts
            term = np.zeros(np.shape(arg)[0])
            term.fill(linear_params[j])
            
            products = arg.copy()
            products[:, :2] = (products[:, :2] - nonlin_params[:2])**orders[j,:2]

            #original pes code
            #compute morse_arguments fix this 
            #products[:, :2] = (1 - np.exp(-1*nonlin_params[3:5]*
            #                         (products[:, :2] - nonlin_params[:2])
            #                             ) 
            #                  )**orders[j,:2]
            #
            ##compute angle values
            products[:, 2] = (np.cos(products[:, 2]) 
                              - np.cos(nonlin_params[2]))**orders[j, 2]
            
            #compute term
            for k in range(3):
                term *= products[:,k]
                
            dms_values += term

        return dms_values  
    return dms
