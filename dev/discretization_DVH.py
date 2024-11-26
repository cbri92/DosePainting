#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 20:32:59 2024

@author: cbri3325
"""
import numpy as np


def find_range(index_store, x, y):
    # Loop through the dictionary and check if x falls within range for any index
    x_index = None
    y_index = None
    
    # Check x
    for index, range_ in index_store.items():
        if range_[0] <= x < range_[1]:
            x_index = index
            break
    
    # Check y
    for index, range_ in index_store.items():
        if range_[0] <= y < range_[1]:
            y_index = index
            break
    
    # Return the result as a tuple
    return (x_index, y_index)

index_store = {0: [0,10], 
               1: [10,20],
               2: [20,30], 
               3: [30,40], 
               4: [40,50], 
               5: [50,60],
               6: [60,70], 
               7: [70,80],
               8: [80,90],
               9 :[90,100],
               10: [100,110], 
               11: [110, 120]}





"""For each patient"""
patient_list = []
for patient in patients: 
    # patient is 80 by 80 array
    dosePs, volPs = patient
    # Generate a 80x1 array where each element is a tuple (row, column)
    coordinate_array = np.array([tuple(idx) for idx in zip(dosePs, volPs)])
    patient_list.append(coordinate_array)
# Reshape the result to be 80x1 (since you want a column vector)
#coordinate_array = coordinate_array.reshape(-1, 1)

Points = np.vstack(patient_list) # 800 x 2 array
print(np.shape(Points)) # check shape


array = np.empty((10, 12, 2), dtype=object)
# Initialize each position with an empty list
for idx in np.ndindex(array.shape):
    array[idx] = []
    
# 80 points for 10 patients, transform to 800 length array with each point as (x,y)

Points = np.array([[1.5, 90], [9.5, 50], [40, 80]])

for point in Points:
    x_index, y_index = find_range(index_store, point[0], point[1])
    array[x_index, y_index, 0].append(point[0])
    array[x_index, y_index, 1].append(point[1])
    
    
print(array)
    


