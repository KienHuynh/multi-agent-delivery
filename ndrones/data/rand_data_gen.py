# Script to generate random data for the package hand-off problem
# 1 source & 1 destination
# v0.1


import numpy as np
import random

import os

def generate(save_path,
        problem_type,
        desginated_point_type,
        ndrone, 
        min_speed, 
        max_speed,
        max_range):
    
    drones = np.zeros((ndrone, 3)).astype(np.int32)
    speeds = np.linspace(min_speed, max_speed, num=ndrone).astype(np.int32)
    
    # Generate s and t 
    s = np.array([0,0])
    t = np.array([0,max_range])
    
    drones[0,0:2] = np.random.rand(2)*max_range
    drones[0,2] = speeds[0]
    time_to_s = np.sqrt(np.sum((drones[0,0:2] - s)**2))/speeds[0]
    #print(drones[0,:])
    scale_factor = 1
    
    for i in range(1, ndrone):
        while True:
            drone_loc = np.random.rand(2)*max_range*scale_factor
            drone_speed = speeds[i]
            time_to_s_i = np.sqrt(np.sum((drone_loc - s)**2))/drone_speed
            if (time_to_s_i < time_to_s):
                scale_factor *= 1.05
            else:
                drones[i,0:2] = drone_loc
                drones[i,2] = drone_speed
                time_to_s = time_to_s_i
                #print(drones[i,:])
                break
    
    with open(save_path, 'w') as f:
        f.write(problem_type + '\n')
        f.write(desginated_point_type + '\n')
        f.write('1\n')
        f.write('%d %d\n' % (s[0], s[1]))
        f.write('1\n')
        f.write('%d %d\n' % (t[0], t[1]))
        f.write(str(ndrone) + '\n')
        
        for i in range(ndrone):
            f.write('%d %d %d\n' % (drones[i,0], drones[i,1], drones[i,2]))
                    
    
if __name__ == '__main__':

    out_dir = './data/'
    
    if (not os.path.isdir(out_dir)):
        os.mkdir(out_dir)
    
    rand_seed = 1311
    random.seed(rand_seed)
    np.random.seed(rand_seed)
    n_data = 1000
    
    problem_type = 'TWODIM EUCLID DISCRETE SINGLE_ID'
    desginated_point_type = 'SINGLE_POINT'
    
    min_ndrone = 5
    max_ndrone = 15
    min_speed = 1
    max_speed = 100
    max_range = 100
    
    for i in range(n_data):
        save_path = out_dir + ('%05d' % i) + '.txt'
        ndrone = np.random.randint(min_ndrone, max_ndrone+1)
        generate(save_path,
            problem_type,
            desginated_point_type,
            ndrone,
            min_speed,
            max_speed,
            max_range)
            
        print(ndrone)