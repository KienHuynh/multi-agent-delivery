# -*- coding: utf-8 -*-
import yaml
import io
import pdb
import os
import numpy as np

# Define data


# Write YAML file
#with io.open('data.yaml', 'w', encoding='utf8') as outfile:
#    yaml.dump(data, outfile, default_flow_style=False, allow_unicode=True)

# Read YAML file
# Read YAML file
sample_step = 80000
max_data_range = 1000000
for i in range(0, 100):
    #pdb.set_trace()
    try:
        file_path = "pho_experiment_1_exsize_100/data_exsize_100_repnum_%d/data_exsize_100_repnum_%d.yml" % (i, i)
        with open(file_path, 'rb') as stream:
            datagen = yaml.load_all(stream, Loader=yaml.UnsafeLoader)
            data = list(datagen)[0]
    except Exception as e:
        print(i)
        print(e)
        continue
    drone_info = data['drone_info']
    package = data['source']
    target = data['target'] 
    odw_straight = data['odw_straight']
    odw_cvx = data['odw_cvx']
    pdb.set_trace()
    write_path = 'pho_experiment_1_exsize_100_myformat/data_exsize_100_repnum_%d/' % (i)
    file_name = 'data_exsize_100_repnum_%d.txt' % (i)
    if (not os.path.isdir(write_path)):
        os.mkdir(write_path)
        
    # Finding min-max coordinates
    arr = np.asarray([]).reshape((0, 2))
    for i in range(len(drone_info)):
        arr = np.concatenate((arr, drone_info[i][0].reshape(1, 2)), 0)
    
    min_x = np.min(arr[:,0]) * max_data_range
    max_x = np.max(arr[:,0]) * max_data_range
    min_y = np.min(arr[:,1]) * max_data_range
    max_y = np.max(arr[:,1]) * max_data_range
    
    with open(write_path + file_name, 'w') as f:
        f.write('0 0\n')
        f.write('%d %d %f\n' % (min_x, max_x, sample_step))
        f.write('%d %d %f\n' % (min_y, max_y, sample_step))
        f.write('0 1\n')
        f.write('%d %d\n' % (int(package[0] * max_data_range), int(package[1] * max_data_range)))
        f.write('0 1\n')
        f.write('%d %d\n' % (int(target[0] * max_data_range), int(target[1] * max_data_range)))
        f.write('%d\n' % len(drone_info))
        for i in range(len(drone_info)):
            info = drone_info[i]
            x = info[0][0] * max_data_range
            y = info[0][1] * max_data_range
            v = info[1] * max_data_range
            
            f.write('%d %d %d\n' % (int(x), int(y), int(v)))

pdb.set_trace()