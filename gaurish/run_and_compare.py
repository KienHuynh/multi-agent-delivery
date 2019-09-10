import numpy as np
import subprocess
import os
import yaml

import pdb


def load_data():
    odw = []
    mine = []
    for i in range(0, 100):
        file_name = 'data_exsize_100_repnum_' + str(i)
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
        
        odw_cvx = data['odw_cvx']
        odw.append(odw_cvx['makespan'])
        
        myfile = 'D:/study/phd_research/main/horsefly/gaurish/output/' + file_name + '/' + file_name + '.txt'
        with open(myfile, 'r') as stream:
            for line in stream:
                mine.append(float(line))
                break
                
    return np.asarray(odw), np.asarray(mine)


def run_experiments():
    file_name = 'data_exsize_100_repnum_'
    input_dir = 'D:/study/phd_research/main/horsefly/gaurish/pho_experiment_1_exsize_100_myformat'
    output_dir = 'D:/study/phd_research/main/horsefly/gaurish/output'
    exe_path = "D:/study/phd_research/main/horsefly/ndrones/Debug/ndrones.exe "
    for i in range(100):
        print('Running exp %d' % i)
        
        file_name_i = file_name + str(i)
        input_path = os.path.join(input_dir, file_name_i, file_name_i + '.txt')
        output_path = os.path.join(output_dir, file_name_i)
        
        if (not os.path.isfile(input_path)):
            continue
        
        if (not os.path.isdir(output_path)):
            os.mkdir(output_path)
            
        output_path = os.path.join(output_path, file_name_i + '.txt')
        
        cmd_opts = "-g 0 -i " + input_path + " -o " + output_path
        subprocess.call(r'D:/study/phd_research/main/horsefly/ndrones/Debug/ndrones.exe ' + cmd_opts)
        


if __name__ == '__main__':
    run_experiments()
    #odw, mine = load_data()
    #print(np.mean(odw))
    #print(np.mean(mine))
    