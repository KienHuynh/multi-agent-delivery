# Script to call the main C++ program and get results


import glob
import os
import subprocess
from shutil import copyfile

import pdb


if __name__ == '__main__':
    data_dir = 'data/'
    tmp_dir = './tmp/'
    out_dir = './result/' + data_dir
    file_list = sorted(glob.glob(data_dir + '*.txt'))
    exe_path = r'../Debug/ndrones.exe '
    
    sampling_methods = ['GRID', 'LOGGRID']
    
    if (not os.path.isdir(out_dir)):
        os.makedirs(out_dir)
        
    if (not os.path.isdir(tmp_dir)):
        os.makedirs(tmp_dir)
    
    for sm in sampling_methods:
        if (not os.path.isdir(os.path.join(out_dir, sm))):
            os.makedirs(os.path.join(out_dir, sm))
    
    #grid_params = ['10 10\n', '20 20\n']
    #loggird_params = ['10 10 1.05\n', '20 20 1.05\n']
    #circular_params = ['10 10\n', '20 20\n']
    
    grid_params = ['10 10\n']
    loggird_params = ['10 10 1.05\n']
    
    for f_path in file_list:
        
        for sm in sampling_methods:
            params = []
            if (sm == 'GRID'):
                params = grid_params
            if (sm == 'LOGGRID'):
                params = loggird_params
            if (sm == 'CIRCULAR'):
                params = circular_params
                
            for i in range(len(params)):
                pdb.set_trace()
                param = params[i]
                fname = os.path.basename(f_path)
                tmp_path = os.path.join(tmp_dir, 'input.txt')
                
                # Copy input file to tmp folder and concatenate the file with details on sampling method
                copyfile(f_path, tmp_path)
                
                with open(tmp_path, 'a') as f:
                    f.write(sm + '\n')
                    f.write(param)
                
                output_path = os.path.join(out_dir, sm, str(i) + '_' + fname)
                cmd_opts = "-g 0 -i " + tmp_path + " -o " + output_path
                print(f)
                subprocess.call(exe_path + cmd_opts)