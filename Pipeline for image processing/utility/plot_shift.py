import pickle
import matplotlib.pyplot as plt
import os
import scipy.io as sio

paths = [
    'Z:/trace_fear_conditioning/train/gzc_rasgrf-ai148d-93/My_V4_Miniscope/',
    'Z:/trace_fear_conditioning/train/gzc_rasgrf-ai148d-94/My_V4_Miniscope/',
    'Z:/trace_fear_conditioning/train/gzc_rasgrf-ai148d-367/My_V4_Miniscope/',
    'Z:/trace_fear_conditioning/train/gzc_rasgrf-ai148d-369/My_V4_Miniscope/',
    'Z:/trace_fear_conditioning/train/gzc_rasgrf-ai148d-370/My_V4_Miniscope/',
    'Z:/trace_fear_conditioning/train/gzc_rasgrf-ai148d-371/My_V4_Miniscope/', 
    ]
for path in paths:
    
    name=(os.path.join(path, 'AMF_despeckle_motion_corrected.pkl'))
    with open(name, 'rb') as file:
        cen = pickle.load(file)
    # plt.figure()  # 创建新的图
    # plt.plot(cen['shifts_rig'])  # % plot rigid shifts
    # plt.legend(['x shifts', 'y shifts'])
    # plt.xlabel('frames')
    # plt.ylabel('pixels')
    
    plt.figure()  # 创建新的图
    x_coords = [coord[0] for coord in cen['shifts_rig']]
    y_coords = [coord[1] for coord in cen['shifts_rig']]
    sio.savemat(os.path.join(path, 'rig_shift.mat'),{"x_coords":x_coords,"y_coords":y_coords})
    plt.plot(y_coords)  # % plot rigid shifts
    plt.plot(x_coords)