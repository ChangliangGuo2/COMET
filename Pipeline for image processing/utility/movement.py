# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 22:30:17 2024

@author: PKU
"""
import pickle
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import maximum_filter1d
import tifffile
import os
import cupy as cp
import cupyx.scipy.ndimage


path='H:/CM2scope_experimental_data/trace_fear_conditioning/train/gzc_rasgrf-ai148d-93/My_V4_Miniscope/'
name=(os.path.join(path, 'shifts_rig.pkl'))
with open(name, 'rb') as file:
    cen = pickle.load(file)
threshold = 2
# plt.figure()  # 创建新的图
# plt.plot(cen['shifts_rig'])  # % plot rigid shifts
# plt.legend(['x shifts', 'y shifts'])
# plt.xlabel('frames')
# plt.ylabel('pixels')

plt.figure()  # 创建新的图
x_coords = [coord[0] for coord in cen['shifts_rig']]
y_coords = [coord[1] for coord in cen['shifts_rig']]
plt.plot(y_coords)  # % plot rigid shifts

y_coords = y_coords[30:len(y_coords)-30]
y_diff = np.diff(y_coords)
y_diff_abs = np.abs(y_diff)
plt.figure()  # 创建新的图
plt.plot(y_diff_abs)  # % plot rigid shifts

y_diff_abs[y_diff_abs<0.4]=0

filtered = maximum_filter1d(y_diff_abs, size=20)
filtered[filtered!=y_diff_abs]=0
plt.figure()  # 创建新的图
plt.plot(filtered)  # % plot rigid shifts

indices=np.nonzero(filtered)[0].tolist()

new_indices = []



image = tifffile.imread(os.path.join(path, 'AMF_MC_denoised_8bit.tif'))
num_slices = image.shape[0]
indices.insert(0, 0)
indices.append(num_slices)

for i in range(len(indices) - 1):
    start, end = indices[i], indices[i+1]
    new_indices.append(start)
    
    # 计算两个元素之间可以插入的元素数量
    num_elements = (end - start) // 200 - 1
    
    # 在两个元素之间插入新的元素
    for j in range(num_elements):
        new_element = start + 200 * (j + 1)
        new_indices.append(new_element)

# 添加最后一个元素
new_indices.append(indices[-1])

print(new_indices)
 
for i in range(len(indices)-1):
    min_val = np.min(image[indices[i]+1:indices[i+1]-1], axis=0)
    image[indices[i]:indices[i+1]] = np.subtract(image[indices[i]:indices[i+1]], min_val, dtype=np.uint8)



def AMF(image,threshold = 2):
    num_slices = image.shape[0]
    print("Dimensions of the image:", image.shape)

    start_index = 0
    while start_index < num_slices:
        end_index = min(start_index + 500, num_slices - 1)
        gpu_image = cp.array(image[start_index:end_index + 1,:, :])

        for i in range(end_index - start_index + 1):
            image_gpu = gpu_image[i, :, :]
            filtered_image = cupyx.scipy.ndimage.median_filter(image_gpu, size=(3, 3))
            mask = image_gpu.astype(cp.float16) > (filtered_image.astype(cp.float16) * threshold)
            image_gpu[mask] = filtered_image[mask]
            gpu_image[i, :, :] = image_gpu

        image_cpu = cp.asnumpy(gpu_image)
        image[start_index:end_index + 1,:, :]=image_cpu
        start_index = end_index + 1
        return image
    
image=AMF(image,threshold)
with tifffile.TiffWriter(os.path.join(path, 'Mins_min.tif'), bigtiff=False, imagej=True) as tif:
    tif.write(image, compression='none')
    

