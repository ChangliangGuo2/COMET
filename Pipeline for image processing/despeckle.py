import cupy as cp
import cupyx.scipy.ndimage
import tifffile
import os
import numpy as np

paths = [
    r'E:\CM2scope\trace_fear_conditioning\recall\gzc_rasgrf-ai148d-93\My_V4_Miniscope'.replace('\\','\\'),
    r'E:\CM2scope\trace_fear_conditioning\recall\gzc_rasgrf-ai148d-94\My_V4_Miniscope'.replace('\\','\\'),
    r'E:\CM2scope\trace_fear_conditioning\recall\gzc_rasgrf-ai148d-367\My_V4_Miniscope'.replace('\\','\\'),
    r'E:\CM2scope\trace_fear_conditioning\recall\gzc_rasgrf-ai148d-369\My_V4_Miniscope'.replace('\\','\\'),
    r'E:\CM2scope\trace_fear_conditioning\recall\gzc_rasgrf-ai148d-370\My_V4_Miniscope'.replace('\\','\\'),
    r'E:\CM2scope\trace_fear_conditioning\recall\gzc_rasgrf-ai148d-371\My_V4_Miniscope'.replace('\\','\\'),
    # Add the rest of the paths here...
]

batch_size=500

for path in paths:
    image = tifffile.imread(os.path.join(path, 'AMF.tif'))
    image_cpu= np.zeros_like(image)
    num_slices = image.shape[0]
    print("Dimensions of the image:", image.shape)

    start_index = 0
    while start_index < num_slices-150:
        end_index = min(start_index + batch_size, num_slices )
        gpu_image = cp.array(image[start_index:end_index ,:, :])
        # filtered_image = cupyx.scipy.ndimage.maximum_filter(gpu_image, size=(5, 1, 1))
        filtered_image = cupyx.scipy.ndimage.uniform_filter1d(gpu_image, size=10,axis= 0)
        filtered_image = cupyx.scipy.ndimage.minimum_filter1d(filtered_image, size=100,axis= 0)

        filtered_image = filtered_image.astype(cp.int16) - cupyx.scipy.ndimage.median_filter(filtered_image.astype(cp.int16), size=(1, 3, 3))
        gpu_image=gpu_image.astype(cp.int16) - filtered_image
        del filtered_image
        gpu_image[gpu_image<0]=0
        gpu_image=gpu_image.astype(cp.uint8)
        # gpu_image=gpu_image-cupyx.scipy.ndimage.median_filter(filtered_image1, size=(3, 1, 1))
        # gpu_image=filtered_image1

        if start_index == 0 & end_index == num_slices - 1:
            image_cpu[0:end_index] = cp.asnumpy(gpu_image[0:end_index+1, :, :])
        elif start_index == 0:
            image_cpu[0:batch_size - 75] = cp.asnumpy(gpu_image[0:batch_size - 75, :, :])
        elif end_index == num_slices:
            image_cpu[start_index + 75:end_index] = cp.asnumpy(gpu_image[75:end_index  - start_index, :, :])
        else:
            image_cpu[start_index + 75:end_index - 75] = cp.asnumpy(gpu_image[75:batch_size - 75, :, :])

        start_index = end_index - 150 

    with tifffile.TiffWriter(os.path.join(path, 'AMF_despeckle.tif'), bigtiff=False, imagej=True) as tif:
        tif.write(image_cpu[:,:,322:322+1944], compression='none')
