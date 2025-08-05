import cupy as cp
import cupyx.scipy.ndimage
import tifffile
import os

paths = [
    'H:/CM2scope_experimental_data/trace_fear_conditioning/train/gzc_rasgrf-ai148d-93/My_V4_Miniscope/',
    # 'H:/CM2scope_experimental_data/trace_fear_conditioning/train/gzc_rasgrf-ai148d-94/My_V4_Miniscope/',
    # 'H:/CM2scope_experimental_data/trace_fear_conditioning/train/gzc_rasgrf-ai148d-367/My_V4_Miniscope/',
    # 'H:/CM2scope_experimental_data/trace_fear_conditioning/train/gzc_rasgrf-ai148d-369/My_V4_Miniscope/',
    # 'H:/CM2scope_experimental_data/trace_fear_conditioning/train/gzc_rasgrf-ai148d-370/My_V4_Miniscope/',
    # 'H:/CM2scope_experimental_data/trace_fear_conditioning/train/gzc_rasgrf-ai148d-371/My_V4_Miniscope/',
    # Add the rest of the paths here...
]

threshold = 1.2

for path in paths:
    image = tifffile.imread(os.path.join(path, '0.tif'))
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
    with tifffile.TiffWriter(os.path.join(path, 'AMF1.2.tif'), bigtiff=False, imagej=True) as tif:
        tif.write(image, compression='none')
