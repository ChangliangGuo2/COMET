import cv2
import os
import tifffile
import numpy as np
import pims
try:
    cv2.setNumThreads(0)
except:
    pass

def avi_to_tif(avi_file_path, output_dir):
    cap = cv2.VideoCapture(avi_file_path)

    length = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    width  = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    height = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))

    dims = [length, height, width]  
    print(dims)             

    ###############################
    # OpenCV codepath
    ###############################

    subindices = [np.r_[range(dims[0])]]
    start_frame = 0
    # Extract the data
    input_arr = np.zeros((dims[0], height, width), dtype=np.uint8)
    counter = 0
    cap.set(1, start_frame)
    current_frame = start_frame
    while counter < dims[0]:
        # Capture frame-by-frame
        if current_frame != subindices[0][counter]:
            current_frame = subindices[0][counter]
            cap.set(1, current_frame)
        ret, frame = cap.read()
        if not ret:
            break
        input_arr[counter] = frame[:, :, 0]
        counter += 1
        current_frame += 1
    # handle spatial subindices
    if len(subindices) > 1:
        input_arr = input_arr[:, subindices[1]]
    if len(subindices) > 2:
        input_arr = input_arr[:, :, subindices[2]]
    # When everything is done, release the capture
    cap.release()
    cv2.destroyAllWindows()
    with tifffile.TiffWriter(output_dir, bigtiff=False, imagej=True) as tif:
        # tif.write(input_arr[:,:,322:322+1944].astype(np.uint8), compression='none')
        tif.write(input_arr[:,:,:].astype(np.uint8), compression='none')

def rgb2gray(rgb):
    # Standard mathematical conversion
    return rgb[...,0]

def avi_to_tif2(avi_file_path, output_dir):
    pims_movie = pims.PyAVReaderTimed(avi_file_path) #   # PyAVReaderIndexed()
    length = len(pims_movie)
    height, width = pims_movie.frame_shape[0:2]    # shape is (h, w, channels)
    dims = [length, height, width]
    if length <= 0 or width <= 0 or height <= 0:
        raise OSError(f"pims fallback failed to handle AVI file {file_name}. Giving up")

    subindices = [np.r_[range(dims[0])]]
    start_frame = 0
        
    # Extract the data (note dims[0] is num frames)
    input_arr = np.zeros((dims[0], height, width), dtype=np.uint8)
    for i, ind in enumerate(subindices[0]):
        input_arr[i] = rgb2gray(pims_movie[ind]).astype(np.uint8)

    # spatial subinds
    if len(subindices) > 1:
        input_arr = input_arr[:, subindices[1]]
    if len(subindices) > 2:
        input_arr = input_arr[:, :, subindices[2]]
    with tifffile.TiffWriter(output_dir, bigtiff=False, imagej=True) as tif:
        tif.write(input_arr[:,:,322:322+1944].astype(np.uint8), compression='none')





directory=r'Z:\labyrinth_longterm\gzc_rasgrf-ai148d-94\My_V4_Miniscope'.replace('\\','\\')
output_dir=r'Z:\labyrinth_longterm\gzc_rasgrf-ai148d-94\crop'.replace('\\','\\')
i = 0
while True:
    avi_file_path = os.path.join(directory+'\\', f"{i}.avi")
    output_path= os.path.join(output_dir, f"{i}.tif")
    if not os.path.exists(avi_file_path):
        break  # If the file does not exist, we assume that we have processed all the files
    print(f"converting {avi_file_path}")   
    avi_to_tif2(avi_file_path, output_path)
    i += 1



