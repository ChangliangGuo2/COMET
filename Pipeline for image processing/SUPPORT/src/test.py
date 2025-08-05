import os
os.environ["CUDA_VISIBLE_DEVICES"]="0"

import numpy as np
import torch
# import skimage.io as skio
import torch.nn as nn
import tifffile
from tifffile import imread
import math
import sys



from tqdm import tqdm
from src.utils.dataset import DatasetSUPPORT_test_stitch
from model.SUPPORT import SUPPORT
from src.utils.util import get_coordinate

def normalize(image):
    """
    Normalize the image to [mean/std]=[0/1]

    Arguments:
        image: image stack (Pytorch Tensor with dimension [T, X, Y])

    Returns:
        image: normalized image stack (Pytorch Tensor with dimension [T, X, Y])
        mean_image: mean of the image stack (np.float)
        std_image: standard deviation of the image stack (np.float)
    """
    mean_image = torch.mean(image)
    std_image = torch.std(image)

    image -= mean_image
    image /= std_image

    return image, mean_image, std_image


if __name__ == '__main__':
    ########## Change it with your data ##############
    data_files = [
    r"H:/CM2scope_experimental_data/trace_fear_conditioning/train/gzc_rasgrf-ai148d-369/My_V4_Miniscope/denoise/test_data/AMF-512_MC.tif",
    # r"Y:/CM2scope/trace_fear_conditioning/train/gzc_rasgrf-ai148d-369/My_V4_Miniscope/AMF_despeckle_MC.tif",
    # r"Y:/CM2scope/trace_fear_conditioning/train/gzc_rasgrf-ai148d-370/My_V4_Miniscope/AMF_despeckle_MC.tif",
    # r"Y:/CM2scope/trace_fear_conditioning/train/gzc_rasgrf-ai148d-371/My_V4_Miniscope/AMF_despeckle_MC.tif",
    # Add the rest of the paths here...
    ]
    for data_file in data_files:
        model_file = "H:/CM2scope_experimental_data/code_github/Basic_pipeline/SUPPORT/pth/model_3.pth" 
        path, filename = os.path.split(data_file)
        base, ext = os.path.splitext(filename)
        new_filename = f"{base}_denoised_8bit{ext}"
        output_file = os.path.join(path, new_filename)
        patch_size = [61, 256, 256]
        patch_interval = [1, 256-32, 256-32]
        batch_size = 16   # lower it if memory exceeds.
        bs_size = 3    # modify if you changed bs_size when training.
        ##################################################
        torch.backends.cudnn.benchmark = True
        model = SUPPORT(in_channels=61, mid_channels=[16, 32, 64, 128, 256], depth=5,\
                blind_conv_channels=64, one_by_one_channels=[32, 16], last_layer_channels=[64, 32, 16], bs_size=bs_size).cuda()
        # 检查可用的GPU数量
        num_gpus = torch.cuda.device_count()

        if num_gpus > 1:
            # 如果有多于一个GPU，使用DataParallel
            model = nn.DataParallel(model, device_ids=range(num_gpus))
            model.module.load_state_dict(torch.load(model_file))
        else:
            # 只有一个或没有GPU，不使用DataParallel
            model.load_state_dict(torch.load(model_file))

        # demo_tif = torch.from_numpy(skio.imread(data_file).astype(np.float32)).type(torch.FloatTensor)
        demo_tif = torch.from_numpy(imread(data_file).astype(np.float16))
        demo_tif, mean_image, std_image = normalize(demo_tif[:, :, :])
        denoised_stack = np.zeros(demo_tif.shape, dtype=np.float16)
        img_size=demo_tif.size()
        whole_s, whole_h, whole_w = img_size
        img_s, img_h, img_w = patch_size
        gap_s, gap_h, gap_w = patch_interval
        num_w = math.ceil((whole_w-img_w+gap_w)/gap_w)
        num_h = math.ceil((whole_h-img_h+gap_h)/gap_h)
        num_s = math.ceil((whole_s-img_s+gap_s)/gap_s)
        coordinate = get_coordinate(img_size, patch_size, patch_interval)
        with torch.no_grad():
            model.eval()
            for  N in tqdm(range(0,num_w*num_h), desc="validate"):     
                
                indices = [coord for coord in coordinate if coord['xy'] == N]
                init_h = indices[0]['init_h']
                end_h = indices[0]['end_h']
                init_w = indices[0]['init_w']
                end_w = indices[0]['end_w']
                demo_tif_cuda=demo_tif[:,init_h:end_h,init_w:end_w].float().cuda()
                testset = DatasetSUPPORT_test_stitch(demo_tif_cuda, indices)
                test_dataloader = torch.utils.data.DataLoader(testset, batch_size=batch_size)
                    # initialize denoised stack to NaN array.
                                    
                # stitching denoised stack
                # insert the results if the stack value was NaN
                # or, half of the output volume
                for count, (noisy_image,single_coordinate) in enumerate(test_dataloader):
                    noisy_image_denoised = model(noisy_image)
                    noisy_image_denoised=noisy_image_denoised.cpu().detach().numpy().astype(np.float16)
                    T = noisy_image.size(1)
                    for bi in range(noisy_image.size(0)): 
                        stack_start_w = int(single_coordinate['stack_start_w'][bi])
                        stack_end_w = int(single_coordinate['stack_end_w'][bi])
                        patch_start_w = int(single_coordinate['patch_start_w'][bi])
                        patch_end_w = int(single_coordinate['patch_end_w'][bi])
        
                        stack_start_h = int(single_coordinate['stack_start_h'][bi])
                        stack_end_h = int(single_coordinate['stack_end_h'][bi])
                        patch_start_h = int(single_coordinate['patch_start_h'][bi])
                        patch_end_h = int(single_coordinate['patch_end_h'][bi])
        
                        stack_start_s = int(single_coordinate['init_s'][bi])
                        
        
                        denoised_stack[stack_start_s+(T//2), stack_start_h:stack_end_h, stack_start_w:stack_end_w] \
                            = noisy_image_denoised[bi,:,patch_start_h:patch_end_h, patch_start_w:patch_end_w]
        
        
                # change nan values to 0 and denormalize
        del demo_tif
        denoised_stack *= std_image.numpy()    # 原地操作，乘以标准差
        denoised_stack += mean_image.numpy()     # 原地操作，加上均值
        denoised_stack[denoised_stack<0]=0
        print(denoised_stack.shape)
        # skio.imsave(output_file, denoised_stack[(model.module.in_channels-1)//2:-(model.module.in_channels-1)//2, : , :], metadata={'axes': 'TYX'})
        with tifffile.TiffWriter(output_file, bigtiff=False, imagej=True) as tif:
            if num_gpus > 1:
                tif.write(denoised_stack[(model.module.in_channels-1)//2:-(model.module.in_channels-1)//2, : , :].astype(np.uint8), compression='none')
            else:
                tif.write(denoised_stack[(model.in_channels-1)//2:-(model.in_channels-1)//2, : , :].astype(np.uint8), compression='none')
        del denoised_stack