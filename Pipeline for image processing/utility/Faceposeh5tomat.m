h5path='H:\CM2scope_experimental_data\running\gzc_rasgrf-ai148d-370\My_WebCam\0_FacemapPose.h5';
bodypart={'eye(back)','eye(bottom)','eye(front)','eye(top)','lowerlip','mouth','nose(bottom)','nose(r)','nose(tip)','nose(top)'};
facepose = struct();
for i=1:length(bodypart)
        facepose(i).bodypart = bodypart{i};
        facepose(i).x = h5read(h5path,['/Facemap/',bodypart{i},'/x/']);
        facepose(i).y = h5read(h5path,['/Facemap/',bodypart{i},'/y/']);;
        facepose(i).likelihood = h5read(h5path,['/Facemap/',bodypart{i},'/likelihood/']);;

end