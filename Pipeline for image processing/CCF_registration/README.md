# Registering each extracted neuron to the CCF atlas
## Parameters：
- front_pos：cranial window Anterior edge to the brgma point _(unit:mm)_
- back_pos：cranial window Posterior edge to the brgma point _(unit:mm)_
- RSPd_MO_point： Displacement adjustment between the standard brain and the FOV of the image;example:[vertical_displacement，horizotal_displacement] _(unit:pixel)_
- angle：rotate angle, clockwise is negtive _(unit:°)_
- instrument_pixel_size：optical instrument pixel size _(unit:mm/pixel)_
- instrument_size：captured image size,example:[1944,1944] _(unit:pixel)_
## How to use cortical_mapping_CM2.m：
Fill in `front_pos`, `back_pos`, `instrument_pixel_size`   `instrument_size` parameters based on the actual surgical conditions and optical equipment parameters.
Then, initially set the endpoint after line 65, iteratively run the program and modify the `RSPd_MO_point` and `angle` parameters until the drawn brain map matches the actual data. 
Once the parameters are set, remove the breakpoint and run the program entirely.
## Result：
<div align=center>
<img src="/imgs/CCF_map2.png" width="500" height="500"/><img src="/imgs/CCF_map1.png" width="500" height="500"/>
</div>
