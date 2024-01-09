SPM.mat
=======
This is the parametric SPM.mat file from the faces example:

Henson, R.N.A., Shallice, T., Gorno-Tempini, M.-L. and Dolan, R.J. (2002)
Face repetition effects in implicit and explicit memory tests as 
measured by fMRI. Cerebral Cortex, 12, 178-186.

On 2nd May 2023 the script face_rep_spm12_batch.m was downloaded from: 
https://www.fil.ion.ucl.ac.uk/spm/data/face_rep/

and used to generate the SPM.mat. 

It has the following regressors:

1. N1 - first presentation of each Non-famous face.
2. N2 - second presentation
3. N2 x lag^1 - parametric effect of lag (# interceding faces)
4. N2 x lag^2 - quadratic effect of lag
5. F1 - first presentation of each Famous face
6. F2 - ...
7. F2 x lag^1 - ...
8. F2 x lag^2 - ...

VOI files
=========
An effects of interest F-contrast was specifed and two VOIs were extracted (left and right Fus), guided by the main effect of task.