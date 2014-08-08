function [t1_norm lesion_norm] = normalize_and_write(t1_path,t1_name,lesion_name)
%
% [t1_norm lesion_norm] = normalize_and_write(t1_path,t1_name,lesion_name)
% 
% function calls SPM8 code to segment T1 and write normalizes T1 and lesion
% area in MNI space
% 
% INPUT: 
% t1_path         % path to data 
% t1_name         % T1 name (.nii or .img)
% lesion_name     % lesion output name from ITKsnap (.nii or .img)
%
% OUTPUT:
% t1_norm: name or normalized T1
% lesion_norm: name or normalized lesion
%
% DH, function calls SPM 8 code to segment and normalize

% input:
t1_tosegment    =[t1_path t1_name];
lesion_towrite    =[t1_path lesion_name];

% SPM function call for segmentation:
% segment options:
options.tpm     = char({'/Users/dorahermes-miller/Documents/m-files/spm8/tpm/grey.nii',...
                    '/Users/dorahermes-miller/Documents/m-files/spm8/tpm/white.nii',...
                    '/Users/dorahermes-miller/Documents/m-files/spm8/tpm/csf.nii'});
options.ngaus   = [2 2 2 4];
options.regtype = 'mni';
options.warpreg = 1;
options.warpco  = 25;
options.biasreg = 1*10.^-4;
options.biasfwhm= 60;
options.samp    = 3;
options.msk     = '';

% do the segmentation:
results = spm_preproc(t1_tosegment,options);

disp([' segmented ' t1_tosegment])

% save parameters from segmentation:
[out(1).sn,out(1).isn]   = spm_prep2sn(results);
[pth,nam]     = spm_fileparts(t1_tosegment);

% save normalization parameters native --> MNI
out(1).snfile = fullfile(pth,[nam '_seg_sn.mat']);
VG=out.sn.VG;
VF=out.sn.VF;
Tr=out.sn.Tr;
Affine=out.sn.Affine;
flags=out.sn.flags;
save(out(1).snfile,'VG','VF','Tr','Affine','flags')
clear VG VF Tr Affine flags

% save normalization parameters MNI --> native
out(1).isnfile = fullfile(pth,[nam '_seg_inv_sn.mat']);
VG=out.isn.VG;
VF=out.isn.VF;
Tr=out.isn.Tr;
Affine=out.isn.Affine;
flags=out.isn.flags;
save(out(1).isnfile,'VG','VF','Tr','Affine','flags')
clear VG VF Tr Affine flags

% save gray and white matter segmentation
job.output.GM=[0 0 1];
job.output.WM=[0 0 1];
job.output.CSF=[0 0 0];
job.output.biascor=1;
job.output.cleanup=0;
spm_preproc_write(out(1).sn,job.output);

disp([' written segmentation output ' t1_tosegment])

% SPM function call for writing normalized images
% write normalized:
flags.preserve  = 0;
flags.bb        = [-100 -130 -70; 100 96 120];
flags.vox       = [1 1 1]; % here is the voxel size
flags.interp    = 0;
flags.wrap      = [0 0 0];
flags.prefix    = 'w';
job.subj.matname{1}=out(1).snfile;
job.subj.resample{1}=t1_tosegment;
job.roptions=flags;
spm_run_normalise_write(job);

disp([' normalized ' t1_tosegment])

job.subj.resample{1}=lesion_towrite;
spm_run_normalise_write(job);

disp([' normalized ' lesion_towrite])

t1_norm         =[t1_path 'w' t1_name];
lesion_norm     =[t1_path 'w' lesion_name];
disp(['written ' t1_norm ' ' lesion_norm])
