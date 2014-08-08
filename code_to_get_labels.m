function [] = code_to_get_labels(t1_path,t1_name,lesion_name)

close all
% just a trick to get all spm functions in the path, spm_defaults does not
% add all we need
spm fmri
close all

% t1_path='../test_patient1MRI/test_folder/';
% t1_name='s507762-0401-00004-000001-01.img';
% lesion_name='test_tumor_native.nii';

% lesion volume 
lesion_info=spm_vol([t1_path lesion_name]); 
[data_l]=spm_read_vols(lesion_info);
temp=spm_read_hdr([t1_path lesion_name]);
native_voxel_volume=prod(abs( temp.dime.pixdim(2:4) ));
native_lesion_volume=native_voxel_volume*length(find(data_l(:)>0));
native_nrvox=length(find(data_l(:)>0));
clear data_l

[t1_norm lesion_norm] = normalize_and_write(t1_path,t1_name,lesion_name);

%%
talairach_label_info='./standard_data/talairach.nii';

lesion_info=spm_vol(lesion_norm); 
[data_l,xyz_lesion]=spm_read_vols(lesion_info);
xyz_lesion=xyz_lesion(:,data_l>1);
clear data_l

tal_info=spm_vol(talairach_label_info); 
[data_tal,xyz_tal]=spm_read_vols(tal_info);
data_tal_lesion=zeros(size(data_tal));
data_tal=data_tal(:);


%% get labels:

% info from: http://www.talairach.org/labels.txt
xyz_lesion_talairach = icbm_spm2tal(xyz_lesion);

% read in labels text file
fid = fopen('./standard_data/talairach_labels2.txt');
C = textscan(fid, '%n %s %s %s %s %s','delimiter', '.');      
label_nr=C{1};
fclose(fid);

%% assign gray matter labels to lesions

lesion_label_all=zeros(size(xyz_lesion,2),1);

for k=1:size(xyz_lesion_talairach,2)
    if mod(k,50)==0, disp(['location ' int2str(k) ' of ' int2str(size(xyz_lesion_talairach,2))]),end
    dist_vect=sqrt(sum([xyz_tal(1,:)-xyz_lesion_talairach(1,k); xyz_tal(2,:)-xyz_lesion_talairach(2,k); xyz_tal(3,:)-xyz_lesion_talairach(3,k)].^2));
    [~,min_dist]=min(dist_vect);
    lesion_label_all(k)=data_tal(min_dist);
end
disp('done assigning gray matter labels to lesion area')

clear xyz_lesion_talairach xyz_tal data_tal

%% assign white matter labels to lesion area 

disp('assigning white matter labels to lesion area')

if ~exist('./standard_data/MNI_JHU_tracts_prob.nii','file')
    gunzip('./standard_data/MNI_JHU_tracts_prob.nii.gz')
    disp('unzipping MNI_JHU_tracts_prob.nii.gz such that SPM can read it')
end
wm_mni_info=spm_vol('./standard_data/MNI_JHU_tracts_prob.nii'); 
[data_wm_mni,xyz_wm_mni]=spm_read_vols(wm_mni_info);

% convert chance to code:
% data_wm_score=zeros(size(data_wm_mni,1),size(data_wm_mni,2),size(data_wm_mni,3));
clear data_wmsc
for k=1:20
    data_wmsc(k).data=data_wm_mni(:,:,:,k);
    data_wmsc(k).data(data_wmsc(k).data<10)=0;
end
clear data_wm_mni

wm_label_all=zeros(size(xyz_lesion,2),1);
for k=1:size(xyz_lesion,2)
    if mod(k,50)==0, disp(['location ' int2str(k) ' of ' int2str(size(xyz_lesion,2))]),end
    dist_vect=sqrt(sum([xyz_wm_mni(1,:)-xyz_lesion(1,k); xyz_wm_mni(2,:)-xyz_lesion(2,k); xyz_wm_mni(3,:)-xyz_lesion(3,k)].^2));
    [d,min_dist]=min(dist_vect);
    wm_prob=0;
    for m=1:20
        if data_wmsc(m).data(min_dist)>wm_prob
            wm_prob=data_wmsc(m).data(min_dist);
            wm_label_all(k,1)=m;
            wm_label_all(k,2)=wm_prob; % in percentage
        end
    end
end

%% create output structure:

% lesion_out (output):
% nr of voxels
% hemisphere
% gray and white matter
% lobes involved 
% brodmann areas involved

% general info: nr of voxels and voxel size
lesion_out.nr_voxels=length(lesion_label_all);
lesion_out.voxel_size=abs(diag(lesion_info.mat(1:3,1:3)));

% hemisphere and gray/white
l_hemi=zeros(length(lesion_label_all),1);
gray_white=zeros(length(lesion_label_all),1);
for k=1:length(lesion_label_all)
    if isequal(C{2}{label_nr==lesion_label_all(k)},'Left Cerebrum')
        l_hemi(k)=1;
    elseif isequal(C{2}{label_nr==lesion_label_all(k)},'Right Cerebrum')
        l_hemi(k)=2;
    elseif isequal(C{2}{label_nr==lesion_label_all(k)},'Inter-Hemispheric')
        l_hemi(k)=3;
    end
    if isequal(C{5}{label_nr==lesion_label_all(k)},'Gray Matter')
        gray_white(k)=1;
    elseif isequal(C{5}{label_nr==lesion_label_all(k)},'White Matter')
        gray_white(k)=2;
    elseif isequal(C{5}{label_nr==lesion_label_all(k)},'Cerebro-Spinal Fluid')
        gray_white(k)=3;
    end
end
lesion_out.hemi.left=length(find(l_hemi==1))/lesion_out.nr_voxels;
lesion_out.hemi.right=length(find(l_hemi==2))/lesion_out.nr_voxels;
lesion_out.hemi.other=length(find(l_hemi==0 | l_hemi==3))/lesion_out.nr_voxels;
lesion_out.gray.gray=length(find(gray_white==1))/lesion_out.nr_voxels;
lesion_out.gray.white=length(find(gray_white==2))/lesion_out.nr_voxels;
lesion_out.gray.other=length(find(gray_white==0 | gray_white==3))/lesion_out.nr_voxels;

% lobe, brodmann area and gyrus
lesion_roi=unique(lesion_label_all);

% lobe:
lobe_label=[];
lobe_nr_v=[];
for k=1:length(lesion_roi)
    lobe_label{k}=C{3}{find(label_nr==lesion_roi(k))};
    lobe_nr_v(k)=length(find(lesion_label_all==lesion_roi(k)));
end
[~,~,a]=unique(lobe_label);
lesion_out.lobe=[];
for k=1:max(a)
    temp_label=lobe_label(a==k);
    lesion_out.lobe{k,1}=temp_label{1};
    lesion_out.lobe{k,2}=sum(lobe_nr_v(a==k));
end
clear lobe_label lobe_nr_v a k 

% gyrus:
gyrus_label=[];
gyrus_nr_v=[];
for k=1:length(lesion_roi)
    gyrus_label{k}=C{4}{find(label_nr==lesion_roi(k))};
    gyrus_nr_v(k)=length(find(lesion_label_all==lesion_roi(k)));
end
[~,~,a]=unique(gyrus_label);
lesion_out.gyrus=[];
for k=1:max(a)
    temp_label=gyrus_label(a==k);
    lesion_out.gyrus{k,1}=temp_label{1};
    lesion_out.gyrus{k,2}=sum(gyrus_nr_v(a==k));
end
clear gyrus_label gyrus_nr_v a k 

% brodmann area:
bm_label=[];
bm_nr_v=[];
for k=1:length(lesion_roi)
    bm_label{k}=C{6}{find(label_nr==lesion_roi(k))};
    bm_nr_v(k)=length(find(lesion_label_all==lesion_roi(k)));
end
[~,~,a]=unique(bm_label);
lesion_out.brodmann=[];
for k=1:max(a)
    temp_label=bm_label(a==k);
    lesion_out.brodmann{k,1}=temp_label{1};
    lesion_out.brodmann{k,2}=sum(bm_nr_v(a==k));
end

fid = fopen('./standard_data/MNI_JHU_tracts_prob.txt');
C = textscan(fid, '%n %s','delimiter', ',');      
label_nr=C{1};
fclose(fid);

% white matter tracts:
wm_label=[];
wm_nr_v=[];

lesion_roi=unique(wm_label_all(:,1));

for k=1:length(lesion_roi)
    wm_label{k}=C{2}{find(label_nr==lesion_roi(k))};
    wm_nr_v(k)=length(find(wm_label_all(:,1)==lesion_roi(k)));
end
[~,~,a]=unique(wm_label);
lesion_out.wmtracts=[];
for k=1:max(a)
    temp_label=wm_label(a==k);
    lesion_out.wmtracts{k,1}=temp_label{1};
    lesion_out.wmtracts{k,2}=sum(wm_nr_v(a==k));
end

%% writing text output:

fid = fopen([t1_path 'lesion_output1.txt'],'w');

fprintf(fid,'\n%s\n',['NATIVE']);
fprintf(fid,'%s\n',['lesion size (nr of voxels): ' int2str(native_nrvox)]);
fprintf(fid,'%s\n',['voxel volume (mm^3): ' num2str(native_voxel_volume)]);
fprintf(fid,'%s\n',['lesion size (mm^3): ' num2str(native_lesion_volume)]);

fprintf(fid,'\n%s\n',['MNI space']);
fprintf(fid,'%s\n',['nr voxels in lesion: ' int2str(lesion_out.nr_voxels)]);
fprintf(fid,'%s\n',['voxel size (mm x mm x mm): ' int2str(lesion_out.voxel_size')]);
fprintf(fid,'%s\n',['lesion size (mm^3): ' int2str(prod(lesion_out.voxel_size) * lesion_out.nr_voxels)]);

%%%% following labels from talairach atlas:
fprintf(fid,'\n%s\n',['following labels from: http://www.talairach.org/labels.txt']);

fprintf(fid,'\n%s\n',['HEMISPHERE (proportion)']);
fprintf(fid,'%s\n',['left hemipshere: ' num2str(lesion_out.hemi.left)]);
fprintf(fid,'%s\n',['right hemipshere: ' num2str(lesion_out.hemi.right)]);
fprintf(fid,'%s\n',['other: ' num2str(lesion_out.hemi.other)]);

fprintf(fid,'\n%s\n',['GRAY OR WHITE MATTER (proportion)']);
fprintf(fid,'%s\n',['gray matter: ' num2str(lesion_out.gray.gray)]);
fprintf(fid,'%s\n',['white matter: ' num2str(lesion_out.gray.white)]);
fprintf(fid,'%s\n',['other: ' num2str(lesion_out.gray.other)]);

fprintf(fid,'\n%s\n',['LOBE (voxels and proportion)']);
for k=1:length(lesion_out.lobe(:,1))
    fprintf(fid,'%s\n',[lesion_out.lobe{k,1} ': ' num2str(lesion_out.lobe{k,2}) ' voxels ' num2str(lesion_out.lobe{k,2}/lesion_out.nr_voxels)]);
end

fprintf(fid,'\n%s\n',['GYRUS (voxels and proportion)']);
for k=1:length(lesion_out.gyrus(:,1))
    fprintf(fid,'%s\n',[lesion_out.gyrus{k,1} ': ' num2str(lesion_out.gyrus{k,2}) ' voxels ' num2str(lesion_out.gyrus{k,2}/lesion_out.nr_voxels)]);
end

fprintf(fid,'\n%s\n',['BRODMANN (voxels and proportion)']);
for k=1:length(lesion_out.brodmann(:,1))
    fprintf(fid,'%s\n',[lesion_out.brodmann{k,1} ': ' num2str(lesion_out.brodmann{k,2}) ' voxels ' num2str(lesion_out.brodmann{k,2}/lesion_out.nr_voxels)]);
end

%%%% following labels from Mori atlas
fprintf(fid,'\n%s\n',['following labels from Zhang et al., 2008, Neuroimage']);

fprintf(fid,'\n%s\n',['WHITE MATTER TRACTS INVOLVED (voxels and proportion)']);
for k=1:length(lesion_out.wmtracts(:,1))
    fprintf(fid,'%s\n',[lesion_out.wmtracts{k,1} ': ' num2str(lesion_out.wmtracts{k,2}) ' voxels ' num2str(lesion_out.wmtracts{k,2}/lesion_out.nr_voxels)]);
end
fprintf(fid,'%s\n',['optic radiation and vertical occipital fasciculus are not assigned']);

fclose(fid);

save([t1_path 'lesion_output1.mat'],'lesion_out');

%% and put out a figure

lesion_info=spm_vol(lesion_norm); 
[data_l]=spm_read_vols(lesion_info);
data_l(data_l>1)=1;

t1_info=spm_vol(t1_norm); 
[data_t1]=spm_read_vols(t1_info);
data_t1=data_t1/max(data_t1(:));

% calculates Centre Of Mass for lesion
clear com
% disp('calculating COM')
[com.sort_x,com.pos_x]=sort(data_l,1);
[com.sort_y,com.pos_y]=sort(data_l,2);
[com.sort_z,com.pos_z]=sort(data_l,3);
%datasorted
com.x_find=round((sum(com.sort_x(:).*com.pos_x(:)))/(sum(com.sort_x(:))));
com.y_find=round((sum(com.sort_y(:).*com.pos_y(:)))/(sum(com.sort_y(:))));
com.z_find=round((sum(com.sort_z(:).*com.pos_z(:)))/(sum(com.sort_z(:))));

%%
figure('Position',[0 0 1000 500])
subplot(2,3,1)
t1_slice=flipud(squeeze(0.8*data_t1(com.x_find,:,:))');
l_slice=flipud(squeeze(data_l(com.x_find,:,:))');
rgbplaatje1=squeeze(cat(3,t1_slice,...
    t1_slice+l_slice,...
    t1_slice));
imshow(rgbplaatje1);

subplot(2,3,2)
t1_slice=flipud(squeeze(0.8*data_t1(:,com.y_find,:))');
l_slice=flipud(squeeze(data_l(:,com.y_find,:))');
rgbplaatje1=squeeze(cat(3,t1_slice,...
    t1_slice+l_slice,...
    t1_slice));
imshow(rgbplaatje1);
title('normalized T1 with lesion area (green)')

subplot(2,3,3)
t1_slice=flipud(squeeze(0.8*data_t1(:,:,com.z_find))');
l_slice=flipud(squeeze(data_l(:,:,com.z_find))');
rgbplaatje1=squeeze(cat(3,t1_slice,...
    t1_slice+l_slice,...
    t1_slice));
imshow(rgbplaatje1);

% and now plot on top of MNI
name_template_coreg='./standard_data/rsingle_subj_T1_1mm.nii';
template_info=spm_vol(name_template_coreg); 
[data_mni]=spm_read_vols(template_info);
data_mni=data_mni/max(data_mni(:));

subplot(2,3,4)
t1_slice=flipud(squeeze(0.8*data_mni(com.x_find,:,:))');
l_slice=flipud(squeeze(data_l(com.x_find,:,:))');
rgbplaatje1=squeeze(cat(3,t1_slice,...
    t1_slice+l_slice,...
    t1_slice));
imshow(rgbplaatje1);

subplot(2,3,5)
t1_slice=flipud(squeeze(0.8*data_mni(:,com.y_find,:))');
l_slice=flipud(squeeze(data_l(:,com.y_find,:))');
rgbplaatje1=squeeze(cat(3,t1_slice,...
    t1_slice+l_slice,...
    t1_slice));
imshow(rgbplaatje1);
title('MNI template with lesion area (green)')

subplot(2,3,6)
t1_slice=flipud(squeeze(0.8*data_mni(:,:,com.z_find))');
l_slice=flipud(squeeze(data_l(:,:,com.z_find))');
rgbplaatje1=squeeze(cat(3,t1_slice,...
    t1_slice+l_slice,...
    t1_slice));
imshow(rgbplaatje1);

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',strcat([t1_path 'figure_normalizedT1_lesion']));

