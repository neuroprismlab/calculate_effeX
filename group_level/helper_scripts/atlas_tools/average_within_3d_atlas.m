function data_avg = average_within_3d_atlas(data,atlas,return_3d)
% average 3D map within predefined regions of atlas
% input:    atlas: must be image from load_nii
%           data: must be same space and orientation as atlas (this script will only check dimensions)
%           return_3d: return 3D image containing average value for each voxel within an area; otherwise will return single value per area
% output: 3D or 2D matrix representing average within each area (see "return_3d"
%
% personal atlases:
%   node-level: /Users/steph/Lab/Misc/Software/data/bioimagesuite/images/shenetal_neuroimage2013/shen_1mm_268_parcellation.nii.gz
%   community-level: /Users/steph/Lab/Misc/Software/data/bioimagesuite/images/shenetal_neuroimage2013/shen_1mm_268_parcellation__in_subnetworks.nii.gz

% check stuff

if ~isequal(size(data) , size(atlas))
    error(['Sizes don''t match. Data size is ',sprintf('%d ',size(data)),', atlas size is ',sprintf('%d ',size(atlas))]);
end

if min(atlas(:)) < 0
    warning('Negative entries found in atlas. Unusual - you may want to check this.')
end

% summarize

unique_atlas_regions=unique(atlas(:));
data_avg=zeros(size(unique_atlas_regions));

for i=unique_atlas_regions
    data_avg(i)=mean(data(atlas==i));
end

if return_3d
    data_avg_3d=atlas;
    for i=unique_atlas_regions
        data_avg_3d(atlas==i)=data_avg(i);
    end
    data_avg=data_avg_3d;
end

%for i = min(atlas(atlas~=0)):max(atlas(:))
%    ind=find(atlas==i);
%    data_atlas(ind)=mean(data(ind));
%end


