function data_avg = average_within_atlas(data,atlas,return_structured)
% average map within predefined regions of atlas
%
% input:    atlas: n x m integer matrix
%           data: n x m numeric matrix (if loaded from nii, ensure this is the same space and orientation as atlas (this script will only check dimensions))
%           return_structured: return matrix structured in the same space as the original atlas (i.e., each entry contains average value for its group); otherwise return single value per group
%
% output: 1 x n_groups or n x m matrix of average within each area


% check stuff
if ~isequal(size(data) , size(atlas))
    error(['Sizes don''t match. Data size is ',sprintf('%d ',size(data)),', atlas size is ',sprintf('%d ',size(atlas))]);
end

if min(atlas(:)) < 0
    warning('Negative entries found in atlas. Unusual - you may want to check this.')
end

% filter out zero values
non_zero_indices = atlas(:) ~= 0;
filtered_atlas = atlas(non_zero_indices);
filtered_data = data(non_zero_indices);

% calculate mean
data_avg = accumarray(filtered_atlas, filtered_data, [], @mean);

% optional: structure data
if return_structured
    data_avg_structured = zeros(size(atlas)); % Initialize with zeros
    unique_atlas_regions = unique(filtered_atlas);
    % Create a mask for all regions at once
    [~, loc] = ismember(atlas, unique_atlas_regions);
    % Assign values in a vectorized manner, ignoring zeros
    non_zero_mask = atlas ~= 0;
    data_avg_structured(non_zero_mask) = data_avg(loc(non_zero_mask));
    data_avg = data_avg_structured;
end



