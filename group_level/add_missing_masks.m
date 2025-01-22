% add missing brain masks from group level datasets

data_dir = '/work/neuroprism/effect_size/data/group_level/';
sub_data_dir = '/work/neuroprism/effect_size/data/subject_level';

% for each group level study, check if results.study_info.mask exists
filenames = dir(data_dir);
filenames = filenames(~[filenames.isdir]); 
studies = {filenames.name}; 

% initialize arrays and counts
n_s_wo_masks = 0;
studies_wo_masks = {};

for i = 1:length(studies)
    study = studies{i};
    fprintf(['Checking study: ',study,'\n'])
    
    data_path = [data_dir, study];

    S = load(data_path,'results');
    
    % check if results.study_info.mask exists
    if ~isfield(S.results.study_info, 'mask')
        n_s_wo_masks = n_s_wo_masks + 1;
        studies_wo_masks{n_s_wo_masks,1} = study;
        studies_wo_masks{n_s_wo_masks,2} = [S.results.study_info.dataset, '_', S.results.study_info.map];
    end
end

datasets_to_load = unique(studies_wo_masks(:,2));
available_datasets = dir(sub_data_dir);
available_datasets = available_datasets(~[available_datasets.isdir]); 
available_datasets = {available_datasets.name};

% load each dataset that contains missing masks
for d = 1:length(datasets_to_load)
    dataset = datasets_to_load{d};
    
    if strcmp(dataset, 'hcp_ep_fc')
        dataset = 'hcpep_fc';
    end
    
    % find matching available dataset
    matches = contains(available_datasets, dataset);
    matching_idx = find(matches);
    % if more than one match, warn user and use the first match
    if size(matching_idx,2) > 1
        warning('More than one matching dataset for this name and map. Using the first one, please check.')
        matching_idx = matching_idx(1);
    end
    matching_dataset = available_datasets(matching_idx);
    
    % load dataset
    D = load([sub_data_dir, '/', matching_dataset{1}]);
    
    % find studies that need this dataset's mask
    idx = contains(studies_wo_masks(:,2), dataset);
    studies_to_do = studies_wo_masks(idx,1);
    
    % loop through studies to do, grab mask, and add
    for a = 1:length(studies_to_do)
        
        study_file = fullfile(data_dir, studies_to_do{a}); % Adjust path if needed
    
        % Open the .mat file
        loaded_data = load(study_file);
        
        results = loaded_data.results;
        
        if isfield(D.study_info, 'mask')
            results.study_info.mask = D.study_info.mask;
            save(study_file, 'results');
        else
            warning(['Mask for', studies_to_do{a}, 'study is not in study_info. Not saved.'])
        end

    end
end
    
    
    