%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Prepares study subject-level data and meta-data for Brain EffeX
%
% Outputs 3 structures:
%
%   - study_info
%       - dataset
%       - test
%       - map
%       - brain_mask
%
%   - brain_data
%       - <name of condition 1>
%           - sub_ids
%           - data -------------- last dim is n_sub long
%           - motion
%       - <name of condition 2>
%           ...
%
%   - outcome
%       - <name of outcome> ----- e.g., "male-female", "EMOTION-REST"
%           - sub_ids
%           - score ------------- continuous, integers (ordinal), binary (for t-test), or strings/factors (categorical)
%                                 currently allowing multiple columns for unique contrast - TODO: consider specifying vs. design mat
%           - contrast ---------- cell array, one contrast per column, for t-test between brain data conditions
%           - category ---------- cognitive, biometric, etc.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Setup - USER-DEFINED

% set paths and strings
data_dir='~/project/HCP_S1200/';
output_dir='~/project/effect_size/';

traits_filter={'PMAT24_A_CR', 'Gender', 'Age', 'ProcSpeed_AgeAdj'};
traits_categories={'cognitive', 'socioeconomic', 'biometric', 'cognitive'};

brain_data_dir=[data_dir,'REST/'];
all_traits_path=[data_dir,'unrestricted_beh_2_4_2022_20_32_4.csv'];
motion_path=[data_dir,'motion/rfMRI_REST1_allsubs_Movement_RelativeRMS_mean.txt'];
motion_subids_path=[data_dir,'motion/subids'];

% study metadata - used for filename and study info
study_info.dataset = 'hcp';
study_info.map = 'fc'; % "fc": functional connectivity, "act": activation map
study_info.test = 'r'; % r, t, or t2

% brain data info
task='REST'; % descriptive


% 2. Load & Prepare Brain Data

% get matrix subject IDs
flist=dir(brain_data_dir);
flist(1:2)=[];
n_subs=length(flist);

% make upper triangle mask
tmp=importdata([brain_data_dir,flist(1).name]);
brain_mask=triu(true(size(tmp,1)),1);

% preallocate
n_edges=sum(+(brain_mask(:)));
d2=zeros(n_edges,n_subs);

% load matrices
for i=1:n_subs

    subids_matrices(i)=str2double(flist(i).name(1:6));
    d=importdata([brain_data_dir,flist(i).name]);
    d2(:,i) = d(brain_mask);

    % progress bar
    if mod(i,50)==0; fprintf('%d/%d\n',i,n_subs); end
end

% add to master variables
study_info.mask=brain_mask;
brain_data.(task).sub_ids=subids; 
brain_data.(task).data=d2;


% 3. Load & Prepare Outcome Data

% load all traits + filter by specified subset
traits_all=readtable(all_traits_path);
traits_filtered = traits_all(:, [{'Subject'}, traits_filter]);

% convert any strings to categorical vars
for col = 2:width(traits_filtered)  % skip first 'Subject' column
    if iscell(traits_filtered{1, col})
        this_trait=traits_filtered.Properties.VariableNames{col};
        traits_filtered.(this_trait) = categorical(traits_filtered.(this_trait));
    end
end

% save separately for each trait
colnames = traits_filtered.Properties.VariableNames;
outcome = struct();
for i = 2:length(colnames)
    colname = colnames{i};
    outcome.(colname).score = traits_filtered.(colname);
    outcome.(colname).sub_ids = traits_filtered.sub_ids;
    outcome.(colname).contrast=NaN; % for contrasts using brain_data fields
    outcome.(colname).category=traits_categories{i-1}; % cognitive, biometric, etc.
end

% 4. Load & Prepare Motion Data

tmp=readmatrix(motion_path);
motion_subids=readmatrix(motion_subids_path);
motion=tmp(:,3);

% align motion to brain data (remove subs without brain data, add NaN if no motion)
% note: prob some error if we have one without the other, esp motion without brain)
motion_filtered=nan(length(subids_matrices), 1);
[common_subids, brain_idx, motion_idx] = intersect(subids_matrices, motion_subids);
motion_filtered(brain_idx) = motion(motion_idx);

% add to master variable
brain_data.(task).motion=traits_filtered(:,size(traits_filtered,2));  


% 5. Save: brain data, traits, motion, subids, brain mask, and study_info

out_path=[output_dir,strjoin({'study',study_info.dataset,study_info.map,study_info.test,outcome_label},'_'),'.mat'];
save(out_path, 'brain_data', 'outcome','study_info')



