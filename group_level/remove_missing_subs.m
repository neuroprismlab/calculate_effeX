% remove subjects that don't have brain_data data, motion, and data2

function [data1, data2, varargout] = remove_missing_subs(data1, data2, S, test_type, test, varargin)
    % varargin: 
    % if test_type is r or t2, provide condition, motion.
    % if test_type is t, provide condition1, condition2, motion1, motion2.
    
    %varargout:
    % if test_type is r or t2, varargout{1} is motion
    % if test_type is t, varargout{1} is motion1, varargout{2} is motion2
   

    if strcmp(test_type, "r")
    
        condition = varargin{1};
        motion = varargin{2};
        
        % filter data to only include subjects with brain, score, and motion data
        subids_intersect_brain_motion = intersect(S.brain_data.(condition).sub_ids, S.brain_data.(condition).sub_ids_motion);
        subids_intersect_all = intersect(subids_intersect_brain_motion, S.outcome.(test).sub_ids);
        
        brain_idx = ismember(S.brain_data.(condition).sub_ids, subids_intersect_all);
        brain_data = data1(brain_idx,:);
        brain_subids = S.brain_data.(condition).sub_ids(brain_idx);
        
        score_idx = ismember(S.outcome.(test).sub_ids, subids_intersect_all);
        score_data = data2(score_idx);
        score_subids = S.outcome.(test).sub_ids(score_idx);
        
        motion_idx = ismember(S.brain_data.(condition).sub_ids_motion, subids_intersect_all);
        motion = motion(motion_idx);
        motion_subids = S.brain_data.(condition).sub_ids_motion(motion_idx);

        % reorder all variables so subjects are in the same order
        [~, brain_idx] = ismember(subids_intersect_all, brain_subids);
        [~, score_idx] = ismember(subids_intersect_all, score_subids);
        [~, motion_idx] = ismember(subids_intersect_all, motion_subids);
        data1 = brain_data(brain_idx,:);
        data2 = score_data(score_idx);
        motion = motion(motion_idx);
        
        varargout{1} = motion;
       

    elseif strcmp(test_type, "t")

        if length(S.outcome.(test).contrast) == 1

            condition = varargin{1};
            motion = varargin{2};
            
            % filter data to only include subjects with brain and motion data
            subids_intersect_all = intersect(S.brain_data.(condition).sub_ids, S.brain_data.(condition).sub_ids_motion);
            
            brain_idx = ismember(S.brain_data.(condition).sub_ids, subids_intersect_all);
            brain_data = data1(brain_idx,:);
            brain_subids = S.brain_data.(condition).sub_ids(brain_idx);
            
            motion_idx = ismember(S.brain_data.(condition).sub_ids_motion, subids_intersect_all);
            motion = motion(motion_idx);
            motion_subids = S.brain_data.(condition).sub_ids_motion(motion_idx);

            % reorder all variables so subjects are in the same order
            [~, brain_idx] = ismember(subids_intersect_all, brain_subids);
            [~, motion_idx] = ismember(subids_intersect_all, motion_subids);
            data1 = brain_data(brain_idx,:);
            motion = motion(motion_idx);
            
            varargout{1} = motion;
        
        else % paired sample
            
            condition1 = varargin{1};
            condition2 = varargin{2};
            motion1 = varargin{3};
            motion2 = varargin{4};
            
            % filter data to only include subjects with brain and motion data in both conditions
            subids_intersect_all = intersect(S.brain_data.(condition1).sub_ids, S.brain_data.(condition2).sub_ids);

            braincond1_idx = ismember(S.brain_data.(condition1).sub_ids, subids_intersect_all);
            braincond1_data = data1(braincond1_idx,:);
            braincond1_subids = S.brain_data.(condition1).sub_ids(braincond1_idx);
            
            braincond2_sub_index = ismember(S.brain_data.(condition2).sub_ids, subids_intersect_all);
            braincond2_data = data2(braincond2_sub_index,:);
            braincond2_subids = S.brain_data.(condition2).sub_ids(braincond2_sub_index);
        
            motion1 = motion1(braincond1_idx,:);
            motion1_subids = S.brain_data.(condition1).sub_ids(braincond1_idx);
            
            motion2 = motion2(braincond2_sub_index,:);
            motion2_subids = S.brain_data.(condition2).sub_ids(braincond2_sub_index);

            % reorder all variables so subjects are in the same order
            [~, braincond1_idx] = ismember(subids_intersect_all, braincond1_subids);
            [~, braincond2_idx] = ismember(subids_intersect_all, braincond2_subids);
            [~, motion1_idx] = ismember(subids_intersect_all, motion1_subids);
            [~, motion2_idx] = ismember(subids_intersect_all, motion2_subids);
            data1 = braincond1_data(braincond1_idx,:);
            data2 = braincond2_data(braincond2_idx,:);
            motion1 = motion1(motion1_idx);
            motion2 = motion2(motion2_idx);
            
            varargout{1} = motion1;
            varargout{2} = motion2;
           
        end
       

    elseif strcmp(test_type, "t2")

        condition = varargin{1};
        motion = varargin{2};
        
        % filter data to only include subjects with brain, group assignments, and motion data
        subids_intersect_brain_motion = intersect(S.brain_data.(condition).sub_ids, S.brain_data.(condition).sub_ids_motion);
        subids_intersect_all = intersect(subids_intersect_brain_motion, S.outcome.(test).sub_ids);
        
        brain_idx = ismember(S.brain_data.(condition).sub_ids, subids_intersect_all);
        brain_data = data1(brain_idx,:);
        brain_subids = S.brain_data.(condition).sub_ids(brain_idx);

        group_idx = ismember(S.outcome.(test).sub_ids, subids_intersect_all);
        group_data = data2(group_idx);
        group_subids = S.outcome.(test).sub_ids(group_idx);

        motion_idx = ismember(S.brain_data.(condition).sub_ids_motion, subids_intersect_all);
        motion = motion(motion_idx);
        motion_subids = S.brain_data.(condition).sub_ids_motion(motion_idx);

        % reorder all variables so subjects are in the same order
        [~, brain_idx] = ismember(subids_intersect_all, brain_subids);
        [~, group_idx] = ismember(subids_intersect_all, group_subids);
        [~, motion_idx] = ismember(subids_intersect_all, motion_subids);
        brain_data = brain_data(brain_idx,:);
        group_data = group_data(group_idx);
        motion = motion(motion_idx);
        
        % SPECIAL - check for duplicate subjects in subids_intersect_all and assign to one of two groups
        [unique_subids, ~, unique_idx] = unique(subids_intersect_all, 'stable');
        duplicate_subids = unique_subids(histc(unique_idx, 1:numel(unique_subids)) > 1);

        if ~isempty(duplicate_subids)
            
            unique_groups = unique(group_data);
            
            % permute and assign half to group 1 and half to group 2 based on original group assignments
            for i = 1:length(duplicate_subids)
                dup_idx = find(subids_intersect_all == duplicate_subids(i));
                permuted_idx = dup_idx(randperm(length(dup_idx)));
                half = floor(length(permuted_idx) / 2);

                % assign subjects to original groups
                group1_idx = group_data(permuted_idx(1:half)) == unique_groups(1);
                group2_idx = group_data(permuted_idx(half+1:end)) == unique_groups(2);

                % filter out subjects that are not in their assigned group
                brain_data(permuted_idx(~group1_idx), :) = [];
                group_data(permuted_idx(~group1_idx)) = [];
                motion(permuted_idx(~group1_idx), :) = [];

                brain_data(permuted_idx(~group2_idx), :) = [];
                group_data(permuted_idx(~group2_idx)) = [];
                motion(permuted_idx(~group2_idx), :) = [];
            end
        end

        data1 = brain_data;
        data2 = group_data;
        
        varargout{1} = motion;

    end
 
 end

