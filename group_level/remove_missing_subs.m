% remove subjects that don't have brain data, motion, and score

function [m, score, varargout] = remove_missing_subs(m, score, S, test_type, test, varargin)
    % varargin: 
    % if test_type is r or t2, provide condition, motion.
    % if test_type is t, provide condition1, condition2, motion1, motion2.
    
    %varargout:
    % if test_type is r or t2, varargout{1} is motion
    % if test_type is t, varargout{1} is motion1, varargout{2} is motion2
   

    if strcmp(test_type, "r")
        condition = varargin{1};
        motion = varargin{2};
        
        overlap_brain_motion = intersect(S.brain_data.(condition).sub_ids, S.brain_data.(condition).sub_ids_motion);
        overlap_all = intersect(overlap_brain_motion, S.outcome.(test).sub_ids);
        brain_sub_index = ismember(S.brain_data.(condition).sub_ids, overlap_all);
        m = m(brain_sub_index,:);
        brain_ids = S.brain_data.(condition).sub_ids(brain_sub_index);
        score_sub_index = ismember(S.outcome.(test).sub_ids, overlap_all);
        score = score(score_sub_index);
        score_ids = S.outcome.(test).sub_ids(score_sub_index);
        motion_sub_index = ismember(S.brain_data.(condition).sub_ids_motion, overlap_all);
        motion = motion(motion_sub_index);
        motion_ids = S.brain_data.(condition).sub_ids_motion(motion_sub_index);

        % make sure brain data, motion, and score are in the same order of
        % subjects
        % create index for reordering each data
        [~, brain_idx] = ismember(overlap_all, brain_ids);
        [~, score_idx] = ismember(overlap_all, score_ids);
        [~, motion_idx] = ismember(overlap_all, motion_ids);
        % reorder
        m = m(brain_idx,:);
        score = score(score_idx);
        motion = motion(motion_idx);
        varargout{1} = motion;
        
    elseif strcmp(test_type, "t")
        if length(S.outcome.(test).contrast) == 1
            condition = varargin{1};
            motion = varargin{2};
            
            overlap_all = intersect(S.brain_data.(condition).sub_ids, S.brain_data.(condition).sub_ids_motion);
            % overlap_all = intersect(overlap_brain_motion, S.outcome.(test).sub_ids);
            brain_sub_index = ismember(S.brain_data.(condition).sub_ids, overlap_all);
            m = m(brain_sub_index,:);
            brain_ids = S.brain_data.(condition).sub_ids(brain_sub_index);
            motion_sub_index = ismember(S.brain_data.(condition).sub_ids_motion, overlap_all);
            motion = motion(motion_sub_index);
            motion_ids = S.brain_data.(condition).sub_ids_motion(motion_sub_index);

            % make sure brain data, motion, and score are in the same order of
            % subjects
            % create index for reordering each data
            [~, brain_idx] = ismember(overlap_all, brain_ids);
            [~, motion_idx] = ismember(overlap_all, motion_ids);
            % reorder
            m = m(brain_idx,:);
            motion = motion(motion_idx);
            varargout{1} = motion;
        
        else
            
            condition1 = varargin{1};
            condition2 = varargin{2};
            motion1 = varargin{3};
            motion2 = varargin{4};

            overlap_all = intersect(S.brain_data.(condition1).sub_ids, S.brain_data.(condition2).sub_ids);
            % overlap_both_motion = intersect(S.brain_data.(condition1).sub_ids_motion, S.brain_data.(condition2).sub_ids_motion);
            %overlap_conds_motion = intersect(overlap_both_conds, overlap_both_motion);
            %overlap_all = intersect(overlap_both_conds, S.outcome.(test).sub_ids);

            % brain data (condition1)
            cond1_sub_index = ismember(S.brain_data.(condition1).sub_ids, overlap_all);
            m = m(cond1_sub_index,:);
            cond1_ids = S.brain_data.(condition1).sub_ids(cond1_sub_index);
            % other brain data (called score for now, TODO: might change)
            cond2_sub_index = ismember(S.brain_data.(condition2).sub_ids, overlap_all);
            score = score(cond2_sub_index,:);
            cond2_ids = S.brain_data.(condition2).sub_ids(cond2_sub_index);
            % motion from first condition 
            % motion1_sub_index = ismember(S.brain_data.(condition1).sub_ids_motion, overlap_all);
            motion1 = motion1(cond1_sub_index,:);
            motion1_ids = S.brain_data.(condition1).sub_ids(cond1_sub_index);
            % motion from second condition (score)
            % motion2_sub_index = ismember(S.brain_data.(condition2).sub_ids_motion, overlap_all);
            motion2 = motion2(cond2_sub_index,:);
            motion2_ids = S.brain_data.(condition2).sub_ids(cond2_sub_index);

            % make sure brain data, motion, and score are in the same order of
            % subjects
            % create index for reordering each data
            [~, cond1_idx] = ismember(overlap_all, cond1_ids);
            [~, cond2_idx] = ismember(overlap_all, cond2_ids);
            [~, motion1_idx] = ismember(overlap_all, motion1_ids);
            [~, motion2_idx] = ismember(overlap_all, motion2_ids);

            % reorder
            m = m(cond1_idx,:);
            score = score(cond2_idx,:);
            motion1 = motion1(motion1_idx);
            motion2 = motion2(motion2_idx);
            varargout{1} = motion1;
            varargout{2} = motion2;
            
           
        end
        
        
    elseif strcmp(test_type, "t2")
        condition = varargin{1};
        motion = varargin{2};
        
        % only keep subjects with brain data, score, and motion
        % only keep subjects that have brain data, motion, and score
        overlap_brain_motion = intersect(S.brain_data.(condition).sub_ids, S.brain_data.(condition).sub_ids_motion);
        overlap_all = intersect(overlap_brain_motion, S.outcome.(test).sub_ids);
        brain_sub_index = ismember(S.brain_data.(condition).sub_ids, overlap_all);
        m = m(brain_sub_index,:);
        brain_ids = S.brain_data.(condition).sub_ids(brain_sub_index);

        score_sub_index = ismember(S.outcome.(test).sub_ids, overlap_all);
        score = score(score_sub_index);
        score_ids = S.outcome.(test).sub_ids(score_sub_index);

        motion_sub_index = ismember(S.brain_data.(condition).sub_ids_motion, overlap_all);
        motion = motion(motion_sub_index);
        motion_ids = S.brain_data.(condition).sub_ids_motion(motion_sub_index);

        % make sure brain data, motion, and score are in the same order of
        % subjects
        % create index for reordering each data
        [~, brain_idx] = ismember(overlap_all, brain_ids);
        [~, score_idx] = ismember(overlap_all, score_ids);
        [~, motion_idx] = ismember(overlap_all, motion_ids);
        % reorder
        m = m(brain_idx,:);
        score = score(score_idx);
        motion = motion(motion_idx);
        varargout{1} = motion;
    end
    