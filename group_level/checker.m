%% Checker for subject level data

% Things to check:
% - data object contains all fields necessary
% - dimensions are correct for each field
% - data type is correct for each field

% starting out with immediate problems:
% - motion should be an array, not table
% - discard any tests that don't fit r, t, or t2 format
%   - keep if: score is continuous, score is categorical with two options,
%     or contrast is provided 

% if it's a cell, convert to array

% load data
data_path = '/work/neuroprism/effect_size/data/subject_level/';
data_filename = 's_abcd_fc_rosenblatt.mat';

S = load([data_path, data_filename]);


% check that motion is an array
conditions = fieldnames(S.brain_data);

for i = 1:length(conditions)
    cond = conditions{i};
    disp(['checking class of motion from ', cond])
    if istable(S.brain_data.(cond).motion)
        S.brain_data.(cond).motion = table2array(S.brain_data.(cond).motion);
        disp(['transformed motion from table to array for ', cond])
    end
    disp(['class of ', cond, ' is ', class(S.brain_data.(cond).motion)])
    
end

        
% discard tests that don't work for r, t, or t2
tests = fieldnames(S.outcome);

for i = 1:length(tests)
    test = tests{i};
    disp(['checking ', test])
    disp([test, ' is class ', class(S.outcome.(test).score)])
    
    % check cell arrays
    if iscell(S.outcome.(test).score)
        
        % remove empty and NaN cells
        emptyCells = cellfun(@isempty, S.outcome.(test).score);
        % remove empty cells from score
        S.outcome.(test).score(emptyCells) = [];
        % remove also from subject ID list
        S.outcome.(test).sub_ids(emptyCells) = [];
        disp(['removed ', num2str(sum(emptyCells)), ' empty cells'])

        % remove NaN cells
        nanCells = cellfun(@(x) any(isnan(x)), S.outcome.(test).score);
        % remove from score
        S.outcome.(test).score(nanCells) = [];
        % remove from sub id list
        S.outcome.(test).sub_ids(nanCells) = [];
        disp(['removed ', num2str(sum(emptyCells)), ' NaN cells'])
        
        % if cell array contains numbers, convert to numeric array
        if all(cellfun(@isnumeric, S.outcome.(test).score))
            S.outcome.(test).score = cell2mat(S.outcome.(test).score);
            disp([test, ' converted from cell to numeric array'])
            
            % TODO: figure out how to account for the case when there could
            % be 0s, 1s, and 2s (3 levels). Would be cell with numbers, but
            % not continuous, so more like an ANOVA and this isn't
            % accounted for currently. Have not seen in the data though.
        
        % if there are only two unique values, change to zeros and ones
        elseif length(unique(S.outcome.(test).score)) == 2
            keys = [0, 1];
            values = unique(rmmissing(S.outcome.(test).score));
            
            level_map = [];
            for idx = 1:length(keys)
                key = keys(idx);
                key_name = ['key_', num2str(key)];
                level_map.(key_name).key = key;
                level_map.(key_name).value = values(idx);
                S.outcome.(test).score(cellfun(@(x) strcmp(x, values{idx}), S.outcome.(test).score)) = {key};
            end
            disp([test, ' has two levels. Converted to binary.'])
            
            % transform from cell array to matrix
            S.outcome.(test).score = cell2mat(S.outcome.(test).score);
            disp([test, ' converted from cell to numeric array'])
            
            % save level_map
            S.outcome.(test).level_map = level_map;
            
        elseif length(unique(S.outcome.(test).score)) > 2
            disp([test, ' has more than 2 unique levels (and not numeric). Removing it'])
            S.outcome = rmfield(S.outcome,test);
        end
       
       
%         
%         try 
%             S.outcome.(test).score = convertCharsToStrings(S.outcome.(test).score);
%             disp([test, ' changed to string array and saved'])
%            
%         catch 
%             S.outcome = rmfield(S.outcome,test);
%             disp(['could not transform ', test, ' to a matrix. ', test, 'reremoved.'])
%         end
        
    elseif isa(S.outcome.(test).score, 'categorical')
        nan_idx = isundefined(S.outcome.(test).score);
        S.outcome.(test).score = S.outcome.(test).score(~nan_idx);
        S.outcome.(test).sub_ids = S.outcome.(test).sub_ids(~nan_idx);
        
        disp(['removed ', num2str(sum(nan_idx)), ' empty cells'])
        
        % change categorical to binary 0s and 1s
        if length(unique(S.outcome.(test).score)) == 2

            keys = [0, 1];
            values = unique(S.outcome.(test).score);

            level_map = [];
            for idx = 1:length(keys)
                key = keys(idx);
                key_name = ['key_', num2str(key)];
                level_map.(key_name).key = key;
                level_map.(key_name).value = values(idx);
            end
            S.outcome.(test).score = double(ismember(S.outcome.(test).score, values(2)));
            disp([test, ' has two levels. Converted to binary.'])
            
            % save level_map
            S.outcome.(test).level_map = level_map;
       
        elseif length(unique(S.outcome.(test).score)) > 2
            disp([test, ' has more than 2 unique levels (and not numeric). Removing it'])
            S.outcome = rmfield(S.outcome,test);
        end
        
        
    else 
        % remove NaN values
        nan_idx = isnan(S.outcome.(test).score);
        S.outcome.(test).score = S.outcome.(test).score(~nan_idx);
        S.outcome.(test).sub_ids = S.outcome.(test).sub_ids(~nan_idx);
        
        disp(['removed ', num2str(sum(nan_idx)), ' empty cells'])
        
    end  
        
end



        
        
            
            
            
