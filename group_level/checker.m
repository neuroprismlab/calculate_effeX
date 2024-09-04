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
% data_path = '/work/neuroprism/effect_size/data/subject_level/';
% data_filename = 's_abcd_fc_rosenblatt.mat';
% 
% S = load([data_path, data_filename]);

% turn checker into a function that takes a subject-level dataset S
% returns a checked and cleaned dataset S

function S = checker(S)

    % Check motion format
    
    disp('Convert motion to array if needed.')
    
    conditions = fieldnames(S.brain_data);
    for i = 1:length(conditions)
        cond = conditions{i};
        if istable(S.brain_data.(cond).motion)
            S.brain_data.(cond).motion = table2array(S.brain_data.(cond).motion);
            disp(['   > converted cond "', cond,'"'])
        end
        
        % check dimnsions of brain data
        if size(S.brain_data.(cond).data,1) ~= length(S.brain_data.(cond).sub_ids)
            S.brain_data.(cond).data = S.brain_data.(cond).data';
        end
        
        % check dimensions of motion
        % want n_subs x 1
        if size(S.brain_data.(cond).motion, 1) == 1
            S.brain_data.(cond).motion = S.brain_data.(cond).motion';
        end
        
        % check dimensions of sub_ids
        % want n_subs x 1
        if size(S.brain_data.(cond).sub_ids, 1) == 1
            S.brain_data.(cond).sub_ids = S.brain_data.(cond).sub_ids';
        end
        
        % check dimensions of sub_ids_motion
        % want n_subs x 1
        if isfield(S.brain_data.(cond), 'sub_ids_motion')
            if size(S.brain_data.(cond).sub_ids_motion, 1) == 1
                S.brain_data.(cond).sub_ids_motion = S.brain_data.(cond).sub_ids_motion';
            end
        end
        
            
    end


    % Discard tests that don't work for r, t, or t2
    
    disp('Checking tests')
    
    tests = fieldnames(S.outcome);
    for i = 1:length(tests)
        test = tests{i};
        disp(['   > test "', test,'" (class ', class(S.outcome.(test).score),')'])

        % check cell arrays
        if iscell(S.outcome.(test).score)

            % remove empty and NaN cells
            emptyCells = cellfun(@isempty, S.outcome.(test).score);
            nanCells = cellfun(@(x) any(isnan(x)), S.outcome.(test).score); % TODO: check
            to_be_removed=[emptyCells,nanCells];
            % remove empty cells from score
            S.outcome.(test).score(to_be_removed) = [];
            % remove also from subject ID list
            S.outcome.(test).sub_ids(to_be_removed) = [];
            disp(['     removed ', num2str(sum(emptyCells)), ' empty cells & ',num2str(sum(emptyCells)),' NaN cells'])
            % TODO: should check whether there are nans (any case) as string

            % if cell array contains numbers, convert to numeric array
            if all(cellfun(@isnumeric, S.outcome.(test).score))
                S.outcome.(test).score = cell2mat(S.outcome.(test).score);
                disp(['     ',test, ' converted from cell to numeric array'])

                % TODO: figure out how to account for the case when there could
                % be 0s, 1s, and 2s (3 levels). Would be cell with numbers, but
                % not continuous, so more like an ANOVA and this isn't
                % accounted for currently. Have not seen in the data though.

            % if there are only two unique values, change to zeros and ones
            elseif length(unique(rmmissing(S.outcome.(test).score))) == 2
                keys = [0, 1];
                values = unique(rmmissing(S.outcome.(test).score));

                level_map = [];
                for idx = 1:length(keys)
                    key = keys(idx);
                    key_name = ['key_', num2str(key)];
                    level_map.(key_name).key = key;
                    level_map.(key_name).value = values(idx);
                    
                    S.outcome.(test).score(cellfun(@(x) strcmp(x, values{idx}), S.outcome.(test).score)) = {key}; % TODO: why assign a cell?
                end
                % TODO: instead of the above and assignment below, consider instead:
                %S.outcome.(test).level_map=table(keys, values);

                disp(['     ', test, ' has two levels. Converted to binary.'])

                % transform from cell array to matrix
                S.outcome.(test).score = cell2mat(S.outcome.(test).score);
                disp(['     ',test, ' converted from cell to numeric array'])

                % save level_map
                S.outcome.(test).level_map = level_map;

            elseif length(unique(S.outcome.(test).score)) > 2
                
                warning(['Test ',test, ' has more than 2 unique levels and is not numeric. No functionality currently exists to run tests accordingly.'])
                S.outcome = rmfield(S.outcome,test);
            end

        elseif isa(S.outcome.(test).score, 'categorical')
            % TODO: much of this can probably be combined with the above
            undefined_idx = isundefined(S.outcome.(test).score);
            S.outcome.(test).score = S.outcome.(test).score(~undefined_idx);
            S.outcome.(test).sub_ids = S.outcome.(test).sub_ids(~undefined_idx);

            disp(['     removed ', num2str(sum(undefined_idx)), ' undefined cells'])

            % change categorical to binary 0s and 1s
            if length(rmmissing(unique(S.outcome.(test).score))) == 2
                % TODO: should check whether there are nans (any case) as string

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
                disp(['     ',test, ' has two levels. Converted to binary.'])

                % save level_map
                S.outcome.(test).level_map = level_map;

            elseif length(unique(S.outcome.(test).score)) > 2
                warning(['Test ',test, ' has more than 2 unique levels and is not numeric. No functionality currently exists to run tests accordingly.'])
                S.outcome = rmfield(S.outcome,test);
            end


        else 
            % remove NaN values
            nan_idx = isnan(S.outcome.(test).score);
            S.outcome.(test).score = S.outcome.(test).score(~nan_idx);
            S.outcome.(test).sub_ids = S.outcome.(test).sub_ids(~nan_idx);

            disp(['     removed ', num2str(sum(nan_idx)), ' empty cells'])

        end  
        
        % check dimensions of score
        % want n_subs x 1
        if isfield(S.outcome, test)
            if size(S.outcome.(test).score,1) == 1
                S.outcome.(test).score = S.outcome.(test).score';
            end
        end
    end
    
end



        
        
            
            
            

