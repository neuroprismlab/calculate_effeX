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
    %if isdouble(S.outcome.(test).score)
        %if length(unique(S.outcome.(test).score)) ==  2
    disp(['checking ', test])
    
    if iscell(S.outcome.(test).score)
        disp([test, ' is a cell array'])
        % if it's character array and there are more than two unique, discard
        %if length(unique(S.outcome.(test).score)) > 2
        
        if length(unique(rmmissing(S.outcome.(test).score))) > 2
            disp([test, ' has more than 2 unique levels. Removing it'])
            S.outcome = rmfield(S.outcome,test);
        end
        
        % TODO: maybe add something so that if it could be a double it is
        
        try 
            S.outcome.(test).score = convertCharsToStrings(S.outcome.(test).score);
            disp([test, ' changed to string array and saved'])
           
        catch 
            S.outcome = rmfield(S.outcome,test);
            disp(['could not transform ', test, ' to a matrix. ', test, 'reremoved.'])
        end
        
    end
end

        
        
            
            
            

