% infer test type from data

function test_type = infer_test_type(S, test)

    % initialize test_type as 'undefined'
    test_type = 'unknown';
    
    % if the score is continuous, test_type is r
    if isa(S.outcome.(test).score, 'double') && length(unique(S.outcome.(test).score)) > 2
   
        test_type = 'r';
    
    % if there are exactly two unique scores, either t or t2
    elseif length(unique(S.outcome.(test).score)) == 2
        % S.outcome.(test).score has one score per subject in sub_ids
        
        % if duplicates in subject ids, then t2, else t
        
        % if there are equal number of subjects in both groups, t2 - TODO: check what the next line accomplishes
        if length(unique(S.outcome.(test).sub_ids)) == length(S.outcome.(test).sub_ids)
            test_type = 't2';
          
        % if unequal number of subs in both groups, t
        else
            test_type = 't';
        end
    
    % if a contrast is provided, test_type is either t or t2
    elseif iscell(S.outcome.(test).contrast)
        if ~isnan(S.outcome.(test).contrast{1})
            if length(S.outcome.(test).contrast)==1
                % if brain_data sub IDs are repeated for conditions 1 and 2, test is t - TODO: check what is meant by this
                condition1 = S.outcome.(test).contrast{1};
                %condition2 = NaN;
                test_type='t';
            elseif length(S.outcome.(test).contrast)==2
                condition1 = S.outcome.(test).contrast{1};
                condition2 = S.outcome.(test).contrast{2};
            
                % if both conditions of the contrast have some of the same sub ids, test is t
                % Note: this will assign "t" even if just one sub ID is repeated, e.g., accidentally
                if ~isempty(intersect(S.brain_data.(condition1).sub_ids, S.brain_data.(condition2).sub_ids))
                    test_type = 't';
                % otherwise, t2
                else
                    condition1 = S.outcome.(test).contrast{1};
                    condition2 = S.outcome.(test).contrast{2};
                    test_type = 't2';
                end
            else 
                error('Contrast provided - expected one or two conditions but more than two conditions given.') 
            end
        end
    else
        error('Could not infer test type. May be categorical.')
    end
    
end


