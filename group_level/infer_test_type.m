% infer test type from data

function test_type = infer_test_type(S, test)

    % initialize test_type as 'undefined'
    test_type = 'unknown';
    
    % if the score is continuous, test_type is r
    if isa(S.outcome.(test).score, 'double') && length(unique(S.outcome.(test).score)) > 2
        test_type = 'r';
    end
    
    % if a contrast is provided, test_type is either t or t2
    if iscell(S.outcome.(test).contrast)
        if ~isnan(S.outcome.(test).contrast{1})
            if length(S.outcome.(test).contrast)==1
                % if brain_data sub IDs are repeated for conditions 1 and 2, test
                % is t
                condition1 = S.outcome.(test).contrast{1};
                %condition2 = NaN;
                test_type='t';
            else
                condition1 = S.outcome.(test).contrast{1};
                condition2 = S.outcome.(test).contrast{2};
            
                % if both conditions of the contrast have the same sub ids, t
                if ~isempty(intersect(S.brain_data.(condition1).sub_ids, S.brain_data.(condition2).sub_ids))
                    test_type = 't';
                    % TODO: currently will say t if even just one sub ID is
                    % repeated across conditions. I think this is okay - steph?

                % otherwise, t2
                else
                    condition1 = S.outcome.(test).contrast{1};
                    condition2 = S.outcome.(test).contrast{2};
                    test_type = 't2';
                end
            end
        end
    end
    
    % if there are exactly two unique scores, either t or t2
    if length(unique(S.outcome.(test).score)) == 2
        % S.outcome.(test).score has one score per subject in sub_ids
        % check if there are duplicates in subject ids, this would suggest
        % that it is t2. if no duplicates, then test is t
        
        % if there are equal number of subjects in both groups, t2
        if length(unique(S.outcome.(test).sub_ids)) == length(S.outcome.(test).sub_ids)
            test_type = 't2';
          
        % if unequal number of subs in both groups, t
        else
            test_type = 't';
        end
    end
    
end
    
    
    