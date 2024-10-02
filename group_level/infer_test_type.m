% infer test type from data

function test_type = infer_test_type(S, test)

    % initialize test_type as undefined
    test_type = 'unknown';

    if length(unique(S.outcome.(test).score)) == 1
        % if score is single number -> t
   
        test_type = 't';
    
    elseif isa(S.outcome.(test).score, 'double') && length(unique(S.outcome.(test).score)) > 2
        % if score is continuous -> r
   
        test_type = 'r';
    
    elseif length(unique(S.outcome.(test).score)) == 2
        % if two unique entries in score -> t (paired) or t2
        
        % if both "conditions" of the score have the same sub ids -> t; otherwise t2
        unique_conditions = unique(S.outcome.(test).score);
        sub_ids_cond1 = S.outcome.(test).sub_ids(S.outcome.(test).score==unique_conditions(1));
        sub_ids_cond2 = S.outcome.(test).sub_ids(S.outcome.(test).score==unique_conditions(2));
        
        if all(ismember(sub_ids_cond1, sub_ids_cond2)) && all(ismember(sub_ids_cond2, sub_ids_cond1))
            test_type = 't';
        elseif isempty(intersect(sub_ids_cond1, sub_ids_cond2))
            test_type = 't2';
        else
            error('Some subjects are duplicated and some are not. Remove all duplicates to run a 2-sample t-test, or remove all non-duplicates to run a paired-sample t-test.')
        end
            
    elseif iscell(S.outcome.(test).contrast)
        % if contrast provided -> t or t2
        
        if ~isnan(S.outcome.(test).contrast{1})

            if length(S.outcome.(test).contrast)==1
                % single condition in contrast -> t

                %condition1 = S.outcome.(test).contrast{1};
                test_type='t';

            elseif length(S.outcome.(test).contrast)==2
                % if two conditions -> t (paired) or t2

               sub_ids_cond1 = S.brain_data.(S.outcome.(test).contrast{1}).sub_ids;
               sub_ids_cond2 = S.brain_data.(S.outcome.(test).contrast{2}).sub_ids;
                
                if all(ismember(sub_ids_cond1, sub_ids_cond2)) && all(ismember(sub_ids_cond2, sub_ids_cond1))
                    test_type = 't';
                elseif isempty(intersect(sub_ids_cond1, sub_ids_cond2))
                    test_type = 't2';
                else
                    error('Some subjects are duplicated and some are not. Remove all duplicates to run a 2-sample t-test, or remove all non-duplicates to run a paired-sample t-test.')
                end
            
            else 
                error('Contrast provided, but more than two conditions given. Provide one or two conditions.') 
            end
        else
            error('Contrast provided but is NaN. Rename contrast to match relevant brain condition.')
        end
    else
        error('Could not infer test type. May be categorical.')
    end
    
end


