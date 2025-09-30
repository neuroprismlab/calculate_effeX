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
        
        n_subs_duplicated = length(intersect(sub_ids_cond1, sub_ids_cond2));
        n_subs_without_duplicate =  length(sub_ids_cond1) + length(sub_ids_cond2) - (2 * n_subs_duplicated); % remove duplicates from both groups

        if n_subs_without_duplicate == 0
            test_type = 't';
        elseif n_subs_duplicated == 0 
            test_type = 't2';
        else
            warning('Some subjects are duplicated and some are not. Selecting paired or 2-sample t-test based on which group size would be larger and will remove some subjects accordingly (i.e., remove all duplicates to run a 2-sample t-test, or remove all non-duplicates to run a paired-sample t-test.)')
            if n_subs_duplicated >= ((n_subs_duplicated + n_subs_without_duplicate) / 2) % compare group sizes - for t2, would assign the duplicated subjects to one of two groups
                test_type = 't';
            else
                test_type = 't2';
            end
            
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
                
                n_subs_duplicated = length(intersect(sub_ids_cond1, sub_ids_cond2));
                n_subs_without_duplicate =  length(sub_ids_cond1) + length(sub_ids_cond2) - (2 * n_subs_duplicated); % remove duplicates from both groups
        
                if n_subs_without_duplicate == 0
                    test_type = 't';
                elseif n_subs_duplicated == 0 
                    test_type = 't2';
                else
                    warning('Some subjects are duplicated and some are not. Selecting paired or 2-sample t-test based on which group size would be larger and will remove some subjects accordingly (i.e., remove all duplicates to run a 2-sample t-test, or remove all non-duplicates to run a paired-sample t-test.)')
                    if n_subs_duplicated >= ((n_subs_duplicated + n_subs_without_duplicate) / 2) % compare group sizes - for t2, would assign the duplicated subjects to one of two groups
                        test_type = 't';
                    else
                        test_type = 't2';
                    end
                        
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


