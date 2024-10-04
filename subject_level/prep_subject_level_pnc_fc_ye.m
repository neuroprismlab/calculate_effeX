% pseudocode documentation

%data_dir='data/subject_level/to_be_organized/s_pnc_fc_ye_corr/';

%brain = [data_dir,'output'];
%score = [data_dir,'behav'];

% original organized data that had all outcome measures
orig_data='/home/s.noble/work/effect_size/data/subject_level/to_be_organized/s_pnc_fc_ye.mat';
new_data='/home/s.noble/work/effect_size/data/subject_level/s_pnc_fc_ye.mat';

load(orig_data)

% count nans
nan_counts = [];
test_names = fieldnames(outcome);
for i = 1:numel(test_names)
    test_name = test_names{i};
    score_data = outcome.(test_name).score;
    score_label{i}=outcome.(test_name).score_label;

    if ~iscell(score_data)
        nan_count = sum(isnan(score_data));
        nan_counts = [nan_counts; nan_count];
    else
        % if score is cell array, append NaN
        nan_counts = [nan_counts; NaN];
    end
end

% find measures with at least 1k subjects (also happens to be 90% of the data)
good_measures = find(nan_counts<298);

%{
% previously used:
pnc_fc_r_rest._age.mat
% 120 "age_at_cnb" ?
pnc_fc_r_rest_lnb_acc.mat
pnc_f_r_rest_lnb_rt.mat
% ?
pnc_f_r_rest_pvrt_acc.mat
pnc_fc_r_rest_pvrt_rt.mat
% ?
pnc_f_r_rest_pfmt_acc.mat
pnc_fc_r_rest_pfmt_rt.mat
% 131-132  {'pfmt_ifac_tot'}    {'pfmt_ifac_rtc'} ? 
pnc_fc_r_rest_pwmt_acc.mat
pnc_fc_r_rest_pwmt_rt.mat
% 151-152 {'pwmt_kiwrd_tot'}    {'pwmt_kiwrd_rtc'}
pnc_fc_r_rest_pmat_rc.mat
% not sure whether this is reaction correct or ratio correct, but we could use {'pmat_pc'} 
pnc_f_r_rest_wrat.mat
% 200 {'wrat_cr_std'} ?
pnc_fc_t2_malerest_femalerest.mat
% 
%}

% for now, just take first unique measures for most "scales" unless reasons otherwise
% added 2 to get "sex" (race is 1 for future ref)
select_measures = [2;good_measures([3,7,11,13,18,67,71,81,95,118, 120, 121, 131,132, 133, 151,152, 153, 174, 175,179,186,189, 200])];
% full list of first two:
%select_measures = good_measures(3,4,5,6,7,8,11,12,13,14,15,16,18,19,25,26,67,68,69,70,71,72,74,75,81,82,84,85,92,93,94,95,96,116,117,118,119,120,121,122, 131,132, 133,134,135,136,135,136, 151,152, 153,154,155,171, 174, 175,176,179,180,186,187,189,190, 200);
%score_label(select_measures)

categories = {'demographic','psychiatric','psychiatric','psychiatric','psychiatric','psychiatric','psychiatric','psychiatric','psychiatric','psychiatric','psychiatric','demographic','cognitive','cognitive','cognitive','cognitive','cognitive','cognitive','cognitive','cognitive','cognitive','cognitive','cognitive','cognitive','cognitive'};
%           - category ---------- demographic, cognitive, biometric, psychiatric, etc.


% replace outcome
outcome_old=outcome;
outcome_new=struct();
for i = 1:numel(select_measures)
    test_index = select_measures(i);
    test_name = test_names{test_index};
    outcome_new.(test_name) = outcome_old.(test_name);
    outcome_new.(test_name).category = categories{i};
end

outcome = outcome_new;

save(new_data,'study_info','brain_data','outcome')




