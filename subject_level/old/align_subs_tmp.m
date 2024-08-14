

% may wish to align motion to brain data (remove subs without brain data, add NaN if no motion)
mot_filtered = nan(length(subs_brain), 1);
[common_subs, brain_idx, mot_idx] = intersect(subs_brain, subs_mot);
mot_filtered(brain_idx) = mot(mot_idx);


