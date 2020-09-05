fs = 1000;
load('speech_response_table_t_highgamma_-0.15_0.59.mat');

speech_response_table= gen_q_tbl(speech_response_table);

% find contact_list oi: left_all, right, left_session1, left_session2, 
lead_left_list = setdiff(find(strcmp(speech_response_table.side,'left')),1:33)';
lead_right_list = setdiff(find(strcmp(speech_response_table.side,'right')),1:33)';
lead_left_session1_list = lead_left_list(find(cell2mat(speech_response_table.ith_session(lead_left_list)) == 1));
lead_left_session2_list = setdiff(lead_left_list, lead_left_session1_list);

SigUp_list = find(speech_response_table.p_120 * height(speech_response_table) < 0.05 & speech_response_table.h_120 == 1);

lead_left_SigUp_list = intersect(SigUp_list, lead_left_list)';
lead_right_SigUp_list = intersect(SigUp_list, lead_right_list)';
lead_left_session1_SigUp_list = intersect(SigUp_list, lead_left_session1_list)';
lead_left_session2_SigUp_list = intersect(SigUp_list, lead_left_session2_list)';


zGrand = [];
epochGrand = table();
numerator = 0;
clearvars word_num_list nonword_num_list

for idx = lead_left_SigUp_list
    numerator = numerator + 1;
    contact_id = speech_response_table.contact_id(idx);
    session_id = speech_response_table.ith_session{idx}; % num
    
    if strcmp('left',speech_response_table.side(idx))
        load(['contact_' num2str(contact_id) '_session' num2str(session_id) '_highgamma_ref.mat']);
    else % right
        load(['contact_' num2str(contact_id) '_highgamma_ref.mat']);
    end
        
    epoch_used = epoch_oi(:,{'id','starts','ends','onset_word','offset_word','onset_vowel','offset_vowel','ITI_starts','stimulus_starts','trial_id'});
    
    zGrand = [zGrand, z_oi]; epochGrand = [epochGrand;epoch_used];
    
    word_num_list(numerator) = sum(mod(epoch_used.trial_id,2)==0 & epoch_used.trial_id<=60); % word num of this unit
    nonword_num_list(numerator) = sum(mod(epoch_used.trial_id,2)==1 & epoch_used.trial_id<=60); % nonword num of this unit       
end

bases_starts = num2cell(round((epochGrand.stimulus_starts - epochGrand.starts - 1) * fs)');
bases_ends = num2cell(round((epochGrand.stimulus_starts - epochGrand.starts) * fs)');
basesGrand = cellfun(@(x,y,z) x(y:z),zGrand,bases_starts,bases_ends,'UniformOutput',0);

roi_starts = num2cell(round((epochGrand.onset_word - epochGrand.starts - 2) * fs)');
roi_ends = num2cell(round((epochGrand.onset_word - epochGrand.starts + 2) * fs)');
trialsGrand = cellfun(@(x,y,z) x(y:z),zGrand,roi_starts,roi_ends,'UniformOutput',0);

word_idxs = find(mod(epochGrand.trial_id,2)==0 & epochGrand.trial_id<=60); % word indexs
nonword_idxs = find(mod(epochGrand.trial_id,2)==1 & epochGrand.trial_id<=60); % nonword indexs

% respectively pick out word trials, nonword trials, word bases and nonword
% bases
word_trials = trialsGrand(word_idxs); nonword_trials = trialsGrand(nonword_idxs);
word_bases = basesGrand(word_idxs); nonword_bases = basesGrand(nonword_idxs);


word_trials_s = cellfun(@(x) smoothdata(x,'movmean',200),word_trials,'UniformOutput',0);
nonword_trials_s = cellfun(@(x) smoothdata(x,'movmean',200),nonword_trials,'UniformOutput',0);

clearvars word_trials_split word_trials_s_split nonword_trials_split nonword_trials_s_split
for parts_id = 1:length(word_num_list)
    
    if parts_id ==1
        start_oi = 1; end_oi = word_num_list(parts_id);
        
        word_trials_split{parts_id} = word_trials(start_oi:end_oi);
        word_trials_s_split{parts_id} = word_trials_s(start_oi:end_oi);
    else
        start_oi = start_oi + word_num_list(parts_id-1);
        end_oi = end_oi + word_num_list(parts_id);
        
        word_trials_split{parts_id} = word_trials(start_oi:end_oi);
        word_trials_s_split{parts_id} = word_trials_s(start_oi:end_oi);
    end
end

for parts_id = 1:length(nonword_num_list)
    
    if parts_id ==1
        start_oi = 1; end_oi = nonword_num_list(parts_id);
        nonword_trials_split{parts_id} = nonword_trials(start_oi:end_oi);
        nonword_trials_s_split{parts_id} = nonword_trials_s(start_oi:end_oi);
    else
        start_oi = start_oi + nonword_num_list(parts_id-1);
        end_oi = end_oi + nonword_num_list(parts_id);
        nonword_trials_split{parts_id} = nonword_trials(start_oi:end_oi);
        nonword_trials_s_split{parts_id} = nonword_trials_s(start_oi:end_oi);
    end
end
avg_word_trial_perunit = cellfun(@(x) mean(cell2mat(x')),word_trials_split,'UniformOutput',0);
avg_word_trial_s_perunit = cellfun(@(x) mean(cell2mat(x')),word_trials_s_split,'UniformOutput',0);
avg_nonword_trial_perunit = cellfun(@(x) mean(cell2mat(x')),nonword_trials_split,'UniformOutput',0);
avg_nonword_trial_s_perunit = cellfun(@(x) mean(cell2mat(x')),nonword_trials_s_split,'UniformOutput',0);

word_bases_mat = cell2mat(word_bases');
word_base_one = mean(word_bases_mat);
word_norm_mean = mean(word_base_one); word_norm_std = std(word_base_one);


nonword_bases_mat = cell2mat(nonword_bases');
nonword_base_one = mean(nonword_bases_mat);
nonword_norm_mean = mean(nonword_base_one); nonword_norm_std = std(nonword_base_one);

avg_word_trial_perunit_global_norm = cellfun(@(x) (x-word_norm_mean)./word_norm_std,avg_word_trial_perunit,'UniformOutput',0);
avg_word_trial_s_perunit_global_norm = cellfun(@(x) (x-word_norm_mean)./word_norm_std,avg_word_trial_s_perunit,'UniformOutput',0);

avg_nonword_trial_perunit_global_norm = cellfun(@(x) (x-nonword_norm_mean)./nonword_norm_std,avg_nonword_trial_perunit,'UniformOutput',0);
avg_nonword_trial_s_perunit_global_norm = cellfun(@(x) (x-nonword_norm_mean)./nonword_norm_std,avg_nonword_trial_s_perunit,'UniformOutput',0);

word_vals = cell2mat(cellfun(@(x) mean(x(2000:2500)),avg_word_trial_perunit_global_norm,'UniformOutput',0));
nonword_vals = cell2mat(cellfun(@(x) mean(x(2000:2500)),avg_nonword_trial_perunit_global_norm,'UniformOutput',0));
[h,p,~,stat] = ttest(word_vals, nonword_vals)

word_vals_s = cell2mat(cellfun(@(x) mean(x(2000:2500)),avg_word_trial_s_perunit_global_norm,'UniformOutput',0));
nonword_vals_s = cell2mat(cellfun(@(x) mean(x(2000:2500)),avg_nonword_trial_s_perunit_global_norm,'UniformOutput',0));
[h,p,~,stat] = ttest(word_vals_s, nonword_vals_s)


word_trials_sn = cellfun(@(x) (x- word_norm_mean)./ word_norm_std, word_trials_s, 'UniformOutput', 0);
nonword_trials_sn = cellfun(@(x) (x- nonword_norm_mean)./ nonword_norm_std, nonword_trials_s, 'UniformOutput', 0);

word_trials_sn_mat = cell2mat(word_trials_sn');
word_tt = mean(word_trials_sn_mat); word_std = std(word_trials_sn_mat); word_se = word_std./sqrt(size(word_trials,2));
word_tt_d = resample(word_tt,4,1000); word_se_d = resample(word_se,4,1000);


nonword_trials_sn_mat = cell2mat(nonword_trials_sn');
nonword_tt = mean(nonword_trials_sn_mat); nonword_std = std(nonword_trials_sn_mat); nonword_se = nonword_std./sqrt(size(nonword_trials,2));
nonword_tt_d = resample(nonword_tt,4,1000); nonword_se_d = resample(nonword_se,4,1000);

Tline = -2:0.25:2;
figure;

sig_line = plot([0.1 0.5], [8.5 8.5],'k','LineWidth',5);

hold on; l_s1 = shadedErrorBar(Tline,word_tt_d,word_se_d,'lineprops','-blue','patchSaturation',0.1);
l_s1.edge(1).Visible = 'off';l_s1.edge(2).Visible = 'off';
l_s1.mainLine.LineWidth = 2; l_s1.mainLine.Color =[0.2,0.4,0.9];l_s1.patch.FaceColor =  [0.2,0.4,0.9];

hold on; r_s = shadedErrorBar(Tline,nonword_tt_d,nonword_se_d,'lineprops','-r','patchSaturation',0.1);
r_s.edge(1).Visible = 'off';r_s.edge(2).Visible = 'off';
r_s.mainLine.LineWidth = 2; r_s.mainLine.Color =[1,0,0];r_s.patch.FaceColor = [1,0,0];


h_ax = gca;
ylim([-1 9])
set(h_ax,'XTick',-2:0.5:2);
set(h_ax,'box','on');
set(h_ax,'TickLength',[0.005,0.005]);
xlabel(h_ax,'Time relative to speech onset (s)');
ylabel(h_ax,'z-scored high gamma power');

hold on; 
cueMean = mean(epochGrand.stimulus_starts(word_idxs) - epochGrand.onset_word(word_idxs));
spoffMean = mean(epochGrand.offset_word(word_idxs) - epochGrand.onset_word(word_idxs));
word_cue_line = plot([cueMean,cueMean], ylim);
word_cue_line.Color = [0.2,0.4,0.9];
word_cue_line.LineStyle = '--';    

word_spoff_line = plot([spoffMean,spoffMean], ylim);
word_spoff_line.Color = [0.2,0.4,0.9];
word_spoff_line.LineStyle = '--';


cueMean = mean(epochGrand.stimulus_starts(nonword_idxs) - epochGrand.onset_word(nonword_idxs));
spoffMean = mean(epochGrand.offset_word(nonword_idxs) - epochGrand.onset_word(nonword_idxs));
nonword_cue_line = plot([cueMean,cueMean], ylim);
nonword_cue_line.Color = [1,0,0];
nonword_cue_line.LineStyle = '--';

nonword_spoff_line = plot([spoffMean,spoffMean], ylim);
nonword_spoff_line.Color = [1,0,0];
nonword_spoff_line.LineStyle = '--';


saveas(gcf, 'Figure6a.fig');
%saveas(gcf, 'Figure6c.fig');


function [tidy_table , two_session_contact_id] = gen_q_tbl(tbl)

unique_contact = unique(tbl.contact_id);
hasi = 0;
discard_idx = [];
for i = 1:length(unique_contact)
    contact_i = unique_contact(i);
    temp = find(tbl.contact_id == contact_i);
    if length(temp) == 1
    else
        hasi = hasi + 1;
    discard_idx(hasi) =  temp(3);
    end
end

two_session_contact_id = tbl.contact_id(discard_idx);
tbl(discard_idx,:) = [];

dbs4039_idx = find(strcmp(tbl.subject_id,'DBS4039'));
tbl(dbs4039_idx,:) = [];
tidy_table = tbl;
end