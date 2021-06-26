% Source code used to generate Figure 6c

% Add shadedErrorBar.m to path

% sampling frequency
fs = 1000;

% load in word trials, nonword trials, and their respective timing data
load('/Users/dengyu/Dropbox (Brain Modulation Lab)/Dengyu/Thal_lexicality_paper/CodeandDataAvailability/Figure5/Figure 5-Source Data 3.mat');

% get each time points' mean value and se
word_mean = mean(word_trials); word_std = std(word_trials); word_se = word_std./sqrt(size(word_trials,1));

nonword_mean = mean(nonword_trials); nonword_std = std(nonword_trials); nonword_se = nonword_std./sqrt(size(nonword_trials,1));

% down-sample for plotting
word_mean_d = resample(word_mean,4,1000); word_se_d = resample(word_se,4,1000);

nonword_mean_d = resample(nonword_mean,4,1000); nonword_se_d = resample(nonword_se,4,1000);

% time axis
Tline = -1:0.25:3;

% plot
figure; hold on;

hold on; word_line = shadedErrorBar(Tline,word_mean_d,word_se_d,'lineprops','-blue','patchSaturation',0.1);
word_line.edge(1).Visible = 'off';word_line.edge(2).Visible = 'off';
word_line.mainLine.LineWidth = 2; word_line.mainLine.Color =[0.2,0.4,0.9];word_line.patch.FaceColor =  [0.2,0.4,0.9];

hold on; nonword_line = shadedErrorBar(Tline,nonword_mean_d,nonword_se_d,'lineprops','-r','patchSaturation',0.1);
nonword_line.edge(1).Visible = 'off';nonword_line.edge(2).Visible = 'off';
nonword_line.mainLine.LineWidth = 2; nonword_line.mainLine.Color =[1,0,0];nonword_line.patch.FaceColor = [1,0,0];

% set figure parameters
h_ax = gca;
ylim([-8 4])
set(h_ax,'XTick',-1:0.5:3);
set(h_ax,'box','on');
set(h_ax,'TickLength',[0.005,0.005]);
xlabel(h_ax,'Time relative to stimulus presentation (s)');
ylabel(h_ax,'z-scored beta power');

% plot timing lines
hold on;

% timing lines for word
sponmean = mean(word_trials_timing.onset_word - word_trials_timing.stimulus_starts);
spoffMean = mean(word_trials_timing.offset_word - word_trials_timing.stimulus_starts);
word_spon_line = plot([sponmean,sponmean], ylim);
word_spon_line.Color = [0.2,0.4,0.9];
word_spon_line.LineStyle = '--';    

word_spoff_line = plot([spoffMean,spoffMean], ylim);
word_spoff_line.Color = [0.2,0.4,0.9];
word_spoff_line.LineStyle = '--';

% timing lines for nonword
sponmean = mean(nonword_trials_timing.onset_word - nonword_trials_timing.stimulus_starts);
spoffMean = mean(nonword_trials_timing.offset_word - nonword_trials_timing.stimulus_starts);
nonword_spon_line = plot([sponmean,sponmean], ylim);
nonword_spon_line.Color = [1,0,0];
nonword_spon_line.LineStyle = '--';

nonword_spoff_line = plot([spoffMean,spoffMean], ylim);
nonword_spoff_line.Color = [1,0,0];
nonword_spoff_line.LineStyle = '--';