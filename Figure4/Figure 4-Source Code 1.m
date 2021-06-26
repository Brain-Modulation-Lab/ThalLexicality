% Source code used to generate Figure 5a

% sampling frequency
fs = 1000;

% load in trial-level LFP data of the recording
load('/Users/dengyu/Dropbox (Brain Modulation Lab)/Dengyu/Thal_lexicality_paper/CodeandDataAvailability/Figure4/Figure 4-Source Data 1.mat', 'z_oi');
% load in timing info of the trials
load('/Users/dengyu/Dropbox (Brain Modulation Lab)/Dengyu/Thal_lexicality_paper/CodeandDataAvailability/Figure4/Figure 4-Source Data 1.mat', 'epoch_oi');
% load in trial order according to reaction time
load('/Users/dengyu/Dropbox (Brain Modulation Lab)/Dengyu/Thal_lexicality_paper/CodeandDataAvailability/Figure4/Figure 4-Source Data 1.mat', 'ephy_timing');

% sort z_oi and epoch_oi according to the trialorder stored in ephy_timing
st_z = z_oi(ephy_timing.TrialOrder); % each trial is 2s precue to 2s post speech offset
st_epoch = epoch_oi(ephy_timing.TrialOrder,:);

% delete trials whose ActOn is missing (because trial_oi not activated)
tr2del = find(isnan(ephy_timing.ActOn)); 

st_z(tr2del) = [];
st_epoch(tr2del,:) = [];
ephy_timing(tr2del,:) = [];

% smooth
st_z_smooth = cellfun(@(x) smooth(x',200)',st_z,'UniformOutput',0);

%baseline (which is 1s pre-cue)
Bt = cellfun(@(x) x(1001:2000), st_z_smooth, 'UniformOutput',0);

% normalize to baseline
zz = cellfun(@(x,y) (x-mean(y))./std(y), st_z_smooth, Bt,'UniformOutput',0);

% define region to be plotted
roi_starts = num2cell((st_epoch.stimulus_starts - st_epoch.starts - 0.5) * fs)'; % 0.5s precue

% default roi_ends is 3s post cue, but if the min length of trials
% is less than 5s, then use the shortest trial length as common
% length
if min(cell2mat(cellfun(@length, zz, 'UniformOutput',0))) < 5000 %  
   length_used = floor(min(cell2mat(cellfun(@length, zz, 'UniformOutput',0)))/100) * 100;

   roi_ends = num2cell(repmat(length_used,[1 length(zz)]));

else
   roi_ends = num2cell((st_epoch.stimulus_starts - st_epoch.starts + 3) * fs)';
end

% prepare data to be plotted
zz_crop = cellfun(@(x,y,z) z(x:y), roi_starts, roi_ends, zz,'UniformOutput',0);
plot_mat = cell2mat(zz_crop');

% time points
tp = linspace(-0.5, length(plot_mat)/1000 -0.5 - 0.001, length(plot_mat));

%plot
imagesc(tp,1:size(plot_mat,1),plot_mat); set(gca, 'YDir', 'Normal');

% set some plotting parameters
hold on; time_spon =plot(ephy_timing.ReactionT,(1:length(ephy_timing.ReactionT)),'-k','LineWidth',3);
hold on; time_cue = plot([0,0],ylim,'k--');
hold on; plot(ephy_timing.ActOn,(1:length(ephy_timing.ActOn)),'r*','MarkerSize',5); % red for beta

hh = gca;
set(hh,'box','on');
set(hh,'TickLength',[0.005,0.005]);
set(hh,'XTickLabel',-0.5:0.5:3);
colormap parula; caxis([-6,6]);
pcolorbar = colorbar;
set(pcolorbar,'TickLength',0.005);
pcolorbar.Label.String = "z-scored power";

xlabel("Time relative to stimulus presentation (s)");
ylabel("Trial number");