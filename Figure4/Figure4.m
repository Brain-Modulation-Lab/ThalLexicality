% data sampling frequency and spectral gradient (fq)
fs = 1000;
fq = 2:2:200;

% load in contact_info
load('contact_info_newMER.mat');

% trial_selection
trial_selection = {'120'};
% region selection
region_selection = {'all'};
% generate idx lists for each region
region_idx_lists{1} = setdiff(1:95, 13:18);
% locking selection
locking_selection = {'sp'};
locking_name_lists = {'Speech Onset'};
% get xlim for different locking selection
xlim_locking = [-2,2];

% loop through trial groups
trial_group_id = 1;
trial_group_name = trial_selection{trial_group_id};
    
% loop through locking selection
locking_id = 1:length(locking_selection);
locking_name = locking_selection{locking_id};
xlim_lower = xlim_locking(locking_id,1);xlim_upper = xlim_locking(locking_id,2);
        
% loop through regions
region_id = 1:length(region_selection);
region_name = region_selection{region_id};
region_idx = region_idx_lists{region_id};
            
clearvars counter z_total timing_total
counter = 0 ; timing_total = [0,0];

% loop through contacts in this region
for contact_id = region_idx
    if length(contact_info(contact_id).session) == 1 % matrix
        counter = counter + 1;
        load([dionysis 'Users/dwang/VIM/datafiles/preprocessed_new/v2/WordVsNonword/contact_spectrogram/'...
        'ref_contact_' num2str(contact_id) '_' trial_group_name '_' locking_name '.mat']);
        timing_total = timing_total + timing_oi;
        z_total(:,:,counter) = data_final;
    else % two sessions
        for session_id = 1:2 % load in one at a time
            counter = counter + 1;
            load([dionysis 'Users/dwang/VIM/datafiles/preprocessed_new/v2/WordVsNonword/contact_spectrogram/'...
            'ref_contact_' num2str(contact_id) '_session' num2str(session_id) '_' trial_group_name '_' locking_name '.mat']);    
            timing_total = timing_total + timing_oi;
            z_total(:,:,counter) = data_final;
        end
    end
end

% get grand averaged z and timing
z_grand_avg = mean(z_total,3);
avg_timing = timing_total./counter;


%default colormap (parula)
% time points
t=linspace(xlim_lower,xlim_upper,size(z_grand_avg,1));
% frequencies
fq;
figure; imagesc(t, fq, z_grand_avg(:,:)');set(gca, 'YDir', 'Normal');
% modify figures
caxis([-3,2]); box on;
set(gca,'TickLength',[0.005,0.005]);

% axis label
xlabel(['Time Relative to ' locking_name_lists{locking_id} ' (s)']); % x-axis label
ylabel('Frequency (Hz)'); % y-axis label

% plot timepoints, black dashed line
hold on;
plot([avg_timing(1),avg_timing(1)],ylim,'k--');
plot([0,0],ylim,'k--'); 
plot([avg_timing(2),avg_timing(2)],ylim,'k--');


% plot colorbar
cb = colorbar;
cb.Label.String = 'z-scored power';
set(cb,'TickLength',0.005);

saveas(gcf,'Figure4.fig');
close all;