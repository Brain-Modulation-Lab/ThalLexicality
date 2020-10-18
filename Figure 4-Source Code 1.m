% Load in time frequency data
load('/Users/dengyu/Dropbox (Brain Modulation Lab)/Dengyu/Thal_lexicality_paper/CodeandDataAvailability/Figure4/Figure 4-Source Data 1.mat','z_grand_avg');
% Load in time points of average stimulus onset and speech offset
load('/Users/dengyu/Dropbox (Brain Modulation Lab)/Dengyu/Thal_lexicality_paper/CodeandDataAvailability/Figure4/Figure 4-Source Data 1.mat','avg_timing');

% data sampling frequency and spectral gradient (fq)
fs = 1000;
fq = 2:2:200;

% Time range of the time-frequency plot
xlim_lower = -2;
xlim_upper = 2;

%default colormap (parula)

% time points
t=linspace(xlim_lower,xlim_upper,size(z_grand_avg,1));

% plot
figure; imagesc(t, fq, z_grand_avg(:,:)');set(gca, 'YDir', 'Normal');

% set figure parameters
caxis([-3,2]); box on;
set(gca,'TickLength',[0.005,0.005]);

% axis label
xlabel('Time relative to speech onset (s)'); % x-axis label
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