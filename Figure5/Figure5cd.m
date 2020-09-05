
load([dionysis 'contact_93_beta_ref.mat'])
% speech onset latency ~ band response interval: p = 0.3779, rho = 0.0896
% speech onset latency ~ band response to speech interval: p = 3.4333e-11,
% rho = 0.6046

Reaction_T = ephy_timing.ReactionT;
ActOn = ephy_timing.ActOn;
ActOn2P = Reaction_T - ActOn;

tr2del = find(isnan(ActOn));

Reaction_T(tr2del) = [];
ActOn(tr2del) = [];
ActOn2P(tr2del) = [];

figure; hold on;
Plot_ActOn = plot(Reaction_T,ActOn, 'b*');
Plot_ActOn2P = plot(Reaction_T,ActOn2P, 'r*');
Plot_ActOn.MarkerSize = 10;
Plot_ActOn2P.MarkerSize = 10;

plt_axis = gca;
set(plt_axis,'XLim',[0.4 2.2],'YLim',[-1 2.5]);
set(plt_axis,'box','on');
set(plt_axis,'TickLength',[0.005,0.005]);
lsline(plt_axis);
xlabel('Speech Onset Latency (s)');
ylabel('Interval Duration (s)');

blue_text.rho = text(1.6,0.65,'\rho = 0.09','Color','blue');
blue_text.p = text(1.78,0.66,'p = 0.38', 'Color','blue');
% delete(blue_text.rho);
% delete(blue_text.p);
red_text.rho = text(1.6,1.0,'\rho = 0.60','Color','red');
red_text.rho.Rotation = 18;
red_text.p = text(1.77,1.16,'p = 3.4 \times 10^{-11}', 'Color','red');
red_text.p.Rotation = 18;
% delete(red_text.rho);
% delete(red_text.p);

% save
saveas(gcf, 'Figure5d.fig');

% turn to high gamma
clearvars
close all;

load('contact_64_session1_highgamma_ref.mat')
% speech onset latency ~ band response interval: p = 1.9287e-15, rho = 0.6664
% speech onset latency ~ band response to speech interval: p = 0.3329,
% rho = -0.0932
Reaction_T = ephy_timing.ReactionT;
ActOn = ephy_timing.ActOn;
ActOn2P = Reaction_T - ActOn;

tr2del = find(isnan(ActOn));

Reaction_T(tr2del) = [];
ActOn(tr2del) = [];
ActOn2P(tr2del) = [];

figure; hold on;
Plot_ActOn = plot(Reaction_T,ActOn, 'b*');
Plot_ActOn2P = plot(Reaction_T,ActOn2P, 'r*');
Plot_ActOn.MarkerSize = 10;
Plot_ActOn2P.MarkerSize = 10;

plt_axis = gca;
set(plt_axis,'XLim',[0.4 2.2],'YLim',[-1 2.5]);
set(plt_axis,'box','on');
set(plt_axis,'TickLength',[0.005,0.005]);
lsline(plt_axis);
xlabel('Speech Onset Latency (s)');
ylabel('Interval Duration (s)');

blue_text.rho = text(1.6,1.66,'\rho = 0.67','Color','blue');
blue_text.rho.Rotation = 21;
blue_text.p = text(1.76,1.85,'p = 1.9 \times 10^{-15}', 'Color','blue');
blue_text.p.Rotation = 21;
% delete(blue_text.rho);
% delete(blue_text.p);
red_text.rho = text(1.6,-0.03,'\rho = -0.09','Color','red');
red_text.rho.Rotation = -4;
red_text.p = text(1.79,-0.06,'p = 0.33', 'Color','red');
red_text.p.Rotation = -3;
% delete(red_text.rho);
% delete(red_text.p);

% save
saveas(gcf, 'Figure5c.fig');