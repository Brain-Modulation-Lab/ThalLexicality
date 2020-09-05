load('SumTimingTable_highgamma_ref.mat');

sp_lock_idx = find(SumTimingTable.ActOnCorrP < 0.05 & abs(SumTimingTable.ActOn2PCorrRho) < abs(SumTimingTable.ActOnCorrRho));
cue_lock_idx = find(SumTimingTable.ActOn2PCorrP < 0.05 & abs(SumTimingTable.ActOn2PCorrRho) > abs(SumTimingTable.ActOnCorrRho));
not_lock_idx = find(SumTimingTable.ActOnCorrP >= 0.05 & SumTimingTable.ActOn2PCorrP >=0.05);

figure; hold on;
xlim([-0.2 1]); ylim([-0.2 1]);
plot(SumTimingTable.ActOn2PCorrRho(sp_lock_idx), SumTimingTable.ActOnCorrRho(sp_lock_idx),'bo','MarkerSize',6,'MarkerEdgeColor','k',...
    'MarkerFaceColor','b');
plot(SumTimingTable.ActOn2PCorrRho(cue_lock_idx), SumTimingTable.ActOnCorrRho(cue_lock_idx),'ro', 'MarkerSize',6,'MarkerEdgeColor','k',...
    'MarkerFaceColor','r');
plot(SumTimingTable.ActOn2PCorrRho(not_lock_idx), SumTimingTable.ActOnCorrRho(not_lock_idx),'ko',  'MarkerSize',6,'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0.5,0.5,0.5]);

plot([-0.2 1], [-0.2 1],'k:');
plt_axis = gca;
set(plt_axis,'box','on');
set(plt_axis,'TickLength',[0.005,0.005]);
legend('speech locked', 'cue locked', 'not locked','Location','northwest')
xlabel('\rho (speech onset latency vs. beta response to speech interval)');
ylabel('\rho (speech onset latency vs. beta response latency)');

saveas(gcf, 'Figure5e.fig');