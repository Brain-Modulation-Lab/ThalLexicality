
load('speech_response_table_t_highgamma_-0.15_0.59.mat');

load('ref_highgamma_-0.15_0.59.mat');

% get combined-session rows id
unique_contact = unique(speech_response_table.contact_id);
hasi = 0;
discard_idx = [];
for i = 1:length(unique_contact)
    contact_i = unique_contact(i);
    temp = find(speech_response_table.contact_id == contact_i);
    if length(temp) == 1
    else
        hasi = hasi + 1;
    discard_idx(hasi) =  temp(3);
    end
end
speech_response_table(discard_idx,:) = [];

Unit_Selectivity_left_combined = Unit_Selectivity_table(discard_idx,:); % 1

Unit_Selectivity_table(discard_idx,:) = [];

speech_response_table(13:18,:) = [];Unit_Selectivity_table(13:18,:) = [];



lead_left_list = setdiff(find(strcmp(speech_response_table.side,'left')),1:33)';
lead_right_list = setdiff(find(strcmp(speech_response_table.side,'right')),1:33)';
lead_left_session1_list = lead_left_list(find(cell2mat(speech_response_table.ith_session(lead_left_list)) == 1));
lead_left_session2_list = setdiff(lead_left_list, lead_left_session1_list);

SigUp_list = find(speech_response_table.p_120 * height(speech_response_table) < 0.05 & speech_response_table.h_120 == 1);

lead_left_SigUp_list = intersect(SigUp_list, lead_left_list)';
lead_left_session1_SigUp_list = intersect(SigUp_list, lead_left_session1_list)';
lead_left_session2_SigUp_list = intersect(SigUp_list, lead_left_session2_list)';
lead_right_SigUp_list = intersect(SigUp_list, lead_right_list)';

%%%% 2 - 9
Unit_Selectivity_left = Unit_Selectivity_table(lead_left_list,:);
Unit_Selectivity_left_s1 = Unit_Selectivity_table(lead_left_session1_list,:);
Unit_Selectivity_left_s2 = Unit_Selectivity_table(lead_left_session2_list,:);

Unit_Selectivity_left_SigUp = Unit_Selectivity_table(lead_left_SigUp_list,:);
Unit_Selectivity_left_s1_SigUp = Unit_Selectivity_table(lead_left_session1_SigUp_list,:);
Unit_Selectivity_left_s2_SigUp = Unit_Selectivity_table(lead_left_session2_SigUp_list,:);

Unit_Selectivity_right = Unit_Selectivity_table(lead_right_list,:);
Unit_Selectivity_right_SigUp = Unit_Selectivity_table(lead_right_SigUp_list,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%% end *

%%%%%%%%%%%% firstly show the mean and median of the data
mean(Unit_Selectivity_left_combined.selectivity), median(Unit_Selectivity_left_combined.selectivity)
mean(Unit_Selectivity_left.selectivity), median(Unit_Selectivity_left.selectivity)
mean(Unit_Selectivity_left_s1.selectivity), median(Unit_Selectivity_left_s1.selectivity)
mean(Unit_Selectivity_left_s2.selectivity), median(Unit_Selectivity_left_s2.selectivity)
mean(Unit_Selectivity_left_SigUp.selectivity), median(Unit_Selectivity_left_SigUp.selectivity)
mean(Unit_Selectivity_left_s1_SigUp.selectivity), median(Unit_Selectivity_left_s1_SigUp.selectivity)
mean(Unit_Selectivity_left_s2_SigUp.selectivity), median(Unit_Selectivity_left_s2_SigUp.selectivity)
mean(Unit_Selectivity_right.selectivity), median(Unit_Selectivity_right.selectivity)
mean(Unit_Selectivity_right_SigUp.selectivity), median(Unit_Selectivity_right_SigUp.selectivity)

%%%%%%%%%%%%%%%%%%%%%%%%% categorical scatter plot with error bar

%%%%%% 1: left s1 vs. left s2 vs. right sigup


%%%% prepare the data
data_1 = - Unit_Selectivity_left_s1_SigUp.selectivity; % the higher the more nonword selective

mean_1 = mean(data_1); sem_1 = std(data_1)/sqrt(length(data_1)) ; % sem = sd/sqrt(n)

data_2 = - Unit_Selectivity_left_s2_SigUp.selectivity; % the higher the more nonword selective

mean_2 = mean(data_2); sem_2 = std(data_2)/sqrt(length(data_2)) ;

data_3 = - Unit_Selectivity_right_SigUp.selectivity; % the higher the more nonword selective

mean_3 = mean(data_3); sem_3 = std(data_3)/sqrt(length(data_3)) ;


jitter = 0.5;
jitScale = jitter * 0.55;

figure (1); hold on
%%%%%%%%%%%%%%%% X = 1, left_s1 sigup
% plot mean horizontal line
h(1).mu=plot([1-jitScale,1+jitScale],[mean_1,mean_1],'-k',...
                'linewidth',1);
            
%%% plot vertical sem line
h(1).sem_vert=plot([1,1],[mean_1-sem_1,mean_1+sem_1],...
'-k', 'linewidth',1);

%%% plot horizontal sem line

h(1).sem_h1 = plot([1 - jitScale/2,1 + jitScale/2],[mean_1+sem_1,mean_1+sem_1],...
'-k', 'linewidth',1);

h(1).sem_h2 = plot([1 - jitScale/2,1 + jitScale/2],[mean_1-sem_1,mean_1-sem_1],...
'-k', 'linewidth',1);


%%%%%% % make jittered X

X = repmat(1,length(data_1),1);
[counts,~,bins] = histcounts(data_1,10);
inds = find(counts~=0);
counts = counts(inds);
Xr = X;
for jj=1:length(inds)
    tWidth = jitter * (1-exp(-0.1 * (counts(jj)-1)));
    xpoints = linspace(-tWidth*0.8, tWidth*0.8, counts(jj));
    Xr(bins==inds(jj)) = xpoints;
end
X = X+Xr;

%%%% plot dot and give color
h(1).data=plot(X, data_1, 'o','Color','k','MarkerFaceColor',[0.2,0.4,0.9]);


%%%%%%%%%%%%%%%% X = 2, left_s2 sigup
% plot mean horizontal line
h(2).mu=plot([2-jitScale,2+jitScale],[mean_2,mean_2],'-k',...
                'linewidth',1);
            
%%% plot vertical sem line
h(2).sem_vert=plot([2,2],[mean_2-sem_2,mean_2+sem_2],...
'-k', 'linewidth',1);

%%% plot horizontal sem line

h(2).sem_h1 = plot([2 - jitScale/2,2 + jitScale/2],[mean_2+sem_2,mean_2+sem_2],...
'-k', 'linewidth',1);

h(2).sem_h2 = plot([2 - jitScale/2,2 + jitScale/2],[mean_2-sem_2,mean_2-sem_2],...
'-k', 'linewidth',1);


%%%%%% % make jittered X

X = repmat(2,length(data_2),1);
[counts,~,bins] = histcounts(data_2,10);
inds = find(counts~=0);
counts = counts(inds);
Xr = X;
for jj=1:length(inds)
    tWidth = jitter * (1-exp(-0.1 * (counts(jj)-1)));
    xpoints = linspace(-tWidth*0.8, tWidth*0.8, counts(jj));
    Xr(bins==inds(jj)) = xpoints;
end
X = X+Xr;

%%%% plot dot and give color
h(2).data=plot(X, data_2, 'o','Color','k','MarkerFaceColor',[0.4,0.8,0.5]);


%%%%%%%%%%%%%%%% X = 3, right sigup
% plot mean horizontal line
h(3).mu=plot([3-jitScale,3+jitScale],[mean_3,mean_3],'-k',...
                'linewidth',1);
            
%%% plot vertical sem line
h(3).sem_vert=plot([3,3],[mean_3-sem_3,mean_3+sem_3],...
'-k', 'linewidth',1);

%%% plot horizontal sem line

h(3).sem_h1 = plot([3 - jitScale/2,3 + jitScale/2],[mean_3+sem_3,mean_3+sem_3],...
'-k', 'linewidth',1);

h(3).sem_h2 = plot([3 - jitScale/2,3 + jitScale/2],[mean_3-sem_3,mean_3-sem_3],...
'-k', 'linewidth',1);


%%%%%% % make jittered X

X = repmat(3,length(data_3),1);
[counts,~,bins] = histcounts(data_3,10);
inds = find(counts~=0);
counts = counts(inds);
Xr = X;
for jj=1:length(inds)
    tWidth = jitter * (1-exp(-0.1 * (counts(jj)-1)));
    xpoints = linspace(-tWidth*0.8, tWidth*0.8, counts(jj));
    Xr(bins==inds(jj)) = xpoints;
end
X = X+Xr;

%%%% plot dot and give color
h(3).data=plot(X, data_3, 'o','Color','k','MarkerFaceColor',[1,0,0]);

%%%%%% tidy ticks and lims
set(gca,'XTick',[1,2,3]);
set(gca,'XTickLabel',[]);
set(gca,'Box','on');
set(gca,'TickLength',[0.005,0.005]);
xlim([0,4])
ylim([-2,4.5]);

%%%% plot statistical line

stat_h(1).line =plot([2,3],[3.3,3.3],...
'-k', 'linewidth',1);

stat_h(1).text =text(2.5,3.4,'**','FontSize',30,'HorizontalAlignment','center');

stat_h(2).line =plot([1,2],[3.5,3.5],...
'-k', 'linewidth',1);

stat_h(2).text =text(1.5,3.75,'ns','FontSize',20,'HorizontalAlignment','center');

stat_h(3).line =plot([1,3],[3.95,3.95],...
'-k', 'linewidth',1);

stat_h(3).text =text(2,4.05,'**','FontSize',30,'HorizontalAlignment','center');

y_label = ylabel(gca, 'High gamma lexical selectivity index');

x_tick1 = text(1,-2.3,'Left session 1','HorizontalAlignment','center');
x_tick2 = text(2,-2.3,'Left session 2','HorizontalAlignment','center');
x_tick3 = text(3,-2.3,'Right','HorizontalAlignment','center');


saveas(gcf, 'Figure6e.fig');