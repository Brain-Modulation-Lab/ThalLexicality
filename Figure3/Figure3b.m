performance_tbl = readtable('/subject_based_behavioral_performance.xlsx');

data_1 = performance_tbl.speech_duration(1:3:end);
SE_list_1 = performance_tbl.speech_duration_SE(1:3:end);
% group level
mean_1 = mean(data_1); SEM_1 = std(data_1)/sqrt(length(data_1));

data_2 = performance_tbl.speech_duration(2:3:end);
SE_list_2 = performance_tbl.speech_duration_SE(2:3:end);
% group level
mean_2 = mean(data_2); SEM_2 = std(data_2)/sqrt(length(data_2));

data_3 = performance_tbl.speech_duration(3:3:end);
SE_list_3 = performance_tbl.speech_duration_SE(3:3:end);
% group level
mean_3 = mean(data_3); SEM_3 = std(data_3)/sqrt(length(data_3));

%%%%%%%%%%%%%%%%%%%%%%%%% categorical scatter plot with error bar

jitter = 0.3;
jitScale = jitter * 0.55;

figure; hold on
%%%%%%%%%%%%%%%% X = 1

%%%%%% % make jittered X

X = repmat(1,length(data_1),1);
[counts,~,bins] = histcounts(data_1,5);
inds = find(counts~=0);
counts = counts(inds);
Xr = X;
for jj=1:length(inds)
    tWidth = jitter * (1-exp(-0.1 * (counts(jj)-1)));
    xpoints = linspace(-tWidth*0.8, tWidth*0.8, counts(jj));
    Xr(bins==inds(jj)) = xpoints;
end
X = X+Xr;
X_1 = X;

for dot = 1:length(data_1)
    h(1).dot_bar(dot).vert = plot([X(dot) X(dot)], [data_1(dot) - SE_list_1(dot), data_1(dot) + SE_list_1(dot)],'-k');
    h(1).dot_bar(dot).sem_up = plot([X(dot) - 0.01, X(dot) + 0.01], [data_1(dot) + SE_list_1(dot), data_1(dot) + SE_list_1(dot)],'-k');
    h(1).dot_bar(dot).sem_low = plot([X(dot) - 0.01, X(dot) + 0.01], [data_1(dot) - SE_list_1(dot), data_1(dot) - SE_list_1(dot)],'-k');
    h(1).dot_bar(dot).dot = plot(X(dot), data_1(dot), 'o','Color','k','MarkerFaceColor',[0.5 0.5 0.5]);
end



% clearvars data_1 SE_list_1 mean_1 SEM_1
%%%%%%%%%%%%%%%% X = 2

%%%%%% % make jittered X

X = repmat(2,length(data_2),1);
[counts,~,bins] = histcounts(data_2,5);
inds = find(counts~=0);
counts = counts(inds);
Xr = X;
for jj=1:length(inds)
    tWidth = jitter * (1-exp(-0.1 * (counts(jj)-1)));
    xpoints = linspace(-tWidth*0.8, tWidth*0.8, counts(jj));
    Xr(bins==inds(jj)) = xpoints;
end
X = X+Xr;
X_2 = X;

for dot = 1:length(data_2)
    h(2).dot_bar(dot).vert = plot([X(dot) X(dot)], [data_2(dot) - SE_list_2(dot), data_2(dot) + SE_list_2(dot)],'-k');
    h(2).dot_bar(dot).sem_up = plot([X(dot) - 0.01, X(dot) + 0.01], [data_2(dot) + SE_list_2(dot), data_2(dot) + SE_list_2(dot)],'-k');
    h(2).dot_bar(dot).sem_low = plot([X(dot) - 0.01, X(dot) + 0.01], [data_2(dot) - SE_list_2(dot), data_2(dot) - SE_list_2(dot)],'-k');
    h(2).dot_bar(dot).dot = plot(X(dot), data_2(dot), 'o','Color','k','MarkerFaceColor',[0.5 0.5 0.5]);
end



% clearvars data_2 SE_list_2 mean_2 SEM_2
%%%%%%%%%%%%%%%% X = 3

%%%%%% % make jittered X

X = repmat(3,length(data_3),1);
[counts,~,bins] = histcounts(data_3,5);
inds = find(counts~=0);
counts = counts(inds);
Xr = X;
for jj=1:length(inds)
    tWidth = jitter * (1-exp(-0.1 * (counts(jj)-1)));
    xpoints = linspace(-tWidth*0.8, tWidth*0.8, counts(jj));
    Xr(bins==inds(jj)) = xpoints;
end
X = X+Xr;
X_3 = X;

for dot = 1:length(data_3)
    h(3).dot_bar(dot).vert = plot([X(dot) X(dot)], [data_3(dot) - SE_list_3(dot), data_3(dot) + SE_list_3(dot)],'-k');
    h(3).dot_bar(dot).sem_up = plot([X(dot) - 0.01, X(dot) + 0.01], [data_3(dot) + SE_list_3(dot), data_3(dot) + SE_list_3(dot)],'-k');
    h(3).dot_bar(dot).sem_low = plot([X(dot) - 0.01, X(dot) + 0.01], [data_3(dot) - SE_list_3(dot), data_3(dot) - SE_list_3(dot)],'-k');
    h(3).dot_bar(dot).dot = plot(X(dot), data_3(dot), 'o','Color','k','MarkerFaceColor',[0.5 0.5 0.5]);
end
% clearvars data_3 SE_list_3 mean_3 SEM_3


for dot_id = 1:length(data_2)
    corr_line(dot_id) = plot([X_2(dot_id) X_3(dot_id)],[data_2(dot_id), data_3(dot_id)],'Color',[0.5,0.5,0.5]);
end


% plot mean horizontal line
h(1).mu=plot([1-jitScale,1+jitScale],[mean_1,mean_1],'-k',...
                'linewidth',2);
            
%%% plot vertical sem line
h(1).sem_vert=plot([1,1],[mean_1-SEM_1,mean_1+SEM_1],...
'-k', 'linewidth',2);

%%% plot horizontal sem line

h(1).sem_h1 = plot([1 - jitScale/2,1 + jitScale/2],[mean_1+SEM_1,mean_1+SEM_1],...
'-k', 'linewidth',2);

h(1).sem_h2 = plot([1 - jitScale/2,1 + jitScale/2],[mean_1-SEM_1,mean_1-SEM_1],...
'-k', 'linewidth',2);

% plot mean horizontal line
h(2).mu=plot([2-jitScale,2+jitScale],[mean_2,mean_2],'-k',...
                'linewidth',2);
            
%%% plot vertical sem line
h(2).sem_vert=plot([2,2],[mean_2-SEM_2,mean_2+SEM_2],...
'-k', 'linewidth',2);

%%% plot horizontal sem line

h(2).sem_h1 = plot([2 - jitScale/2,2 + jitScale/2],[mean_2+SEM_2,mean_2+SEM_2],...
'-k', 'linewidth',2);

h(2).sem_h2 = plot([2 - jitScale/2,2 + jitScale/2],[mean_2-SEM_2,mean_2-SEM_2],...
'-k', 'linewidth',2);

% plot mean horizontal line
h(3).mu=plot([3-jitScale,3+jitScale],[mean_3,mean_3],'-k',...
                'linewidth',2);
            
%%% plot vertical sem line
h(3).sem_vert=plot([3,3],[mean_3-SEM_3,mean_3+SEM_3],...
'-k', 'linewidth',2);

%%% plot horizontal sem line

h(3).sem_h1 = plot([3 - jitScale/2,3 + jitScale/2],[mean_3+SEM_3,mean_3+SEM_3],...
'-k', 'linewidth',2);

h(3).sem_h2 = plot([3 - jitScale/2,3 + jitScale/2],[mean_3-SEM_3,mean_3-SEM_3],...
'-k', 'linewidth',2);

%%%%%% tidy ticks and lims
set(gca,'XTick',[1,2,3]);
set(gca,'XTickLabel',[]);
set(gca,'Box','on');
set(gca,'TickLength',[0.005,0.005]);
xlim([0,4])
ylim([0.3,1.2]);
set(gca,'YTick',0.3:0.2:1.2);

y_label = ylabel(gca, 'Speech duration (s)');
x_tick1 = text(1,0.3-0.036,'All trials','HorizontalAlignment','center');
x_tick2 = text(2,0.3-0.036,'Word trials','HorizontalAlignment','center');
x_tick3 = text(3,0.3-0.036,'Nonword trials','HorizontalAlignment','center');

% plot statistical line
stat_h.line =plot([2,3],[1.11,1.11],...
'-k', 'linewidth',1);

stat_h.text =text(2.5,1.15,'ns','FontSize',20,'HorizontalAlignment','center');

saveas(gcf, 'Figure3b.fig');