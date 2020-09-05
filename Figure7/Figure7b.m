
%%% * load in files needed

load('speech_response_table_t_highgamma_-0.15_0.59.mat','speech_response_table');

load('ref_highgamma_-0.15_0.59.mat','Unit_Selectivity_table');

load('contact_info_newMER.mat', 'contact_info');
%%% done *


speech_response_table = gen_q_tbl(speech_response_table);
Unit_Selectivity_table = gen_q_tbl(Unit_Selectivity_table);

%%%%%%%%%%%%%%%% ***: generate index of contact list oi
lead_left_list = setdiff(find(strcmp(speech_response_table.side,'left')),1:33)';
lead_left_session1_list = lead_left_list(find(cell2mat(speech_response_table.ith_session(lead_left_list)) == 1));
lead_left_session2_list = setdiff(lead_left_list, lead_left_session1_list);
lead_right_list = setdiff(find(strcmp(speech_response_table.side,'right')),1:33)';
SigUp_list = find(speech_response_table.p_120 * height(speech_response_table) < 0.05 & ...
    speech_response_table.h_120 == 1);
lead_left_SigUp_list = intersect(SigUp_list, lead_left_list)';
lead_left_session1_SigUp_list = intersect(SigUp_list, lead_left_session1_list)';
lead_left_session2_SigUp_list = intersect(SigUp_list, lead_left_session2_list)';
lead_right_SigUp_list = intersect(SigUp_list, lead_right_list)';
%%%%%% done ***

tbl_left_sigup = Unit_Selectivity_table(lead_left_SigUp_list,:);

tbl_oi = DW_add_loc_info(tbl_left_sigup);

VLp_idx = find(strcmp('VLp',tbl_oi.Group_Used));
VLAA_idx = setdiff(1:height(tbl_oi), VLp_idx);
[h, p, ~, stat] = ttest2(- tbl_oi.selectivity(VLAA_idx),  - tbl_oi.selectivity(VLp_idx));

figure; hold on;
hh = gca;

VLAA_plt = scatter(tbl_oi.coord_y(VLAA_idx), -tbl_oi.selectivity(VLAA_idx), [], ...
    'MarkerFaceColor', [0.433,0.738,0.621], 'MarkerEdgeColor',[0.433,0.738,0.621]);
VLp_plt = scatter(tbl_oi.coord_y(VLp_idx), -tbl_oi.selectivity(VLp_idx),[],...
    'MarkerFaceColor', [0.125,0.324,0.830],'MarkerEdgeColor',[0.125,0.324,0.830]);

set(hh,'XLim',[-22 -8],'YLim',[-2 3]);
set(hh,'box','on');
set(hh,'TickLength',[0.005,0.005]);

% plot regression line
P = polyfit(tbl_oi.coord_y,-tbl_oi.selectivity,1);
x0 = -22 ; x1 = -8;
xi = linspace(x0,x1) ;
yi = P(1)*xi+P(2);
regression_line = plot(xi,yi,'k');

% add spearman stat result
stat_text.rho = text(-11.5, 2.7,'Pearson''s \rho = 0.49');
stat_text.p = text(-11.5, 2.5,'p = 4.1 \times 10^{-4}');
%

xlabel('posterior-anterior (mm)');
ylabel('High gamma lexical selectivity index');
hh_legend = legend('VA+VLa', 'VLp','Location','northwest');
saveas(gcf, 'Figure7b.fig');
close all


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