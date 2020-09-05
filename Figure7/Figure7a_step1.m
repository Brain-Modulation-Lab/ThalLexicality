%%% * load in files needed
load('speech_response_table_t_highgamma_-0.15_0.59.mat', 'speech_response_table');

load('ref_highgamma_-0.15_0.59.mat', 'Unit_Selectivity_table');

load('contact_info_newMER.mat', 'contact_info');


speech_response_table = gen_q_tbl(speech_response_table);
Unit_Selectivity_table = gen_q_tbl(Unit_Selectivity_table);


%%%%%%%%%%%%%%%% ***: generate index of contact list oi
lead_left_list = setdiff(find(strcmp(speech_response_table.side,'left')),1:33)';
lead_left_session1_list = lead_left_list(find(cell2mat(speech_response_table.ith_session(lead_left_list)) == 1));
lead_left_session2_list = setdiff(lead_left_list, lead_left_session1_list);
lead_right_list = setdiff(find(strcmp(speech_response_table.side,'right')),1:33)';
SigUp_list = find(speech_response_table.p_120 * height(speech_response_table) < 0.05 & ...
    speech_response_table.h_120 == 1);
lead_left_session1_SigUp_list = intersect(SigUp_list, lead_left_session1_list)';
lead_left_session2_SigUp_list = intersect(SigUp_list, lead_left_session2_list)';
lead_right_SigUp_list = intersect(SigUp_list, lead_right_list)';
%%%%%% done ***

%%%%%%% **** prepare union merged table
% ls1_sigup
tbl_ls1_sigup = Unit_Selectivity_table(lead_left_session1_SigUp_list,:);
contact_id_oi = tbl_ls1_sigup.contact_id;
 % get all contacts' coordinates
coords_mess = extractfield(contact_info,'mni_coords'); 
coords_all = reshape(coords_mess,[3,length(coords_mess)/3])';
coords_oi = coords_all(contact_id_oi,:);
tbl_ls1_sigup.X = coords_oi(:,1);
tbl_ls1_sigup.Y = coords_oi(:,2);
tbl_ls1_sigup.Z = coords_oi(:,3);
tbl_ls1_sigup.selectivity = - tbl_ls1_sigup.selectivity;

% ls2_sigup
tbl_ls2_sigup = Unit_Selectivity_table(lead_left_session2_SigUp_list,:);
contact_id_oi = tbl_ls2_sigup.contact_id;
 % get all contacts' coordinates
coords_mess = extractfield(contact_info,'mni_coords'); 
coords_all = reshape(coords_mess,[3,length(coords_mess)/3])';
coords_oi = coords_all(contact_id_oi,:);
tbl_ls2_sigup.X = coords_oi(:,1);
tbl_ls2_sigup.Y = coords_oi(:,2);
tbl_ls2_sigup.Z = coords_oi(:,3);
tbl_ls2_sigup.selectivity = - tbl_ls2_sigup.selectivity;

% merge activity of left s1 and left s2 
contact_intersect = intersect(tbl_ls1_sigup.contact_id, tbl_ls2_sigup.contact_id);
contact_union = union(tbl_ls2_sigup.contact_id, tbl_ls1_sigup.contact_id);

tbl_left_union = table();
stack_tbl{1} = tbl_ls1_sigup; stack_tbl{2} = tbl_ls2_sigup;
for i_order = 1:length(contact_union)
    contact_id = contact_union(i_order);
    tbl_left_union.contact_id(i_order) = contact_id;
    session_oi = find([ismember(contact_id,tbl_ls1_sigup.contact_id), ismember(contact_id, tbl_ls2_sigup.contact_id)]);
    tbl_left_union.ith_session{i_order} = session_oi; % cell
    
    subject_id = stack_tbl{session_oi(1)}.subject_id{find(stack_tbl{session_oi(1)}.contact_id == contact_id)};
    side = stack_tbl{session_oi(1)}.side{find(stack_tbl{session_oi(1)}.contact_id == contact_id)};
    XX = stack_tbl{session_oi(1)}.X(find(stack_tbl{session_oi(1)}.contact_id == contact_id));
    YY = stack_tbl{session_oi(1)}.Y(find(stack_tbl{session_oi(1)}.contact_id == contact_id));
    ZZ = stack_tbl{session_oi(1)}.Z(find(stack_tbl{session_oi(1)}.contact_id == contact_id));
    
    tbl_left_union.subject_id{i_order} = subject_id;
    tbl_left_union.side{i_order} = side;
    tbl_left_union.X(i_order) = XX;
    tbl_left_union.Y(i_order) = YY;
    tbl_left_union.Z(i_order) = ZZ;
    
    if length(session_oi) == 2
        selectivity_s1 = stack_tbl{session_oi(1)}.selectivity(find(stack_tbl{session_oi(1)}.contact_id == contact_id));
        selectivity_s2 = stack_tbl{session_oi(2)}.selectivity(find(stack_tbl{session_oi(2)}.contact_id == contact_id));
        
        selectivity = mean([selectivity_s1,selectivity_s2]);
        
    elseif session_oi == 1
        selectivity_s1 = stack_tbl{1}.selectivity(find(stack_tbl{1}.contact_id == contact_id));
        selectivity_s2 = NaN;
        selectivity = selectivity_s1;
    else % session_oi == 2
        selectivity_s1 = NaN;
        selectivity_s2 = stack_tbl{2}.selectivity(find(stack_tbl{2}.contact_id == contact_id));
        selectivity = selectivity_s2;
    end
    tbl_left_union.selectivity_s1(i_order) = selectivity_s1;
    tbl_left_union.selectivity_s2(i_order) = selectivity_s2;
    tbl_left_union.selectivity(i_order) = selectivity;
end

lead_left_contact_id = speech_response_table.contact_id(lead_left_session1_list);
%% left union merged + right
figure; hold on;
hh = gca;
coords_mess = extractfield(contact_info,'mni_coords'); 
coords_all = reshape(coords_mess,[3,length(coords_mess)/3])';
coords_oi = coords_all(lead_left_contact_id, :);

sig_idx = find(tbl_left_union.selectivity_s1 > 1.645 | tbl_left_union.selectivity_s2 > 1.645);

for i = 1:length(lead_left_contact_id)
    if ismember(lead_left_contact_id(i), tbl_left_union.contact_id)
        idx_oi = find(tbl_left_union.contact_id == lead_left_contact_id(i));
        if ismember(idx_oi, sig_idx)

            scatter3(tbl_left_union.X(idx_oi), tbl_left_union.Y(idx_oi), ...
            tbl_left_union.Z(idx_oi), 80, 'k','filled','o')
        else
        
            scatter3(tbl_left_union.X(idx_oi), tbl_left_union.Y(idx_oi), ...
                tbl_left_union.Z(idx_oi), 80, 'k','o')
        end
    else
        scatter3(coords_oi(i,1), coords_oi(i,2), coords_oi(i,3), 80, 'k','o')
    end
end

grid on;

% label X, Y and Z axis
xlabel('lateral-medial (mm)');
ylabel('anterior-posterior (mm)');
zlabel('ventral-dorsal (mm)');

cam_used = campos;
save('Figure7a_campos.mat', 'cam_used');


campos(cam_used);
% save figure
saveas(gcf, 'Figure7a_step1.fig');
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