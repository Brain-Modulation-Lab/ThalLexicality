% Color list
% 1 red (255, 0, 0)
% 2 lime (0, 255, 0)
% 3 gray (128, 128, 128)
% 4 yellow (255, 255, 0)
% 5 aqua (0, 255, 255)
% 6 purple (128, 0, 128)
% 7 navy (0,0,128)
% 8 green (0, 128, 0)
% 9 pink (255,105,180)
% 10 orange (255, 165, 0)
% 11 blue (0, 0, 255)

% mni coords: left - right; posterior - anterior; bottom - top
DW_machine;
load('contact_info_newMER.mat');

contact_info(13:18) = [];

%% This chunk is to get the camera position needed
% find right contact indexes
side_info = extractfield(contact_info,'side');
right_idx = find(strcmp('right',side_info));

coords_mess = extractfield(contact_info,'mni_coords');
coords_all = reshape(coords_mess,[3,length(coords_mess)/3])'; % get all contacts' coordinates

% flip right contacts to left so as to plot left and right contacts
% together
coords_all(right_idx,1) = -coords_all(right_idx,1);

figure; plot3(coords_all(:,1),coords_all(:,2),coords_all(:,3),'o')
cam_pos_left = campos;
close all
%%

subject_list = unique(extractfield(contact_info,'subject_id'))';
color_list = {[255, 0, 0]/255;[0, 255, 0]/255;[128, 128, 128]/255;[255, 255, 0]/255;[0, 255, 255]/255;[128, 0, 128]/255;[0,0,128]/255;[0, 128, 0]/255;[255,105,180]/255; [255, 165, 0]/255; [0, 0, 255]/255};

figure (1); hold on;
for i = 1:length(contact_info)
    if strcmp(contact_info(i).side,'left') % left 
    hh(i) = plot3(contact_info(i).mni_coords(1),contact_info(i).mni_coords(2),contact_info(i).mni_coords(3),'o','MarkerSize',12,'MarkerEdgeColor','k',...
            'MarkerFaceColor',cell2mat(color_list(strcmp(contact_info(i).subject_id,subject_list))));
    else % right, flip to left
    hh(i) = plot3( - contact_info(i).mni_coords(1),contact_info(i).mni_coords(2),contact_info(i).mni_coords(3),'^','MarkerSize',12,'MarkerEdgeColor','k',...
            'MarkerFaceColor',cell2mat(color_list(strcmp(contact_info(i).subject_id,subject_list))));        
    end
    campos(cam_pos_left);
end

set(gca, 'XTick',[-20 -15 -10 ]);
ylim([-25 -5]); zlim([-5 15]);
set(gca, 'YTick',[-25:5:-5]); set(gca, 'ZTick', [-5:5:15]);
set(gca, 'YTick',[-25:5:-5]); set(gca, 'ZTick',[-5:5:15]);
grid on;

% legend([hh(1) hh(13) hh(19) hh(25) hh(31) hh(40) hh(48) hh(56) hh(64) hh(72) hh(80) hh(88)],'Subject 1', 'Subject 2','Subject 3', 'Subject 4', 'Subject 5', 'Subject 6', 'Subject 7',...
%     'Subject 8', 'Subject 9', 'Subject 10', 'Subject 11', 'Subject 12');

xlabel('lateral-medial (mm)');
ylabel('anterior-posterior (mm)');
zlabel('ventral-dorsal (mm)');

saveas(gcf, 'Figure2b.fig');