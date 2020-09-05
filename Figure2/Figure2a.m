
% load in thalamus atlas
load('subcort_Ewert_v2.1.7_Medium.mat');
% rename this one as atlases_thal
atlases_thal = atlases;

% load in sub-structure atlas
load('subcort_Ewert_v2.1.7.mat');

% atlas index of VA, VLa, VLp
group{1} = [36 41 67 68];
group{2} = [65 66 99 100 101];
group{3} = [64 69 70 71 72 82 89 90 92 98];

figure; hold on;
thal = patch('vertices',atlases_thal.fv{16,2}.vertices,'faces',atlases_thal.fv{16,2}.faces,...
'facecolor',atlases_thal.colormap(16,:),'edgecolor','none',...
'facelighting', 'gouraud', 'specularstrength', .50, 'facealpha', 0.1);

% display VLp
for i = 1:length(group{3})
    VLp{i} = patch('vertices',atlases.fv{group{3}(i),2}.vertices,'faces',atlases.fv{group{3}(i),2}.faces,...
'facecolor',[0.125,0.324,0.830],'edgecolor','none',...
'facelighting', 'gouraud', 'specularstrength', .50, 'facealpha', 0.2);
end

% display VA
for i = 1:length(group{1})
    VA{i} = patch('vertices',atlases.fv{group{1}(i),2}.vertices,'faces',atlases.fv{group{1}(i),2}.faces,...
'facecolor',[0.433,0.738,0.621],'edgecolor','none',...
'facelighting', 'gouraud', 'specularstrength', .50, 'facealpha', 0.2);
end

% display VLa
for i = 1:length(group{2})
    VLa{i} = patch('vertices',atlases.fv{group{2}(i),2}.vertices,'faces',atlases.fv{group{2}(i),2}.faces,...
'facecolor',[0.433,0.738,0.621],'edgecolor','none',...
'facelighting', 'gouraud', 'specularstrength', .50, 'facealpha', 0.2);
end

vim = patch('vertices',atlases_thal.fv{18,2}.vertices,'faces',atlases_thal.fv{18,2}.faces,...
'facecolor',[255, 99,71]/255,'edgecolor','none',...
'facelighting', 'gouraud', 'specularstrength', .50, 'facealpha', 0.5);
axis off; axis equal;

% load in contact_info_newMER
load('contact_info_newMER.mat');

% find right contact indexes
side_info = extractfield(contact_info,'side');
right_idx = find(strcmp('right',side_info));

coords_mess = extractfield(contact_info,'mni_coords');
coords_all = reshape(coords_mess,[3,length(coords_mess)/3])'; % get all contacts' coordinates

% flip right contacts to left so as to plot left and right contacts
% together
coords_all(right_idx,1) = -coords_all(right_idx,1);

coords_left = coords_all(setdiff(1:length(coords_all), right_idx),:);
coords_right = coords_all(right_idx,:);

hold on;
contacts_new_left = plot3(coords_left(:,1), coords_left(:,2), coords_left(:,3),'ko', 'Color', [0,0,0], 'markersize', 6,'MarkerFaceColor','k');
contacts_new_right = plot3(coords_right(:,1), coords_right(:,2), coords_right(:,3),'k^', 'Color', [0,0,0], 'markersize', 6,'MarkerFaceColor','k');

% plot mri

nii = load_nii('MNI_ICBM_2009b_NLIN_ASYM/t2.nii');
nii_reorder = permute(nii.img,[2 1 3]);

load mri
hold on;
colormap(map)

[X,Y,Z] = meshgrid(-98:0.5:98.5, -134:0.5:98.5, -72:0.5:116.5);

h = slice(X,Y,Z, nii_reorder,5,[],[])
set(h, 'EdgeColor','none')
% xlabel('x')
% ylabel('y')
% zlabel('z')

% adjust the camera view

%save this figure
saveas(gcf, 'Figure2a.fig');