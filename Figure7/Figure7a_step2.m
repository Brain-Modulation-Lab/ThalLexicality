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

% open pre-created figure and hold on
open('Figure7a_step1.fig');
hold on;

% display entire thalamus
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

axis equal

xlabel('');
ylabel('');
zlabel('');
colorbar off

grid off; axis off;

load('Figure7a_campos.mat', 'cam_used');
campos(cam_used);

saveas(gcf, 'Figure7a_step2.fig');