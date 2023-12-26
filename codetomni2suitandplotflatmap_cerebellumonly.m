%% codetomni2suitandplotflatmap.m contains the following:
% clc; clear all;
addpath /n04dat04/ypwang/software/spm12/
addpath /n04dat04/ypwang/software/spm12/compat
addpath /n04dat04/ypwang/software/spm12/toolbox/DARTEL
addpath /n04dat04/ypwang/software/spm12/toolbox/suit % suit-3.5
% spm fmri

suit_mni2suit('image_cerebellumonly_nifti.nii');

figure
set(gca, 'Visible', 'off')
Data = suit_map2surf('Wm2s_image_cerebellumonly_nifti.nii','space','SUIT', 'stats',@mode);
suit_plotflatmap(Data, 'cmap', colormap(jet(256)), 'cscale', [min(Data), max(Data)]);
caxis([0,1]);
colorbar('horiz','position',[0.4 0.14 0.23 0.03]);
savefig('figure')

fig = openfig('figure.fig');
filename = 'image_cerebellumonly.png';
% print(fig ,'-dpng','-r300',filename)
saveas(fig, filename)
clearvars
