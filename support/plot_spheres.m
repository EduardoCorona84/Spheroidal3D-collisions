function plot_spheres(Ct, r, lims,titleStr)
if ~isa(Ct,'cell')
    Ct = {Ct};
end
if ~exist('r','var')
    r = ones(size(Ct{1},1),1);
end
if ~exist('lims', 'var') || isempty(lims)
    R = computeMaxRad(Ct{1});
    lims = [-2*R 2*R];
end
if ~exist('titleStr','var') || isempty(titleStr)
    titleStr = ['Vesicle with ' num2str(size(Ct{1},1)) ' particles'];
end
%%
fig = uifigure('Name', titleStr, 'Position', [100 100 700 500]);
ax = uiaxes(fig, 'Position', [50 100 600 370]);
xlim(ax,lims);
ylim(ax, lims);
zlim(ax, lims);
xlabel(ax, 'X');
ylabel(ax, 'Y');
zlabel(ax, 'Z');
title(ax, titleStr);
view(ax, 45,45)
grid(ax, 'on')
% Add slider
sld = uislider(fig, ...
    'Position', [100, 50, 500, 3], ...
    'MajorTicks', 10:10:numel(Ct), ...
    "Limits",[1, numel(Ct)+1], ...
    "Value",1);
% Change the plots when the slider moves
sld.ValueChangingFcn = @(src,event) updateSurf(ax, Ct, r, src.Value);
% Initialize plot
updateSurf(ax, Ct, r, 1)

function updateSurf(ax, Ct, r, i)
persistent Sx Sy Sz
if isempty(Sx)
    [Sx, Sy, Sz] = sphere(40);
end
% Clear old plots
cla(ax);

i = round(i);
C = Ct{i};
n3 = size(C,1);
for k=1:n3
    Sx_=r(k)*Sx;
    Sy_=r(k)*Sy;
    Sz_=r(k)*Sz;
    Sc_= C(k,:);surf(ax, Sx_+Sc_(1),Sy_+Sc_(2),Sz_+Sc_(3));
    hold(ax,'on')
end
hold(ax,'off')
drawnow

function maxRad = computeMaxRad(C)
n = size(C,1);
maxRad = max(arrayfun(@(i) norm(C(i,:)), 1:n));