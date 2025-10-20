function plot_spheres(Ct, lims)
if ~exist('lims', 'var') || isempty(lims)
    lims = [-3 3 -3 3 -1 1];
end
%%
set(0,'DefaultFigureWindowStyle','docked')
fig = figure;
axis(lims)
grid on
view(0,90)

n3 = size(Ct{1},1);
nt = numel(Ct);
if ~exist('r','var')
    r = ones(n3,1);
end
% F(nt) = struct('cdata',[],'colormap',[]);
hold on 
hndls = updateSurf(cell(n3,1), Ct, r, 1);
hold off


title('Simple Amphi');
xlabel('X'); ylabel('Y'); zlabel('Z');

% Add slider
hSlider = uicontrol('Style', 'slider',...
    'Min', 1, 'Max', nt, 'Value', 1,...
    'Units', 'normalized',...
    'Position', [0.2 0.02 0.6 0.05],...
    'Callback', @(src, event) updateSurf(hndls, Ct, r, src.Value));

% Add a text label to display slider value
uicontrol('Style', 'text', 'Units', 'normalized', ...
    'Position', [0.82 0.02 0.1 0.05], ...
    'String', 'ixTime = ');
end

function hndls = updateSurf(hndls, Ct, r, i)
    i = round(i);
    n3 = numel(hndls);
    [Sx, Sy, Sz]=sphere(40);
    C = Ct{i};
    for k=1:n3
        Sx_=r(k)*Sx;
        Sy_=r(k)*Sy;
        Sz_=r(k)*Sz;
        Sc_= C(k,:);
        if i ==1 
            hndls{k} = surf(Sx_+Sc_(1),Sy_+Sc_(2),Sz_+Sc_(3));
        else 
            hndls{k}.XData=Sx_+Sc_(1);
            hndls{k}.YData=Sy_+Sc_(2);
            hndls{k}.ZData=Sz_+Sc_(3);
        end
    end
end

% 
% %%
% fig = figure();
% axis([-3 3 -3 3 -1 1])
% grid on
% view(0,90)
% movie(fig,F,2)

