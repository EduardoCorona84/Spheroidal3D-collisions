% torus surface with facelighting, borrowed from Alex's BIE3D
%
% 05/28/24

clear
addpath('./utils/')

type = 'torus';
so.a = 1.0;   % major radius
b = 0.5;    % mean poloidal radius
wc = 0.1;  % surf modulation ampl
wm = 3;   % # wobbles in minor, poloidal (note swapped from 2013)
wn = 5;   % # wobbles in toroidal, major
f = @(t,p) b + wc*cos(wm*t+wn*p);         % wobble func, then its partials...
ft = @(t,p) -wc*wm*sin(wm*t+wn*p); fp = @(t,p) -wc*wn*sin(wm*t+wn*p);
so.b = {f,ft,fp};     % pass in instead of b param
so.p = 6;
o = [];
[pan N] = create_panels(type,so,o);
[x,~,w] = getallnodes(pan);

figure(1),clf
h=showsurffunc(pan, w, struct('nofig',true));
% cmp = getPyPlot_cMap('rainbow', [], [], '"/Users/hzhu/.pyenv/versions/3.9.13/bin/python"');
cmp = getPyPlot_cMap('rainbow', [], [], '"python3"');
colormap(cmp)
lightangle(45,0); % comment out if want to turn off lighting
% plot3 if want to plot 3d line