function [p,ns,centers,u0s,as,oblates,Rmats] = readSystem(fileName)
% Read stored set up info on a system of spheroid suspension from file.
format shortg
fid = fopen(fileName,'r');

% p
line1=fgetl(fid);
line1split=split(line1,'=');
p=str2double(line1split(end));

% ns total
line1=fgetl(fid);
line1split=split(line1,'=');
ns=str2double(line1split(end));

% centers
fgetl(fid);
centers=zeros(ns,3);
for ii=1:ns
    posii=fgetl(fid);
    posiisplit=split(posii,', ');
    centers(ii,:)=str2num(char(posiisplit))';
end
% fgetl(fid);
if fgetl(fid)~=-1
    fprintf("\n PAUSE: not done with centers yet!")
end

% u0's
fgetl(fid);
u0s=zeros(1,ns);
for ii=1:ns
    u0ii=fgetl(fid);
    u0s(ii)=str2num(u0ii);
end
% fgetl(fid);
if fgetl(fid)~=-1
    fprintf("\n PAUSE: not done with u0's yet!")
end

% a's
fgetl(fid);
as=zeros(1,ns);
for ii=1:ns
    aii=fgetl(fid);
    as(ii)=str2num(aii);
end
% fgetl(fid);
if fgetl(fid)~=-1
    fprintf("\n PAUSE: not done with a's yet!")
end

% oblate boolean list
fgetl(fid);
oblates=zeros(1,ns);
for ii=1:ns
    obii=fgetl(fid);
    oblates(ii)=str2double(obii);
end
% fgetl(fid);
if fgetl(fid)~=-1
    fprintf("\n PAUSE: not done with oblates yet!")
end

% Rotation matrices
fgetl(fid);
Rmats=zeros(3,3,ns);
for ii=1:ns
    for jj=1:3
        R1=fgetl(fid);
        R1split=split(R1,', ');
        Rmats(jj,:,ii)=str2num(char(R1split))';
    end
    fgetl(fid);
end

% check file is now empty
if fgetl(fid)~=-1
    fprintf("\n File is not empty when readfile ends. Check readSystem.m.")
end

fclose(fid);

end

