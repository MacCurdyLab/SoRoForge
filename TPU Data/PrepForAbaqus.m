%% compute nominal stress strain vector

clear; clc; close all

A = 6.13*1.25;
L = 33;

nData = 14;

load allRaw.mat

x = allRaw{nData,1}(:,2);
f = allRaw{nData,1}(:,3);

s = x/L;
S = f/A;

maxs = 2;

Si = find(s>maxs,1);

nSkip = 25;

Data = [S(1:nSkip:Si) s(1:nSkip:Si)];