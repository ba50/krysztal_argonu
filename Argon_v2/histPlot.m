close all;
clear;

p = dlmread('p.dat');

figure
hist(p(:,1),50);
figure
hist(p(:,2),50);
figure
hist(p(:,3),50);