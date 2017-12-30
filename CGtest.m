% Strong isospin test example
%
% mikael.mieskolainen@cern.ch, 2017
clear;
close all;
addpath ./src

% SU(2) isospin states
proton  = [1/2 1/2];
neutron = [1/2 -1/2];

piplus  = [1  1];
pizero  = [1  0];
piminus = [1 -1];

% Calculate Glesbch-Gordan ratios (based on Wigner-Eckart theorem)
CGR = @(p1,p2,k) clebschgordan(p1(1),p2(1), p1(2),p2(2), k(1),k(2));

ratio = CGR(neutron,piplus, proton) / CGR(proton,pizero, proton);
abs(ratio)^2
