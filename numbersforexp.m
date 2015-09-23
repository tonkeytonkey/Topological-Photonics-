close all; 
clear all; 
clc; 


%we need to check that calculations agree with the paper 

%old rods

%numbers from data sheet match nature paper, other batch of rods is just
%shorter and should be equally useful 
eps=14.63; 
Ms=1912*10^3/(4*pi);
loss=13*10^3/(4*pi); %convert from Oe to A/m


%{
%new rods
eps=14.87; 
Ms=1827;
loss=10; 
%}