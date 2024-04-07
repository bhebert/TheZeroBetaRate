function [ZinputR] = RandomFeatures(Zinput,c,gms)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

rng(12302018);

L = size(Zinput,1);
T = size(Zinput,2);

P = floor(c/2 * T)+1;

ws = mvnrnd(zeros(L,1),eye(L),P);

gammas = randsample(gms,T,true);

vals = gammas.*(ws * Zinput);

cosvec = cos(vals);
sinvec = sin(vals);

ZinputR = vertcat(cosvec,sinvec) / sqrt(P);

end