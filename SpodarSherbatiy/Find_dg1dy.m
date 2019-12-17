function [ dg1dy ] = Find_dg1dy( y )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

dg1dy=[(0*ones(1,size(y,2)))', 2*y(2,:)'-1];

end

