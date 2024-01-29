function [DL_UL_indicators] = TTI_indicator(  Frame_structure)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
   
if strcmp(Frame_structure,'TDD_1')
    DL_UL_indicators=[1,1,0,0,1,1,1,0,0,1];% DL 1, UL 0.
end

if strcmp(Frame_structure,'FDD')
    DL_UL_indicators=[1,0,1,1,1,1,0,1,1,1];% DL 1, UL 0.
end


