function [ theta_matrix_periodogram, Averaged_periodogram ] = averaging( PSD, interval )
%UNTITLED2 Summary of this function goes here
%  Detailed explanation goes here
Averaged_periodogram=[];
first = 1;
last=floor(length(PSD)/(interval+1));
for i = 1:interval+1
Averaged_periodogram = [Averaged_periodogram, mean(PSD(first:last))]; 
first = floor(first+length(PSD)/(1+interval));
last = floor(last+length(PSD)/(1+interval));
end
theta_matrix_periodogram=0:(1/(interval)):1;
end