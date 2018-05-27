clear all; close all; clc;

tic
y = zeros(1,50000);
for i = 1:50000
   
    y(i) = 5*i;
    
end
toc