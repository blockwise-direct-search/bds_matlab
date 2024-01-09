function [fv,ifail,icount]=imfil_matcutest(x)
% f_matcutest
% Simple example of using imfil.m
%
global fun_imfil
x = double(x(:));
fv = fun_imfil(x);
%
% This function never fails to return a value 
%
ifail=0;
%
% and every call to the function has the same cost.
%
icount=1;
