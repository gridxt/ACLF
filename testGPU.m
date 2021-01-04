clc;
clear all;

load('system.mat');

% tic;
% %A = dissect(J);
% x = J\F;
% cpu = toc;

% f = @()gpuArray(gather(J)\gather(F));
% gpu = gputimeit(f);

tic;
J_gpu = gpuArray(J);
F_gpu = gpuArray(F);
Mat_gpu = J_gpu*F_gpu;
%mm = gather(Mat_gpu);
%m = chol
Gtime = toc;

tic;
Mat_cpu = J*F;
Ctime = toc;