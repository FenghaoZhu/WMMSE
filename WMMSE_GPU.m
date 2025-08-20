% WMMSE_GPU.m - GPU-optimized version of WMMSE algorithm
% This version is specifically optimized for GPU acceleration and requires
% MATLAB's Parallel Computing Toolbox

clc;clear;

% Force GPU usage (will error if not available)
if gpuDeviceCount == 0
    error('No GPU detected. Please use WMMSE.m for CPU version or ensure GPU is available.');
end

gpu = gpuDevice();
fprintf('Using GPU: %s (Compute Capability: %.1f)\n', gpu.Name, gpu.ComputeCapability);
fprintf('GPU Memory: %.1f GB\n', gpu.TotalMemory/1e9);

K = 1; % 基站个数，此版本固定为1
T = 128; % 发射天线个数
R = 4; % 接收天线个数
epsilon = 1e-3; % 收敛条件
sigma2 = 1; % 噪声功率
snr = 10; % 信噪比
P = db2pow(snr)*sigma2; % 发射功率

I = 16; % 用户个数
alpha1 = gpuArray(ones(I,K)); % 权重系数，直接创建为GPU数组

d = 4; % 假设每个用户都有d条路独立的数据流

max_iter = 100;

% Initialize all matrices directly on GPU
fprintf('Initializing matrices on GPU...\n');
tic;

% 初始化信道向量直接在GPU上
H = zeros(R,T,I,'gpuArray'); 
for i=1:I
    H(: , :, i) = sqrt(1/2)*(randn(R,T,'gpuArray')+1i*randn(R,T,'gpuArray')); 
end

rate = []; % 初始化一个空向量记录rate

% 初始化W和U矩阵直接在GPU上
U = randn(R,d,I,'gpuArray') + 1j*randn(R,d,I,'gpuArray');
W = zeros(d,d,I,'gpuArray');
for i=1:I
    W(:,:,i) = eye(d,d,'gpuArray');
end

% 初始化波束赋形矩阵直接在GPU上
V = zeros(T,d,I,'gpuArray'); 
for i=1:I
    v = sqrt(1/2)*(randn(T,d,'gpuArray')+1i*randn(T,d,'gpuArray'));
    V(:,:, i) = sqrt(P/(I*trace(v*v')))*v;
end 

fprintf('GPU initialization completed.\n');

% 求初始化发射波束V后求系统和速率
rate_old = sum_rate(H,V,sigma2,R,I,alpha1,true);
rate = [rate gather(rate_old)];

fprintf('Starting WMMSE iterations on GPU...\n');
iter1 = 1;
while(1)
    U = find_U(H,V,sigma2, P, R,I,d,true); 
    W = find_W(U,H,V, R , I,d,true); 
    V = find_V(alpha1,H,sigma2,U,W,T, R, I,d ,P,true); 
    rate_new = sum_rate(H,V,sigma2,R,I,alpha1,true);
    rate = [rate gather(rate_new)];
    iter1 = iter1 + 1;
    if abs(gather(rate_new)-gather(rate_old)) / gather(rate_old) < epsilon || iter1 > max_iter
        break;
    end
    rate_old = rate_new;
end

toc;
fprintf('WMMSE algorithm completed in %d iterations.\n', iter1-1);

% Plot results
plot(0:iter1-1,rate,'r-o')
grid on
xlabel('Iterations')
ylabel('Sum rate (bits per channel use)')
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1)
title('WMMSE (GPU Accelerated), K=1, T=128, R=4, d=4, 10dB, \epsilon=1e-3','Interpreter','tex')

% Display GPU memory usage
fprintf('GPU Memory Used: %.1f GB\n', (gpu.TotalMemory - gpu.FreeMemory)/1e9);

% Clean up GPU memory
clear H U V W alpha1;
gpuDevice([]); % Reset GPU
fprintf('GPU memory cleared.\n');