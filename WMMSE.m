clc;clear;
K = 1; % 基站个数，此版本固定为1
T = 128; % 发射天线个数
R = 4; % 接收天线个数
epsilon = 1e-3; % 收敛条件
sigma2 = 1; % 噪声功率
snr = 10; % 信噪比
P = db2pow(snr)*sigma2; % 发射功率

I = 16; % 用户个数
alpha1 = ones(I,K); % 权重系数，都假设相同

d = 4; % 假设每个用户都有d条路独立的数据流

max_iter = 100;
tic;
% 初始化信道向量
H = zeros(R,T,I); % 信道系数 用户数量*每个用户天线数量*基站天线数量
for i=1:I
    H(: , :, i)=sqrt(1/2)*(randn(R,T)+1i*randn(R,T)); % 圆复高斯信道
end

rate = []; % 初始化一个空向量记录rate

% 初始化W和U矩阵
U =randn(R,d,I) + 1j*randn(R,d,I);
W = zeros(d,d,I);
for i=1:I
    W(:,:,i)=eye(d,d);
end

% 初始化波束赋形矩阵
V = zeros(T,d,I); % 每个基站的波束赋形矩阵
for i=1:I
        v = sqrt(1/2)*(randn(T,d)+1i*randn(T,d));
        V(:,:, i) = sqrt(P/(I*trace(v*v')))*v;
end 

% 求初始化发射波束V后求系统和速率
rate_old = sum_rate(H,V,sigma2,R,I,alpha1);
rate = [rate rate_old];


iter1 = 1;
while(1)
    U = find_U(H,V,sigma2, P, R,I,d); 
    W = find_W(U,H,V, R , I,d); 
    V = find_V(alpha1,H,sigma2,U,W,T, R, I,d ,P); 
    rate_new = sum_rate(H,V,sigma2,R,I,alpha1); % 计算和速率
    rate = [rate rate_new];
    iter1 = iter1 + 1;
    if abs(rate_new-rate_old) / rate_old < epsilon || iter1 > max_iter
        break;
    end
    rate_old = rate_new;
end
toc;
plot(0:iter1-1,rate,'r-o')
grid on
xlabel('Iterations')
ylabel('Sum rate (bits per channel use)')
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1)
title('WMMSE, K=1, T=128, R=4, d=4, 10dB, \epsilon=1e-3','Interpreter','tex')