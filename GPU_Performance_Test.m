% GPU_Performance_Test.m - Performance comparison between CPU and GPU versions
% This script compares the performance of CPU vs GPU implementations

clc; clear;

fprintf('=== WMMSE GPU Performance Comparison ===\n\n');

% Check GPU availability
useGPU = false;
if gpuDeviceCount > 0
    gpu = gpuDevice();
    fprintf('GPU detected: %s\n', gpu.Name);
    useGPU = true;
else
    fprintf('No GPU detected. Only CPU version will be tested.\n');
end

% Test parameters
K = 1; 
T = 128; 
R = 4; 
epsilon = 1e-3; 
sigma2 = 1; 
snr = 10; 
P = db2pow(snr)*sigma2; 
I = 16; 
alpha1 = ones(I,K); 
d = 4; 
max_iter = 100;

% Number of test runs for averaging
num_runs = 3;

fprintf('\nTest Parameters:\n');
fprintf('Transmit antennas (T): %d\n', T);
fprintf('Receive antennas (R): %d\n', R);
fprintf('Users (I): %d\n', I);
fprintf('Data streams per user (d): %d\n', d);
fprintf('Number of test runs: %d\n\n', num_runs);

%% CPU Performance Test
fprintf('Running CPU performance test...\n');
cpu_times = zeros(num_runs, 1);
cpu_rates = cell(num_runs, 1);

for run = 1:num_runs
    fprintf('  CPU Run %d/%d...', run, num_runs);
    
    % Initialize matrices
    H = zeros(R,T,I);
    for i=1:I
        H(:,:,i) = sqrt(1/2)*(randn(R,T)+1i*randn(R,T));
    end
    
    U = randn(R,d,I) + 1j*randn(R,d,I);
    W = zeros(d,d,I);
    for i=1:I
        W(:,:,i) = eye(d,d);
    end
    
    V = zeros(T,d,I);
    for i=1:I
        v = sqrt(1/2)*(randn(T,d)+1i*randn(T,d));
        V(:,:,i) = sqrt(P/(I*trace(v*v')))*v;
    end
    
    % Run algorithm
    tic;
    rate_old = sum_rate(H,V,sigma2,R,I,alpha1,false);
    rate = rate_old;
    
    iter1 = 1;
    while(1)
        U = find_U(H,V,sigma2, P, R,I,d,false); 
        W = find_W(U,H,V, R , I,d,false); 
        V = find_V(alpha1,H,sigma2,U,W,T, R, I,d ,P,false); 
        rate_new = sum_rate(H,V,sigma2,R,I,alpha1,false);
        rate = [rate rate_new];
        iter1 = iter1 + 1;
        if abs(rate_new-rate_old) / rate_old < epsilon || iter1 > max_iter
            break;
        end
        rate_old = rate_new;
    end
    
    cpu_times(run) = toc;
    cpu_rates{run} = rate;
    fprintf(' %.2f seconds\n', cpu_times(run));
end

%% GPU Performance Test
if useGPU
    fprintf('\nRunning GPU performance test...\n');
    gpu_times = zeros(num_runs, 1);
    gpu_rates = cell(num_runs, 1);
    
    for run = 1:num_runs
        fprintf('  GPU Run %d/%d...', run, num_runs);
        
        % Initialize matrices on GPU
        H = zeros(R,T,I,'gpuArray');
        for i=1:I
            H(:,:,i) = sqrt(1/2)*(randn(R,T,'gpuArray')+1i*randn(R,T,'gpuArray'));
        end
        
        alpha1_gpu = gpuArray(alpha1);
        U = randn(R,d,I,'gpuArray') + 1j*randn(R,d,I,'gpuArray');
        W = zeros(d,d,I,'gpuArray');
        for i=1:I
            W(:,:,i) = eye(d,d,'gpuArray');
        end
        
        V = zeros(T,d,I,'gpuArray');
        for i=1:I
            v = sqrt(1/2)*(randn(T,d,'gpuArray')+1i*randn(T,d,'gpuArray'));
            V(:,:,i) = sqrt(P/(I*trace(v*v')))*v;
        end
        
        % Run algorithm
        tic;
        rate_old = sum_rate(H,V,sigma2,R,I,alpha1_gpu,true);
        rate = gather(rate_old);
        
        iter1 = 1;
        while(1)
            U = find_U(H,V,sigma2, P, R,I,d,true); 
            W = find_W(U,H,V, R , I,d,true); 
            V = find_V(alpha1_gpu,H,sigma2,U,W,T, R, I,d ,P,true); 
            rate_new = sum_rate(H,V,sigma2,R,I,alpha1_gpu,true);
            rate = [rate gather(rate_new)];
            iter1 = iter1 + 1;
            if abs(gather(rate_new)-gather(rate_old)) / gather(rate_old) < epsilon || iter1 > max_iter
                break;
            end
            rate_old = rate_new;
        end
        
        gpu_times(run) = toc;
        gpu_rates{run} = rate;
        fprintf(' %.2f seconds\n', gpu_times(run));
        
        % Clear GPU memory for next run
        clear H U V W alpha1_gpu;
    end
    
    % Reset GPU
    gpuDevice([]);
end

%% Results Analysis
fprintf('\n=== Performance Results ===\n');
fprintf('CPU Performance:\n');
fprintf('  Average time: %.3f ± %.3f seconds\n', mean(cpu_times), std(cpu_times));
fprintf('  Min time: %.3f seconds\n', min(cpu_times));
fprintf('  Max time: %.3f seconds\n', max(cpu_times));

if useGPU
    fprintf('\nGPU Performance:\n');
    fprintf('  Average time: %.3f ± %.3f seconds\n', mean(gpu_times), std(gpu_times));
    fprintf('  Min time: %.3f seconds\n', min(gpu_times));
    fprintf('  Max time: %.3f seconds\n', max(gpu_times));
    
    speedup = mean(cpu_times) / mean(gpu_times);
    fprintf('\nSpeedup: %.2fx faster on GPU\n', speedup);
    
    % Verify convergence consistency
    cpu_final_rate = cpu_rates{1}(end);
    gpu_final_rate = gpu_rates{1}(end);
    rate_diff = abs(cpu_final_rate - gpu_final_rate) / cpu_final_rate * 100;
    fprintf('Convergence difference: %.4f%% (should be < 0.1%%)\n', rate_diff);
    
    % Plot comparison
    figure;
    subplot(1,2,1);
    plot(0:length(cpu_rates{1})-1, cpu_rates{1}, 'b-o', 'LineWidth', 1.5);
    title('CPU Performance');
    xlabel('Iterations');
    ylabel('Sum Rate');
    grid on;
    
    subplot(1,2,2);
    plot(0:length(gpu_rates{1})-1, gpu_rates{1}, 'r-o', 'LineWidth', 1.5);
    title('GPU Performance');
    xlabel('Iterations');
    ylabel('Sum Rate');
    grid on;
    
    sgtitle(sprintf('WMMSE Algorithm: CPU vs GPU (%.2fx speedup)', speedup));
end

fprintf('\nPerformance test completed.\n');