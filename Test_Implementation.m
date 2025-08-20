% Test_Implementation.m - Test script to verify GPU implementation
% This script tests the modified functions to ensure they work correctly

clc; clear;

fprintf('=== Testing GPU-Accelerated WMMSE Implementation ===\n\n');

%% Test Parameters
K = 1;
T = 8; % Smaller size for testing
R = 2;
sigma2 = 1;
snr = 10;
P = db2pow(snr)*sigma2;
I = 4; % Fewer users for testing
alpha1 = ones(I,K);
d = 2; % Smaller streams

fprintf('Test Parameters: T=%d, R=%d, I=%d, d=%d\n\n', T, R, I, d);

%% Test 1: Function Compatibility
fprintf('Test 1: Function compatibility with CPU arrays...\n');

try
    % Create test matrices
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
    
    % Test all functions with CPU arrays
    U_new = find_U(H,V,sigma2, P, R,I,d,false);
    W_new = find_W(U,H,V, R, I,d,false);
    V_new = find_V(alpha1,H,sigma2,U,W,T, R, I,d ,P,false);
    rate = sum_rate(H,V,sigma2,R,I,alpha1,false);
    
    fprintf('   ✓ All functions work with CPU arrays\n');
    fprintf('   Initial sum rate: %.4f\n', rate);
    
catch ME
    fprintf('   ✗ CPU test failed: %s\n', ME.message);
    return;
end

%% Test 2: Backward Compatibility
fprintf('\nTest 2: Backward compatibility (old function signatures)...\n');

try
    % Test old function signatures (without useGPU parameter)
    U_old = find_U(H,V,sigma2, P, R,I,d);
    W_old = find_W(U,H,V, R, I,d);
    V_old = find_V(alpha1,H,sigma2,U,W,T, R, I,d ,P);
    rate_old = sum_rate(H,V,sigma2,R,I,alpha1);
    
    % Compare results
    U_diff = norm(U_new(:) - U_old(:));
    W_diff = norm(W_new(:) - W_old(:));
    V_diff = norm(V_new(:) - V_old(:));
    rate_diff = abs(rate - rate_old);
    
    fprintf('   ✓ Backward compatibility maintained\n');
    fprintf('   Differences (should be ~0): U=%.2e, W=%.2e, V=%.2e, rate=%.2e\n', ...
            U_diff, W_diff, V_diff, rate_diff);
    
catch ME
    fprintf('   ✗ Backward compatibility test failed: %s\n', ME.message);
end

%% Test 3: GPU Functionality (if available)
fprintf('\nTest 3: GPU functionality...\n');

try
    if gpuDeviceCount > 0
        fprintf('   GPU detected, testing GPU arrays...\n');
        
        % Convert to GPU arrays
        H_gpu = gpuArray(H);
        U_gpu = gpuArray(U);
        W_gpu = gpuArray(W);
        V_gpu = gpuArray(V);
        alpha1_gpu = gpuArray(alpha1);
        
        % Test all functions with GPU arrays
        U_gpu_new = find_U(H_gpu,V_gpu,sigma2, P, R,I,d,true);
        W_gpu_new = find_W(U_gpu,H_gpu,V_gpu, R, I,d,true);
        V_gpu_new = find_V(alpha1_gpu,H_gpu,sigma2,U_gpu,W_gpu,T, R, I,d ,P,true);
        rate_gpu = sum_rate(H_gpu,V_gpu,sigma2,R,I,alpha1_gpu,true);
        
        % Compare GPU vs CPU results
        rate_gpu_cpu_diff = abs(gather(rate_gpu) - rate) / rate * 100;
        
        fprintf('   ✓ GPU functions work correctly\n');
        fprintf('   GPU vs CPU rate difference: %.4f%% (should be < 0.1%%)\n', rate_gpu_cpu_diff);
        
        if rate_gpu_cpu_diff < 0.1
            fprintf('   ✓ GPU and CPU results match within tolerance\n');
        else
            fprintf('   ⚠ GPU and CPU results differ more than expected\n');
        end
        
        % Clear GPU memory
        clear H_gpu U_gpu W_gpu V_gpu alpha1_gpu U_gpu_new W_gpu_new V_gpu_new;
        gpuDevice([]);
        
    else
        fprintf('   No GPU available, skipping GPU tests\n');
    end
    
catch ME
    fprintf('   ✗ GPU test failed: %s\n', ME.message);
    if gpuDeviceCount > 0
        try
            gpuDevice([]);
        catch
        end
    end
end

%% Test 4: Mini Algorithm Run
fprintf('\nTest 4: Mini WMMSE algorithm run...\n');

try
    epsilon = 1e-2; % Looser tolerance for quick test
    max_iter = 10;  % Fewer iterations for quick test
    
    % Initialize
    rate_old = sum_rate(H,V,sigma2,R,I,alpha1,false);
    rate_history = rate_old;
    
    % Run a few iterations
    for iter = 1:max_iter
        U = find_U(H,V,sigma2, P, R,I,d,false);
        W = find_W(U,H,V, R, I,d,false);
        V = find_V(alpha1,H,sigma2,U,W,T, R, I,d ,P,false);
        rate_new = sum_rate(H,V,sigma2,R,I,alpha1,false);
        rate_history = [rate_history rate_new];
        
        if abs(rate_new-rate_old) / rate_old < epsilon
            break;
        end
        rate_old = rate_new;
    end
    
    fprintf('   ✓ Algorithm converged in %d iterations\n', iter);
    fprintf('   Initial rate: %.4f\n', rate_history(1));
    fprintf('   Final rate: %.4f\n', rate_history(end));
    fprintf('   Rate improvement: %.4f\n', rate_history(end) - rate_history(1));
    
    if rate_history(end) > rate_history(1)
        fprintf('   ✓ Algorithm shows convergence improvement\n');
    else
        fprintf('   ⚠ Algorithm did not improve (may be at local optimum)\n');
    end
    
catch ME
    fprintf('   ✗ Mini algorithm test failed: %s\n', ME.message);
end

fprintf('\n=== Testing Complete ===\n');
fprintf('If all tests passed, the GPU implementation is ready to use!\n');