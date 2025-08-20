% GPU_Setup_Guide.m - Setup and troubleshooting guide for GPU acceleration
% This script helps diagnose GPU setup and provides troubleshooting tips

clc; clear;

fprintf('=== WMMSE GPU Setup and Diagnostic Guide ===\n\n');

%% GPU Detection and Setup
fprintf('1. Checking GPU Availability...\n');
try
    gpu_count = gpuDeviceCount;
    if gpu_count > 0
        fprintf('   ✓ %d GPU(s) detected\n', gpu_count);
        
        for i = 1:gpu_count
            gpu = gpuDevice(i);
            fprintf('   GPU %d: %s\n', i, gpu.Name);
            fprintf('          Compute Capability: %.1f\n', gpu.ComputeCapability);
            fprintf('          Total Memory: %.1f GB\n', gpu.TotalMemory/1e9);
            fprintf('          Free Memory: %.1f GB\n', gpu.FreeMemory/1e9);
        end
        
        % Test basic GPU operations
        fprintf('\n2. Testing Basic GPU Operations...\n');
        try
            A = randn(1000,'gpuArray');
            B = randn(1000,'gpuArray');
            C = A * B;
            clear A B C;
            fprintf('   ✓ Basic GPU matrix operations working\n');
        catch ME
            fprintf('   ✗ GPU operations failed: %s\n', ME.message);
        end
        
        % Test complex GPU operations
        fprintf('\n3. Testing Complex GPU Operations...\n');
        try
            A = randn(100) + 1j*randn(100,'gpuArray');
            B = inv(A);
            C = det(A);
            clear A B C;
            fprintf('   ✓ Complex GPU operations working\n');
        catch ME
            fprintf('   ✗ Complex GPU operations failed: %s\n', ME.message);
        end
        
    else
        fprintf('   ✗ No GPUs detected\n');
    end
    
catch ME
    fprintf('   ✗ GPU detection failed: %s\n', ME.message);
    fprintf('   This likely means Parallel Computing Toolbox is not installed\n');
end

%% Parallel Computing Toolbox Check
fprintf('\n4. Checking Parallel Computing Toolbox...\n');
try
    license_available = license('test', 'Distrib_Computing_Toolbox');
    if license_available
        fprintf('   ✓ Parallel Computing Toolbox license available\n');
        
        % Check version
        v = ver('parallel');
        if ~isempty(v)
            fprintf('   Version: %s\n', v.Version);
        end
    else
        fprintf('   ✗ Parallel Computing Toolbox license not available\n');
    end
catch ME
    fprintf('   ✗ Parallel Computing Toolbox check failed: %s\n', ME.message);
end

%% Memory Requirements Estimation
fprintf('\n5. Memory Requirements for WMMSE Algorithm...\n');
T = 128; R = 4; I = 16; d = 4;

% Estimate memory usage
complex_double_size = 16; % bytes (8 for real + 8 for imaginary)
H_memory = R * T * I * complex_double_size;
V_memory = T * d * I * complex_double_size;
U_memory = R * d * I * complex_double_size;
W_memory = d * d * I * complex_double_size;

total_memory = (H_memory + V_memory + U_memory + W_memory) * 2; % Factor of 2 for intermediate calculations
total_memory_GB = total_memory / 1e9;

fprintf('   Estimated GPU memory needed: %.2f GB\n', total_memory_GB);
fprintf('   For T=%d, R=%d, I=%d, d=%d\n', T, R, I, d);

if gpu_count > 0
    gpu = gpuDevice(1);
    if gpu.FreeMemory/1e9 > total_memory_GB
        fprintf('   ✓ Sufficient GPU memory available\n');
    else
        fprintf('   ⚠ May not have sufficient GPU memory\n');
        fprintf('   Consider reducing problem size or using CPU version\n');
    end
end

%% Performance Recommendations
fprintf('\n6. Performance Recommendations...\n');
fprintf('   • For T < 64, R < 8, I < 8: CPU may be faster due to overhead\n');
fprintf('   • For T ≥ 128, R ≥ 4, I ≥ 16: GPU acceleration recommended\n');
fprintf('   • For very large problems (T > 256): Monitor GPU memory usage\n');
fprintf('   • Use WMMSE.m for automatic GPU/CPU selection\n');
fprintf('   • Use WMMSE_GPU.m only when GPU is confirmed working\n');

%% Troubleshooting
fprintf('\n7. Common Issues and Solutions...\n');
fprintf('   Issue: "Undefined function gpuArray"\n');
fprintf('   Solution: Install Parallel Computing Toolbox\n\n');
fprintf('   Issue: "No supported GPU found"\n');
fprintf('   Solution: Update GPU drivers, check CUDA compatibility\n\n');
fprintf('   Issue: "Out of memory on device"\n');
fprintf('   Solution: Reduce problem size (T, R, I, d) or use CPU version\n\n');
fprintf('   Issue: "GPU operations slower than CPU"\n');
fprintf('   Solution: Problem size may be too small for GPU benefit\n\n');

fprintf('=== Diagnostic Complete ===\n');
fprintf('If all checks pass, you can use GPU acceleration with confidence!\n');