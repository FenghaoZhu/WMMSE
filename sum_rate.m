function system_rate = sum_rate(H,V,sigma2,R,I,alpha1,useGPU)
% GPU-accelerated version of sum_rate function
% useGPU: boolean flag indicating whether to use GPU acceleration

    if nargin < 7
        useGPU = false; % Default to CPU if not specified
    end
    
    if useGPU
        % GPU-accelerated implementation
        rate = zeros(I,1,'like',H);
        for i=1:I
            denominator = zeros(R,R,'like',H);
            H_i = H(:,:,i);
            
            % Vectorized computation for denominator
            for l=1:I
                V_l = V(:,:,l);
                denominator = denominator + H_i * V_l * V_l' * H_i';
            end
            
            V_i = V(:,:,i);
            numerator = H_i * V_i * V_i' * H_i';
            denominator = denominator - numerator + sigma2*eye(R,'like',H);

            % Use GPU-optimized determinant and matrix operations
            temp_matrix = eye(R,'like',H) + numerator / denominator;
            rate(i) = log2(det(temp_matrix));
        end
        system_rate = real(sum(rate.*alpha1(:,1),'all'));
    else
        % Original CPU implementation
        rate = zeros(I,1);
        for i=1:I
                denominator = zeros(R,R);
                for l=1:I
                    denominator = denominator + H(:,:,i)*V(:,:,l)*V(:,:,l)'*H(:,:,i)';
                end
                numerator = H(:,:,i)*V(:,:,i)*V(:,:,i)'*H(:,:,i)';
                denominator = denominator - numerator + sigma2*eye(R);

                rate(i) = log2(det(eye(R)+numerator / denominator));
        end
        system_rate = real(sum(rate.*alpha1,'all'));
    end
end