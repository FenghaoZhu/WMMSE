function U = find_U(H,V,sigma2, P, R,I,d,useGPU)
% GPU-accelerated version of find_U function
% useGPU: boolean flag indicating whether to use GPU acceleration

    if nargin < 8
        useGPU = false; % Default to CPU if not specified
    end
    
    if useGPU
        % GPU-accelerated implementation
        J = zeros(R,R,I,'like',H);  % Create GPU array with same type as H
        U = zeros(R,d,I,'like',H);
        
        for i=1:I
            J_temp = zeros(R,R,'like',H);
            for l=1:I
                H_i = H(:,:,i);
                V_l = V(:,:,l);
                J_temp = J_temp + H_i * V_l * V_l' * H_i';
            end
            J(:,:,i) = J_temp;
            
            % Use GPU-optimized linear solver
            U(:,:,i) = (J(:,:,i) + sigma2*eye(R,R,'like',H)) \ (H(:,:,i)*V(:,:,i));
        end
    else
        % Original CPU implementation
        J = zeros(R,R,I);  %计算不含噪声项的矩阵
        U = zeros(R,d,I);

        for i=1:I
                for l=1:I
                        J(:,:,i) = J(:,:,i) + H(:,:,i)*V(:,:,l)*(V(:,:,l)')*(H(:,:,i)'); 
                end
      
                U(:,:,i) = (J(:,:,i) + sigma2*eye(R,R)) \ (H(:,:,i)*V(:,:,i)); 
        end
    end
end