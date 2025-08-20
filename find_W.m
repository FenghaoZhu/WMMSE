function W = find_W(U,H,V, R, I,d,useGPU)
% GPU-accelerated version of find_W function
% useGPU: boolean flag indicating whether to use GPU acceleration

    if nargin < 7
        useGPU = false; % Default to CPU if not specified
    end
    
    if useGPU
        % GPU-accelerated implementation
        W = zeros(d,d,I,'like',U);
        for i=1:I
            % Use GPU-optimized matrix inversion
            temp_matrix = eye(d,'like',U) - U(:,:,i)' * H(:,:,i) * V(:,:,i);
            W(:,:,i) = inv(temp_matrix);
        end
    else
        % Original CPU implementation
        W = zeros(d,d,I);
        for i=1:I
                W(:,:,i) = inv(eye(d)-U(:,:,i)'*H(:,:,i)*V(:,:,i)); 
        end
    end
end