function V = find_V(alpha1, H, sigma2, U, W, T , R ,I ,d ,P, useGPU)
% GPU-accelerated version of find_V function
% useGPU: boolean flag indicating whether to use GPU acceleration

    if nargin < 11
        useGPU = false; % Default to CPU if not specified
    end
    
    if useGPU
        % GPU-accelerated implementation
        J = zeros(T, T,'like',H);
        V = zeros(T,d, I,'like',H);
       
        % Compute J matrix on GPU
        for l=1:I
            H_l = H(:,:,l);
            U_l = U(:,:,l);
            W_l = W(:,:,l);
            J = J + alpha1(l, 1) * H_l' * U_l * W_l * U_l' * H_l;   
        end

        max_iter = 100; % 二分法查找最优对偶变量\mu
        mu = zeros(1,1,'like',H);
        mu_min = 0;
        mu_max = 10;
        iter = 0;
        
        while(1)
            mu1 = (mu_max+mu_min) / 2;
            P_tem = 0;

            for i=1:I % 计算功率和
                H_i = H(:,:,i);
                U_i = U(:,:,i);
                W_i = W(:,:,i);
                V_tem = (J + mu1*eye(T,'like',H)) \ (alpha1(i,1) * H_i' * U_i * W_i); 
                P_tem = P_tem + real(trace(V_tem * V_tem'));
            end

            if gather(P_tem) > P  % gather for CPU comparison
                mu_min = mu1;
            else
                mu_max = mu1;
            end
            iter = iter + 1;

            if abs(mu_max - mu_min) < 1e-5 || iter > max_iter
                break
            end
        end

        mu = mu1;

        for l=1:I
            H_l = H(:,:,l);
            U_l = U(:,:,l);
            W_l = W(:,:,l);
            V(:,:,l) = (J + mu*eye(T, T,'like',H)) \ (alpha1(l, 1) * H_l' * U_l * W_l); 
        end
    else
        % Original CPU implementation
        J=zeros(T, T);
        V=zeros(T,d, I);
       
                for l=1:I
                    J = J + alpha1(l, 1) * H(:,:,l)'*U(:,:,l)*W(:,:,l)*(U(:,:,l)')*(H(:,:,l));   
                end

                max_iter = 100; % 二分法查找最优对偶变量\mu
                mu = zeros(1,1);
                mu_min = 0;
                mu_max = 10;
                iter = 0;
                while(1)
                    mu1 = (mu_max+mu_min) / 2;
                    P_tem = 0;

                    for i=1:I % 计算功率和
                        V_tem = ((J+mu1*eye(T))) \ (alpha1(i,1)*(H(:,:,i)'*U(:,:,i)*W(:,:,i))); 
                        P_tem = P_tem + real(trace(V_tem*V_tem'));
                    end

                    if P_tem > P
                        mu_min = mu1;
                    else
                        mu_max = mu1;
                    end
                    iter = iter + 1;

                    if abs(mu_max - mu_min) < 1e-5 || iter > max_iter
                        break
                    end
                end

                mu = mu1;

                for l=1:I
                    V(:,:,l) = (J+mu*eye(T, T)) \ (alpha1(l, 1) * (H(:,:,l)'*U(:,:,l)*W(:,:,l))); 
                end
    end
end