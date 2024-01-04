function U = find_U(H,V,sigma2, P, R,I,d)

    J = zeros(R,R,I);  %计算不含噪声项的矩阵
    U = zeros(R,d,I);

    for i=1:I
            for l=1:I
                    J(:,:,i) = J(:,:,i) + H(:,:,i)*V(:,:,l)*(V(:,:,l)')*(H(:,:,i)'); 
            end
  
            U(:,:,i) = (J(:,:,i) + sigma2*eye(R,R)) \ (H(:,:,i)*V(:,:,i)); 
    end     
end