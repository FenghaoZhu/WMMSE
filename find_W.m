function W = find_W(U,H,V, R, I,d)
    W = zeros(d,d,I);
    for i=1:I
            W(:,:,i) = inv(eye(d)-U(:,:,i)'*H(:,:,i)*V(:,:,i)); 
    end
end