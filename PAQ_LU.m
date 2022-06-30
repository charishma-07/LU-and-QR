function [L, U, P,Q] = PAQ_LU(A)
    [row ,column] = size(A);
   
    U = A;
    L = eye(row);       
    P = eye(row);
    Q = eye(row);
    
    
    for i=1:column - 1
        [v , k] = max(abs(U(i:column,i:row)));   % v&k are  row vector
        [~ ,pos] = max(abs(v));                     % ~ is unused  in remaining part of code, so we used this symbol
        j_max = k(pos)+i-1;                       % row position of Max element
        c_max = pos+i-1;                          % coloumn position of Max element
        
        if j_max ~= i || c_max ~= i
            if(c_max ~= i)
            U(:,[i,c_max]) = U(:,[c_max,i]);
            Q(:,[i,c_max]) = Q(:,[c_max,i]);
            end
            if j_max ~= i
            U([i,j_max],:) = U([j_max,i],:);
            P([i,j_max],:) = P([j_max,i],:);
            L([i,j_max],1:i-1) = L([j_max,i],1:i-1);
            end
        end
    
        for j=i+1 : row
            L(j,i) = U(j,i)/U(i,i);
            U(j,:) = U(j,:) - L(j,i)*U(i,:);
            if U(j,:)== 0
                U(j,:) = 0;
            end
        end
    end
end