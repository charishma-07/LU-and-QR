% QR factorization
%TATA CHARISHMA
%2020MCB1252
clc
clear all
A = input('Enter the matrix A\n');
[m n] = size(A);
Q = eye(n);         
R = eye(n);         
X = zeros(n,1)       
if m~=n
    print('The matrix A is not a square matrix\n');
    return
end
% A = [a1 , a2, ......... , an]
% Q  = [q1, q2, ..........., qn]
% R = [r1, r2,.........., rn]
b1 = norm(A(1:m,1),2) ;     % Assign b1 to ||a1||  
Q(1:m,1) = A(1:m,1)/b1 ;        % q1 = a1/b1
R(1,1) = b1;
for i = 2: n
    g =zeros(n,1);
    for j = 1:i-1
      g = g +(dot(A(1:m,i),Q(1:m,j)))*Q(1:m,j)  % Summation of all projections of ai on to q1, q2,....,qi-1
      R(j,i)= dot(A(1:m,i),Q(1:m,j))    
    end
    X = A(1:m,i) - g ;      % Here, we obtain the perpendicular projections to Q
    Q(1:m, i) = X/norm(X,2) ;  % qi = Per(ai)/||ai||
    R(i,i) = norm(X,2);
end
Q 
R