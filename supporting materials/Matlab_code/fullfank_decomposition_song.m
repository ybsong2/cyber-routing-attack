function [P,Q]=fullfank_decomposition_song(A)
% 本demo实现对任意矩阵的满秩分解
matrix_A=rank(A);
[m,n]=size(A);
B_row=rref([A,eye(m)]);%首先对A初等行变换
B_row_P=inv(B_row(:,n+1:n+m));
B_column=rref([B_row(:,1:n);eye(n)]');%对行变换后的A再进行列变换
B_column_T=B_column';
B_row_Q=inv(B_column_T(m+1:m+n,:));
P=B_row_P(:,1:matrix_A);
Q=B_row_Q(1:matrix_A,:);