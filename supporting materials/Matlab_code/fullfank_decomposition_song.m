function [P,Q]=fullfank_decomposition_song(A)
% ��demoʵ�ֶ������������ȷֽ�
matrix_A=rank(A);
[m,n]=size(A);
B_row=rref([A,eye(m)]);%���ȶ�A�����б任
B_row_P=inv(B_row(:,n+1:n+m));
B_column=rref([B_row(:,1:n);eye(n)]');%���б任���A�ٽ����б任
B_column_T=B_column';
B_row_Q=inv(B_column_T(m+1:m+n,:));
P=B_row_P(:,1:matrix_A);
Q=B_row_Q(1:matrix_A,:);