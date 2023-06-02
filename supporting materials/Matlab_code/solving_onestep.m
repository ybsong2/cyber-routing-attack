%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���������ڼ��㹥��Ŀ��One step objective��OO��������ȡ���ŵĹ�������
%�������ִ��˳��System_model-->������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clc,clear
x_dim=size(Ad,1);
time=40;
%ele_row1=[1 0 0 0];ele_row2=[0 1 0 0];ele_row3=[0 0 1 0];ele_row4=[0 0 0 1];
%% �������еĿ��н⣬����attack_mati
quanpailie=perms([1 2 3 4 5]);
[m,n]=size(quanpailie);
attack_mati=cell(1,m);          %�洢���еĿ��н�
for i=1:m
    attack_mati{i}=zeros(5,5);%�����еĹ������󸳳�ֵΪ0
end
for i=1:m
    for j=1:n
       attack_mati{i}(j,quanpailie(i,j))=1;
    end
end
%% ���L_0
L0=cell(1,time);%�洢Ŀ�꺯���Ĺ��̱���\mathcal{L}_0
for i=1:time
    middle_var=0;
    for j=0:i-1
        middle_var=middle_var+Ad^j*Bd*U(:,i-j)*U(:,i-j)'*Bd'*Ad'^j;
    end
    L0{i}=Ad^i*P0*Ad'^i+middle_var;
end

%% ������ֵ������OOΪĿ�꺯��������Ż�����P1
bar_P=P0;
I=eye(x_dim);
K=bar_P*C'*inv(C*bar_P*C'+R);
P_corrupt=cell(1,m); %�ڵ�ǰʱ�̴洢���в����µĹ�������
optimal_stratagy=cell(1,time); %�洢���Ź�������
P_corrupt{1}=bar_P;
P_corrupt_tihuan=bar_P;     %�洢��ǰʱ���µ����Ŀ�꺯��ֵ
P_corrupt_trace=zeros(1,m);
P_max=trace(bar_P);
P_max_change=zeros(1,time);     %�洢���й��������ڵ�Ŀ�꺯��ֵ
%xxxx_OO=zeros(1:time);
for k=1:time
    for i=1:m
        Att_matrix=attack_mati{i};
        P_corrupt{i+1}=(I-K*Att_matrix*C)*(Ad*P_corrupt_tihuan*Ad'+Q)*(I-K*Att_matrix*C).'...
            +(K*C+K*Att_matrix*C)*L0{k}*(K*C+K*Att_matrix*C).'+K*Att_matrix*R*Att_matrix'*K.';
        P_corrupt_trace(i)=trace(P_corrupt{i+1});
        if P_max<P_corrupt_trace(i)
            P_max=P_corrupt_trace(i);
            optimal_stra=i;
        end
    end
    xxxx_OO(k)=optimal_stra;  %չʾ���ŵĹ�������
    P_corrupt_tihuan=P_corrupt{optimal_stra+1};%�洢�����\tilde{P}_{k-1}��������һ��ѭ��ʹ��
    P_max_change(k)=P_max;      %������ŵ�Ŀ�꺯��ֵ\tilde{P}_k
    optimal_stratagy{k}=attack_mati{optimal_stra};      %���ÿ��ʱ�̵�����routing��������
end



  

