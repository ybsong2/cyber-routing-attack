%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%本程序用于计算攻击目标One step objective（OO），并获取最优的攻击策略
%本程序的执行顺序：System_model-->本程序
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clc,clear
x_dim=size(Ad,1);
time=40;
%ele_row1=[1 0 0 0];ele_row2=[0 1 0 0];ele_row3=[0 0 1 0];ele_row4=[0 0 0 1];
%% 生成所有的可行解，存入attack_mati
quanpailie=perms([1 2 3 4 5]);
[m,n]=size(quanpailie);
attack_mati=cell(1,m);          %存储所有的可行解
for i=1:m
    attack_mati{i}=zeros(5,5);%给所有的攻击矩阵赋初值为0
end
for i=1:m
    for j=1:n
       attack_mati{i}(j,quanpailie(i,j))=1;
    end
end
%% 求解L_0
L0=cell(1,time);%存储目标函数的过程变量\mathcal{L}_0
for i=1:time
    middle_var=0;
    for j=0:i-1
        middle_var=middle_var+Ad^j*Bd*U(:,i-j)*U(:,i-j)'*Bd'*Ad'^j;
    end
    L0{i}=Ad^i*P0*Ad'^i+middle_var;
end

%% 求得最大值，即以OO为目标函数，求解优化问题P1
bar_P=P0;
I=eye(x_dim);
K=bar_P*C'*inv(C*bar_P*C'+R);
P_corrupt=cell(1,m); %在当前时刻存储所有策略下的攻击矩阵
optimal_stratagy=cell(1,time); %存储最优攻击策略
P_corrupt{1}=bar_P;
P_corrupt_tihuan=bar_P;     %存储当前时刻下的最大目标函数值
P_corrupt_trace=zeros(1,m);
P_max=trace(bar_P);
P_max_change=zeros(1,time);     %存储所有攻击区间内的目标函数值
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
    xxxx_OO(k)=optimal_stra;  %展示最优的攻击策略
    P_corrupt_tihuan=P_corrupt{optimal_stra+1};%存储受损的\tilde{P}_{k-1}，用于下一次循环使用
    P_max_change(k)=P_max;      %存放最优的目标函数值\tilde{P}_k
    optimal_stratagy{k}=attack_mati{optimal_stra};      %存放每个时刻的最优routing攻击策略
end



  

