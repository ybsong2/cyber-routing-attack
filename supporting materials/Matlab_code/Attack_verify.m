%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%本程序用于计算最优的攻击策略，攻击起始时刻是由变量m0给定的，攻击时长为time-m0
%本程序的执行顺序：System_model-->本程序
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clc,clear
x_dim=size(Ad,1);
time=41;
%ele_row1=[1 0 0 0];ele_row2=[0 1 0 0];ele_row3=[0 0 1 0];ele_row4=[0 0 0 1];
%% 生成所有的可行解，存入attack_mati
quanpailie=perms([1 2 3 4 5]);
[m,n]=size(quanpailie);
attack_mati=cell(1,m);%存储所有的可行解
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
        middle_var=middle_var+Ad^j*Bd*U*U'*Bd'*Ad'^j;
    end
    L0{i}=Ad^i*P0*Ad'^i+middle_var;
end

%% 求解最优化问题
m0=1;
bar_P=P0;
I=eye(x_dim);
K=bar_P*C'*inv(C*bar_P*C'+R);
P_max_change_TO=zeros(1,time);
optimal_stratagy_TO=cell(1,time);
for k=m0+1:time
    J_T_trace=zeros(1,m);
    P_max_TO=trace(bar_P);
    for i=1:m
        Att_matrix=attack_mati{i};
        middle_var1=zeros(x_dim);
        for j=0:k-1-m0
            middle_var1=middle_var1+(Ad-K*Att_matrix*C*Ad)^j*(I-K*Att_matrix*C)*...
                            Q*(I-K*Att_matrix*C).'*((Ad-K*Att_matrix*C*Ad).')^j;
        end
        middle_var2=zeros(x_dim);
        for j=0:k-1-m0              
            middle_var2=middle_var2+(Ad-K*Att_matrix*C*Ad)^j*(K*Att_matrix*R*Att_matrix'*K')*((Ad-K*Att_matrix*C*Ad).')^j;
        end
        middle_var4=zeros(x_dim);
        for j=1:k-1-m0              %这里少算了一项，用middle_var3来给它添加上
            middle_var4=middle_var4+(Ad-K*Att_matrix*C*Ad)^j*(K*C+K*Att_matrix*C)*L0{j}*(K*C+K*Att_matrix*C)'...
                *((Ad-K*Att_matrix*C*Ad).')^j;
        end
        middle_var3=(K*C+K*Att_matrix*C)*bar_P*(K*C+K*Att_matrix*C).';
        J_T_trace(i)=trace((Ad-K*Att_matrix*C*Ad)^(k-m0)*bar_P*...                           
                        ((Ad-K*Att_matrix*C*Ad).')^(k-m0)+middle_var1+middle_var2+middle_var3+middle_var4);
        if P_max_TO<J_T_trace(i)
            P_max_TO=J_T_trace(i);
            optimal_stra_TO=i;
        end
    end
    xxxx_TO(k)=optimal_stra_TO;
    P_max_change_TO(k)=P_max_TO;      %存放最优的目标函数值\tilde{P}_k
    optimal_stratagy_TO{k}=attack_mati{optimal_stra_TO};
end
plot(1:time,log10(P_max_change_TO));


P_max_change_HO=zeros(1,time);
optimal_stratagy_HO=cell(1,time);
for k=m0+1:time
    J_H_trace=zeros(1,m);
    P_max_HO=trace(bar_P);
    for i=1:m
        Att_matrix=attack_mati{i};
        middle_var1_HO=zeros(x_dim);
        for j=0:k-m0
            middle_var1_HO=middle_var1_HO+(Ad-K*Att_matrix*C*Ad)^j*bar_P*((Ad-K*Att_matrix*C*Ad).')^j;
        end
        middle_var2_HO=zeros(x_dim);
        for l=0:k-1-m0
            for j=0:l
            middle_var2_HO=middle_var2_HO+(Ad-K*Att_matrix*C*Ad)^j*(I-K*Att_matrix*C)*Q*(I-K*Att_matrix*C).'*((Ad-K*Att_matrix*C*Ad).')^j;
            end
        end
        middle_var4_HO=zeros(x_dim);
        for l=0:k-1-m0
            for j=1:l             %这里少算了一项，用middle_var3来给它添加上
            middle_var4_HO=middle_var4_HO+(Ad-K*Att_matrix*C*Ad)^j*(K*C+K*Att_matrix*C)*L0{j}*(K*C+K*Att_matrix*C)'...
                *((Ad-K*Att_matrix*C*Ad).')^j;
            end
        end
        middle_var3_HO=(k-m0)*(K*C+K*Att_matrix*C)*bar_P*(K*C+K*Att_matrix*C).';
        middle_var5_HO=zeros(x_dim);
        for l=0:k-m0-1
            for j=0:l 
                middle_var5_HO=middle_var5_HO+(Ad-K*Att_matrix*C*Ad)^j*K*Att_matrix*R*Att_matrix'*K'*((Ad-K*Att_matrix*C*Ad).')^j;
            end
        end        
        J_H_trace(i)=trace(middle_var1_HO+middle_var2_HO+middle_var3_HO+middle_var4_HO+middle_var5_HO);
        if P_max_HO<J_H_trace(i)
            P_max_HO=J_H_trace(i);
            optimal_stra_HO=i;
        end
    end
    xxxx_HO(k)=optimal_stra_HO;
    P_max_change_HO(k)=P_max_HO;      %存放最优的目标函数值\tilde{P}_k
    optimal_stratagy_HO{k}=attack_mati{optimal_stra_HO};
end
figure
plot(1:time,P_max_change_HO);



