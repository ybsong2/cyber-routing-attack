%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���������ڼ��㹥��Ŀ��Historlic objective��HO��������ȡ���ŵĹ�������
%�������ִ��˳��System_model-->sloving_onestep-->solving_terminalobjective-->������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clc,clear
x_dim=size(Ad,1);
time=40;
%ele_row1=[1 0 0 0];ele_row2=[0 1 0 0];ele_row3=[0 0 1 0];ele_row4=[0 0 0 1];
%% �������еĿ��н⣬����attack_mati
quanpailie=perms([1 2 3 4 5]);
[m,n]=size(quanpailie);
attack_mati=cell(1,m);  %�洢���еĿ��н�
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

%% ������Ż�����

bar_P=P0;
I=eye(x_dim);
K=bar_P*C'*inv(C*bar_P*C'+R);
P_max_change_HO=zeros(1,time);
optimal_stratagy_HO=cell(1,time);
%xxxx_HO=zeros(1:time);
for k=1:time
    J_H_trace=zeros(1,m);
    P_max_HO=trace(bar_P);
    for i=1:m
        Att_matrix=attack_mati{i};
        middle_var1_HO=zeros(x_dim);
        for j=0:k
            middle_var1_HO=middle_var1_HO+(Ad-K*Att_matrix*C*Ad)^j*((I-K*Att_matrix*C)*bar_P*(I-K*Att_matrix*C).'+...
                            ((K*C+K*Att_matrix*C)*bar_P*(K*C+K*Att_matrix*C).'))*((Ad-K*Att_matrix*C*Ad).')^j;
        end
        middle_var2_HO=zeros(x_dim);
        for l=0:k-1
            for j=0:l
            middle_var2_HO=middle_var2_HO+(Ad-K*Att_matrix*C*Ad)^j*(I-K*Att_matrix*C)*Q*(I-K*Att_matrix*C).'*((Ad-K*Att_matrix*C*Ad).')^j;
            end
        end
        middle_var4_HO=zeros(x_dim);
        for l=0:k-1
            for j=1:l             %����������һ���middle_var3�����������
            middle_var4_HO=middle_var4_HO+(Ad-K*Att_matrix*C*Ad)^j*(K*C+K*Att_matrix*C)*L0{j}*(K*C+K*Att_matrix*C)'...
                *((Ad-K*Att_matrix*C*Ad).')^j;
            end
        end
        middle_var3_HO=k*(K*C+K*Att_matrix*C)*bar_P*(K*C+K*Att_matrix*C).';
        middle_var5_HO=zeros(x_dim);
        for l=0:k
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
    P_max_change_HO(k)=P_max_HO;      %������ŵ�Ŀ�꺯��ֵ\tilde{P}_k
    optimal_stratagy_HO{k}=attack_mati{optimal_stra_HO};
end
plot(1:time,xxxx_OO,'--',1:time,xxxx_TO,'-.',1:time,xxxx_HO,'-','linewidth',2);
legend({'One-step objective','Terminal objective','Holistic objective',},  'Interpreter','latex', 'FontSize', 18, 'location', 'northwest');
xlabel('Sampling instant' ,'Interpreter','latex','FontSize',18);%(a) $T=2N=24, \ \alpha=0.8$'
ylabel({'Optimal attack strtegy'}, 'Interpreter','latex','FontSize',18);
set(gca,'FontSize',18);%�����������С
   



