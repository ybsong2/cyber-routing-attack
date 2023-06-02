
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%本程序用于计算攻击目标，并比较不同攻击目标下的攻击效果(OO,TO,HO,without attack)
%本程序的执行顺序：System_model-->本程序
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% times=2;
% bar_P=P0;
% x_dim=size(Ad,1);    %给定输出的维度
% syms a1 a2 a3 a4 b1 b2 b3 b4 c1 c2 c3 c4 d1 d2 d3 d4
% Att_matrix=[a1,a2,a3,a4;b1,b2,b3,b4;c1,c2,c3,c4;d1,d2,d3,d4];%假设攻击矩阵R^C如左边所示
% %% 计算中间变量\mathcal{L}_0(k)
% L0=cell(1,times);
% for i=1:times
%     middle_var=0;
%     for j=0:i-1
%         middle_var=middle_var+Ad^j*Bd*U*U'*Bd'*Ad'^j;
%     end
%     L0{i}=Ad^i*P0*Ad'^i+middle_var;
% end
% %% 计算误差协方差的迭代式
% I=eye(x_dim);
% K=bar_P*C'*inv(C*bar_P*C'+R);
% P_corrupt=cell(1,times);
% huajian_P_corrupt=cell(1,times);
% %huajian_P_corrupt=(1,times);
% P_corrupt{1}=bar_P;
% for i=1:times
%     P_corrupt{i+1}=(I-K*Att_matrix*C)*(A*P_corrupt{i}*A'+Q)*(I-K*Att_matrix*C).'...
%         +(K*C+K*Att_matrix*C)*L0{i}*(K*C+K*Att_matrix*C).'+K*Att_matrix*R*Att_matrix'*K.';
%     P_corrupt_trace=trace(P_corrupt{i+1});
%     huajian_P_corrupt{i}=vpa(P_corrupt_trace,2)
%     %huajian_P_corrupt=vpa()
% end


%% 实验：攻击发起于任意时刻
x_dim=size(Ad,1);
time=100;
%ele_row1=[1 0 0 0];ele_row2=[0 1 0 0];ele_row3=[0 0 1 0];ele_row4=[0 0 0 1];
%% 生成所有的可行解，存入attack_mati
quanpailie=perms([1 2 3 4 5]);
[m,n]=size(quanpailie);
attack_mati=cell(1,m);  %存储所有的可行解
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

%% 给定攻击存在的时间区间
k_upperbound=20;
k_interval1=10;
k_interval2=40;
k_interval3=70;
k_interval=[10,40,70];

%% 求解最优化问题
times_bound=100;
bar_P=P0;
I=eye(x_dim);
K=bar_P*C'*inv(C*bar_P*C'+R);
P_max_change_TO=zeros(1,time);
P_max_change_HO=zeros(1,time);
P_max_change=zeros(1,time);     %存储所有攻击区间内的目标函数值
optimal_stratagy_TO=cell(1,time);
xxxx=zeros(1,time);%存放OO的攻击行为
xxxx_TO=zeros(1,time);%存放TO的攻击行为
xxxx_HO=zeros(1,time);%存放HO的攻击行为
% detection=zeros(1:time);%存放OO的检测情况
% detection_TO=zeros(1:time);
% detection_HO=zeros(1:time);

%% 无攻击区间
for times=1:times_bound
    if times<=k_interval1
        P_max_change_TO(times)=trace(bar_P);
        P_max_change_HO(times)=trace(bar_P);
        P_max_change(times)=trace(bar_P);
        xxxx(times)=NaN;
        xxxx_HO(times)=NaN;
        xxxx_TO(times)=NaN;
        detection(times)=NaN;
        detection_TO(times)=NaN;
        detection_HO(times)=NaN;
    elseif times>k_interval1+k_upperbound&&times<=k_interval2
        P_max_change_TO(times)=trace(bar_P);
        P_max_change_HO(times)=trace(bar_P);
        P_max_change(times)=trace(bar_P);
        xxxx(times)=NaN;
        xxxx_HO(times)=NaN;
        xxxx_TO(times)=NaN;
        detection(times)=NaN;
        detection_TO(times)=NaN;
        detection_HO(times)=NaN;
    else times>k_interval2+k_upperbound&&times<=k_interval3
        P_max_change_TO(times)=trace(bar_P);
        P_max_change_HO(times)=trace(bar_P);
        P_max_change(times)=trace(bar_P);
        xxxx(times)=NaN;
        xxxx_HO(times)=NaN;
        xxxx_TO(times)=NaN;
        detection(times)=NaN;
        detection_TO(times)=NaN;
        detection_HO(times)=NaN;
    end
end

% TO攻击区间
%       for n=1:size(k_interval,2)
%         for k=k_interval(n)+1:k_interval(n)+k_upperbound
%             J_T_trace=zeros(1,m);
%             P_max_TO=trace(bar_P);
%             for i=1:m
%                 Att_matrix=attack_mati{i};
%                 middle_var1=zeros(x_dim);
%                 for j=k_interval(n)+1:k-1
%                     middle_var1=middle_var1+(Ad-K*Att_matrix*C*Ad)^j*(I-K*Att_matrix*C)*...
%                         Q*(I-K*Att_matrix*C).'*((Ad-K*Att_matrix*C*Ad).')^j;
%                 end
%                 middle_var2=zeros(x_dim);
%                 for j=k_interval(n)+1:k-1
%                     middle_var2=middle_var2+(Ad-K*Att_matrix*C*Ad)^j*(K*Att_matrix*R*Att_matrix'*K')*((Ad-K*Att_matrix*C*Ad).')^j;
%                 end
%                 middle_var4=zeros(x_dim);
%                 for j=k_interval(n)+1:k-1              %这里少算了一项，用middle_var3来给它添加上
%                     middle_var4=middle_var4+(Ad-K*Att_matrix*C*Ad)^j*(K*C+K*Att_matrix*C)*L0{j}*(K*C+K*Att_matrix*C)'...
%                         *((Ad-K*Att_matrix*C*Ad).')^j;
%                     %trace(middle_var4)
%                 end
%                 %middle_var3=(K*C+K*Att_matrix*C)*bar_P*(K*C+K*Att_matrix*C).';
%                 %         J_T_trace(i)=trace((Ad-K*Att_matrix*C*Ad)^(k-chenk_k)*((I-K*Att_matrix*C)*bar_P*(I-K*Att_matrix*C).'+...
%                 %                             ((K*C+K*Att_matrix*C)*bar_P*(K*C+K*Att_matrix*C).'))*...
%                 %                         ((Ad-K*Att_matrix*C*Ad).')^(k-chenk_k)+middle_var1+middle_var2+middle_var4);
%                 J_T_trace(i)=trace((Ad-K*Att_matrix*C*Ad)^(k-k_interval(n))*bar_P*...
%                     ((Ad-K*Att_matrix*C*Ad).')^(k-k_interval(n))+middle_var1+middle_var2+middle_var4);
%                 if P_max_TO<J_T_trace(i)
%                     P_max_TO=J_T_trace(i);
%                     optimal_stra_TO=i;
%                 end
%             end
%             xxxx_TO=optimal_stra_TO
%             P_max_change_TO(k)=P_max_TO;      %存放最优的目标函数值\tilde{P}_k
%             optimal_stratagy_TO{k}=attack_mati{optimal_stra_TO};
%         end
%       end

% 仅按照最优策略计算TO
%Att_matrix=attack_mati{8};
%optimal_stra=[1 1 1 1 1 14 14 14 8 8 8 8 8 8 8 8 8 8 8 8];
optimal_stra=[28 28 28 80 80 38 38 38 38 38 38 38 38 38 38 38 38 38 38 38];
for n=1:size(k_interval,2)
    P_max_change_TO_matrix=bar_P;
    middle_cunchuliang=1;
        for k=k_interval(n)+1:k_interval(n)+k_upperbound
            xxxx_TO(k)=38;
            %Att_matrix=attack_mati{optimal_stra(middle_cunchuliang)};
            Att_matrix=attack_mati{38};
%              P_max_change_TO_matrix=(I-K*Att_matrix*C)*(Ad*P_max_change_TO_matrix*Ad'+Q)*(I-K*Att_matrix*C).'...
%                 +(K*C+K*Att_matrix*C)*L0{k}*(K*C+K*Att_matrix*C).'+K*Att_matrix*R*Att_matrix'*K.';
        middle_var1=zeros(x_dim);
        for j=0:k-1-k_interval(n)
            middle_var1=middle_var1+(Ad-K*Att_matrix*C*Ad)^j*(I-K*Att_matrix*C)*...
                            Q*(I-K*Att_matrix*C).'*((Ad-K*Att_matrix*C*Ad).')^j;
        end
        middle_var2=zeros(x_dim);
        for j=0:k-1-k_interval(n)              
            middle_var2=middle_var2+(Ad-K*Att_matrix*C*Ad)^j*(K*Att_matrix*R*Att_matrix'*K')*((Ad-K*Att_matrix*C*Ad).')^j;
        end
        middle_var4=zeros(x_dim);
        for j=1:k-1-k_interval(n)             %这里少算了一项，用middle_var3来给它添加上
            middle_var4=middle_var4+(Ad-K*Att_matrix*C*Ad)^j*(K*C+K*Att_matrix*C)*L0{j}*(K*C+K*Att_matrix*C)'...
                *((Ad-K*Att_matrix*C*Ad).')^j;
        end
        middle_var3=(K*C+K*Att_matrix*C)*bar_P*(K*C+K*Att_matrix*C).';
        P_max_change_TO(k)=trace((Ad-K*Att_matrix*C*Ad)^(k-k_interval(n))*bar_P*...                           
                        ((Ad-K*Att_matrix*C*Ad).')^(k-k_interval(n))+middle_var1+middle_var2+middle_var3+middle_var4);
            middle_cunchuliang=middle_cunchuliang+1;
        end
end

%% 平均误差
%optimal_stra1=[1 1 1 1 1 1 14 14 8 8 8 8 8 8 8 8 8 8 8 8];
optimal_stral=[28 28 28 80 80 38 38 38 38 38 38 38 38 38 38 38 38 38 38 38];
for n=1:size(k_interval,2)
    P_max_change_HO_matrix=bar_P;
    middle_cunchuliang=1;
        for k=k_interval(n)+1:k_interval(n)+k_upperbound
            xxxx_HO(k)=38;            
           % Att_matrix=attack_mati{optimal_stra(middle_cunchuliang)}; 
            Att_matrix=attack_mati{38};
            middle_var1_HO=zeros(x_dim);
        for j=0:k-k_interval(n)
            middle_var1_HO=middle_var1_HO+(Ad-K*Att_matrix*C*Ad)^j*bar_P*((Ad-K*Att_matrix*C*Ad).')^j;
        end
        middle_var2_HO=zeros(x_dim);
        for l=0:k-1-k_interval(n)
            for j=0:l
            middle_var2_HO=middle_var2_HO+(Ad-K*Att_matrix*C*Ad)^j*(I-K*Att_matrix*C)*Q*(I-K*Att_matrix*C).'*((Ad-K*Att_matrix*C*Ad).')^j;
            end
        end
        middle_var4_HO=zeros(x_dim);
        for l=0:k-1-k_interval(n)
            for j=1:l             %这里少算了一项，用middle_var3来给它添加上
            middle_var4_HO=middle_var4_HO+(Ad-K*Att_matrix*C*Ad)^j*(K*C+K*Att_matrix*C)*L0{j}*(K*C+K*Att_matrix*C)'...
                *((Ad-K*Att_matrix*C*Ad).')^j;
            end
        end
        middle_var3_HO=(k-k_interval(n))*(K*C+K*Att_matrix*C)*bar_P*(K*C+K*Att_matrix*C).';
        middle_var5_HO=zeros(x_dim);
        for l=0:k-k_interval(n)-1
            for j=0:l 
                middle_var5_HO=middle_var5_HO+(Ad-K*Att_matrix*C*Ad)^j*K*Att_matrix*R*Att_matrix'*K'*((Ad-K*Att_matrix*C*Ad).')^j;
            end
        end        
        P_max_change_HO(k)=trace(middle_var1_HO+middle_var2_HO+middle_var3_HO+middle_var4_HO+middle_var5_HO)/(k-k_interval(n));%求解攻击下的平均误差协方差
         middle_cunchuliang=middle_cunchuliang+1;
        end
end
      
%% one-step攻击区间
P_corrupt=cell(1,m); %在当前时刻存储所有策略下的攻击矩阵
optimal_stratagy=cell(1,time); %存储最优攻击策略
P_corrupt{1}=bar_P;
P_corrupt_trace=zeros(1,m);
P_max=trace(bar_P);

for n=1:size(k_interval,2)
    P_corrupt_tihuan=bar_P;     %存储当前时刻下的最大目标函数值
    for k=k_interval(n)+1:k_interval(n)+k_upperbound
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
        xxxx(k)=optimal_stra;  %展示最优的攻击策略
        P_corrupt_tihuan=P_corrupt{optimal_stra+1};%存储受损的\tilde{P}_{k-1}，用于下一次循环使用
        P_max_change(k)=P_max;      %存放最优的目标函数值\tilde{P}_k
        optimal_stratagy{k}=attack_mati{optimal_stra};      %存放每个时刻的最优routing攻击策略\
        detection(k)=3;
    end
end

for n=2:size(k_interval,2)
    for k=k_interval(n)+1:k_interval(n)+k_upperbound
        detection_HO(k)=2;
        detection_TO(k)=1;
    end
end

P_without_attack=zeros(1,times);
for i=1:times
    P_without_attack(i)=trace(bar_P);
end
plot(1:times,log10(P_without_attack),'--',1:times,log10(P_max_change),'-',1:times,log10(P_max_change_TO),':',1:times,log10(P_max_change_HO),'-.','linewidth',2);
legend({'Without attack','One-step objective','Terminal objective','Holistic objective',},  'Interpreter','latex', 'FontSize', 18, 'location', 'northwest');
xlabel('Sampling instant' ,'Interpreter','latex','FontSize',18);%(a) $T=2N=24, \ \alpha=0.8$'
ylabel({'$\log(J(R^{C\star}))$'}, 'Interpreter','latex','FontSize',18);
set(gca,'FontSize',18);%坐标轴字体大小
figure
subplot(3,1,1)
% scatter(1:times,xxxx,'diamond');
% hold on
% scatter(1:times,xxxx_TO,'.');
% hold on
% scatter(1:times,xxxx_HO,'o');
plot(1:times,xxxx,'diamond',1:times,xxxx_TO,'.',1:times,xxxx_HO,'o','linewidth',2)
xlim([0,100]);
ylim([0,10]);
legend({'One-step objective','Terminal objective','Holistic objective',},  'Interpreter','latex', 'FontSize', 18, 'location', 'southeast');
xlabel('Sampling instant' ,'Interpreter','latex','FontSize',18);%(a) $T=2N=24, \ \alpha=0.8$'
ylabel({'Strategy$(R^{C\star})$'}, 'Interpreter','latex','FontSize',18);
set(gca,'FontSize',18);%坐标轴字体大小
subplot(3,1,2)
% scatter(1:times,detection,'pentagrm');
% hold on
% scatter(1:times,detection_TO,'diamond');
% hold on
% scatter(1:times,detection_HO,'o');
plot(1:times,detection,'pentagrm',1:times,detection_TO,'diamond',1:times,detection_HO,'o','linewidth',2);
set(gca,'YTick',0:1:4);
set(gca,'YTicklabel',{'','HO','TO','OO',''});
%legend({'One-step objective','Terminal objective','Holistic objective',},  'Interpreter','latex', 'FontSize', 18, 'location', 'northwest');
xlabel('Sampling instant' ,'Interpreter','latex','FontSize',18);%(a) $T=2N=24, \ \alpha=0.8$'
ylabel({'Alarm'}, 'Interpreter','latex','FontSize',18);
set(gca,'FontSize',18);%坐标轴字体大小
xlim([0,100]);
ylim([0,4]);
defense_action=xxxx;
defense_action(11)=NaN;
defense_action(41)=NaN;
defense_action(71)=NaN;
subplot(3,1,3)
plot(1:times,defense_action,'pentagrm','linewidth',2);
xlim([0,100]);
ylim([0,10]);
xlabel('Sampling instant' ,'Interpreter','latex','FontSize',18);%(a) $T=2N=24, \ \alpha=0.8$'
ylabel({'Strategy$(R^{D\star})$'}, 'Interpreter','latex','FontSize',18);
set(gca,'FontSize',18);%坐标轴字体大小

%% 有防御机制时，不同攻击策略下误差协方差的变化
P_under_defense=P_without_attack;
P_under_defense(11)=P_max_change(11);
P_under_defense(41)=P_max_change_TO(41);
P_under_defense(71)=P_max_change_HO(71);
P_under_attack=P_without_attack;

    for k=k_interval(1)+1:k_interval(1)+k_upperbound
        P_under_attack(k)=P_max_change(k);
        %P_under_attack(k)=P_max_change_TO(k);
        %P_under_attack(k)=P_max_change_HO(k);
    end
for k=k_interval(2)+1:k_interval(2)+k_upperbound
        %P_under_attack(k)=P_max_change(k);
        P_under_attack(k)=P_max_change_TO(k);
        %P_under_attack(k)=P_max_change_HO(k);
end
for k=k_interval(3)+1:k_interval(3)+k_upperbound
        %P_under_attack(k)=P_max_change(k);
        %P_under_attack(k)=P_max_change_TO(k);
        P_under_attack(k)=P_max_change_HO(k);
end
    
figure
plot(1:times,log10(P_without_attack),'--',1:times,log10(P_under_defense),'-',1:times,log10(P_under_attack),':','linewidth',2);
xlabel('Sampling instant' ,'Interpreter','latex','FontSize',18);%(a) $T=2N=24, \ \alpha=0.8$'
ylabel({'$\log(J(R^{C\star},R^{D\star}))$'}, 'Interpreter','latex','FontSize',18);
set(gca,'FontSize',18);%坐标轴字体大小
ylim([0.6,2]);
legend({'Without attack','With defense','With attack',},  'Interpreter','latex', 'FontSize', 18, 'location', 'northwest');
