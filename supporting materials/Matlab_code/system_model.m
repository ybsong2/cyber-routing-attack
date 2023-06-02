%%% 
%本demo用于生成5-Generator&8-Line的离散化系统模型，其中调用fullfank_decompositon_song函数用于测量矩阵的满秩分解

clc,clear
delta_t=0.001;%1e-05;%给定采样周期
state_num=25;%系统状态的维数
%由当前的4个传感器扩展到n个需要修改的参量包括：R1=eye(6)*2e-02; R2=eye(7)*4e-02; R3=eye(5)*6e-02; sensor_num=3; C_original{1}=……
                  %135-136行：sensor_num赋值、R_original{i}赋值；161行：添加观测矩阵；SEB_scheduling_Attack的demo：4行触发参数beta


%% 初始化物理参数——电机参数
Gen_num=5;%同步电机的个数
Line_num=8;%电力传输线的条数
Line_num_connect=7;%节点连接线条数
M=[3.04 3.67 3.35 3.98 3.4];
D=[4.50 4.65 4.43 3.94 4.90];
X_d1=[0.0329 0.0329 0.3290 0.0239 0.2390];%对应变量X'_d
X_d=[0.1016 0.1016 1.016 0.1016 1.016];
T_do1=[5.57 5.57 5.57 5.57 5.57];
b_0=[656.00 656.00 656.00 656.00 656.00];
b_1=[1332.00 1332.00 1332.00 1332.00 1332.00];
c_0=[3.23 3.23 3.23 3.23 3.23];
c_1=[32.3 32.3 32.3 32.3 32.3];
V=[1.03 1.03 1.023 1.03 1.023];
theta=[0 0.1041 0.0933 0.0351 0.0607];
%初始化物理参数——电力传输线参数
% R=[0.00335 0.00313 0.03004 0.00514 0.00701 0.04200 0.01714];
% X=[0.01057 0.00368 0.05242 0.01074 0.02234 0.12685 0.04143];
% G_line=zeros(1,Line_num_connect);
% B_Line=zeros(1,Line_num_connect);
% for i=1:Line_num_connect
%     G_line(i)=R(i)/(R(i)^2+X(i)^2);
%     B_Line(i)=-X(i)/(R(i)^2+X(i)^2);
% end
%G1_G2_Z=R+Xj=0.03317+0.05610j;G1_G2_Y=G+Bj=7.80938-13.20790j;
%G1_G3_Z=R+Xj=0.07361+0.18957j;G1_G3_Y=G+Bj=1.77994-4.58394j;
%G1_G4_Z=R+Xj=0.09254+0.23144j;G1_G4_Y=G+Bj=1.48950-3.72520j;
%G1_G5_Z=R+Xj=0.03705+0.07374j;G1_G5_Y=G+Bj=4.70894-11.44387j;
%G2_G3_Z=R+Xj=0.04670+0.14110j;G2_G3_Y=G+Bj=2.11407-6.38748j;
%G2_G4_Z=R+Xj=0.06563+0.18270j;G2_G4_Y=G+Bj=1.74147-4.84788j;
%G2_G5_Z=R+Xj=0.01014+0.02599j;G2_G5_Y=G+Bj=13.02840-33.39331j;
%G3_G4_Z=R+Xj=0.02563+0.06274j;G3_G4_Y=G+Bj=5.57998-13.65931j;
%G3_G5_Z=R+Xj=0.05058+0.15973j;G3_G5_Y=G+Bj=1.80179-5.69001j;
%G4_G5_Z=R+Xj=0.07129+0.20136j;G4_G5_Y=G+Bj=1.56241-4.41307j;
line_G=[0,23.4596840433188,6.96307528613382,6.62605901832314,19.1995414627126;...
    23.4596840433188,0,11.1845799622121,8.56455045895540,59.5429682736801;...
    6.96307528613382,11.1845799622121,0,24.4280885355771,9.90669243080572;...
    6.62605901832314,8.56455045895540,24.4280885355771,0,7.81412810575724;...
    19.1995414627126,59.5429682736801,9.90669243080572,7.81412810575724,0];
line_B=[0,-12.5259451004693,-5.44008007003978,-5.13285297414408,-12.2304478629469;...
    -12.5259451004693,0,-10.2766825042292,-7.32803565988164,-48.2361123030678;...
    -5.44008007003978,-10.2766825042292,0,-18.8778647801975,-9.54258303856521;...
    -5.13285297414408,-7.32803565988164,-18.8778647801975,0,-6.96777150002868;...
    -12.2304478629469,-48.2361123030678,-9.54258303856521,-6.96777150002868,0];
I_d_com_delta_i=zeros(1,Gen_num);%对应的变量是alpha-I_{di}/alpha-delta_i
I_d_com_V_i=zeros(1,Gen_num);%对应的变量是alpha-I_{di}/alpha-E'_{qi}
I_d_com_delta_j=zeros(Gen_num,Gen_num);%alpha-I_{di}/alpha-delta_j
I_d_com_V_j=zeros(Gen_num,Gen_num);%alpha-I_{di}/alpha-E'_{qj}
alpha_Pei_com_delta_i=zeros(1,Gen_num);%对应的变量是alpha_Pei/alpha-delta_i
alpha_Pei_com_V_i=zeros(1,Gen_num);%对应的变量是alpha_Pei/alpha-E'_{qi}
alpha_Pei_com_delta_j=zeros(Gen_num,Gen_num);%alpha_Pei/alpha-delta_j
alpha_Pei_com_V_j=zeros(Gen_num,Gen_num);
for i=1:Gen_num
    for j=1:Gen_num
        I_d_com_delta_i(i)=I_d_com_delta_i(i)-V(i)*(line_B(i,j)*sin(theta(i)-theta(j))+line_G(i,j)*cos(theta(i)-theta(j)));
        I_d_com_V_i(i)=I_d_com_V_i(i)+line_B(i,j)*cos(theta(i)-theta(j))-line_G(i,j)*sin(theta(i)-theta(j));
        alpha_Pei_com_delta_i(i)=alpha_Pei_com_delta_i(i)+V(i)*(line_B(i,j)*cos(theta(i)-theta(j))-line_G(i,j)*sin(theta(i)-theta(j)))*V(j);
        alpha_Pei_com_V_i(i)=alpha_Pei_com_V_i(i)+(line_B(i,j)*sin(theta(i)-theta(j))+line_G(i,j)*cos(theta(i)-theta(j)))*V(j);
        if i~=j
            I_d_com_delta_j(i,j)=V(i)*(line_B(i,j)*sin(theta(i)-theta(j))+line_G(i,j)*cos(theta(i)-theta(j)));
            alpha_Pei_com_delta_j(i,j)=V(i)*(-line_B(i,j)*cos(theta(i)-theta(j))+line_G(i,j)*sin(theta(i)-theta(j)))*V(j);
            alpha_Pei_com_V_j(i,j)=V(i)*(line_B(i,j)*sin(theta(i)-theta(j))+line_G(i,j)*cos(theta(i)-theta(j)));
        end
    end  
end

%% 计算系统矩阵
Bsub=cell(Gen_num,Gen_num);
Asub_ij=cell(Gen_num,Gen_num);
psi=zeros(1,Gen_num);
for i=1:Gen_num %循环生成系统矩阵F和输入矩阵G
    psi(i)=(X_d(i)-X_d1(i))/T_do1(i);
    Bsub{i}=[0 0 0 1 0]';
    for j=1:Gen_num
        if i==j
            Asub_ij{i,i}=[0 1 0 0 0;-alpha_Pei_com_delta_i(i)/M(i) -D(i)/M(i) -alpha_Pei_com_V_i(i)/M(i) 0 0;...
             psi(i)*I_d_com_delta_i(i) 0 -1/T_do1(i)+psi(i)*I_d_com_V_i(i) b_1(i)/T_do1(i) b_0(i)/T_do1(i);...
             0 0 0 -c_1(i) -c_0(i); 0 0 0 1 0];
            Bsub{i,i}=[0 0 0 1 0]';
        else
            Asub_ij{i,j}=[0 0 0 0 0;-alpha_Pei_com_delta_j(i,j)/M(i) 0 -alpha_Pei_com_V_j(i,j)/M(i) 0 0;...
             psi(i)*I_d_com_delta_j(i,j) 0 -1/T_do1(i)+psi(i)*I_d_com_V_j(i,j) b_1(i)/T_do1(i) b_0(i)/T_do1(i);...
             0 0 0 0 0; 0 0 0 0 0];
         Bsub{i,j}=[0 0 0 0 0]';
        end
    end 
end
A_con=cell2mat(Asub_ij);%将元组子矩阵转化为完整的大A系统矩阵的形式
B_con=cell2mat(Bsub);%将元组子矩阵转化为完整的大B系统矩阵的形式
% Ad=eye(state_num)+A_con*delta_t;
% Bd=B_con*delta_t; %离散化系统参数A,B
[Ad,Bd]=c2d(A_con,B_con,delta_t);  %离散化系统参数A,B
Q=eye(state_num)*1e-03;

% %观测矩阵C是按照如下的方法随机生成的C=round(rand(5,25));rank(obsv(A,C))

%% 求解稳态\bar{p}
sensor_num=4;%给定传感器的个数
R_original=cell(1,sensor_num);
%R_original{1}=eye(7)*1e-05;R_original{2}=eye(7)*4e-02;R_original{3}=eye(5)*6e-01;R_original{4}=eye(6)*2e-02;%此行需要更改
R_original{1}=eye(7)*1e-05;R_original{2}=eye(7)*4e-04;R_original{3}=eye(5)*6e-02;R_original{4}=eye(6)*2e-02;
% R_original{5}=eye(7)*6e-03;R_original{6}=eye(6)*4e-03;R_original{7}=eye(5)*6e-03;R_original{8}=eye(5)*2e-03;
% R_original{9}=eye(6)*6e-03;R_original{10}=eye(6)*4e-03;R_original{11}=eye(6)*6e-03;R_original{12}=eye(5)*2e-03;

% C_original=cell(1,sensor_num);
C_original{1}=[0,1,1,0,1,1,1,0,1,1,1,0,1,0,0,1,1,0,0,0,1,0,1,1,1;...
    1,0,0,0,0,0,1,0,1,1,0,1,1,1,0,1,0,0,0,0,0,0,1,0,0;...
    0,1,1,0,1,1,1,0,0,1,0,1,0,0,0,0,0,1,1,0,0,0,0,1,0;...
    1,0,0,0,0,0,1,0,1,1,0,1,1,1,0,1,0,0,0,0,0,0,1,0,0;...
    0,1,1,0,1,1,1,0,0,1,0,1,0,0,0,0,0,1,1,0,0,0,0,1,0;...
    1,1,1,0,0,0,1,0,0,1,1,1,1,1,0,1,1,0,1,0,1,0,0,1,1;...
    0,1,1,1,1,1,1,0,1,0,1,1,1,0,1,1,1,1,1,0,1,1,0,0,1];
C_original{2}=[1,0,0,1,1,0,0,0,0,0,1,1,1,0,1,1,0,0,1,0,1,1,1,0,0;...
    0,0,1,0,0,1,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1;...
    1,0,0,0,1,1,0,0,0,1,1,1,0,1,0,0,1,1,1,1,0,0,0,1,0;...
    1,1,0,1,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,0,0,1,1;...
    1,1,0,0,1,1,1,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,0,1;...
    1,0,0,0,1,1,0,0,0,1,1,1,0,1,0,0,1,1,1,1,0,0,0,1,0;...
    1,1,0,1,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,0,0,1,1;];
 %C_original{3}=[1,0,0,1,0,1,1,1,1,1,1,1,0,1,1,0,1,0,0,0,0,1,0,1,1;0,1,0,1,1,1,1,1,0,0,0,0,0,1,0,1,0,0,0,1,1,1,1,1,1;0,0,0,1,1,1,0,1,1,0,0,0,0,1,1,1,1,0,0,1,0,1,1,1,0;0,1,0,0,0,0,1,1,1,0,0,0,1,0,0,0,0,1,1,1,0,0,1,0,0;0,0,0,0,0,0,0,0,1,1,0,1,1,1,1,1,1,0,0,0,0,1,0,0,1];
 C_original{3}=[0,0,1,1,0,1,1,1,0,1,1,0,1,1,0,0,1,1,0,1,0,1,1,1,0;1,1,0,1,1,1,0,1,1,1,1,0,0,1,0,0,1,1,1,0,1,0,0,1,0;0,0,1,1,0,0,1,1,1,0,1,0,0,0,0,1,0,1,1,1,1,1,1,0,1;0,1,0,0,0,1,0,1,0,1,0,1,0,1,0,1,1,1,0,1,1,1,1,1,1;1,0,1,0,0,0,0,0,1,1,1,0,1,1,1,0,1,1,1,0,1,0,1,1,1];
C_original{4}=[1,0,1,1,1,1,1,0,1,0,1,1,1,1,0,1,0,1,0,0,0,0,0,0,1;...
    0,1,1,0,1,0,1,1,1,1,1,0,1,0,0,0,0,1,1,0,0,1,0,0,0;...
    0,0,0,1,0,1,0,1,1,0,0,1,0,1,0,0,0,1,0,0,1,0,1,1,1;...
    1,0,0,0,1,1,0,1,0,0,1,1,0,1,0,0,1,0,0,1,0,1,0,1,1;...
    1,0,0,0,1,1,0,1,0,0,1,1,0,1,0,0,1,0,0,1,0,1,0,1,1;...
    0,0,0,0,0,1,1,1,0,0,1,0,1,1,0,0,0,0,1,0,1,0,1,1,1];
 %C_original{1}=round(rand(5,25));
%C_original{2}=round(rand(5,25));
%C_original{3}=round(rand(5,25));
%C_original{4}=round(rand(5,25));

R1=cell(1,sensor_num);
C=cell(1,sensor_num);
M=cell(1,sensor_num);
O=cell(1,sensor_num);
rank_C=zeros(1,sensor_num);%用于存放满秩分解后得到的C的秩
for i=1:sensor_num
    [M{i},O{i}]=fullfank_decomposition_song(C_original{i});
    C{i}=(M{i}'*inv(R_original{i})*M{i})^(1/2)*O{i};
    rank_C(i)=min(size(C{i}));
    R1{i}=eye(rank_C(i));
end

err=1e-5;%定义计算精度
error=10;
% Last_P=P0;   %%注释掉的此部分可用于求解里卡蒂方程的半正定近似解
% j=1;
% while (error>=err)
% barP_1=A*Last_P*A'+Q-A*Last_P*C_1'*inv(C_1*Last_P*C_1'+R_1)*C_1*Last_P*A';
%   error=norm((barP_1-Last_P),'inf');
%   Last_P=barP_1;
%   j=j+1;
% end
Last_P=cell(1,sensor_num);
for i=1:sensor_num
    Last_P{i}=eye(state_num);%给定状态估计误差的初值
end
j=1;
while (error>=err)
    error_mid=zeros(1,sensor_num);
    for i=1:sensor_num
        barP_1=Ad*Last_P{i}*Ad'+Q-Ad*Last_P{i}*C{i}'*inv(C{i}*Last_P{i}*C{i}'+R1{i})*C{i}*Last_P{i}*Ad';
        error_mid(i)=norm((barP_1-Last_P{i}),'inf');
        Last_P{i}=barP_1;
    end
    [error,seq]=max(error_mid);
    j=j+1;
end



%% 系统运行
%给定系统初始参数
N=100; %仿真时间，时间序列总数
W=sqrt(Q)*randn(state_num,N);%既然Q为0，即W=0
V=cell(1,sensor_num);
Z=cell(1,sensor_num);%给定输出
Xkf=cell(1,sensor_num);
P0_1=cell(1,sensor_num);
P0_attack=cell(1,sensor_num);
U=zeros(state_num/5,N);
U(:,1)=[0.06;0.02;0.07;0.05;0.05];%初始化系统输入
X=zeros(state_num,N);%物体真实状态
%X(:,1)=[1;2;1;0;1;1;1;2;1;1;2;1;1;2;1;0;1;1;1;2;1;1;2;1;1];%初始化状态变量
X(:,1)=[0.001;0.002;0.001;0;0.001;0.001;0.001;0.002;0.001;0.001;0.002;0.001;0.001;0.002;0.001;0.000;0.001;0.001;0.001;0.002;0.001;0.001;0.002;0.001;0.001];%初始化状态变量
for i=1:sensor_num
    V{i}=sqrt(R1{i})*randn(rank_C(i),N);
    Z{i}=zeros(rank_C(i),N);
    Z{i}(:,1)=C{i}*X(:,1);
    Xkf{i}=zeros(state_num,N);%卡尔曼估计状态初始化
    Xkf{i}(:,1)=X(:,1);
    P0_1{i}=Last_P{i};
    trace(P0_1{i})
    P0_attack{i}=Last_P{i};
end
I=eye(state_num); %25维系统


%% 用于求解反馈增益L
setlmis([]);
P=lmivar(1,[state_num,1]);
K=lmivar(2,[Gen_num,state_num]);

lmiterm([-1 1 1 P],1,1);

lmiterm([2 1 1 P],-1,1);
lmiterm([2 2 1 K],Bd,1);
lmiterm([2 2 1 P],Ad,1);
lmiterm([2 2 2 P],-1,1);



lmisys=getlmis;
[tmin,xfeas]=feasp(lmisys);
K=dec2mat(lmisys,xfeas,K);
P=dec2mat(lmisys,xfeas,P);

L=K*inv(P);%反馈增益L

%% 重新补充变量
P0=P0_1{1};
C=C{1};
%C=C_original{3};
R=eye(5);

for k=2:N
    %底层物理系统的运行过程
    X(:,k)=Ad*X(:,k-1)+Bd*U(:,k-1)+W(:,k);         %物体下落，受状态方程的驱动
    i=1;  
    Z{i}(:,k)=C*X(:,k)+V{i}(:,k);               %位移传感器对目标进行观测
    
        X_pre=Ad*Xkf{i}(:,k-1)+Bd*U(:,k-1);             %状态预测 
        P_pre=Ad*P0_1{i}*Ad'+Q;                    %协方差预测
        Kg=P_pre*C'*inv(C*P_pre*C'+R1{i});      %计算卡尔曼增益
        %接下来需要判断传输的信息值是否受到攻击
        Xkf{i}(:,k)=X_pre+Kg*(Z{i}(:,k)-C*X_pre);   %状态更新
        P0_1{i}=(I-Kg*C)*P_pre;               %方差更新
            
    U(:,k-1)=L*Xkf{1}(:,k);
end
p_delay=cell(1,N);
%% 延迟验证
for k=1:N
    if k>10&&k<=30
        p_delay{k}=Ad*p_delay{k-1}*Ad'+Q;
    elseif k>40&&k<=60
        p_delay{k}=Ad*p_delay{k-1}*Ad'+Q;
    elseif k>70&&k<=90
        p_delay{k}=Ad*p_delay{k-1}*Ad'+Q;
    else
        p_delay{k}=P0;
    end
    p_delay_trace(k)=trace(p_delay{k});
end
figure
plot(1:N,log10(p_delay_trace),'--','linewidth',2);
legend({'delay',},  'Interpreter','latex', 'FontSize', 18, 'location', 'northwest');
xlabel('Sampling instant' ,'Interpreter','latex','FontSize',18);%(a) $T=2N=24, \ \alpha=0.8$'
ylabel({'$\log({Tr}(P_k))$'}, 'Interpreter','latex','FontSize',18);
set(gca,'FontSize',18);%坐标轴字体大小


figure
subplot(2,2,1)
plot(1:N,Xkf{1}(1,:),1:N,X(1,:));
subplot(2,2,2)
plot(1:N,Xkf{1}(6,:),1:N,X(6,:));
subplot(2,2,3)
plot(1:N,Xkf{1}(11,:),1:N,X(11,:));
subplot(2,2,4)
plot(1:N,Xkf{1}(16,:),1:N,X(16,:));
%  plot(1:N,Xkf{1}(15,:),1:N,Xkf{2}(15,:),1:N,Xkf{3}(15,:),1:N,X(15,:),1:N,Xkf{1}(25,:),1:N,Xkf{2}(25,:),1:N,Xkf{3}(25,:)),1:N,X(25,:);
figure
subplot(2,2,1)
plot(1:N,Xkf{1}(1,:),'-.',1:N,X(1,:),'-','linewidth',1.5);
xlabel('Sampling instant' ,'Interpreter','latex','FontSize',18);%(a) $T=2N=24, \ \alpha=0.8$'
ylabel({'p.u.'}, 'Interpreter','latex','FontSize',18);
set(gca,'FontSize',18);%坐标轴字体大小
legend({'$\Delta \hat{\delta}_1$','$\Delta \delta_1$',},  'Interpreter','latex', 'FontSize', 18, 'location', 'northwest');
title('(a)','Fontsize',18);
ylim([-1,2]);
subplot(2,2,2)
plot(1:N,Xkf{1}(6,:),'-.',1:N,X(6,:),'-','linewidth',1.5);
xlabel('Sampling instant' ,'Interpreter','latex','FontSize',18);%(a) $T=2N=24, \ \alpha=0.8$'
ylabel({'p.u.'}, 'Interpreter','latex','FontSize',18);
set(gca,'FontSize',18);%坐标轴字体大小
legend({'$\Delta \hat{\delta}_2$','$\Delta \delta_2$',},  'Interpreter','latex', 'FontSize', 18, 'location', 'northwest');
title('(b)','Fontsize',18);
ylim([-4,4]);
subplot(2,2,3)
plot(1:N,Xkf{1}(11,:),'-.',1:N,X(11,:),'-','linewidth',1.5);
xlabel('Sampling instant' ,'Interpreter','latex','FontSize',18);%(a) $T=2N=24, \ \alpha=0.8$'
ylabel({'p.u.'}, 'Interpreter','latex','FontSize',18);
set(gca,'FontSize',18);%坐标轴字体大小
legend({'$\Delta \hat{\delta}_3$','$\Delta \delta_3$',},  'Interpreter','latex', 'FontSize', 18, 'location', 'northwest');
title('(c)','Fontsize',18);
ylim([-4,2]);
subplot(2,2,4)
plot(1:N,Xkf{1}(16,:),'-.',1:N,X(16,:),'-','linewidth',1.5);
xlabel('Sampling instant' ,'Interpreter','latex','FontSize',18);%(a) $T=2N=24, \ \alpha=0.8$'
ylabel({'p.u.'}, 'Interpreter','latex','FontSize',18);
set(gca,'FontSize',18);%坐标轴字体大小
legend({'$\Delta \hat{\delta}_4$','$\Delta \delta_4$',},  'Interpreter','latex', 'FontSize', 18, 'location', 'northwest');
title('(d)','Fontsize',18);
ylim([-1,0.5]);
