
clear all
close all
clc


%固定參數
L0=225;     %可能要刪掉
L1=250;     %upper arm
L2=250;   %forearm
L3=150;      %length of end effector

x_base_R=0;   %基準點
y_base_R=0;
z_base_R=0;

x_base_L=0;   %基準點
y_base_L=0;
z_base_L=0;


DEF_DESCRETE_POINT=90;

%{
P7=[71.3397;-5.0000;0]
P8IK=[77.2464;3.069;0]
P8=[80;0;0]
V_r_hIK=[P8IK-P7]
norm(V_r_hIK)

V_r=[P8-P7]
theta1=acos(V_r_hIK'*V_r/(norm(V_r_hIK)*norm(V_r)))
theta1=theta1*(180)/pi
%}
 
 %產生正方形O(60,0,0)  Q(80,-20,0)  R(80,-20,-20)  S(60,0,-20)
 %把此路徑分成90份
% O_R=[300 -50 0];
% Q_R=[300 -200 0];
% R_R=[300 -200 -220];
% S_R=[300 -50 -220];
% 
% O_L=[300 50 0];
% Q_L=[300 200 0];
% R_L=[300 200 -220];
% S_L=[300 50 -220];
 
Path_R=zeros(DEF_DESCRETE_POINT,3);%規畫的路徑點
PathPoint_R=zeros(DEF_DESCRETE_POINT,3);%記錄實際上的點，畫圖使用
PathTheta_R=zeros(DEF_DESCRETE_POINT,7);%記錄每軸角度，畫圖使用
 
Path_L=zeros(DEF_DESCRETE_POINT,3);%規畫的路徑點
PathPoint_L=zeros(DEF_DESCRETE_POINT,3);%記錄實際上的點，畫圖使用
PathTheta_L=zeros(DEF_DESCRETE_POINT,7);%記錄每軸角度，畫圖使用

%畫正方形做IK FK測試
% for t=1:1:DEF_DESCRETE_POINT
%     if t<=25
%         Path_R(t,1:3)=O_R+(Q_R-O_R)*t/25;
%         Path_L(t,1:3)=O_L+(Q_L-O_L)*t/25;
%     elseif t<=50
%         Path_R(t,1:3)=Q_R+(R_R-Q_R)*(t-25)/25;
%         Path_L(t,1:3)=Q_L+(R_L-Q_L)*(t-25)/25;
%     elseif t<=75
%         Path_R(t,1:3)=R_R+(S_R-R_R)*(t-50)/25;
%         Path_L(t,1:3)=R_L+(S_L-R_L)*(t-50)/25;
%     else 
%         Path_R(t,1:3)=S_R+(O_R-S_R)*(t-75)/15;
%         Path_L(t,1:3)=S_L+(O_L-S_L)*(t-75)/15;
%     end
% end

%% 由起始點往前進100進行右邊線的縫紉   布料大小200x200  邊緣10
R_point1 = [300 -10 0];
L_point1 = [300 90 0];

R_point2 = [500 -10 0]; %往前200
L_point2 = [500 90 0]; %往前200

%%進行旋轉90度
R_point3 = [300 -10 0]; %右手鬆開後 x往後退200 
L_point3 = [500 90 0]; %左手不動

R_point4 = [500 -10 0]; %右手x 圓周往前200   
L_point4 = [300 90 0];  %左手x 圓周往後200  


R_point5 = [300 -10 0]; %右手鬆開x往後200
L_point5 = [300 90 0]; %左手不動

%%往前200進行縫紉
R_point6 = [500 -10 0]; %右手x往前200
L_point6 = [500 90 0];  %左手x往前200


for t=1:1:DEF_DESCRETE_POINT
    if t<=DEF_DESCRETE_POINT*0.2
        Path_R(t,1:3)=R_point1+(R_point2-R_point1)*t/(DEF_DESCRETE_POINT*0.2);%往前200
        Path_L(t,1:3)=L_point1+(L_point2-L_point1)*t/(DEF_DESCRETE_POINT*0.2);%往前200
    elseif t<=DEF_DESCRETE_POINT*0.4
        Path_R(t,1:3)=R_point2+(R_point3-R_point2)*(t-DEF_DESCRETE_POINT*0.2)/(DEF_DESCRETE_POINT*0.2);%右手鬆開後 x往後退200
        Path_L(t,1:3)=L_point2+(L_point3-L_point2)*(t-DEF_DESCRETE_POINT*0.2)/(DEF_DESCRETE_POINT*0.2);%左手不動
    elseif t<=DEF_DESCRETE_POINT*0.6
        xcR=(500+300)*0.5;
        ycR=-10;
        zcR=0;
        rR=500-xcR;
        
        xcL=(500+300)*0.5;
        ycL=90;
        zcL=0;
        rL=500-xcL;
        
        
        Path_R(t,1:3)=[xcR ycR zcR]+rR*[cos( pi*(t-DEF_DESCRETE_POINT*0.4)/(DEF_DESCRETE_POINT*0.2) + pi) sin(pi*(t-DEF_DESCRETE_POINT*0.4)/(DEF_DESCRETE_POINT*0.2) + pi) 0]; %右手下到上弧形
        Path_L(t,1:3)=[xcL ycL zcL]+rL*[cos( pi*(t-DEF_DESCRETE_POINT*0.4)/(DEF_DESCRETE_POINT*0.2)) sin(pi*(t-DEF_DESCRETE_POINT*0.4)/(DEF_DESCRETE_POINT*0.2)) 0]; %左手上到下弧形
        
        %Path_R(t,1:3)=R_point3+(R_point4-R_point3)*(t-DEF_DESCRETE_POINT*0.4)/(DEF_DESCRETE_POINT*0.2);
        %Path_L(t,1:3)=L_point3+(L_point4-L_point3)*(t-DEF_DESCRETE_POINT*0.4)/(DEF_DESCRETE_POINT*0.2);
    elseif t<=DEF_DESCRETE_POINT*0.8
        Path_R(t,1:3)=R_point4+(R_point5-R_point4)*(t-DEF_DESCRETE_POINT*0.6)/(DEF_DESCRETE_POINT*0.2);
        Path_L(t,1:3)=L_point4+(L_point5-L_point4)*(t-DEF_DESCRETE_POINT*0.6)/(DEF_DESCRETE_POINT*0.2);
    elseif t<=DEF_DESCRETE_POINT
        Path_R(t,1:3)=R_point5+(R_point6-R_point5)*(t-DEF_DESCRETE_POINT*0.8)/(DEF_DESCRETE_POINT*0.2);
        Path_L(t,1:3)=L_point5+(L_point6-L_point5)*(t-DEF_DESCRETE_POINT*0.8)/(DEF_DESCRETE_POINT*0.2);
    end
    
end


%畫直線
%  O=[0 0 -(L1+L2+L3)]; %初始點
%  Q=[0 -20 -180]; %第一軸不會有到倒轉情況的末點路徑
% 
%  for t=1:1:DEF_DESCRETE_POINT
%         Path(t,1:3)=O+(Q-O)*(t-1)/DEF_DESCRETE_POINT;
%  end

for t=1:1:DEF_DESCRETE_POINT
 
    %輸入參數
    in_x_end_R=Path_R(t,1);
    in_y_end_R=Path_R(t,2);
    in_z_end_R=Path_R(t,3);
    
    in_x_end_L=Path_L(t,1);
    in_y_end_L=Path_L(t,2);
    in_z_end_L=Path_L(t,3);
   
    in_alpha_R=60*(pi/180);
    in_beta_R=0*(t/DEF_DESCRETE_POINT)*(pi/180);
    in_gamma_R=0*(t/DEF_DESCRETE_POINT)*(pi/180);
    
    in_alpha_L=-60*(pi/180);
    in_beta_L=0*(t/DEF_DESCRETE_POINT)*(pi/180);
    in_gamma_L=0*(t/DEF_DESCRETE_POINT)*(pi/180);

    Rednt_alpha_R=-90*(pi/180);
    Rednt_alpha_L=90*(pi/180);
    %輸出參數 initial
    %theta=zeros(1,7);
    
    %
    %末點位置in==>IK==>theta==>FK==>末點位置out
    %
    %inverse kinematic
    y_base_R=-L0;%header0 座標系偏移到shoulder0 座標系 差Y方向的L0
    theta_R=IK_7DOF(L1,L2,L3,x_base_R,y_base_R,z_base_R,in_x_end_R,in_y_end_R,in_z_end_R,in_alpha_R,in_beta_R,in_gamma_R,Rednt_alpha_R);
  	y_base_L=L0;%header0 座標系偏移到shoulder0 座標系 差Y方向的L0
    theta_L=IK_7DOF(L1,L2,L3,x_base_L,y_base_L,z_base_L,in_x_end_L,in_y_end_L,in_z_end_L,in_alpha_L,in_beta_L,in_gamma_L,Rednt_alpha_L);
   
    %forward kinematic
    %theta=[0 0 0 0 0 0 0];
    [out_x_end_R,out_y_end_R,out_z_end_R,out_alpha_R,out_beta_R,out_gamma_R,P_R,RotationM_R] = FK_7DOF(-L0,L1,L2,L3,x_base_R,y_base_R,z_base_R,theta_R);

    [out_x_end_L,out_y_end_L,out_z_end_L,out_alpha_L,out_beta_L,out_gamma_L,P_L,RotationM_L] = FK_7DOF(L0,L1,L2,L3,x_base_L,y_base_L,z_base_L,theta_L);
    %記錄路徑上的點
    PathPoint_R(t,1:3)=[out_x_end_R out_y_end_R out_z_end_R];
    PathPoint_L(t,1:3)=[out_x_end_L out_y_end_L out_z_end_L];
    
    %畫關節點圖
    Draw_7DOF_dual_point(P_R,RotationM_R,PathPoint_R,P_L,RotationM_L,PathPoint_L);
   
    %記錄每軸角度變化
    PathTheta_R(t,1:7)=theta_R*(180/pi);
    PathTheta_L(t,1:7)=theta_L*(180/pi);
    
    In_R=[in_x_end_R in_y_end_R in_z_end_R in_alpha_R in_beta_R in_gamma_R]
    Out_R=[out_x_end_R out_y_end_R out_z_end_R out_alpha_R out_beta_R out_gamma_R]
    
    In_L=[in_x_end_L in_y_end_L in_z_end_L in_alpha_L in_beta_L in_gamma_L]
    Out_L=[out_x_end_L out_y_end_L out_z_end_L out_alpha_L out_beta_L out_gamma_L]
    
    %確認FK 和IK誤差
%     if(out_x_end-in_x_end)>1e-5 || (out_y_end-in_y_end)>1e-5 || (out_z_end-in_z_end)>1e-5 || (out_alpha-in_alpha)>1e-5 || (out_beta-in_beta)>1e-5 || (out_gamma-in_gamma)>1e-5 
%         display('===============')
%         display('IK FK not match')
%         i
%         In=[in_x_end in_y_end in_z_end in_alpha*(180/pi) in_beta*(180/pi) in_gamma*(180/pi)]
%         Out=[out_x_end out_y_end out_z_end out_alpha*(180/pi) out_beta*(180/pi) out_gamma*(180/pi)]
%         
%         break;
%     end
    
    pause(0.1);
end

 %畫JointAngle
%  Draw_7DOF_JointAnglePath(PathTheta);