
clear all
close all
clc


%�T�w�Ѽ�
% L0=0;     %�i��n�R��
% L1=250;     %upper arm
% L2=250;   %forearm
% L3=150;      %length of end effector

x_base_R=0;   %����I
y_base_R=0;
z_base_R=0;

x_base_L=0;   %����I
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
 
 %���ͥ����O(60,0,0)  Q(80,-20,0)  R(80,-20,-20)  S(60,0,-20)
 %�⦹���|����90��
O_R=[500 -50 0];
Q_R=[500 -200 0];
R_R=[500 -200 -220];
S_R=[500 -50 -220];

% O_L=[500 50 0];
% Q_L=[500 200 0];
% R_L=[500 200 -220];
% S_L=[500 50 -220];

O_L=[500 -50 0];
Q_L=[500 -200 0];
R_L=[500 -200 -220];
S_L=[500 -50 -220];
 
Path_R=zeros(DEF_DESCRETE_POINT,3);%�W�e�����|�I
PathPoint_R=zeros(DEF_DESCRETE_POINT,3);%�O����ڤW���I�A�e�Ϩϥ�
PathTheta_R=zeros(DEF_DESCRETE_POINT,7);%�O���C�b���סA�e�Ϩϥ�
 
Path_L=zeros(DEF_DESCRETE_POINT,3);%�W�e�����|�I
PathPoint_L=zeros(DEF_DESCRETE_POINT,3);%�O����ڤW���I�A�e�Ϩϥ�
PathTheta_L=zeros(DEF_DESCRETE_POINT,7);%�O���C�b���סA�e�Ϩϥ�

%�e����ΰ�IK FK����
for t=1:1:DEF_DESCRETE_POINT
    if t<=25
        Path_R(t,1:3)=O_R+(Q_R-O_R)*t/25;
        Path_L(t,1:3)=O_L+(Q_L-O_L)*t/25;
    elseif t<=50
        Path_R(t,1:3)=Q_R+(R_R-Q_R)*(t-25)/25;
        Path_L(t,1:3)=Q_L+(R_L-Q_L)*(t-25)/25;
    elseif t<=75
        Path_R(t,1:3)=R_R+(S_R-R_R)*(t-50)/25;
        Path_L(t,1:3)=R_L+(S_L-R_L)*(t-50)/25;
    else 
        Path_R(t,1:3)=S_R+(O_R-S_R)*(t-75)/15;
        Path_L(t,1:3)=S_L+(O_L-S_L)*(t-75)/15;
    end
end

%�e���u
%  O=[0 0 -(L1+L2+L3)]; %��l�I
%  Q=[0 -20 -180]; %�Ĥ@�b���|������ౡ�p�����I���|
% 
%  for t=1:1:DEF_DESCRETE_POINT
%         Path(t,1:3)=O+(Q-O)*(t-1)/DEF_DESCRETE_POINT;
%  end

%for t=1:1:DEF_DESCRETE_POINT
 
    %��J�Ѽ�
    in_x_end_R=Path_R(t,1);
    in_y_end_R=Path_R(t,2);
    in_z_end_R=Path_R(t,3);
    
    in_x_end_L=Path_L(t,1);
    in_y_end_L=Path_L(t,2);
    in_z_end_L=Path_L(t,3);
   
    in_alpha_R=0*(pi/180);
    in_beta_R=0*(t/DEF_DESCRETE_POINT)*(pi/180);
    in_gamma_R=0*(t/DEF_DESCRETE_POINT)*(pi/180);
    
    in_alpha_L=0*(pi/180);
    in_beta_L=0*(t/DEF_DESCRETE_POINT)*(pi/180);
    in_gamma_L=0*(t/DEF_DESCRETE_POINT)*(pi/180);

    Rednt_alpha_R=(-90)*(pi/180);
    Rednt_alpha_L=(-90)*(pi/180);
    %��X�Ѽ� initial
    %theta=zeros(1,7);
    
    %
    %���I��min==>IK==>theta==>FK==>���I��mout
    %
    %inverse kinematic
    %y_base_R=-L0;%�k�� original
    
    L0=0;     %�Y��ӻH
    L1=250;   %L�� ����
    L2=50;    %L�� �u��
    L3=50;    %L�� �u��
    L4=250;   %L�� ���� 
    L5=150;   %��end-effector
   
    org_upper=(L1^2+L2^2)^0.5;
    org_fore=(L3^2+L4^2)^0.5;
    theta_R=IK_7DOF(org_upper,org_fore,L5,x_base_R,y_base_R,z_base_R,in_x_end_R,in_y_end_R,in_z_end_R,in_alpha_R,in_beta_R,in_gamma_R,Rednt_alpha_R);
  	
    %y_base_L=-L0;%���� fullbend
    theta_L=IK_7DOF_FullBend(L0,L1,L2,L3,L4,L5,x_base_L,y_base_L,z_base_L,in_x_end_L,in_y_end_L,in_z_end_L,in_alpha_L,in_beta_L,in_gamma_L,Rednt_alpha_L);
  
   
    %forward kinematic
    %theta=[0 0 0 0 0 0 0];
    [out_x_end_R,out_y_end_R,out_z_end_R,out_alpha_R,out_beta_R,out_gamma_R,P_R,RotationM_R] = FK_7DOF(L0,org_upper,org_fore,L5,x_base_R,y_base_R,z_base_R,theta_R);
    [out_x_end_L,out_y_end_L,out_z_end_L,out_alpha_L,out_beta_L,out_gamma_L,P_L,RotationM_L] = FK_7DOF_FullBend(L0,L1,L2,L3,L4,L5,x_base_L,y_base_L,z_base_L,theta_L);
    
    %�O�����|�W���I
    PathPoint_R(t,1:3)=[out_x_end_R out_y_end_R out_z_end_R];
    PathPoint_L(t,1:3)=[out_x_end_L out_y_end_L out_z_end_L];
    
    %�e���`�I��
    Draw_7DOF_compare_point(P_R,RotationM_R,PathPoint_R,P_L,RotationM_L,PathPoint_L);
   
    %�O���C�b�����ܤ�
    PathTheta_R(t,1:7)=theta_R*(180/pi);
    PathTheta_L(t,1:7)=theta_L*(180/pi);
    
    In_R=[in_x_end_R in_y_end_R in_z_end_R in_alpha_R in_beta_R in_gamma_R]
    Out_R=[out_x_end_R out_y_end_R out_z_end_R out_alpha_R out_beta_R out_gamma_R]
    
    In_L=[in_x_end_L in_y_end_L in_z_end_L in_alpha_L in_beta_L in_gamma_L]
    Out_L=[out_x_end_L out_y_end_L out_z_end_L out_alpha_L out_beta_L out_gamma_L]
    
    %�T�{FK �MIK�~�t
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
%end

 %�eJointAngle
%  Draw_7DOF_JointAnglePath(PathTheta);