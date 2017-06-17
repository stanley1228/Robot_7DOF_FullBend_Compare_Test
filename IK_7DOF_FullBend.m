%{
%�T�w�Ѽ�
L1=20; %upper arm
L2=20; %forearm
L3=5;  %length of end effector
x_base=0;   %����I
y_base=0;
z_base=0;
%��J�Ѽ�

x_end=10;
y_end=0;
z_end=0;
alpha=0*pi/180;
beta=0;
gamma=0;
Rednt_alpha=pi/2;

%��X�Ѽ�
theta=zeros(1,7);
 %}


%�Q�έ쥻IK_7DOF �p���1,7�b���צ���
%IK_7DOF ���4�b���ɭԤ@�w�n�O�T���ΡA�p�G��y�ܦ������A�Ѻ�|���~
%�]�N�O�s��theat4�@�w�n�j��180-atan(L1/L2)-atan(L1/L2)
function theta = IK_7DOF_FullBend( L0,L1,L2,L3,L4,L5,x_base,y_base,z_base,x_end,y_end,z_end,alpha,beta,gamma,Rednt_alpha)

  DEF_RAD_2_DEG =(180/pi);
 
 theta = IK_7DOF( (L1^2+L2^2)^0.5,(L3^2+L4^2)^0.5,L5,x_base,y_base,z_base,x_end,y_end,z_end,alpha,beta,gamma,Rednt_alpha);
 
 %���հ��F��V�ۤϰ��D ���Ftheta6 ��Lz�b��V���M�l���Ǫ���V�ۤ�
 theta=-theta;
 theta(6)=-theta(6);
 
 %theta(1)=theta(1)-atan(L2/L1);
 %theta(7)=theta(7)-atan(L3/L4);
 
 theta(4)=2*pi-atan(L1/L2)-atan(L4/L3)-(pi-theta(4));%theta(4)�ѥX�Ӫ����Ӭ��t�ȥ[�t���ܥ���
 %theta(4)=2*pi-atan(L1/L2)-atan(L4/L3)-(pi+theta(4));%theta(4)�ѥX�Ӫ����Ӭ��t��
 %oldIK_theta_DEG=theta*DEF_RAD_2_DEG
%  aaaa=atan(L2/L1)*DEF_RAD_2_DEG



%  theta(1)=-theta(1)-atan(L2/L1);
%  theta(7)=-theta(7)-atan(L3/L4);
%  
%  theta(4)=2*pi-atan(L1/L2)-atan(L4/L3)-(pi+theta(4));%theta(4)�ѥX�Ӫ����Ӭ��t��
end
