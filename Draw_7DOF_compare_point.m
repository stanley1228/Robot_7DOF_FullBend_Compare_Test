function r= Draw_7DOF_compare_point(P_R,RotationM_R,PathPoint_R,P_L,RotationM_L,PathPoint_L,Vproj_end_ru_rf,V_rf_extend,Vn_u_f,V_rf_l4,Vn_rfl4_nuf,Vproj_end_rfl4_nuf)
%DRAW_7DOF_POINT Summary of this function goes here
%   Detailed explanation goes here\
%�e�Dfullbend ��fullbend ���ϰ���� 
figure(1)
cla reset


%%���դ��u�V�� ��������
%  AZ=0;
%  EL=0;
% 
% xlim([-100 100]) % ���� X �b�d�� 
% ylim([-100 50]) % ���� Y �b�d�� 
% zlim([-250 20]) % ���� Z �b�d�� 


%Path���ե�  L0=0; L1=100; L2=100; L3=10; 

%  AZ=-50;
%  EL=40;
%  
% xlim([-100 100]) % ���� X �b�d�� 
% ylim([-100 50]) % ���� Y �b�d�� 
% zlim([-60 20]) % ���� Z �b�d�� 


%Path���ե� L0=255; L1=250; L2=250; L3=150;
%  AZ=-180;
%  EL=0;
 AZ=-63;
 EL=28;

%  
% xlim([-200 550]) % ���� X �b�d�� 
% ylim([-500 500]) % ���� Y �b�d�� 
% zlim([-400 200]) % ���� Z �b�d�� 

view(AZ,EL);


xlabel('x');
ylabel('y');

hold on;    grid on;    box on; rotate3d on ;

%% ========�e�s��======== %%
%Right Arm
for i=1:1:8
    if i==1
        plot3([0,P_R(i,1)],[0,P_R(i,2)],[0,P_R(i,3)],'-r','LineWidth',2); %��y�e��Joint1
    else
        plot3([P_R(i-1,1),P_R(i,1)],[P_R(i-1,2),P_R(i,2)],[P_R(i-1,3),P_R(i,3)],'-r','LineWidth',2);
    end
end

%Left Arm
for i=1:1:10
    if i==1
        plot3([0,P_L(i,1)],[0,P_L(i,2)],[0,P_L(i,3)],'-k','LineWidth',2); %��y�e��Joint1
    else
        plot3([P_L(i-1,1),P_L(i,1)],[P_L(i-1,2),P_L(i,2)],[P_L(i-1,3),P_L(i,3)],'-k','LineWidth',2);
    end
end

plot3([P_L(5,1),P_L(7,1)],[P_L(5,2),P_L(7,2)],[P_L(5,3),P_L(7,3)],'--k','LineWidth',2); %��y�e��Joint1  %test

plot3([P_L(3,1),P_L(5,1)],[P_L(3,2),P_L(5,2)],[P_L(3,3),P_L(5,3)],'--k','LineWidth',2); %��y�e��Joint1  %test

%% ========�e���I======== %%
%�e���I
plot3(0,0,0,'ro','MarkerSize',10,'Linewidth',4);text(0,0,0,'Org')

%% ========�e�C�b���`�I======== %%
%Right Arm
for i=1:1:7
  plot3(P_R(i,1),P_R(i,2),P_R(i,3),'bo','MarkerSize',5,'Linewidth',4);
  
  %�Х�
  if i==2
       text(P_R(i,1),P_R(i,2),P_R(i,3),'shoulder');
  elseif i==4
      text(P_R(i,1),P_R(i,2),P_R(i,3),'Elbow');
  elseif i==6
      text(P_R(i,1),P_R(i,2),P_R(i,3),'Wst');   
  end
end

%Left Arm
for i=1:1:10
  plot3(P_L(i,1),P_L(i,2),P_L(i,3),'bo','MarkerSize',5,'Linewidth',4);
  
  %�Х�
  if i==2
       text(P_L(i,1),P_L(i,2),P_L(i,3),'shoulder');
   elseif i==5
      text(P_L(i,1),P_L(i,2),P_L(i,3),'Elbow');
  elseif i==7
      text(P_L(i,1),P_L(i,2),P_L(i,3),'Wst');   
  end
end

%%==�eru rf���W��v�I==%%
%plot3([P_L(7,1),Vproj_end_ru_rf(1)],[P_L(7,2),Vproj_end_ru_rf(2)],[P_L(7,3),Vproj_end_ru_rf(3)],'--g','MarkerSize',2,'Linewidth',1);

%%==�erfl4 nuf���W��v�I==%%
plot3([P_L(7,1),Vproj_end_rfl4_nuf(1)],[P_L(7,2),Vproj_end_rfl4_nuf(2)],[P_L(7,3),Vproj_end_rfl4_nuf(3)],'--g','MarkerSize',2,'Linewidth',1);

%%==�erf�����I==%%
%plot3([P_L(7,1),P_L(7,1)+V_rf_extend(1)],[P_L(7,2),P_L(7,2)+V_rf_extend(2)],[P_L(7,3),P_L(7,3)+V_rf_extend(3)],'-b','MarkerSize',2,'Linewidth',1);

%%==�eV_rf_l4�����I==%%
V_rf_l4_Length=V_rf_l4*100;
plot3([P_L(7,1),P_L(7,1)+V_rf_l4_Length(1)],[P_L(7,2),P_L(7,2)+V_rf_l4_Length(2)],[P_L(7,3),P_L(7,3)+V_rf_l4_Length(3)],'--b','MarkerSize',2,'Linewidth',1);

%%==�eVn_rf  ��Vr_l11�������k�V�q==%%
Vn_rfl4_nuf_Length=Vn_rfl4_nuf*200;
plot3([P_L(7,1),P_L(7,1)+Vn_rfl4_nuf_Length(1)],[P_L(7,2),P_L(7,2)+Vn_rfl4_nuf_Length(2)],[P_L(7,3),P_L(7,3)+Vn_rfl4_nuf_Length(3)],'--b','MarkerSize',2,'Linewidth',1);

%%==�eVn_u_f
Vn_u_f=Vn_u_f*100/norm(Vn_u_f);
plot3([P_L(5,1),P_L(5,1)+Vn_u_f(1)],[P_L(5,2),P_L(5,2)+Vn_u_f(2)],[P_L(5,3),P_L(5,3)+Vn_u_f(3)],'--b','MarkerSize',2,'Linewidth',1);
plot3([P_L(7,1),P_L(7,1)+Vn_u_f(1)],[P_L(7,2),P_L(7,2)+Vn_u_f(2)],[P_L(7,3),P_L(7,3)+Vn_u_f(3)],'--b','MarkerSize',2,'Linewidth',1);

%%==�e���I�M�k�V�q Vn_u_f
%plot3([P_R(8,1),P_R(8,1)+Vn_u_f(1)],[P_R(8,2),P_R(8,2)+Vn_u_f(2)],[P_R(8,3),P_R(8,3)+Vn_u_f(3)],'--m','MarkerSize',2,'Linewidth',1);

%%==�e���I�M�k�V�q Vn_rfl4_nuf
Vn_rfl4_nuf_Length=Vn_rfl4_nuf*200;
plot3([P_R(8,1),P_R(8,1)-Vn_rfl4_nuf_Length(1)],[P_R(8,2),P_R(8,2)-Vn_rfl4_nuf_Length(2)],[P_R(8,3),P_R(8,3)-Vn_rfl4_nuf_Length(3)],'--m','MarkerSize',2,'Linewidth',1);

%%  ========�e���|�W���I======== %%
%Right Arm
plot3(PathPoint_R(:,1),PathPoint_R(:,2),PathPoint_R(:,3),'mo','MarkerSize',2,'Linewidth',1);

%Left Arm
plot3(PathPoint_L(:,1),PathPoint_L(:,2),PathPoint_L(:,3),'mo','MarkerSize',2,'Linewidth',1);

%% ========End effector======== %%
%Right Arm
plot3(P_R(8,1),P_R(8,2),P_R(8,3),'go','MarkerSize',10,'Linewidth',4);text(P_R(8,1),P_R(8,2),P_R(8,3),'R_End');
%Left Arm
plot3(P_L(10,1),P_L(10,2),P_L(10,3),'go','MarkerSize',10,'Linewidth',4);text(P_L(10,1),P_L(10,2),P_L(10,3),'L_End');

%% ========���I���A�y�жb�Х�  orientation V_H_hat_x V_H_hat_y V_H_hat_z ========%%
%Right Arm
V_H_HAT_UNIT_LEN=10;
RotationM_R=RotationM_R*V_H_HAT_UNIT_LEN;
V_H_hat_x=RotationM_R(1:3,1);
V_H_hat_y=RotationM_R(1:3,2);
V_H_hat_z=RotationM_R(1:3,3);
plot3([P_R(8,1),P_R(8,1)+V_H_hat_x(1,1)],[P_R(8,2),P_R(8,2)+V_H_hat_x(2,1)],[P_R(8,3),P_R(8,3)+V_H_hat_x(3,1)],'-m','LineWidth',2); text(P_R(8,1)+V_H_hat_x(1,1),P_R(8,2)+V_H_hat_x(2,1),P_R(8,3)+V_H_hat_x(3,1),'X')
plot3([P_R(8,1),P_R(8,1)+V_H_hat_y(1,1)],[P_R(8,2),P_R(8,2)+V_H_hat_y(2,1)],[P_R(8,3),P_R(8,3)+V_H_hat_y(3,1)],'-g','LineWidth',2); text(P_R(8,1)+V_H_hat_y(1,1),P_R(8,2)+V_H_hat_y(2,1),P_R(8,3)+V_H_hat_y(3,1),'Y')
plot3([P_R(8,1),P_R(8,1)+V_H_hat_z(1,1)],[P_R(8,2),P_R(8,2)+V_H_hat_z(2,1)],[P_R(8,3),P_R(8,3)+V_H_hat_z(3,1)],'-b','LineWidth',2); text(P_R(8,1)+V_H_hat_z(1,1),P_R(8,2)+V_H_hat_z(2,1),P_R(8,3)+V_H_hat_z(3,1),'Z')

%Left Arm
V_H_HAT_UNIT_LEN=10;
RotationM_L=RotationM_L*V_H_HAT_UNIT_LEN;
V_H_hat_x=RotationM_L(1:3,1);
V_H_hat_y=RotationM_L(1:3,2);
V_H_hat_z=RotationM_L(1:3,3);
plot3([P_L(10,1),P_L(10,1)+V_H_hat_x(1,1)],[P_L(10,2),P_L(10,2)+V_H_hat_x(2,1)],[P_L(10,3),P_L(10,3)+V_H_hat_x(3,1)],'-m','LineWidth',2); text(P_L(10,1)+V_H_hat_x(1,1),P_L(10,2)+V_H_hat_x(2,1),P_L(10,3)+V_H_hat_x(3,1),'X')
plot3([P_L(10,1),P_L(10,1)+V_H_hat_y(1,1)],[P_L(10,2),P_L(10,2)+V_H_hat_y(2,1)],[P_L(10,3),P_L(10,3)+V_H_hat_y(3,1)],'-g','LineWidth',2); text(P_L(10,1)+V_H_hat_y(1,1),P_L(10,2)+V_H_hat_y(2,1),P_L(10,3)+V_H_hat_y(3,1),'Y')
plot3([P_L(10,1),P_L(10,1)+V_H_hat_z(1,1)],[P_L(10,2),P_L(10,2)+V_H_hat_z(2,1)],[P_L(10,3),P_L(10,3)+V_H_hat_z(3,1)],'-b','LineWidth',2); text(P_L(10,1)+V_H_hat_z(1,1),P_L(10,2)+V_H_hat_z(2,1),P_L(10,3)+V_H_hat_z(3,1),'Z')


r=0;


axis('equal');%ekXYZ�C�@�檺���Z�۵�
 
xlim([-200 650]) % ���� X �b�d�� 
ylim([-500 500]) % ���� Y �b�d�� 
zlim([-400 200]) % ���� Z �b�d�� 


end

