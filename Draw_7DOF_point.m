function r= Draw_7DOF_point(P,RotationM,PathPoint)
%DRAW_7DOF_POINT Summary of this function goes here
%   Detailed explanation goes here\

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
 AZ=-50;
 EL=40;
 
xlim([-200 550]) % ���� X �b�d�� 
ylim([-500 500]) % ���� Y �b�d�� 
zlim([-400 200]) % ���� Z �b�d�� 

view(AZ,EL);


xlabel('x');
ylabel('y');

hold on;    grid on;    box on; rotate3d on ;

%�e�s��
for i=1:1:8
    if i==1
        plot3([0,P(i,1)],[0,P(i,2)],[0,P(i,3)],'-r','LineWidth',2); %��y�e��Joint1
    else
        plot3([P(i-1,1),P(i,1)],[P(i-1,2),P(i,2)],[P(i-1,3),P(i,3)],'-r','LineWidth',2);
    end
end

%�e���I
plot3(0,0,0,'ro','MarkerSize',10,'Linewidth',4);text(0,0,0,'shoulder')

%�e�C�b���`�I
for i=1:1:7
  plot3(P(i,1),P(i,2),P(i,3),'bo','MarkerSize',5,'Linewidth',4);
  
  %�Х�
  if i==4
      text(P(i,1),P(i,2),P(i,3),'Elbow');
  elseif i==6
      text(P(i,1),P(i,2),P(i,3),'Wst');   
  end
end

%�e���|�W���I
plot3(PathPoint(:,1),PathPoint(:,2),PathPoint(:,3),'mo','MarkerSize',2,'Linewidth',1);

%End effector
plot3(P(8,1),P(8,2),P(8,3),'go','MarkerSize',10,'Linewidth',4);text(P(8,1),P(8,2),P(8,3),'End');

%���I���A�y�жb�Х�  orientation V_H_hat_x V_H_hat_y V_H_hat_z 
V_H_HAT_UNIT_LEN=10;
RotationM=RotationM*V_H_HAT_UNIT_LEN;
V_H_hat_x=RotationM(1:3,1);
V_H_hat_y=RotationM(1:3,2);
V_H_hat_z=RotationM(1:3,3);
plot3([P(8,1),P(8,1)+V_H_hat_x(1,1)],[P(8,2),P(8,2)+V_H_hat_x(2,1)],[P(8,3),P(8,3)+V_H_hat_x(3,1)],'-m','LineWidth',2); text(P(8,1)+V_H_hat_x(1,1),P(8,2)+V_H_hat_x(2,1),P(8,3)+V_H_hat_x(3,1),'X')
plot3([P(8,1),P(8,1)+V_H_hat_y(1,1)],[P(8,2),P(8,2)+V_H_hat_y(2,1)],[P(8,3),P(8,3)+V_H_hat_y(3,1)],'-g','LineWidth',2); text(P(8,1)+V_H_hat_y(1,1),P(8,2)+V_H_hat_y(2,1),P(8,3)+V_H_hat_y(3,1),'Y')
plot3([P(8,1),P(8,1)+V_H_hat_z(1,1)],[P(8,2),P(8,2)+V_H_hat_z(2,1)],[P(8,3),P(8,3)+V_H_hat_z(3,1)],'-b','LineWidth',2); text(P(8,1)+V_H_hat_z(1,1),P(8,2)+V_H_hat_z(2,1),P(8,3)+V_H_hat_z(3,1),'Z')

r=0;


% axis('equal');%ekXYZ�C�@�檺���Z�۵�
% xlim([-150 250]) % ���� X �b�d�� 
% ylim([-250 100]) % ���� y �b�d�� 
% zlim([-300 30]) % ���� z �b�d�� 


end
