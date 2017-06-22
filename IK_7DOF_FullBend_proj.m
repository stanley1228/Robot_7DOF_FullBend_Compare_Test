
%L1 �W�uL������
%L2 �W�uL���u��
%L3 �W�uL���u��
%L4 �W�uL������
%L5 end effector
function theta = IK_7DOF_FullBend_proj(L0,L1,L2,L3,L4,L5,x_base,y_base,z_base,x_end,y_end,z_end,alpha,beta,gamma,Rednt_alpha)
%��X�Ѽ�
theta=zeros(1,7);

%�D�XH_hat_x
%R=R_z1x2z3(alpha,beta,gamma);
R=R_z1x2y3(alpha,beta,gamma);
V_H_hat_x=R(1:3,1);%���X�کԨ��ഫ������x�}�A���X��1�欰X�b�����V�q
V_H_hat_x=V_H_hat_x/norm(V_H_hat_x);
V_H_hat_y=R(1:3,2);%���X�کԨ��ഫ������x�}�A���X��2�欰Y�b�����V�q

%V_H_hat_y=V_H_hat_y/norm(V_H_hat_y);
V_r_end=[x_end-x_base;
         y_end-y_base;
         z_end-z_base];
V_r_h=L5*V_H_hat_x;
V_r_wst=V_r_end-V_r_h;

ru_norm=(L1^2+L2^2)^0.5; %L�����������
rf_norm=(L3^2+L4^2)^0.5;

theta_tmp=acos((ru_norm^2 + rf_norm^2- norm(V_r_wst)^2) / (2*ru_norm*rf_norm));
theta(4)=2*pi-atan(L1/L2)-atan(L4/L3)-theta_tmp;
theta4_org=pi-theta_tmp;

V_r_m=(ru_norm^2-rf_norm^2+norm(V_r_wst)^2)/(2*norm(V_r_wst)^2)*V_r_wst;

%Redundant circle �b�|R
Rednt_cir_R=ru_norm^2-((ru_norm^2-rf_norm^2+norm(V_r_wst)^2)/(2*norm(V_r_wst)))^2;
Rednt_cir_R=Rednt_cir_R^0.5;

%�ꤤ���I��Elbow�V�q V_r_u
V_shx=[1;0;0];
V_shy=[0;1;0];
V_shz=[0;0;1];

V_alpha_hat=cross(V_r_wst,V_shz)/norm(cross(V_r_wst,V_shz));
V_beta_hat=cross(V_r_wst,V_alpha_hat)/norm(cross(V_r_wst,V_alpha_hat));

temp=Rogridues(Rednt_alpha,V_r_wst/norm(V_r_wst))*[Rednt_cir_R*V_beta_hat;1];  %Rednt_alpha����V�M�פ�W����V�ʬۤ�
V_R_u=temp(1:3,1);
V_r_u=V_r_m+V_R_u;

%���� V_r_u  ��V_ru_l1
V_r_f=V_r_wst-V_r_u;
Vn_u_f=cross(V_r_u,V_r_f)/norm(cross(V_r_u,V_r_f)); %ru �� rf���k�V�q
theat_upoff=atan(L2/L1);
temp=Rogridues(-theat_upoff,Vn_u_f)*[V_r_u;1];  %���� V_r_u  ��V_ru_l1
V_ru_l1=temp(1:3,1);
V_ru_l1=V_ru_l1*L1/norm(V_ru_l1); %�վ㦨L1����

theta(1)=atan2(V_ru_l1(1),-V_ru_l1(3));

if theta(1) ~= 0
    theta(2)=atan2(-V_ru_l1(2),V_ru_l1(1)/sin(theta(1)));
else
    theta(2)=atan2(-V_ru_l1(2),-V_ru_l1(3));
end   

%���3�b
%��shy(V_r_u,V_r_f���k�V�q)�g�L1,2�b�����  �PV_r_u,V_r_f �ݭn��3�b��h��
V_n_yrot12=Ry(-theta(1))*Rx(-theta(2))*[-V_shy;1];  %��V�w�q�����Y �]���|�ttheta1 theta2�h�t��
V_n_yrot12=V_n_yrot12(1:3,1);
theta(3)=acos(V_n_yrot12'*Vn_u_f/(norm(V_n_yrot12)*norm(Vn_u_f))); 
theta(3)=-theta(3);%��V�w�q�����Y �]���|�t�t��

theta(5)=0;

temp=Ry(-theta(1))*Rx(-theta(2))*Rz(-theta(3))*Ry(-theta4_org)*[V_shx;1]; 
V_n_x_rot1234=temp(1:3,1);

Vproj_rh_xy = proj_on_plan( V_r_h,V_shx,V_shy);
Vproj_rf_xy = proj_on_plan( V_r_f,V_shx,V_shy);
theta(6)=0.5*pi-acos(Vproj_rh_xy'*Vproj_rf_xy/(norm(Vproj_rh_xy)*norm(Vproj_rf_xy))); 

Vproj_rh_xz = proj_on_plan( V_r_h,V_shx,V_shz);
Vproj_rf_xz = proj_on_plan( V_r_f,V_shx,V_shz);
theta(7)=0.5*pi-acos(Vproj_rh_xz'*Vproj_rf_xz/(norm(Vproj_rh_xz)*norm(Vproj_rf_xz))); 
theta(6)=0;
theta(7)=0;
%theta(7)=gamma;




end
