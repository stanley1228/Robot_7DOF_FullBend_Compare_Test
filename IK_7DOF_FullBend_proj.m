
%L1 上臂L型長邊
%L2 上臂L型短邊
%L3 上臂L型長邊
%L4 上臂L型短邊
%L5 end effector
function theta = IK_7DOF_FullBend_proj(L0,L1,L2,L3,L4,L5,x_base,y_base,z_base,x_end,y_end,z_end,alpha,beta,gamma,Rednt_alpha)
%輸出參數
theta=zeros(1,7);

%求出H_hat_x
%R=R_z1x2z3(alpha,beta,gamma);
R=R_z1x2y3(alpha,beta,gamma);
V_H_hat_x=R(1:3,1);%取出歐拉角轉換的旋轉矩陣，取出第1行為X軸旋轉後向量
V_H_hat_x=V_H_hat_x/norm(V_H_hat_x);
V_H_hat_y=R(1:3,2);%取出歐拉角轉換的旋轉矩陣，取出第2行為Y軸旋轉後向量

V_H_hat_y=V_H_hat_y/norm(V_H_hat_y);
V_r_end=[x_end-x_base;
         y_end-y_base;
         z_end-z_base];
V_r_h=L3*V_H_hat_x;
V_r_wst=V_r_end-V_r_h;

ru_norm=(L1^2+L2^2)^0.5; %L型的鞋邊長度
rf_norm=(L3^2+L4^2)^0.5;

theta_tmp=acos(ru_norm^2 + rf_norm^2- norm(V_r_wst)^2)/(2*ru_norm*rf_norm);
theta(4)=2*pi-atan(L1/L2)-atan(L4/L3)-theta_tmp;

V_r_m=(ru_norm^2-rf_norm^2+norm(V_r_wst)^2)/(2*norm(V_r_wst)^2)*V_r_wst;

%Redundant circle 半徑R
Rednt_cir_R=ru_norm^2-((ru_norm^2-rf_norm^2+norm(V_r_wst)^2)/(2*norm(V_r_wst)))^2;
Rednt_cir_R=Rednt_cir_R^0.5;

%圓中心點到Elbow向量 V_r_u
V_shz=[0;0;1];
V_alpha_hat=cross(V_r_wst,V_shz)/norm(cross(V_r_wst,V_shz));
V_beta_hat=cross(V_r_wst,V_alpha_hat)/norm(cross(V_r_wst,V_alpha_hat));

temp=Rogridues(Rednt_alpha,V_r_wst/norm(V_r_wst))*[Rednt_cir_R*V_beta_hat;1];  %Rednt_alpha的方向和論文上的方向性相反
V_R_u=temp(1:3,1);
V_r_u=V_r_m+V_R_u;


V_r_f=V_r_wst-V_r_u;

Vn_u_f=cross(V_r_u,V_r_f)/norm(cross(V_r_u,V_r_f)); %ru 及 rf的法向量

theat_upoff=atan(L2/L1);


temp=Rogridues(theat_upoff,Vn_u_f)*[V_r_u;1];  %Rednt_alpha的方向和論文上的方向性相反
V_ru_l1=temp(1:3,1);

V_shx=[1;0;0];
V_shy=[0;1;0];
%proj_shx_V_ru_l=V_ru_l'*V_shx/(V_shx'*V_shx)
%proj_shz_V_ru_l=V_ru_l'*V_shz/(V_shz'*V_shz)

Vproj_xz_ru_l1 =proj_on_plan(V_ru_l1,V_shx,V_shz); %vector向量對 ort_plan_axis1 及ort_plan_axis1 兩個正交向量平面做投影的向量
theta(1)=acos( Vproj_xz_ru_l1'*(-V_shz) / (norm(Vproj_xz_ru_l1)*norm(V_shz)) );
  

Vproj_yz_ru_l1=proj_on_plan(V_ru_l1,V_shy,V_shz);
theta(2)=acos( Vproj_yz_ru_l1'*(-V_shz) / (norm(Vproj_yz_ru_l1)*norm(V_shz)) );

V_n_rot12=Ry(-theta(1)) * Rx(theta(2))* [-V_shy;1];
V_n_rot12=V_n_rot12(1:3,1);
theta(3)=acos(V_n_rot12'*Vn_u_f/(norm(V_n_rot12)*norm(Vn_u_f)));


theta(5)=alpha;
theta(6)=beta;
theta(7)=gamma;

end
