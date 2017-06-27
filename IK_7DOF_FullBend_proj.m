
%L1 上臂L型長邊
%L2 上臂L型短邊
%L3 上臂L型短邊
%L4 上臂L型長邊
%L5 end effector
%function theta = IK_7DOF_FullBend_proj(L0,L1,L2,L3,L4,L5,x_base,y_base,z_base,x_end,y_end,z_end,alpha,beta,gamma,Rednt_alpha)
%輸出參數
theta=zeros(1,7);

%求出H_hat_x
%R=R_z1x2z3(alpha,beta,gamma);
R=R_z1x2y3(alpha,beta,gamma);
V_H_hat_x=R(1:3,1);%取出歐拉角轉換的旋轉矩陣，取出第1行為X軸旋轉後向量
V_H_hat_x=V_H_hat_x/norm(V_H_hat_x);
V_H_hat_y=R(1:3,2);%取出歐拉角轉換的旋轉矩陣，取出第2行為Y軸旋轉後向量

%V_H_hat_y=V_H_hat_y/norm(V_H_hat_y);
V_r_end=[x_end-x_base;
         y_end-y_base;
         z_end-z_base];
V_r_h=L5*V_H_hat_x;
V_r_wst=V_r_end-V_r_h;

ru_norm=(L1^2+L2^2)^0.5; %L型的斜邊長度
rf_norm=(L3^2+L4^2)^0.5;

theta_tmp=acos((ru_norm^2 + rf_norm^2- norm(V_r_wst)^2) / (2*ru_norm*rf_norm));
theta(4)=2*pi-atan2(L1,L2)-atan2(L4,L3)-theta_tmp;
theta4_org=pi-theta_tmp;

V_r_m=(ru_norm^2-rf_norm^2+norm(V_r_wst)^2)/(2*norm(V_r_wst)^2)*V_r_wst;

%Redundant circle 半徑R
Rednt_cir_R=ru_norm^2-((ru_norm^2-rf_norm^2+norm(V_r_wst)^2)/(2*norm(V_r_wst)))^2;
Rednt_cir_R=Rednt_cir_R^0.5;

%圓中心點到Elbow向量 V_r_u
V_shx=[1;0;0];
V_shy=[0;1;0];
V_shz=[0;0;1];

V_alpha_hat=cross(V_r_wst,V_shz)/norm(cross(V_r_wst,V_shz));
V_beta_hat=cross(V_r_wst,V_alpha_hat)/norm(cross(V_r_wst,V_alpha_hat));

temp=Rogridues(Rednt_alpha,V_r_wst/norm(V_r_wst))*[Rednt_cir_R*V_beta_hat;1];  %Rednt_alpha的方向和論文上的方向性相反
V_R_u=temp(1:3,1);
V_r_u=V_r_m+V_R_u;

%旋轉 V_r_u  到V_ru_l1
V_r_f=V_r_wst-V_r_u;
Vn_u_f=cross(V_r_u,V_r_f)/norm(cross(V_r_u,V_r_f)); %ru 及 rf的法向量
theat_upoff=atan(L2/L1);
temp=Rogridues(-theat_upoff,Vn_u_f)*[V_r_u;1];  %旋轉 V_r_u  到V_ru_l1
V_ru_l1=temp(1:3,1);
V_ru_l1=V_ru_l1*L1/norm(V_ru_l1); %調整成L1長度

theta(1)=atan2(V_ru_l1(1),-V_ru_l1(3));%org

if theta(1) ~= 0
    theta(2)=atan2(-V_ru_l1(2),V_ru_l1(1)/sin(theta(1)));
else
    theta(2)=atan2(V_ru_l1(2),-V_ru_l1(3));   
end   


%% ==算第3軸
%看shy(V_r_u,V_r_f的法向量)經過1,2軸旋轉後  與V_r_u,V_r_f 需要第3軸轉多少

V_n_yrot12=Ry(-theta(1))*Rx(-theta(2))*[-V_shy;1];  %方向定義的關係 因此會差theta1 theta2多負號
V_n_yrot12=V_n_yrot12(1:3,1);

Vn_nuf_nyrot12=cross(Vn_u_f,V_n_yrot12);
Vn_nuf_nyrot12=Vn_nuf_nyrot12/norm(Vn_nuf_nyrot12);

temp=V_n_yrot12'*Vn_u_f/norm(V_n_yrot12)/norm(Vn_u_f); 

%Vn_u_f 和 V_n_yrot12的法向量   與 V_ru_l1同方向 theta(3)需要加負號
if norm(Vn_nuf_nyrot12 - V_ru_l1/norm(V_ru_l1)) < 1.e-7
    theta(3)=-acos(temp)
else
    theta(3)=acos(temp)
end
%theta(3)=-theta(3);%方向定義的關係 因此會差負號

%% ==算第5軸
theta(5)=0;
ct=(Vn_u_f'*V_r_wst-Vn_u_f'*V_r_end)/(norm(Vn_u_f)^2); %在V_r_u,V_r_f平面上，且經過V_r_end點的直線參數式的t

Vproj_end_ru_rf=V_r_end+ct*Vn_u_f;%V_r_end 沿著V_r_u,V_r_f平面法向量投影在平面上的點

V_wst_to_projend=Vproj_end_ru_rf-V_r_wst;
V_wst_to_end=V_r_end-V_r_wst;

V_rf_extend=L5*V_r_f/norm(V_r_f); 

%theta(7)=-acos(V_r_f'*V_wst_to_projend/(norm(V_r_f)*norm(V_wst_to_projend)))-(0.5*pi-atan(L4/L3)); 
%theta(6)=acos(abs(V_wst_to_projend'*V_wst_to_end)/norm(V_wst_to_projend)/norm(V_wst_to_end));


%旋轉V_r_f 到 V_rf_l4
theat_lowoff=atan(L3/L4);
temp=Rogridues(theat_lowoff,Vn_u_f)*[V_r_f;1];  %旋轉 V_r_f  V_rf_l4
V_rf_l4=temp(1:3,1);
V_rf_l4=V_rf_l4*L4/norm(V_rf_l4); %調整成L4長度

%V_n_rfl4 及V_n_rf形成的平面 的法向量
Vn_rfl4_nuf=cross(V_rf_l4,Vn_u_f)/norm(cross(V_rf_l4,Vn_u_f));
t_rfl4_nuf=(Vn_rfl4_nuf'*V_r_wst-Vn_rfl4_nuf'*V_r_end)/(norm(Vn_rfl4_nuf)^2); %V_n_rf,V_n_rfl4平面上，且經過V_r_end點的直線參數式的t 為rfl4_nuf
Vproj_end_rfl4_nuf=V_r_end+t_rfl4_nuf*Vn_rfl4_nuf;%V_r_end 沿著V_n_rfl4,V_n_rf平面法向量投影在平面上的點
V_wst_to_projend_rfl4_nuf=Vproj_end_rfl4_nuf-V_r_wst;

%% 算第6軸
%防止在acos(1.000000.....)的時候會出現虛部的情況
temp=V_rf_l4'*V_wst_to_projend_rfl4_nuf/norm(V_rf_l4)/norm(V_wst_to_projend_rfl4_nuf);
if abs(temp-1)<1.e-7 
   if temp >0
       temp=1;
   else
       temp=-1;
   end
end

%平面法向量 和 Vn_rfl4_nuf  同邊要加負號  判斷theta6要往上或往下轉
Vn_rfl4_WstToProjEndRfl4Nuf=cross(V_rf_l4/norm(V_rf_l4),V_wst_to_projend_rfl4_nuf/norm(V_wst_to_projend_rfl4_nuf));
Vn_rfl4_WstToProjEndRfl4Nuf=Vn_rfl4_WstToProjEndRfl4Nuf/norm(Vn_rfl4_WstToProjEndRfl4Nuf);
if norm(Vn_rfl4_WstToProjEndRfl4Nuf - Vn_rfl4_nuf) < 1.e-7
    theta(6)=-acos(temp); 
else
    theta(6)=acos(temp); 
end


%% 算第7軸
temp=Rogridues(-theta(6),Vn_rfl4_nuf)*[Vn_u_f;1]; 
Vn_nuf_rotx6_along_NRfl4Nuf=temp(1:3,1);%nuf 沿著 Vn_rfl4_nuf 旋轉第6軸角度得到投影點與目標點平面的法向量
Vn_nuf_rotx6_along_NRfl4Nuf=Vn_nuf_rotx6_along_NRfl4Nuf/norm(Vn_nuf_rotx6_along_NRfl4Nuf);
Vn_WstToEnd_WstToProjEndRfl4Nuf=cross(V_wst_to_end,V_wst_to_projend_rfl4_nuf);%V_wst_to_projend 和 V_wst_to_end的法向量
Vn_WstToEnd_WstToProjEndRfl4Nuf=Vn_WstToEnd_WstToProjEndRfl4Nuf/norm(Vn_WstToEnd_WstToProjEndRfl4Nuf);

%利用法向量方向 判斷theta7旋轉方向
temp=V_wst_to_projend_rfl4_nuf'*V_wst_to_end/norm(V_wst_to_projend_rfl4_nuf)/norm(V_wst_to_end);
if norm(Vn_WstToEnd_WstToProjEndRfl4Nuf - Vn_nuf_rotx6_along_NRfl4Nuf) < 1.e-7
    theta(7)=-acos(temp); 
else
    theta(7)=acos(temp); 
end
    

%% 測試redundant左右手角度，這兩軸要相反
%theta(3)=-theta(3);
%theta(6)=-theta(6);
%theta(6)=0
%theta(7)=0
%theta(7)=-(0.5*pi-atan(L4/L3));
%theta(6)=acos(abs(V_wst_to_projend'*V_wst_to_end)/norm(V_wst_to_projend)/norm(V_wst_to_end));

%theta(7)=0
%theta(7)=0;
%theta(7)=-acos(( L5^2+norm(V_wst_to_projend)^2-norm(V_rf_to_end)^2 ) /(%2*L5*norm(V_wst_to_projend)) ); %餘弦定理

%theta(6)=acos(norm(V_wst_to_projend) /L5);
% 
% V_r_f'*V_n_x_rot1234%test
% cross(V_r_f,V_n_x_rot1234)%test
% 
% figure(2)%test
% V_r_f_unit=V_r_f/norm(V_r_f);%test
% V_n_x_rot1234_unit=V_n_x_rot1234/norm(V_n_x_rot1234);%test
% 
% ddd=acos(V_r_f_unit'*V_n_x_rot1234_unit/(norm(V_r_f_unit)*norm(V_n_x_rot1234_unit)))
% plot([0,V_r_f_unit(1)],[0,V_r_f_unit(2)]);%test
% hold on;
% plot([0,V_n_x_rot1234_unit(1)],[0,V_n_x_rot1234_unit(2)]);%test
% axis('equal');%ekXYZ每一格的間距相等

% Vproj_rh_xy = proj_on_plan( V_r_h,V_shx,V_shy);
% Vproj_rf_xy = proj_on_plan( V_r_f,V_shx,V_shy);
% theta(7)=-acos(Vproj_rh_xy'*Vproj_rf_xy/(norm(Vproj_rh_xy)*norm(Vproj_rf_xy))); 
% 
% Vproj_rh_xz = proj_on_plan( V_r_h,V_shx,V_shz);
% Vproj_rf_xz = proj_on_plan( V_r_f,V_shx,V_shz);
% theta(6)=acos(Vproj_rh_xz'*Vproj_rf_xz/(norm(Vproj_rh_xz)*norm(Vproj_rf_xz))); 
% theta(6)=0;
% theta(7)=0;
%theta(7)=gamma;




%end
