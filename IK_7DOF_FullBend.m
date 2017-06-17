%{
%固定參數
L1=20; %upper arm
L2=20; %forearm
L3=5;  %length of end effector
x_base=0;   %基準點
y_base=0;
z_base=0;
%輸入參數

x_end=10;
y_end=0;
z_end=0;
alpha=0*pi/180;
beta=0;
gamma=0;
Rednt_alpha=pi/2;

%輸出參數
theta=zeros(1,7);
 %}


%利用原本IK_7DOF 計算後1,7軸角度扣除
%IK_7DOF 算第4軸的時候一定要是三角形，如果手肘變成往後折，解算會錯誤
%也就是新的theat4一定要大於180-atan(L1/L2)-atan(L1/L2)
function theta = IK_7DOF_FullBend( L0,L1,L2,L3,L4,L5,x_base,y_base,z_base,x_end,y_end,z_end,alpha,beta,gamma,Rednt_alpha)

  DEF_RAD_2_DEG =(180/pi);
 
 theta = IK_7DOF( (L1^2+L2^2)^0.5,(L3^2+L4^2)^0.5,L5,x_base,y_base,z_base,x_end,y_end,z_end,alpha,beta,gamma,Rednt_alpha);
 
 %測試馬達方向相反問題 除了theta6 其他z軸方向都和衍祥學長方向相反
 theta=-theta;
 theta(6)=-theta(6);
 
 %theta(1)=theta(1)-atan(L2/L1);
 %theta(7)=theta(7)-atan(L3/L4);
 
 theta(4)=2*pi-atan(L1/L2)-atan(L4/L3)-(pi-theta(4));%theta(4)解出來的應該為負值加負號變正值
 %theta(4)=2*pi-atan(L1/L2)-atan(L4/L3)-(pi+theta(4));%theta(4)解出來的應該為負值
 %oldIK_theta_DEG=theta*DEF_RAD_2_DEG
%  aaaa=atan(L2/L1)*DEF_RAD_2_DEG



%  theta(1)=-theta(1)-atan(L2/L1);
%  theta(7)=-theta(7)-atan(L3/L4);
%  
%  theta(4)=2*pi-atan(L1/L2)-atan(L4/L3)-(pi+theta(4));%theta(4)解出來的應該為負值
end
