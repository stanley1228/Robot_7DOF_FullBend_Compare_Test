 %vector�V�q�� ort_plan_axis1 ��ort_plan_axis1 ��ӥ���V�q��������v���V�q
function ans_vec = proj_on_plan( vecor,ort_plan_axis1,ort_plan_axis2 )

ans_vec=vecor'*ort_plan_axis1/(ort_plan_axis1'*ort_plan_axis1)*ort_plan_axis1+ ...
        vecor'*ort_plan_axis2/(ort_plan_axis2'*ort_plan_axis2)*ort_plan_axis2;
    

end

