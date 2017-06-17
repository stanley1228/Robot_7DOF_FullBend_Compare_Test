 %vector向量對 ort_plan_axis1 及ort_plan_axis1 兩個正交向量平面做投影的向量
function ans_vec = proj_on_plan( vecor,ort_plan_axis1,ort_plan_axis2 )

ans_vec=vecor'*ort_plan_axis1/(ort_plan_axis1'*ort_plan_axis1)*ort_plan_axis1+ ...
        vecor'*ort_plan_axis2/(ort_plan_axis2'*ort_plan_axis2)*ort_plan_axis2;
    

end

