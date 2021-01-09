function [c, ceq]=eval_zono_highway_cons(AH,K,A_con,b_con,s_con )
%we need max(Ax -b)> 0, since fmincon requires nonlin<=0
%we specify constraint as min(b-Ax)<=0
 k1c = AH.v_array(vd_idx);k1g = AH.zono_full.kvdg;
            k2c = AH.y_array(y_des_idx);k2g = AH.zono_full.kyg;
            c_k = [k1c; k2c];
            g_k = [k1g; k2g];
    epsilon = 1e-3;
    ceq = [];

    c = A_con*x + b_con;
    c = reshape(c, 6, []);
    c = min(c, [], 1) + epsilon;
end