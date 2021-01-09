function dz = closed_loop_highwaySS(tdummy,z,udummy)
%             x y h vx  delta vx_des y_des
% options.x0 = [0;0;0;kvc;0;    kvc;   kyc];
            % extract the states
            y = z(2);
            psi =  z(3)+z(8);
            vx  =  z(4)+z(9);
            delta = z(5)+z(10);
            
            vx_nom = z(6); % this is the output from network,network output discrete, converted to num by planner
            y_nom =  z(7); % output of network as well
            

            Tc = 2;
            m = 1558;
            lf = 1.462884;
            lr = 1.405516;
            l = lf + lr;

            Cf = 1.432394487827058e+05;
            Cr = 2.214094963126969e+05;
            g = 9.80655;
            % get feedback control inputs
            
            ua = 5*(vx_nom-vx);
            delta_cmd = 1*(y_nom-y)-0.1*psi-0.1*delta;
            delta_dot = 5*(delta_cmd-delta);


            kus = (m*g*lr/(l*Cf)-m*g*lf/(l*Cr));

            yr = tan(delta)*vx/(l+kus*vx^2/g);

            vy = yr*(lr-m*vx^2*lr/(Cr*l));
      

            % saturate the inputs % no saturation, todo: warning for
            % saturation
%             delta_dot = bound_values(delta_dot,A.max_del_dot) ;
%             ua = bound_values(ua,A.max_accel) ;
%             
            % calculate the derivatives
            dz = [vx*cos(psi)-vy*sin(psi);...
                     vx*sin(psi)+vy*cos(psi);...
                     yr;...
                     -Tc/m*vx+Tc*ua;delta_dot;0;0;0;0;0;1];
  %            x y h vx  delta vx_des y_des
%options.x0 = [0;0;0;kvc;0;    kvc;   kyc];

end