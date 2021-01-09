function q = rotm2quat_sym(R)
    R11 = R(1,1) ;
    R12 = R(1,2) ;
    R13 = R(1,3) ;
    R21 = R(2,1) ;
    R22 = R(2,2) ;
    R23 = R(2,3) ;
    R31 = R(3,1) ;
    R32 = R(3,2) ;
    R33 = R(3,3) ;
    
    qw = sqrt(1 + R11 + R22 + R33)/2 ;
    qx = (R31 - R13)/(4*qw) ;
    qy = (R13 - R31)/(4*qw) ;
    qz = (R21 - R12)/(4*qw) ;
    
    q = [qw qx qy qz] ;
end