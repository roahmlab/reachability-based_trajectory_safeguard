% offline processing to save the zonotope object in a more compact form
%1. Use a dictionary to refer to the actions that should be used
%2. Combine ERS with PRS so there is only one object to refer to during
%online
%UNDER CONSTRUCTION
clear
load zono_full_7.13_1spd.mat
M = containers.Map;
key = zeros(2,0);
% Build key table
% zono_peak_peak = AH.zono_full.res{v_ini_idx ,y_des_idx,h_ini_idx,del_idx,vd_idx,1};
%                                    2 options  always 1
for u0 = v_range
    for hi = h_range
        for delta = del_range
            M(char("u0="+num2str(u0)+"hi="+num2str(hi)+"delta="+num2str(delta)+"_tb"))= [];
            %     M(char("u0="+num2str(u0)+"_Aytb"))= [];
            M(char("u0="+num2str(u0)+"hi="+num2str(hi)+"delta="+num2str(delta)+"_zono"))= cell(0);
            for Au = v_des_range
                key = [key [u0; hi; delta; Au]];
                %2    2  7    2
            end
        end
    end
end
%%
for i = 1:length(v_range)
    for j = 1:length(h_range)
        for k = 1:length(del_range)
            for l = 1:length(v_des_range)
                %[~,y_des_idx] =min (abs(AH.y_array - (K_user(2)))); %always 1
                zono_peak =  res{i,1,j,k,l,1};
                zono_mid  =  res{i,1,j,k,l,2};
                zono_stop =  res{i,1,j,k,l,3};
                zono_err  = error_table{i ,1,j,k,l};
                zono_peak_mid_stop = [zono_peak; zono_mid;zono_stop];
                n = length(zono_peak_mid_stop);

                c = center(zono_peak{end}{end});
                dx = c(1);
                dy = c(2);
                M(char("u0="+num2str(v_range(i))+"hi="+num2str(h_range(j))+"delta="+num2str(del_range(k))+"_tb")) = [ M(char("u0="+num2str(v_range(i))+"hi="+num2str(h_range(j))+"delta="+num2str(del_range(k))+"_tb")) [v_des_range(l);dx;dy]];
                
                FRS =  cell(0);
                %create empty FRS obj
                for t_idx = 1: n
                    
                    zono_one = zono_peak_mid_stop{t_idx}{1};
                    center_bb = zono_err(1:2,t_idx);
                    h =  zono_err(3,t_idx);
                    gen = zono_err(4:5,t_idx);
                    len = gen(1);
                    width = gen(2);
                    ego_gen = [[cos(h)*len; sin(h)*len], [sin(-h)*width; cos(-h)*width]];
                    gen_err = zeros(7,2);gen_err(1:2,1:2) = ego_gen;
                    err_zono  = zonotope([[center_bb(1);center_bb(2);0;0;0;0;0], gen_err]);
                    %add err and prs together
                    FRS{end+1} = deleteAligned(zono_one + err_zono);
                end
                MAu = M(char("u0="+num2str(v_range(i))+"hi="+num2str(h_range(j))+"delta="+num2str(del_range(k))+"_zono"));
                MAu{end+1} = FRS;
                M(char("u0="+num2str(v_range(i))+"hi="+num2str(h_range(j))+"delta="+num2str(del_range(k))+"_zono")) = MAu;
            end
        end
    end
end
%%
save('zono_full_7.13_1spd_cleaned.mat','M')
% peak_idx = 65; %this is end of 3.25 sec before braking;
%     for j = 1:length(FRS) %Add footprint, in this case probably all 0
%         ft = deleteAligned(project(FRS{j},5));
%         headingc = center(ft);
%         gen = [2.4 1.1];
%         h =  headingc;
%         len = gen(1);
%         width = gen(2);
%         ego_gen = [[cos(h)*len; sin(h)*len], [sin(-h)*width; cos(-h)*width]];
%         gen = zeros(14,2);gen(1:2,1:2) = ego_gen;
%         err_zono  = zonotope([zeros(14,1), gen]);
%         FRS{j}= FRS{j}+err_zono;
%     end
    
    
% end
% 
