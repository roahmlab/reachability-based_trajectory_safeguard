num_v_des = 6;
connections = cell(size(res,1),size(res,2),size(res,3),num_v_des);
for i = 1:size(res,1)
    for j = 1:size(res,2)
        [i j h vd_idx]
        for h = 1:size(res,3)
            for vd_idx = 1:num_v_des
                
                %each res{i,j}, is the end of this smaller than the beginning of
                %res{ii,jj}
                genendvd = v_range(vd_idx) %+ 2*(vd_idx-2);
%                 % 22 + 2*(1- 2);
%                 if genendvd < 22 || genendvd > 32
%                     continue;
%                 end
                genend=generators(box( res{i,j,h,vd_idx}{end}{1}));
                genendv=genend(4,4);
                genendh=genend(3,3);
                genendd=genend(5,5);
                cend = center( res{i,j,h,vd_idx}{end}{1});
                v_minend = cend(4) - genendv;
                v_maxend= cend(4) +genendv;
                h_minend = cend(3) - genendh;
                h_maxend= cend(3) +genendh;
                d_minend = cend(5) - genendd;
                d_maxend= cend(5) +genendd;
                
                con = [];
                %         figure(1);clf; hold on;
                for ii = 1:size(res,1)
                    for jj = 1:size(res,2)
                        for hh = 1:size(res,3)
                            for vdd = 1:num_v_des
                                 genbeginvd = v_range(ii) + 2*(vdd-2);
                                % 22 + 2*(1- 2);
                                if genbeginvd < 22 || genbeginvd > 32
                                    continue;
                                end
                                %[0.3 0.8 1.3] --> [0.5 0.8 1]
                                %   vmin1   vmax1  vminend  vmaxend
                                genbegin = generators(box(res{ii,jj,hh,vdd}{1}{1}));
                                genbeginv = genbegin(4,4);
                                genbeginh = genbegin(3,3);
                                genbegind = genbegin(5,5);
                                cbegin = center( res{ii,jj,hh,vdd}{1}{1});
                                v_min1 = cbegin(4) - genbeginv;
                                v_max1 = cbegin(4) + genbeginv;
                                h_min1 = cbegin(3) - genbeginh;
                                h_max1 = cbegin(3) + genbeginh;
                                d_min1 = cbegin(5) - genbegind;
                                d_max1 = cbegin(5) + genbegind;
                                if v_minend > v_min1 && v_maxend < v_max1 && h_minend > h_min1&& h_maxend < h_max1%don't worry about delta lol &&  d_minend > d_min1&& d_maxend < d_max1
                                    con = [con; ii jj hh vdd];
                                end
                            end
                        end
                    end
                end
                connections{i,j,h,vd_idx} = con;
            end
        end
    end
end

save(['connections.mat'], 'connections','-v7.3');
%%
load Zonotope_simple_3sec.mat
num_v_des = 6;
zono_poly = cell(size(res_simple,1),size(res_simple,2),size(res_simple,3),num_v_des,10);
for ii = 1:size(res_simple,1)
    for jj = 1:size(res_simple,2)
        for hh = 1:size(res_simple,3)
            for vdd = 1:num_v_des
%             zono_poly{ii,jj,hh,vdd}  = polyshape();
            for kk = 1:length(res_simple{ii,jj,hh,vdd})
                zono_poly{ii,jj,hh,vdd,kk}  = polyshape(res_simple{ii,jj,hh,vdd}{kk}(1,:),res_simple{ii,jj,hh,vdd}{kk}(2,:));
            end
%             plot(zono_poly{ii,jj,hh,vdd})
            end
        end
    end
end

save(['zonopolygons4.7.mat'], 'zono_poly','-v7.3');
% for ii = 1:size(res_simple,1)
%     for jj = 1:size(res_simple,2)
%         zono_poly{ii,jj,hh,vdd}  = zono_poly{ii,jj,hh,vdd}.Vertices;
%     end
% end

% i = 1; j = 1
% b=box(res{i,j}{1}{1})
% b2=box(res{i,j}{end}{1})
% center(res{i,j}{end}{1})
%containsPoint(Z,p)
%         p=plot(res{i,j}{1}{1},[4,3],'g')
%         p=plot(res{i,j}{end}{1},[4,3],'r')

% connections{1,1} = [1 2; 3 5; 6 1]

% res_simple = cell(size(res,1),size(res,2));
% for i = 1:size(res,1)
%     for j = 1:size(res,2)
%         res_simple(i,j) = box(res{i,j}{end}{1});
%     end
% end
%UNDER CONSTRCTION



