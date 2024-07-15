
%第三步： 该步骤是对于曲线型特征分割点云进行点云规则化处理
%输入的是第二步的Curve_sort_segment，及规则化后的点云距离分辨率resolution，多项式的项数TM  如二次函数X=a(1)T^2+b(1)T+c(1);Y=a(2)T^2+b(2)T+c(2),... TM=2
% 输出：规则化后的点云 Curve_feature_segment
function  [Curve_feature_segment] = CFPR(Curve_sort_segment,resolution,TM)  % Curve Feature point cloud regularization(CFPR)
for i=1:length(Curve_sort_segment)
    [range_resol] = range_resolut(Curve_sort_segment{i});
    [parameter,~] = space_curve_LS(Curve_sort_segment{i},TM);
    pnts=Curve_sort_segment{i};
    n=size(pnts,1);    
    interval=resolution/range_resol;
    T=0:interval:(n-1)*1;
    ZX=polyval(parameter(1,:),T); %曲线拟合后的规则化x方向坐标(1xn)
    ZY=polyval(parameter(2,:),T);%曲线拟合后的规则化y方向坐标(1xn)
    ZZ=polyval(parameter(3,:),T);%曲线拟合后的规则化z方向坐标(1xn)
    curve_fit_pnts=[ZX;ZY;ZZ]';
    neighbor=size(Curve_sort_segment{i},1);
    [idx,dist]=knnsearch(Curve_sort_segment{i}, Curve_sort_segment{i}, 'k', neighbor);
    max_dist=max(dist(:,neighbor));
    fat_idx_start=idx(dist(:,neighbor)==max_dist,1);
    fat_idx_end=idx(dist(:,neighbor)==max_dist,neighbor);
    fat_idx=[fat_idx_start,fat_idx_end];%第一列是初始点索引，第二列是终点索引
    fat_idx=transpose(sort(fat_idx')); %第一列是初始点索引，第二列是终点索引
    fat_idx=unique(fat_idx,'rows');
    fps_pnt=Curve_sort_segment{i}(fat_idx(1,:),:); 
    [idx1,dist1]=knnsearch(curve_fit_pnts, fps_pnt, 'k', 1);
    Curve_feature_segment{i}=curve_fit_pnts(min(idx1):max(idx1),:);
end

