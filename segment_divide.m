%第二步： 该步骤是对于第一步得到的segment进行曲线和直线分离处理，并且保留大于等于number的分割点云
%输入的是第一步的segment，多项式的项数TM，及保留的点数number  如二次函数X=a(1)T^2+b(1)T+c(1);Y=a(2)T^2+b(2)T+c(2),... TM=2
%number  保留的分割的大于等于该数据量的点云
function  [line_sort_segment,Curve_sort_segment] = segment_divide(segment,TM,number)  % Curve Feature point cloud regularization(FPR)


for i=1:length(segment)
    segment_n(i,:)=[size(segment{i},1),i];
end
    sort_segment_n=sortrows(segment_n,1,'descend');
    
 for i=1:length(segment)
     sort_segment{i}=segment{sort_segment_n(i,2)};
 end

    L=1;C=1;
for i=1:length(sort_segment)

    pnts=sort_segment{i};
        if size(pnts,1)>=number
    [line_vector1,mean_pnt1] = space_line_TLS(pnts);
    [PL_dis1] = PL_distance_TLS(pnts, mean_pnt1, line_vector1);
    mean_PLdis=mean(PL_dis1);
    
    [parameter1,~] = space_curve_LS(pnts,TM);  %
    [PC_dis1] = PC_distance_LS(pnts,parameter1); 
    mean_PCdis=mean(PC_dis1); %

    if mean_PLdis<=mean_PCdis
        line_sort_segment{L}=sort_segment{i};
        L=L+1;
    else
        Curve_sort_segment{C}=sort_segment{i};
        C=C+1;
    end
        end
end

    
    