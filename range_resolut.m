%���ڼ�����ƾ���ֱ��ʣ���ȷ�����ڵ�֮��ľ���
function [range_resol] = range_resolut(input_pnts)
neighbor_idx=knnsearch(input_pnts,input_pnts,'k',2);
dis=sqrt(sum((input_pnts-input_pnts(neighbor_idx(:,2),:)).^2,2));
range_resol=mean(dis);