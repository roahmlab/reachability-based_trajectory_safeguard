function [out_data,shift,scale]=shift_and_scale_data(data)
shift=mean(data,2);

shifted_data=data-shift;

scale=max(abs(shifted_data),[],2);
scale(scale==0)=1;
out_data=shifted_data./repmat(scale,[1,size(data,2)]);
end