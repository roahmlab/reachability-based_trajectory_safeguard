function data=filter_array(data,datafield,range)
L=datafield>=range(1)&datafield<=range(2);


data=data(L);


end