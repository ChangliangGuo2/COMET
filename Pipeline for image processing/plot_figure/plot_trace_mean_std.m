[~,idx]=sort(Cells_list_P(:,2),'descend');
Cells_list_P=Cells_list_P(idx,:);
id=A_neuron_good_idx(Cells_list_P(1:500,1))+1;
trace=normalized_C_trace(id,:);
trace_mean=mean(trace);
trace_SEM=std(trace)/sqrt(500);
% plot(trace_mean);
% hold on
a=area([(trace_mean-trace_SEM)',(trace_SEM)',(trace_SEM)']);
a(2).FaceAlpha = 0.3;a(3).FaceAlpha = 0.3;
a(1).LineStyle = 'none';a(3).LineStyle = 'none';
newcolors = [1 1 1; 1 0 0; 1 0 0];
colororder(newcolors)
