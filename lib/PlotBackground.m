function f3 = PlotBackground(GL,Ice,Firn)
            for k = 1:length(GL)
                ind = isnan(GL(k).X+GL(k).Y);
                GL(k).Y(ind)=[];
                GL(k).X(ind)=[];
                patch(GL(k).X,GL(k).Y,'k')
            end
            for k = 1:length(Ice)
                ind = isnan(Ice(k).X+Ice(k).Y);
                Ice(k).Y(ind)=[];
                Ice(k).X(ind)=[];
                patch(Ice(k).X,Ice(k).Y,RGB('light light gray'))
            end
for kk = 1: length(Firn)
    f3 = fill(Firn(kk).X,Firn(kk).Y,'g');
                f3.EdgeColor = RGB('gray');            
            f3.FaceColor = 'none';
end

end