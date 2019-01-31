function [h, g, f3] = PlotBackground(GL,Ice,Firn)
            for k = 1:length(GL)
                ind = isnan(GL(k).X+GL(k).Y);
                GL(k).Y(ind)=[];
                GL(k).X(ind)=[];
                g(k) = patch(GL(k).X,GL(k).Y,RGB('black'));
%                 g(k).EdgeColor = RGB('gray'); 
                g(k).LineWidth = 0.2; 
            end
            for k = 1:length(Ice)
                ind = isnan(Ice(k).X+Ice(k).Y);
                Ice(k).Y(ind)=[];
                Ice(k).X(ind)=[];
                h = patch(Ice(k).X,Ice(k).Y,RGB('light light gray'));
            end
for kk = 1: length(Firn)
    f3 = fill(Firn(kk).X,Firn(kk).Y,'g');
                f3.EdgeColor = RGB('light light gray');            
            f3.FaceColor = 'none';
end

end