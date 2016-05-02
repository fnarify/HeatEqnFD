function visualise(x,t,u,mask,opts)

if nargin<5 || strcmp(opts,'')
    
    opts='';
    
end

%-----------------------------------
%
%-----------------------------------

umax=max(max(u));
umin=min(min(u));

if length(strfind(opts,'x'))
    
else
    
    for k=1:length(mask)
        
        figure(k)
        plot(x,u(:,mask(k)))
        xlabel('x')
        ylabel('u')
        title(sprintf('t=%0.1e',t(mask(k))))
        ax=axis();
        ax(3:4)=[umin umax];
        axis(ax);
        
    end
    
end