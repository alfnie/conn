function font_offset=conn_font_init(varargin)
global CONN_gui;
hfig=[];
try
    h0=get(0,'screensize'); h0=h0(1,3:4)-h0(1,1:2)+1; h0(1)=min(h0(1),2*h0(2)); %h0=h0/max(1,max(abs(h0))/2000);
    minheight=500;
    hfig=figure('units','pixels','position',[0*72+1,h0(2)-max(minheight,.5*h0(1))-48,h0(1)-0*72-1,max(minheight,.5*h0(1))],'visible','off');
    hax=axes('units','norm','position',[0 0 1 1],'parent',hfig);
    h=text(0,-2,'Initializing. Please wait','fontunits','norm','fontsize',1/60,'horizontalalignment','center','verticalalignment','bottom','color',.75*[1 1 1],'parent',hax);
    set(hax,'units','norm','position',[0 0 1 1],'xlim',[-2 2],'ylim',[-2.5 2]);
    drawnow;
    set(h,'fontunits','points');
    tfontsize=get(h,'fontsize');
    font_offset=max(-4,round(tfontsize)-8);
    if ~nargout,
        CONN_gui.font_offset=font_offset;
        set(0,{'defaultuicontrolfontsize','defaulttextfontsize','defaultaxesfontsize'},repmat({8+CONN_gui.font_offset},1,3));
    end
    delete(hfig);
catch
    if ishandle(hfig), delete(hfig); end
    if ~nargout, CONN_gui.font_offset=0;
    else font_offset=0;
    end
end
end
