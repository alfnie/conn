function fh=conn_table_display(data,varargin)
% CONN_TABLE_DISPLAY
% displays contingency table data
%
% conn_table_display((data [,'rlabel', labels_rows, 'clabel', labels_cols, ...]);
%       data        : (N1 x N2) table of counts/values
%       labels_rows : (N1) labels of first dimension / rows
%       labels_cols : (N2) labels of seconddimension / columns
%
% fh = conn_table_display(...); returns function handle for additional options
%
% fh('colormap' [,cmap, dim]);       : colormap (cmap = [#colors x 3]; dim=1 color lines based on rows; dim=2 color lines based on columns)
% fh('fontsize' [, v1, v2, v3]);     : changes fontsize  (v1=labels-rows, v2=labels-columns, v3=values)
% fh('background' [, rgb]);          : changes background color
% fh('showtril' [, state]);          : show lower-triangular portion -for square tables only- ('on' 'off')
% fh('showval' [,state]);            : display numeric values ('on' 'off')
% fh('print' [,filename]);           : prints display to high-resolution .jpg file
%

options=struct(...
    'newfigure',true,...            % create new figure
    'rfontsize',12,...              % fontsize of left labels
    'cfontsize',12,...              % fontsize of right labels
    'vfontsize',10,...              % fontsize of right labels
    'dfontsize',2,...               % fonsize proportional to values
    'color',jet(256),...            % colormap
    'colorsign',[],...              % color based on neg/pos sign
    'colordim',1,...                % color based on left (1) or right (2) order
    'background',[.95 .95 .9],...
    'fontcolor',.25*[1 1 1],...     
    'facealpha',.5,...
    'separation',.25,...            % separation between levels of each factor
    'showval',[],...                % also show values within lines
    'showtril',true,...             % also show lower triangular
    'showtril_idx',[],...
    'showalllabels',false,...       % also show labels with 0 values
    'rlabel',{{}},...               % left labels
    'clabel',{{}} );                % right labels

for n1=1:2:numel(varargin)-1
    assert(isfield(options,varargin{n1}),'unrecognized option %s',varargin{n1});
    options.(varargin{n1})=varargin{n1+1};
end

hstruct=struct();
data(isnan(data))=0;
issquare=size(data,1)==size(data,2);
signdata=sign(data);
if isempty(options.colorsign), options.colorsign=any(data(:)<0)&any(data(:)>0); end
if isempty(options.showtril), if isequal(data,data'), [options.showtril_idx,nill]=dmperm((data'*data)~=0); [i,j]=find(triu(data(options.showtril_idx,options.showtril_idx))); options.showtril=~isempty(intersect(i,j)); else; options.showtril=true; end; end
if isempty(options.showval), options.showval=numel(unique(data(data~=0)))>1; end
data=abs(data);
N1=size(data,1);
N2=size(data,2);
assert(numel(options.rlabel)==0|numel(options.rlabel)==N1,'incorrect number of elements in rlabel (found %d expected %d)',numel(options.rlabel),N1);
assert(numel(options.clabel)==0|numel(options.clabel)==N2,'incorrect number of elements in clabel (found %d expected %d)',numel(options.clabel),N2);
if isempty(options.rlabel), options.rlabel=arrayfun(@(n)sprintf('Row%d',n),1:N1,'uni',0); end
if isempty(options.clabel), options.clabel=arrayfun(@(n)sprintf('Col%d',n),1:N2,'uni',0); end

if options.newfigure, 
    hmsg=conn_msgbox('Initializing. Please wait...','',-1);
    hstruct.hfig=figure('numbertitle','off','menubar','none','name','connection display','units','norm','position',[.3 .3 .3 .6],'color',options.background); 
    hstruct.hax=axes('units','norm','position',[.25 .05 .5 .9],'color',options.background);
    hc1=uimenu(hstruct.hfig,'Label','Effects');
    if N1==N2, hc2=uimenu(hc1,'Label','hide/show lower triangular part','callback',{@conn_table_display_refresh,'showtril'}); end
    hc2=uimenu(hc1,'Label','hide/show values','callback',{@conn_table_display_refresh,'showval'});
    hc2=uimenu(hc1,'Label','line colors','callback',{@conn_table_display_refresh,'colormap'});
    hc2=uimenu(hc1,'Label','fontsize');
    uimenu(hc2,'Label','increase labels fontsize','callback',{@conn_table_display_refresh,'fontsize','+'});
    uimenu(hc2,'Label','decrease labels fontsize','callback',{@conn_table_display_refresh,'fontsize','-'});
    uimenu(hc2,'Label','set labels fontsize','callback',{@conn_table_display_refresh,'fontsize'});
    hc2=uimenu(hc1,'Label','background');
    uimenu(hc2,'Label','white background','callback',{@conn_table_display_refresh,'background',[1 1 1]});
    uimenu(hc2,'Label','light background','callback',{@conn_table_display_refresh,'background',[.95 .95 .9]});
    uimenu(hc2,'Label','dark background','callback',{@conn_table_display_refresh,'background',[.11 .11 .11]});
    uimenu(hc2,'Label','black background','callback',{@conn_table_display_refresh,'background',[0 0 0]});
    hc1=uimenu(hstruct.hfig,'Label','Print');
    uimenu(hc1,'Label','current view','callback',{@conn_table_display_refresh,'print'});
else
    hstruct.hfig=gcf; 
    hstruct.hax=gca;
end

fh=@(varargin)conn_table_display_refresh([],[],varargin{:});
conn_table_display_refresh([],[],'update');

if options.newfigure, 
    hstruct.fh=fh;
    set(hstruct.hfig,'userdata',hstruct); 
    if ishandle(hmsg), delete(hmsg); end
end

    function conn_table_display_refresh(hObject,eventdata,option,varargin)
        switch(option)
            case 'update'
                hstruct.hlines=[];
                hstruct.vtext=[];
                cla(hstruct.hax);
                if options.showtril, 
                    thisdata=data;
                    thisrlabel=options.rlabel;
                    thisclabel=options.clabel;
                else
                    if isempty(options.showtril_idx), [options.showtril_idx,nill]=dmperm((data'*data)~=0); end
                    thisdata=triu(data(options.showtril_idx,options.showtril_idx));
                    thisrlabel=options.rlabel(options.showtril_idx);
                    thisclabel=options.clabel(options.showtril_idx);
                end
                sx1=sum(thisdata,2)';
                sx2=sum(thisdata,1);
                sx=sum(thisdata(:));
                mx=max(thisdata(:));
                x1=1-[0, cumsum(sx1)]/sx;
                x2=1-[0, cumsum(sx2)]/sx;
                wx=[linspace(0,1,256),linspace(1,0,256)];
                wy=.5+.5*tanh(linspace(-3,3,256));
                for n1=1:N1,
                    for n2=1:N2,
                        if thisdata(n1,n2)>0,
                            k1=sum(thisdata(n1,1:n2-1),2)/sx1(n1);
                            k2=k1+thisdata(n1,n2)/sx1(n1);
                            a1=x1(n1)+(x1(n1+1)-x1(n1))*(options.separation/2+(1-options.separation)*k1);
                            a2=x1(n1)+(x1(n1+1)-x1(n1))*(options.separation/2+(1-options.separation)*k2);
                            
                            k1=sum(thisdata(1:n1-1,n2),1)/sx2(n2);
                            k2=k1+thisdata(n1,n2)/sx2(n2);
                            b1=x2(n2)+(x2(n2+1)-x2(n2))*(options.separation/2+(1-options.separation)*k1);
                            b2=x2(n2)+(x2(n2+1)-x2(n2))*(options.separation/2+(1-options.separation)*k2);
                            
                            y=[a1+(b1-a1)*wy, b2+(a2-b2)*wy];
                            if options.colorsign,       color=options.color(round(1+(size(options.color,1)-1)*(signdata(n1,n2)>0)),:);
                            elseif options.colordim==1, color=options.color(round(1+(size(options.color,1)-1)*(x1(n1)+x1(n1+1))/2),:);
                            else,                       color=options.color(round(1+(size(options.color,1)-1)*(x2(n2)+x2(n2+1))/2),:);
                            end
                            hstruct.hlines(end+1)=patch(wx,y,'k','edgecolor','none','facecolor',color,'facealpha',options.facealpha,'parent',hstruct.hax);
                            if issquare&&n2<n1, set(hstruct.hlines(end),'tag','symmetric'); if ~options.showtril, set(hstruct.hlines(end),'visible','off'); end; end
                            hold(hstruct.hax,'on'); 
                            plot([wx(1) wx(1)],[y(1) y(end)],'k-','linewidth',3,'parent',hstruct.hax); 
                            plot([wx(256) wx(256)],[y(256) y(257)],'k-','linewidth',3,'parent',hstruct.hax); 
                            if round(thisdata(n1,n2))>=1e3, ldata=num2str(signdata(n1,n2)*round(thisdata(n1,n2))); else ldata=mat2str(signdata(n1,n2)*thisdata(n1,n2),3); end
                            hstruct.vtext(end+1)=text(wx(1)+.02,(y(1)+y(end))/2,ldata,'horizontalalignment','left','fontsize',max(1,options.vfontsize+options.dfontsize*(-1+2*thisdata(n1,n2)/mx)),'color',options.fontcolor,'tag','values','parent',hstruct.hax); 
                            hstruct.vtext(end+1)=text(wx(256)-.02,(y(256)+y(257))/2,ldata,'horizontalalignment','right','fontsize',max(1,options.vfontsize+options.dfontsize*(-1+2*thisdata(n1,n2)/mx)),'color',options.fontcolor,'tag','values','parent',hstruct.hax); 
                            hold(hstruct.hax,'off');
                            if ~options.showval, set(hstruct.vtext(end-1:end),'visible','off'); end
                            if issquare&&n2<n1, set(hstruct.vtext(end-1:end),'tag','symmetric'); if ~options.showtril, set(hstruct.vtext(end-1:end),'visible','off'); end; end
                            
                        end
                    end
                end
                hstruct.xtext=[];
                hstruct.ytext=[];
                if ~isempty(thisrlabel)
                    for n1=1:N1,
                        if options.showalllabels||sx1(n1)>0
                            hstruct.xtext(end+1)=text(-.025,.5*x1(n1)+.5*x1(n1+1),thisrlabel{n1},'horizontalalignment','right','fontsize',max(1,options.rfontsize+options.dfontsize*(-1+2*sx1(n1)/max(sx1))),'color',options.fontcolor,'parent',hstruct.hax);
                        end
                    end
                end
                if ~isempty(thisclabel)
                    for n2=1:N2,
                        if options.showalllabels||sx2(n2)>0
                            hstruct.ytext(end+1)=text(1.025,.5*x2(n2)+.5*x2(n2+1),thisclabel{n2},'horizontalalignment','left','fontsize',max(1,options.cfontsize+options.dfontsize*(-1+2*sx2(n2)/max(sx2))),'color',options.fontcolor,'parent',hstruct.hax);
                        end
                    end
                end
                axis(hstruct.hax,'off');
                
            case 'colormap'
                cmap=[];
                colordim=[];
                if numel(varargin)>0&&~isempty(varargin{1}), cmap=varargin{1}; if ischar(cmap), cmap=str2num(cmap); end; end
                if numel(varargin)>1&&~isempty(varargin{2}), colordim=varargin{2}; if ischar(colordim), colordim=str2num(colordim); end; end
                if isempty(cmap)||isempty(colordim)
                    if isempty(cmap), cmap=options.color; end
                    if isempty(colordim), if options.colorsign, colordim=0; else colordim=options.colordim; end; end
                    s=conn_menu_inputdlg({'Enter colors: (Nx3 rgb)','Colors represent: (1:rows, 2:columns; 0:sign)'},'conn_table_display',1,{mat2str(cmap),mat2str(colordim)});
                    if isempty(s), return; end
                    cmap=str2num(s{1});
                    colordim=str2num(s{2});
                end
                if size(cmap,2)==3&~isempty(colordim),
                    options.color=cmap;
                    if colordim==0, 
                        options.colorsign=true;
                    else
                        options.colorsign=false;
                        options.colordim=colordim;
                    end
                    conn_table_display_refresh([],[],'update');
                end
                    
            case 'showtril'
                if numel(varargin)>0&&isequal(varargin{1},'on'), options.showtril=true;
                elseif numel(varargin)>0&&isequal(varargin{1},'off'), options.showtril=false;
                else options.showtril=~options.showtril;
                end
                conn_table_display_refresh([],[],'update');
                %h=findobj(hstruct.hfig,'tag','symmetric'); 
                %if options.showtril, set(h,'visible','on'); 
                %else set(h,'visible','off'); 
                %end
                %if ~options.showval, set(hstruct.vtext,'visible','off'); end
                
            case 'showval'
                if numel(varargin)>0&&isequal(varargin{1},'on'), options.showval=true;
                elseif numel(varargin)>0&&isequal(varargin{1},'off'), options.showval=false;
                else options.showval=~options.showval;
                end
                h=findobj(hstruct.hfig,'tag','symmetric'); 
                if options.showval, set(hstruct.vtext,'visible','on'); 
                else set(hstruct.vtext,'visible','off'); 
                end
                if ~options.showtril, set(h,'visible','off'); end
                
            case 'fontsize'
                if numel(varargin)>0&&(isequal(varargin{1},'+')||isequal(varargin{1},'-'))
                    if strcmp(varargin{1},'+'), ds=1; else ds=-1; end
                    h=findobj(hstruct.hfig,'type','text');
                    s=max(1,cell2mat(get(h,'fontsize'))+ds);
                    for ss=reshape(unique(s),1,[]),
                        set(h(s==ss),'fontsize',ss);
                    end
                    options.rfontsize=options.rfontsize+ds;
                    options.cfontsize=options.cfontsize+ds;
                    options.vfontsize=options.vfontsize+ds;
                else
                    if numel(varargin)>0&&~isempty(varargin{1}), 
                        s=varargin{1}; 
                        if ischar(s), s=str2num(s); end
                    else
                        s=conn_menu_inputdlg('Enter fontsize (labels1, labels2, values)','conn_table_display',1,{num2str([options.rfontsize options.cfontsize options.vfontsize])});
                        if ~isempty(s), s=str2num(s{1}); end
                    end
                    if ~isempty(s), 
                        if numel(s)<3, s=s(min(numel(s),1:3)); end
                        set(hstruct.xtext,'fontsize',max(1,s(1))); 
                        set(hstruct.ytext,'fontsize',max(1,s(2))); 
                        set(hstruct.vtext,'fontsize',max(1,s(3))); 
                        options.rfontsize=s(1);
                        options.cfontsize=s(2);
                        options.vfontsize=s(3);
                    end
                end
            case 'background'
                if numel(varargin)>0&&~isempty(varargin{1}), tcolor=varargin{1};
                else tcolor=uisetcolor(options.background,'Select color'); if isempty(tcolor)||isequal(tcolor,0), return; end; 
                end
                if ischar(tcolor), tcolor=str2num(tcolor); end
                h=findobj(hstruct.hfig,'type','text'); 
                nc0=get(hstruct.hfig,'color'); 
                set([hstruct.hfig hstruct.hax],'color',tcolor); 
                hc=cell2mat(get(h,'color')); 
                for nc=unique(hc,'rows')', 
                    idx=all(bsxfun(@eq,nc',hc),2); 
                    set(h(idx),'color',max(0,min(1,tcolor-sign(mean(tcolor-.5))*abs(nc'-nc0)))); 
                end
                options.background=tcolor;
                
            case 'print'
                if numel(varargin)>0, pfilename=varargin{1}; options=varargin(2:end);
                else pfilename=fullfile(pwd,'print01.jpg'); options={};
                end
                conn_print(hstruct.hfig,pfilename,options{:});
        end
    end
end

