function conn_contrasthelp(handle,title,str,values,types)
% gui help for likely contrast vectors/matrices

n=numel(values);
hc1=get(handle,'uicontextmenu');
if ishandle(hc1), delete(hc1); end; 
names={}; contrasts={}; rnames={};
if isempty(str), return; end
if isempty(types),types=ones(n,1); end
utypes=unique(types);
for n0=utypes(:)'
    idx=find(types==n0);
    if numel(idx)>1
        c=zeros(1,n);
        c(idx)=1/numel(idx);
        contrasts{end+1}=c;
        names{end+1}=['<HTML>(T-contrast) Main effect of <i>',strjoinstr(str(values(idx)),' & '),'</i>',' (',mat2str(contrasts{end},3),')'];
        if numel(idx)==numel(types),rnames{end+1}='Average';
        else rnames{end+1}=['<HTML>Average of <i>',strjoinstr(str(values(idx)),' & '),'</i>'];
        end
    end
    if  numel(idx)>1&&numel(idx)<4
        for n1=idx(:)',
            for n2=idx(:)'
                if n1~=n2,
                    contrasts{end+1}=full(sparse([1 1],[n1 n2],[1 -1],1,n));
                    names{end+1}=['<HTML>(T-contrast) <i>',str{values(n1)},'</i> > <i>',str{values(n2)},'</i>',' (',mat2str(contrasts{end},3),')'];
                    rnames{end+1}=['<HTML>Difference <i>',str{values(n1)},'</i> > <i>',str{values(n2)},'</i>'];
                end
            end
        end
    end
    if numel(idx)>1&&numel(utypes)>1
        c=zeros(numel(idx),n);
        c(:,idx)=eye(numel(idx));
        contrasts{end+1}=c;
        names{end+1}=['<HTML>(F-contrast) Any effect among <i>',strjoinstr(str(values(idx)),' or '),'</i>',' (',mat2str(contrasts{end},3),')'];
        if numel(idx)==numel(types), rnames{end+1}='Any effects (F-test)';
        else rnames{end+1}=['<HTML>Any effect among <i>',strjoinstr(str(values(idx)),' or '),'</i>'];
        end
    end
    if numel(idx)>2
        c=zeros(numel(idx)-1,n);
        c(:,idx)=convn(eye(numel(idx)),[1;-1],'valid');
        contrasts{end+1}=c;
        names{end+1}=['<HTML>(F-contrast) Any difference between <i>',strjoinstr(str(values(idx)),' & '),'</i>',' (',mat2str(contrasts{end},3),')'];
        if numel(idx)==numel(types), rnames{end+1}='Any differences (F-test)';
        else rnames{end+1}=['<HTML>Any difference between <i>',strjoinstr(str(values(idx)),' & '),'</i>'];
        end
    end
end
if n==1
    contrasts{end+1}=1;
    names{end+1}=['<HTML>(T-contrast) Effect of <i>',str{values(1)},'</i>',' (',mat2str(contrasts{end},3),')'];
    rnames{end+1}=['<HTML>Effect of <i>',str{values(1)},'</i>'];
else
    for n1=1:n,
        contrasts{end+1}=full(sparse(1,n1,1,1,n));
        names{end+1}=['<HTML>(T-contrast) Simple main effect of <i>',str{values(n1)},'</i>',' (',mat2str(contrasts{end},3),')'];
        rnames{end+1}=['<HTML>Effect of <i>',str{values(n1)},'</i>'];
    end
end
if n>1,%&&numel(utypes)>1
    contrasts{end+1}=eye(n);
    names{end+1}=['(F-contrast) Any effect of interest',' (',mat2str(contrasts{end},3),')'];
    rnames{end+1}='Any effects (F-test)'; %['Any effect of interest'];
end
hc0=get(handle,'userdata');
if isempty(hc0), 
    hc0=conn_menu('popup',get(handle,'position')*[1 0 0 0;0 1 0 0;0 0 1 0;0 -1 0 1],'',[{'<HTML><i>user-defined contrast</i></HTML>'},rnames],['List of suggested possible contrasts for the selected set of ',title],@(varargin)conn_contrasthelpcallback(handle,contrasts));
    set(handle,'value',1,'userdata',hc0);
    set(hc0,'fontsize',get(hc0,'fontsize'));
else set(hc0,'string',[{'<HTML><i>user-defined contrast</i></HTML>'},rnames],'value',1,'callback',@(varargin)conn_contrasthelpcallback(handle,contrasts));
end
contrast_value=full(str2num(get(handle,'string')));
findval=find(cellfun(@(x)isequal(full(x),contrast_value),contrasts));
if numel(findval)==1, set(hc0,'value',1+findval); end
% hc1=uicontextmenu('parent',handle);
% for n1=1:numel(names),
%     uimenu(hc1,'Label',names{n1},'callback',@(varargin)conn_contrasthelpcallback(handle,contrasts,n1));
% end
% set(handle,'uicontextmenu',hc1);
end

function conn_contrasthelpcallback(handle,contrasts,n1)
if nargin<3, n1=get(gcbo,'value')-1; end
if n1
    set(handle,'string',mat2str(contrasts{n1}));
    h=get(handle,'callback');
    if ~isempty(h)
        try
            if isa(h,'function_handle'), feval(h);
            else eval(h);
            end
        end
    end
end
end

function str=strjoinstr(str1,str2)
str=[str1(:)';repmat({str2},1,length(str1))];
str=reshape(str(1:end-1),1,numel(str)-1);
str=[str{:}];
end

