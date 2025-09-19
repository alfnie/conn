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
        names{end+1}=['(T-contrast) Main effect of ',strjoinstr(str(values(idx)),' & '),'',' (',mat2str(contrasts{end},3),')'];
        if numel(idx)==numel(types),rnames{end+1}='Average';
        else rnames{end+1}=['Average of ',strjoinstr(str(values(idx)),' & '),''];
        end
    end
    if  numel(idx)>1&&numel(idx)<4
        for n1=idx(:)',
            for n2=idx(:)'
                if n1~=n2,
                    contrasts{end+1}=full(sparse([1 1],[n1 n2],[1 -1],1,n));
                    names{end+1}=['(T-contrast) ',str{values(n1)},' > ',str{values(n2)},'',' (',mat2str(contrasts{end},3),')'];
                    rnames{end+1}=['Difference ',str{values(n1)},' > ',str{values(n2)},''];
                end
            end
        end
    end
    if numel(idx)>1&&numel(utypes)>1
        c=zeros(numel(idx),n);
        c(:,idx)=eye(numel(idx));
        contrasts{end+1}=c;
        names{end+1}=['(F-contrast) Any effect among ',strjoinstr(str(values(idx)),' or '),'',' (',mat2str(contrasts{end},3),')'];
        if numel(idx)==numel(types), rnames{end+1}='Any effects (F-test)';
        else rnames{end+1}=['Any effect among ',strjoinstr(str(values(idx)),' or '),''];
        end
    end
    if numel(idx)>2
        c=zeros(numel(idx)-1,n);
        c(:,idx)=convn(eye(numel(idx)),[1;-1],'valid');
        contrasts{end+1}=c;
        names{end+1}=['(F-contrast) Any difference between ',strjoinstr(str(values(idx)),' & '),'',' (',mat2str(contrasts{end},3),')'];
        if numel(idx)==numel(types), rnames{end+1}='Any differences (F-test)';
        else rnames{end+1}=['Any difference between ',strjoinstr(str(values(idx)),' & '),''];
        end
    end
end
if n==1
    contrasts{end+1}=1;
    names{end+1}=['(T-contrast) Effect of ',str{values(1)},'',' (',mat2str(contrasts{end},3),')'];
    rnames{end+1}=['Effect of ',str{values(1)},''];
else
    for n1=1:n,
        contrasts{end+1}=full(sparse(1,n1,1,1,n));
        names{end+1}=['(T-contrast) Simple main effect of ',str{values(n1)},'',' (',mat2str(contrasts{end},3),')'];
        rnames{end+1}=['Effect of ',str{values(n1)},''];
    end
end
if n>1,%&&numel(utypes)>1
    contrasts{end+1}=eye(n);
    names{end+1}=['(F-contrast) Any effect of interest',' (',mat2str(contrasts{end},3),')'];
    rnames{end+1}='Any effects (F-test)'; %['Any effect of interest'];
end
hc0=get(handle,'userdata');
if isempty(hc0), 
    hc0=conn_menu('popup',get(handle,'position')*[1 0 0 0;0 1 0 0;0 0 1 0;0 -1 0 1],'',[{'custom contrast'},rnames],['List of suggested possible contrasts for the selected set of ',title],@(varargin)conn_contrasthelpcallback(handle,contrasts));
    set(handle,'value',1,'userdata',hc0);
    set(hc0,'fontsize',get(hc0,'fontsize'));
else set(hc0,'string',[{'custom contrast'},rnames],'value',1,'callback',@(varargin)conn_contrasthelpcallback(handle,contrasts));
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

