function groups = group2(m,nsteps)

if nargin<2||isempty(nsteps), nsteps=10; end % set to 0 for spectral bisection initialization, 1 or higher for additional random-initialization iterations
N=size(m,1);
W=abs(m)/max(abs(m(:))); 
W(1:N+1:end)=0; 
W=(W+W')/2; 
sW=sum(W(:));

bestS=nan;
bestgroups=[];
for iter=0:nsteps
    if iter==0, % start with 2nd smallest eigenvector of laplacian
        D=diag(sum(W,2));
        L=eye(size(D,1))-diag(1./sqrt(diag(D)))*W*diag(1./sqrt(diag(D)));
        [q,d]=eig(L);
        [nill,idx]=sort(q(:,2));
        groups=reshape(idx<=N/2,1,[]);
    else
        groups=randperm(N)<=N/2;
    end
    for n=1:N^2
        S0 = groupscore(W,groups); %fprintf('%.1f ',S0);
        if 1            
            g=2*groups-1;
            k=2*(W'*g' - g*W) - 4*W; % deltaS if we switch (i,j) (assuming symmetric and zero-diagonal)
            k=k.*(groups==1 & groups'==0);
            maxk=max(k(:));
            if maxk<=0, break; end
            [i,j]=find(k==maxk,1);
            groups([i,j])=groups([j,i]);
            %imagesc(groups); drawnow;
        else
            if 1
                maxS=S0; maxidx=[];
                for i=1:N-1,
                    for j=i+1:N
                        if groups(i)~=groups(j)
                            tgroups=groups; tgroups([i,j])=tgroups([j,i]);
                            S = groupscore(W,tgroups);
                            if S>maxS, maxidx=tgroups; maxS=S; end
                        end
                    end
                end
                if isempty(maxidx), break; end
                groups=maxidx;
            else
                maxS=S0; maxidx=[];
                for idx=find(groups(:)'==0),
                    tgroups=groups; tgroups(idx)=1-tgroups(idx);
                    S = groupscore(W,tgroups);
                    if S>maxS, maxidx=idx; maxS=S; end
                end
                if isempty(maxidx), break; end
                maxidx1=maxidx;

                maxS=S0; maxidx=[];
                for idx=find(groups(:)'==1),
                    tgroups=groups; tgroups(idx)=1-tgroups(idx);
                    S = groupscore(W,tgroups);
                    if S>maxS, maxidx=idx; maxS=S; end
                end
                if isempty(maxidx), break; end
                maxidx2=maxidx;
                groups(maxidx1)=1-groups(maxidx1);
                groups(maxidx2)=1-groups(maxidx2);
            end
        end
    end
    if n==N^2, disp('maximum iteration reached'); end
    S = groupscore(W,groups);
    if isnan(bestS)||S>bestS, bestgroups=groups; bestS=S; fprintf('\n%.1f ',2*bestS-sW); end
end
groups=bestgroups;
end

function score = groupscore(W,groups)
within  = groups==groups';
score = W(:)'*within(:);
%scorereal=2*score-sum(sum(W(:)));
%between = groups~=groups';
%score = W(:)'*(within(:)-between(:));
end