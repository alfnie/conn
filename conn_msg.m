function varargout = conn_msg(N)
% 'please wait' message
% conn_msg(N)
%
persistent idx msgs nmsgs

if nargin<1||isempty(N), N=1; end
if isempty(msgs)
    msgs={{'loading','saving','swapping','trying','drafting','worrying over','searching for','looking for','examining','updating','preparing','crafting', 'mining', 'initializing', 'downloading', 'running', 'building', 'discovering', 'rediscovering','sending','receiving','holding','finding','computing','fixing','waiting for','following','ordering','installing','testing','evaluating','considering','disentangling','composing','reciting','checking','distributing','formatting','delineating','describing','analyzing','rearranging','discretizing','summarizing','printing','starting','generating','promoting','expecting','upgrading','restarting','deconstructing','randomizing','plotting','improving','exploring','surveying','recording','simulating','requesting','resolving','accepting','overinterpreting','researching','assembling','enhancing','supporting','dealing with','setting up','developing','imagining','drawing','visualizing','remembering','reconfiguring','harmonizing','handling','identifying','characterizing','restoring','securing','glossing over','introducing','welcoming','selecting','confirming','recommending','maintaining','implementing','familiarizing with'},...
          {'the world','space-time','plans','minions','more minions','free energy','prime numbers','a hadamard product','a kronecker product','eigenvalues','spectra','data','patience','gifts','bits','clusters','aliens','a melody','society','gradients','laplacians','stuff','work','jobs','warp drive','something else','life','real life','unicorns','the universe','another universe','optimism','a flux capacitor','42','something','everything','something unexpected','hope','an alternative','a solution','success','satellites','Massachusetts','a loading screen','webpage','social media','humorous message','humor','infinite','Java','Matlab','Unix','energy sources','electromagnetism','wave-particle duality','Maxwell''s equations','a wave function','quantum leap','a holodeck','Spock','robots','a transporter','a grant','optimizer','buried treasure','a dynamic language','a happy song','the internet','pi','a coffee','a constant','circuits','equations','connections','happiness','Bob','diversity','gravitational constant','vector fields','gravitational fields','random fields','a white rabbit','side project','all numbers','differences','confidence','waves','joy','convergence','available resources','an end game'}};
    nmsgs=cellfun('length',msgs);
end

for n=1:N
    if isempty(idx), idx=ceil(rand(size(msgs)).*nmsgs);
    else idx=1+mod(idx+floor(rand(size(msgs)).*(nmsgs-1)),nmsgs);
    end
    tmsgs=arrayfun(@(i)msgs{i}{idx(i)},1:numel(idx),'uni',0);
    tstr=sprintf(' %s',tmsgs{:});
    if n==1, str=tstr;
    elseif n==N, str=[str,' while',tstr];
    else str=[str,' and',tstr];
    end
end
str=[str,', please wait...']; 
if ~nargout, disp(str); varargout={}; 
else varargout={str};
end
