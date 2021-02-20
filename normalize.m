function xs=normalize(xs,M);
%****************************************************************************
%%%% normalizes energy

ncc=size(xs);
nc=ncc(2);


for comp=1:nc
norme=sqrt(xs(:,comp)'*M*xs(:,comp));
U(:,comp)=xs(:,comp)/norme;
end;
