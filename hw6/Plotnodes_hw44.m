clear all;close all

% File names
nodf='hw44.nod';
belf='hw44.bel';
dndf='hw44.dnd';

nod=load(nodf);
nod=nod(:,2:3); % All nodes
bel=load(belf); 
bel=bel(:,2:3); % Boundary element incidence list
bnod=unique(bel(:)); % Boundary nodes
dnd=load(dndf); 
dnd=dnd(:,2); % Type 1 BC nodes (ground)

% Plot nodes
plot(nod(:,1),nod(:,2),'b.')
title(nodf)
axis equal


hold on
% Plot Boundary nodes
plot(nod(bnod,1),nod(bnod,2),'g.')
% Plot type 1 BC nodes
plot(nod(dnd,1),nod(dnd,2),'bo')
% Plot Current sink
plot(nod(492:493,1),nod(492:493,2),'ro')
% Plot Current source
plot(nod(503,1),nod(503,2),'gx','markersize',18,'linewidth',2)
legend('Nodes','Boundary Nodes','Dirichlet nodes','Boundary current sink','current source')
