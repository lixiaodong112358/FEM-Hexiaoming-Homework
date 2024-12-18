function [boundarynodes]=generate_boundarynodes(Pb,boundary_type)
if boundary_type=='DD'
boundarynodes=[0,0;1,max(size(Pb))];
elseif boundary_type=='DN'
    boundarynodes=[0,1;1,max(size(Pb))];
end
end