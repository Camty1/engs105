pos = [ -0.0611, -0.0526;
    -0.0542, -0.0550;
    -0.0486, -0.0533;
    -0.0597, -0.0501];

x = pos(:,1);
y = pos(:,2);

kappa = 0.642;
m = 0;

function val = phi(xi, nu, basis)
if basis == 1
    val = (1 - xi) * (1 - nu) / 4;
    
elseif basis == 2
    val = (1 + xi) * (1 - nu) / 4;
    
elseif basis == 3
    val = (1 + xi) * (1 + nu) / 4;
    
else
    val = (1 - xi) * (1 + nu) / 4;
    
end
end
