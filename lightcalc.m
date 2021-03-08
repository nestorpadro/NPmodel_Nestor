
%%%%%%%%%%%%%%%%%%%%%%%%%%%   LIGHT FUNCTION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I = lightcalc(z,C,deltaZ,kp,kbg,I0)

int = zeros;
n = length(z);
for i = 2:n
    int(i) = int(i-1)+kp*C(i)*deltaZ;
end

I(1) = I0 * (exp(-kbg * z(1)) * exp(-kp * C(1) * deltaZ));
I = I0 .* (exp(-kbg * z)) .* (exp(-int));

end