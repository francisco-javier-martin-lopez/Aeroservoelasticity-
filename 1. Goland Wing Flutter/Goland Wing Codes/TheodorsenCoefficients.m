

function CoefficientsVector = ...
    TheodorsenCoefficients(c,e)

mu = acos(c);
T1 = -1/3 * sqrt(1 - c^2) * (2 + c^2) + c*mu;
T2 = c*(1-c^2) - sqrt(1-c^2) * (1 + c^2)*mu + c*mu^2;
T3 = -(1/8 + c^2)*mu^2 + 1/4 * c * sqrt(1-c^2) * mu * (7 + 2*c^2) ...
    -1/8*(1-c^2)*(5*c^2+4);
T4 = -mu * c*sqrt(1-c^2);
T5 = -(1-c^2) - mu^2 + 2*c*sqrt(1-c^2)*mu;
T7 = -(1/8+c^2)*mu + 1/8*c*sqrt(1-c^2)*(7+2*c^2);
T8 = -1/3*sqrt(1-c^2)*(2*c^2+1)+c*mu;
T9 = (1/2)*(1/3*(sqrt(1-c^2))^3+e*T4);
T10 = sqrt(1-c^2)+mu;
T11 = mu*(1-2*c)+sqrt(1-c^2)*(2-c);
T12 = sqrt(1-c^2)*(2+c) - mu*(2*c+1);
T13 = -1/2 * (T7 + (c-e)*T1);

CoefficientsVector = [T1;T2;T3;T4;T5;0;T7;T8;T9;T10;T11;T12;T13];

end