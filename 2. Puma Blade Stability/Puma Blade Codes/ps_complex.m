function Res = ps_complex(X, Y)
    S1 = 0;
    S2 = 0;
    S3 = 0;
    S4 = 0;
    for jjj=1:length(X)
        S1 = S1 + real(X(jjj))*real(Y(jjj)) + imag(X(jjj))*imag(Y(jjj));
        S2 = S2 + real(X(jjj))*imag(Y(jjj)) - imag(X(jjj))*real(Y(jjj));
        S3 = norm(X(jjj))^2 + S3;
        S4 = norm(Y(jjj))^2 + S4;
    end
    Res = sqrt(S1^2 + S2^2)/sqrt(S3*S4);
end