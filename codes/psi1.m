function psi1=psi1(s,xi)
global a alpha alphad betae cd eta G v h kpr kptr ktpr ktr kpz kptz ktpz ktz M md T0; 
psi1=-((M.*s.*(alphad.*cd.*(betae.*s + (-kptr + kptz).*xi.^2) + alpha.*md.*(cd.*s + (ktr - ktz).*xi.^2)))./...
   (-2.*betae.^2.*cd.*eta.*G.*M.*s.^2 + 2.*betae.*eta.*G.*M.*(cd.*(kptr - kptz) + (ktpr - ktpz).*md).*s.*xi.^2 + alphad.^2.*cd.*s.*(s + (kpr - kpz).*M.*xi.^2) + ...
    alpha.*alphad.*M.*s.*(2.*betae.*cd.*s + (cd.*(-kptr + kptz) + (-ktpr + ktpz).*md).*xi.^2) + ...
    md.*(cd.*s.*(alpha.^2.*M.*s + 2.*eta.*G.*(s + (kpr - kpz).*M.*xi.^2)) + xi.^2.*(alpha.^2.*(ktr - ktz).*M.*s + 2.*eta.*G.*((-(kptr - kptz)).*(ktpr - ktpz).*M.*xi.^2 + ...
          ktr.*(s + (kpr - kpz).*M.*xi.^2) - ktz.*(s + (kpr - kpz).*M.*xi.^2))))));
end