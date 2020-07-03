function s = dec2ari(dec)

global max_arity

if dec == 0
    s = '0';
else
    s = '';
    while (dec ~= 0)
        s = [s num2str(rem(dec,max_arity))];
        dec = floor(dec/max_arity);
    end
    s = s(length(s):-1:1);
end

end