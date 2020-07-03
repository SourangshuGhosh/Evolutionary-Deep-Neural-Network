function dec = ari2dec(s)

global max_arity

s = double(s)-48;
j = length(s);
dec = 0;
for i = j:-1:1
    dec = dec + s(j-i+1)*max_arity^(i-1);
end

end