function j = wrapIndex(i,L)
    wrapped = mod(i,L);
    j = (wrapped == 0).*L + (wrapped ~= 0).*wrapped;
end