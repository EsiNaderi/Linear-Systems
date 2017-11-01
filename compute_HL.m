function HL=compute_HL(sys,L)

H0=sys.d;
O0=sys.c;
s1=size(sys.d,1);
s2=size(H0,2);

if L==0
    HL=H0;
end
if L==-1
    HL=0;    
end
for i=1:L
    HL=[sys.d zeros(s1,s2);O0*sys.b H0];
    H0=HL;
    s2=size(H0,2);
    O0=compute_OL(sys,i);
end
