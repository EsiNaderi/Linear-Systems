function OL=compute_OL(sys,L)


O0=sys.c;

if L==0;
    OL=O0;
end

for i=1:L
    TM=O0*sys.a;
    OL=[sys.c;TM];
    O0=OL;
end
    