function dh = liczeniestanuciag(t,h,F1_0,FD,C1,A2,alfa1,alfa2)

     tau=150;
    if t>=tau
        F1 = F1_0;
    else
        F1=0;
    end
    dh = zeros(2, 1);
    dh(1)=nthroot(((F1+FD-alfa1*sqrt(h(1)))/C1),3);
    dh(2)=(alfa1*(h(1)^1/2)-alfa2*sqrt(h(2)))/A2;
end