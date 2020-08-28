clear all
close all
clc
% Decomposition Analysis
Colombia2017;%calldatabase;
fy=1990-1989;
ty=2017-1989;
%LDMI I - Additive approach
k=1; 
    for j=fy:ty
D1a(k,j)=C1(k,j)/Et(1,j);  
D1b(k,j)=C2(k,j)/Et(1,j);  
D1c(k,j)=C3(k,j)/Et(1,j);
e1(k,j)=C1(k,j)/Wf1(k,j);%e1j        
e2(k,j)=C2(k,j)/Wf2(k,j);%e2j
e3(k,j)=C3(k,j)/Wf3(k,j);%e3j
u1(k,j)=Wf1(k,j)/Ef1(k,j); %u1j 
u2(k,j)=Wf2(k,j)/Ef2(k,j); %u2j 
u3(k,j)=Wf3(k,j)/Ef3(k,j); %u3j 
m1(k,j)=Ef1(k,j)/Ef(k,j); %m1j
m2(k,j)=Ef2(k,j)/Ef(k,j); %m2j 
m3(k,j)=Ef3(k,j)/Ef(k,j); %m3j 
p(k,j)=Ef(k,j)/Et(k,j); % pj
    end
L1c(k)=(D1c(k,ty)-D1c(k,fy))/(log(D1c(k,ty))-log(D1c(k,fy)));
L1a(k)=(D1a(k,ty)-D1a(k,fy))/(log(D1a(k,ty))-log(D1a(k,fy)));
L1b(k)=(D1b(k,ty)-D1b(k,fy))/(log(D1b(k,ty))-log(D1b(k,fy)));
AD2(k,1)=L1a(k)*log(e1(k,ty)/e1(k,fy));% 
AD2(k,2)=L1b(k)*log(e2(k,ty)/e2(k,fy));% 
AD2(k,3)=L1c(k)*log(e3(k,ty)/e3(k,fy));% 
AD2(k,4)=L1a(k)*log(u1(k,ty)/u1(k,fy));% 
AD2(k,5)=L1b(k)*log(u2(k,ty)/u2(k,fy));% 
AD2(k,6)=L1c(k)*log(u3(k,ty)/u3(k,fy));% 
AD2(k,7)=L1a(k)*log(m1(k,ty)/m1(k,fy));% 
AD2(k,8)=L1b(k)*log(m2(k,ty)/m2(k,fy));% 
AD2(k,9)=L1c(k)*log(m3(k,ty)/m3(k,fy));% 
AD2(k,10)=L1a(k)*log(p(k,ty)/p(k,fy));% 
AD2(k,11)=L1b(k)*log(p(k,ty)/p(k,fy));% 
AD2(k,12)=L1c(k)*log(p(k,ty)/p(k,fy));% 
sum(AD2)-((C1(k,ty)+C2(k,ty)+C3(k,ty))/Et(1,ty)-(C1(k,fy)+C2(k,fy)+C3(k,fy))/Et(1,fy)) 

