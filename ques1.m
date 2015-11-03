function [error] = ques1(n )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
p=@(x)(-2*x/(1+x^2));                                           %p(x)
q=@(x)(1);                                                      %q(x)
r=@(x)((2/(1+x^2))-log(1+x^2));                               %r(x)
aa=0;bb=1;y0=0;yn=log(2);Y=@(x)log(1+x^2);                  %a,b,y(a),y(b),Y(X)
fprintf('\n');

%h=1/2^(n-1);
h=(bb-aa)/n;
y=zeros(n+1,1);
y(1)=y0;
y(n+1)=yn;
a=zeros(n+1,1);
b=zeros(n+1,1);
c=zeros(n+1,1);
blank='  ';
disp(['finite difference method with h= ',num2str(h)]);
fprintf('\n');
%compute tridiagonal matrix
%a(1)=0;
for i=1:n-1   %y(2) is correponding to x(1) so i:1 to n-1
    x=aa+i*h;
    if(i~=1)
        a(i)=-1-((h/2)*feval(p,x));      %a2,a3,....,a(n-1)
    end
    b(i)=2+((h^2)*feval(q,x));              %b1,b2...,b(n-1)
    if(i~=n-1)
        c(i)=-1+((h/2)*feval(p,x));         %c1,c2,...,c(n-1)
    end
end
%compute diagonal and subdiagonal entries of tridiagonal matrix

disp('the subdiagonal of A= ');
disp(a);
disp('the main diagonal of A= ');
disp(b);
disp('the superdiagonal of A= ');
disp(c);
%compute coefficient matrix B
d(1)=-((h^2)*feval(r,aa+h))+y0*(1+((h/2)*feval(p,aa+h)));
d(n-1)=-((h^2)*feval(r,bb-h))+yn*(1-((h/2)*feval(p,bb-h)));
for i=2:n-2
    x=aa+i*h;
    d(i)=-((h^2)*feval(r,x));
end
fprintf('\n');
disp('the coefficinet matrix B= ');
disp(d);
disp('the solution of BVP= ');
fprintf('\n');
disp('xi        yi          Yi          error');
fprintf('\n');
y(2:n)=tridiag(a,b,c,d,n-1);
%sum=0;
%ratio=1;
for i=1:(n+1)
    x=aa+(i-1)*h;
    Yappr(i,:)=y(i);
    Xsol(i,:)=x;
    Ysol(i,:)=feval(Y,x);
    err(i)=abs(y(i)-feval(Y,x));
    %ratio=0;
    %if(i>2)
    %    ratio=err(i)/err(i-1);
    %end    
    fprintf('%f  %f  %f  %f\n',x,y(i),feval(Y,x),err(i));
    %plot(x,y(i),'r')
end
error=max(err);
%sum=sum/n;
figure;plot(Xsol,Yappr,'r')
hold on;
plot(Xsol,Ysol,'b');
hold off;
title('graph for ques1')
legend('red=approx yi','blue=exact Yi')
xlabel(['x(i) at h=' num2str(h)]),ylabel('y(x(i)) or Y(x(i))'),title('QUES1')



%figure;
%xlabel('x(m)'),ylabel('y(m)'),title('deflection vs position'),grid;


end

