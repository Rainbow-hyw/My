function [error] = ques2(n )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
p=@(x)(0);                                           %p(x)
q=@(x)(1);                                                      %q(x)
r=@(x)(6*x-x^3);                               %r(x)
aa=0;bb=1;y0=0;yn=3;Y=@(x)x^3;                  %a,b,y(a),y(b),Y(X)
fprintf('\n');

%h=1/2^(n-1);
h=(bb-aa)/n;
y=zeros(n+1,1);
y(1)=y0;
%y(n+1)=yn;
blank='  ';
disp(['finite difference method with h= ',num2str(h)]);
fprintf('\n');
%compute tridiagonal matrix
%a(1)=0;
for i=1:n   %y(2) is correponding to x(1) so i:1 to n-1
    x=aa+i*h;
    if(i~=1)
        a(i)=-1-((h/2)*feval(p,x));      %a2,a3,....,a(n-1)
    end
    b(i)=2+((h^2)*feval(q,x));              %b1,b2...,b(n-1)
    if(i~=n)
        c(i)=-1+((h/2)*feval(p,x));         %c1,c2,...,c(n-1)
    end
end
a(n)=-2;
%compute diagonal and subdiagonal entries of tridiagonal matrix

disp('the subdiagonal of A= ');
disp(a);
disp('the main diagonal of A= ');
disp(b);
disp('the superdiagonal of A= ');
disp(c);
%compute coefficient matrix B
d(1)=-((h^2)*feval(r,aa+h))+y0*(1+((h/2)*feval(p,aa+h)));
d(n)=-((h^2)*feval(r,bb-h))+yn*(2*h-h*h*feval(p,bb-h));
for i=2:n-1
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
y(2:n+1)=tridiag(a,b,c,d,n);
sum=0;
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
figure;plot(Xsol,Yappr,'r')
hold on;
plot(Xsol,Ysol,'b');
hold off;
title('graph for ques1 ')
legend('red=approx yi','blue=exact Yi')
xlabel(['x(i) at h=' num2str(h)]),ylabel('y(x(i)) or Y(x(i))'),title('QUES2')

end



