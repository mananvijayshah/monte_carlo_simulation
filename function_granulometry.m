function [report, d32] = function_granulometry(h1,clases,h1max)
%GRANULOMETRY gives the granulometric analysis of the particulate system.
%
%   Z = GRANULOMETRY(H1,CLASES) returns the granulometric report of h1 in a
%   given number of clases "CLASES". H1 is a vector.
%   Korina Terrazas Velarde 12.04.2007
    entidades=max(size(h1));
    interval=h1max/clases;
    liminf=[0:interval:(clases-1)*interval]';
    limsup=[interval:interval:clases*interval]';
    diampromedio=(limsup+liminf)/2;
    eachclass=zeros(clases,1);
    for i=1:clases
        encontrados=find(h1>liminf(i)&h1<=limsup(i));
        cuantos=size(encontrados);
        eachclass(i)=cuantos(1,1);
    end
    porcentaje=(eachclass/entidades)*100;
    d32=sum((porcentaje/100).*diampromedio);
    %d32=sum((porcentaje/100).*limsup);
    cumulativa=cumsum(porcentaje);
    q3=porcentaje/interval;
    report=cat(2,liminf,limsup,diampromedio,eachclass,porcentaje,cumulativa,q3);
    %lala=sum(report(:,4));
    
end


