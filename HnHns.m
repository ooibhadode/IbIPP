function [Hn,Hns]=HnHns(nelx,nely,rnmin)
inH = ones((nelx+1)*(nely+1)*(2*(ceil(rnmin)+1))^2,1);
jnH = ones(size(inH)); snH = zeros(size(inH)); k =0;
[elex,eley] = meshgrid(1.5:nelx+0.5,1.5:nely+0.5);
for in1 = 1:nelx+1
    for jn1 = 1:nely+1
        en1 = (in1-1)*(nely+1)+jn1;
        for in2 = max(in1-ceil(rnmin),1):min(in1+ceil(rnmin)-1,nelx)
            for jn2 = max(jn1-ceil(rnmin),1):min(jn1+ceil(rnmin)-1,nely)
                en2 = (in2-1)*nely+jn2; k = k+1; inH(k) = en1;
                jnH(k) = en2;
                snH(k) = max(0,rnmin-sqrt((in1-elex(jn2,in2))^2+(jn1-eley(jn2,in2))^2));
            end
        end
    end
end
Hn=sparse(inH,jnH,snH); Hns = sum(Hn,2);