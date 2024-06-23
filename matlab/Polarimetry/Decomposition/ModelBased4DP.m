function [ms,mv,alpha,delta] = ModelBased4DP(C2, transmitPol)

arguments
    C2 (:,:,2,2)
    transmitPol string
end

switch transmitPol
    case 'H'
        [ms,mv,alpha,delta] = impl_H(C2);
    case 'V'
        [ms,mv,alpha,delta] = impl_V(C2);
    otherwise
        error('错误的发射极化');
end

end


function [ms,mv,alpha,delta] = impl_H(C2)

[height,width,~,~] = size(C2);
s1 = squeeze(abs(C2(:,:,1,1))+abs(C2(:,:,2,2)));
s2 = squeeze(abs(C2(:,:,1,1))-abs(C2(:,:,2,2)));
s3 = squeeze(2*real(C2(:,:,1,2)));
s4 = squeeze(2*imag(C2(:,:,1,2)));
a = 0.75;
b = -2*s1+0.5*s2;
c = s1.^2-s2.^2-s3.^2-s4.^2;

x1 = (-b+sqrt(b.^2-4*a.*c))./(2*a);
x2 = (-b-sqrt(b.^2-4*a.*c))./(2*a);
mv = zeros(height, width);
for j=1:width
    for i=1:height
        if x1(i,j)<s1(i,j)
            mv(i,j) = x1(i,j);
        else
            if x2(i,j)<0
                mv(i,j) = 0;
            else
                mv(i,j) = x2(i,j);
            end
        end
    end
end
assert(class(mv)=="double");
assert(all(mv>=0, 'all') && all(mv<=s1, "all"));

ms = s1-mv;
s2p = (s2-0.5*mv)./ms;
s3p = s3./ms;
s4p = s4./ms;
alpha = 90*acos(s2p)/pi;
alpha(isnan(alpha)) = 0;
alpha = real(alpha);
delta =  180*angle(s3p+1i*s4p)/pi;

end


function [ms,mv,alpha,delta] = impl_V(C2)

[height,width,~,~] = size(C2);
s1 = squeeze(abs(C2(:,:,1,1))+abs(C2(:,:,2,2)));
s2 = squeeze(abs(C2(:,:,1,1))-abs(C2(:,:,2,2)));
s3 = squeeze(2*real(C2(:,:,1,2)));
s4 = squeeze(2*imag(C2(:,:,1,2)));
a = 0.75;
b = -2*s1-0.5*s2;
c = s1.^2-s2.^2-s3.^2-s4.^2;

x1 = (-b+sqrt(b.^2-4*a.*c))./(2*a);
x2 = (-b-sqrt(b.^2-4*a.*c))./(2*a);
mv = zeros(height, width);
for j=1:width
    for i=1:height
        if x1(i,j)<s1(i,j)
            mv(i,j) = x1(i,j);
        else
            if x2(i,j)<0
                mv(i,j) = 0;
            else
                mv(i,j) = x2(i,j);
            end
        end
    end
end
assert(class(mv)=="double");
assert(all(mv>=0, 'all') && all(mv<=s1, "all"));

ms = s1-mv;
s2p = (s2+0.5*mv)./ms;
s3p = s3./ms;
s4p = s4./ms;
alpha = 90*acos(-s2p)/pi;
alpha(isnan(alpha)) = 0;
alpha = real(alpha);
delta =  180*angle(s3p+1i*s4p)/pi;

end