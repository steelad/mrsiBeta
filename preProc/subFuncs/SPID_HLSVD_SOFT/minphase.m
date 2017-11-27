function [hn] = minphase(h)

%changed by JB Poullet 15 December 2006     %%%
r = cplxpair(roots(h));
%function [hn] = minphase(h,r)
%r = cplxpair((r));                       %%%


% Find roots inside & on unit circle

ra = abs(r);

iru  = find(abs(ra-1) <= 1e-8);	% indices for roots on uc
irin = find(ra-1 < -1e-8);	% indices for roots inside uc

% Delete roots outside unit circle, double those inside:

cr = [r(iru) ; r(irin) ; r(irin)];

% Combine complex-conj roots: 

j=1;
k=1;
while(j<=length(cr))
	if(imag(cr(j))~=0)
		cc(k,:) = conv([1 -cr(j)],[1 -cr(j+1)]);
		j=j+2;
		k=k+1;
	else
		cc(k,:) = [0 1 -cr(j)];
		j = j+1;
		k=k+1;
	end;
end;
if(rem(length(cc),2)~=0) cc(k,:) = [ 0 0 1 ]; end;

% Get reference spectrum:

Hi = freqz(h,1,512);
K1 = abs(Hi);
K1 = 10*log10(K1/max(K1));
for i=1:512 
	if(K1(i)<-100) K1(i)=-100; end;
end;

% Expand polynomial:

c=1;
while(size(cc)>0)
	k=0;
	l = size(cc);
	for j=1:l
		ct = conv(c,cc(j,:));
		H = freqz(ct,1,512);

		K2 = abs(H);
		K2 = 10*log10(K2/max(K2));
		for i=1:512 
			if(K1(i)<-100) K2(i)=-100; end;
		end;
		k(j) = sum(abs(K1-K2));
	end;

	[minv,mini] = min(k);

	cq = cc(mini,:);

	if(mini==1)
		cc = cc(2:l,:);
	elseif(mini==l)
		cc = cc(1:l-1,:);
	else
		cc = cc([1:mini-1 mini+1:l],:); 
	end;

	k=0;
	l = size(cc);
	for j=1:l
		ct = conv(cq,cc(j,:));
		H = freqz(ct,1,512);
		k(j) = max(abs(H));
	end;

	[minv,mini] = min(k);

	cq  = conv(cq,cc(mini,:));
	c  = conv(c,cq);

	if(mini==1)
		cc = cc(2:l,:);
	elseif(mini==l)
		cc = cc(1:l-1,:);
	else
		cc = cc([1:mini-1 mini+1:l],:); 
	end;
end;

% Strip off extraneous leading zeros & normalize power:

s = find(c~=0);
hn = real(c(s(1):length(c)));
hn = sqrt(sum(h.*h)/sum(hn.*hn))*hn;

if(0)
	rhn = roots(hn);
	hold off;
	plot(r,'m+');
	hold on;

	%zplane(h,1);
	plot(sin((1:1000)*2*pi/1000)+cos((1:1000)*2*pi/1000)*sqrt(-1),'y-');
	plot(rhn,'w+');
	hold off;
	error
end;

return;
