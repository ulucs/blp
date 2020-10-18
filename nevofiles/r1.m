	% This program computes the random coefficeints discrete choice model
%described
% the "A Research Assistant's Guide to Discrete Choice Models of
%Demand," NBER technical
% paper #221, and "Measuring Market Power in the Ready-to-Eat Cereal
% Industry," NBER WP
% #6387.

	% Written by Aviv Nevo, May 1998.  See http://www.faculty.econ.northwestern.edu/faculty/nevo/supplements/rc_dc_code.htm

	%Eric Rasmusen, Dec. 31, 2005
 %I have added a lot of commands to display matrix size
%I have updated some of the obsolete minimization functions
%I added a diary, a log, with the date you run the program.
 % I wrote this in Matlab 6.
 % I added the separator lines


	%  I also use Bronwyn Hall's, April 2005  modifications for Matlab 7
  %compatibility. See 
%http://emlab.berkeley.edu/users/bhhall/e220c/rc_dc_code.htm

%The inputs to this code are the data matrices ps2.mat and iv.mat

%The output from this program will be results.txt and  mydiary.txt.  . 

disp('                                                                    ')
disp('**********************************************************************')
disp('                                                                    ')


delete mydiary.txt %we don't want any old mydiary.txt file to interfere.
diary mydiary.txt

date = datestr(clock)

	%The "date" command prints out the date, useful for keeping track of the diary.

global invA ns x1 x2 s_jt IV vfull dfull theta1 theti thetj cdid cdindex

	%global ns nmkt nbrn n_inst

	% load data. see description in readme.txt

load ps2.mat
load iv.mat

	 %Use the two files ps2.mat and iv.mat
% http://www.faculty.econ.northwestern.edu/faculty/nevo/supplements/readme.html
% ps2.mat contains the matrices v, demogr, x1, x2, s_jt, id_demo



disp('                                                                    ')
disp('**********************************************************************')
disp('                                                                    ')



ns = 20;       % number of simulated "indviduals" per market %
nmkt = 94;     % number of markets = (# of cities)*(# of quarters)  %
nbrn = 24;     % number of brands per market. if the numebr differs by
%market this requires some "accounting" vector %
 % Default is ns=20, nmkt=94, nbrn=24. Thus we have 24*94=2256 observations
%Set this smaller if you want to have more a dataset you can eyeball better

n_inst = 20   %Number of instruments for price


disp('                                                                    ')
disp('**********************************************************************')
disp('                                                                    ')



IV = [iv(:,2:(n_inst+1)) x1(:,2:(nbrn+1))];

	% The previous line creates a matrix IV from the instruments and the x's.
%20 columns of instruments for price by default, and nbrn brand dummies


IV= IV(1:nmkt*nbrn, :);
s_jt= s_jt(1:nmkt*nbrn, :);
x1=x1(1:nmkt*nbrn, :);
x2= x2(1:nmkt*nbrn, :);

% The previous line  reduces the data matrix, e.g.,  x1, to its first  nmkt*nbrn observations;

disp(['The dimensions of object IV are ' num2str(size(IV)) ])





clear iv


disp('                                                                    ')
disp('**********************************************************************')
disp('                                                                    ')



	% The  vector  below relates each observation to the market it is in %

cdid = kron([1:nmkt]',ones(nbrn,1));

disp(['Object cdid dimensions: ' num2str(size(cdid)) ])


	% the vector below provides for each index the of the last observation   in the data  used here all brands appear in all markets. if this   is not the case the two  vectors, cdid and cdindex, have to be     created in a different fashion but  the rest of the program works fine. 

cdindex = [nbrn:nbrn:nbrn*nmkt]';

disp(['The dimensions of object cdindex are      ' num2str(size(cdindex)) ])

disp('                                                                    ')
disp('**********************************************************************')
disp('                                                                   ')

	% starting values. zero elements in the following matrix correspond to
% coeff that will not be max over,i.e are fixed at zero.  The defaults are in %comments here

%theta2w=    [0.3302   5.4819         0    0.2037         0;
%             2.4526  15.8935    -1.2000        0    2.6342;
%             0.0163  -0.2506         0    0.0511         0;
%             0.2441   1.2650         0   -0.8091         0];


 theta2w=    [0.3302   5.4819         0    0.2037         0;
              2.4526  15.8935    -1.2000        0    2.6342;
             0.0163  -0.2506         0    0.0511         0;
             0.2441   1.2650         0   -0.8091         0];

	%There are 13 parameters to be estimated, and 7 set to zero by assumption
%These are 4 characteristic parameters and 9 interaction  parameters.
 %Some interactions are set to zero (e.g. Sugar and Incomsquared)


disp('                                                                    ')
disp(['The dimensions of object theta2w are      ' num2str(size(theta2w)) ])
 disp('                                                                    ')

	% create a vector of the non-zero elements in the above matrix, and the
% corresponding row and column indices. this facilitates passing values %
% to the functions below. %

[theti, thetj, theta2]=find(theta2w);

horz=['    mean       sigma      income   income^2    age    child'];
vert=['constant  ';
      'price     ';
          'sugar     ';
          'mushy     '];

disp('                                                                    ')
disp('**********************************************************************')
disp('                                                                   ')

	% create weight matrix

invA = inv([IV'*IV]);

disp(['The dimensions of object invA  are ' num2str(size(invA)) ])




disp('                                                                    ')
disp('**********************************************************************')
disp('                                                                   ')

	% Logit results and save the mean utility as initial values for the
%search below.

	% compute the outside good market share by market:

temp = cumsum(s_jt);
sum1 = temp(cdindex,:);
sum1(2:size(sum1,1),:) = diff(sum1);
outshr = 1.0 - sum1(cdid,:);

disp(['Object temp dimensions: ' num2str(size(temp)) ])
disp(['Object sum1 dimensions: ' num2str(size(sum1)) ])
disp(['Object outshr dimensions: ' num2str(size(outshr)) ])




y = log(s_jt) - log(outshr);
mid = x1'*IV*invA*IV';
t = inv(mid*x1)*mid*y;
mvalold = x1*t;
oldt2 = zeros(size(theta2));
mvalold = exp(mvalold);

 	%s_jt and outshr and x1 are  data that Nevo provides

disp(['Object y dimensions: ' num2str(size(y)) ])
disp(['Object mid dimensions: ' num2str(size(mid)) ])
disp(['Object t dimensions: ' num2str(size(t)) ])
disp(['Object mvalold dimensions: ' num2str(size(mvalold)) ])
disp(['Object oldt2 dimensions: ' num2str(size(oldt2)) ])


	%the next command creates a new file, mvalold.mat, with mvaoldold and oldt2 in it, and then clears out the old mvalold from memory. 

save mvalold mvalold oldt2
clear mid y outshr t oldt2 mvalold temp sum1

vfull = v(cdid,:);
dfull = demogr(cdid,:);

	% the matrix v has 80 iid normal random numbers for each of the 94 observations
 %the matrix demogr has random draws from the CPS for 20 individuals per obs.

disp(['Object vfull dimensions: ' num2str(size(vfull)) ])
disp(['Object dfull dimensions: ' num2str(size(dfull)) ])


disp('                                                                    ')
disp('**********************************************************************')
disp('                                                                   ')

	%options = foptions; obsolete matlab
	%options(2) = 0.1;obsolete matlab
	%options(3) = 0.01;obsolete matlab
	%foptions is obsolete in Matlab 6. optimset is its replacement.

	%Below,   set the maximum number of iterations for the main optimization command. Note .that the number of iterations for the meanval.m calculation is separate, and unconstrained. 

	%I can't get the Jacobian-Gradobj to work, so I turned it off. It generates a gigantic value for df in gmmobjg.m, and sends the next iteration into outer space. The problem is probably in the gmmobjg.m file, since Bronwyn Hall did not change the jacob.m file. Or, it might be in the meanval.m file somehow. 

options = optimset('GradObj','off','MaxIter',150,'Display','iter','TolFun',0.1,'TolX',0.01);


tic

	% tic starts the clock to measure  the time it takes for optimization

	%Use either (1) or (2) below, but not both, to do the minimization. Use % to comment out the command you are not using. 

	% (1) The following line computes the estimates using a Quasi-Newton method  with an   *analytic* gradient   The old Matlab command fminu is now fminunc, in versions 6 and  7. What is in  Nevo is the obsolete command:

  % OBSOLETE: % [theta2,  options] = fminu('gmmobj',theta2,options,'gradobj')

         %   The new command is below.    Fminunc is  from the  "optimization toolbox", not in the basic  matlab. What you  should use now for the Quasi-Newton method is the next line: 

  %      [theta2,fval,exitflag,output] = fminunc('gmmobjg',theta2,options);

        %gmmobjg.m  is a  separate matlab file created by Bronwyn Hall to %replace  Nevo's gmmobj.m and gradobj.m with one file for Matlab 6 or 7.
 %gmmobjg calls up meanval.m

	 % (2) The following command computes the estimates using a simplex search method.
	%The old Matlab command fmins is now fminsearch, in version 7.

    [theta2,fval,exitflag,output]  = fminsearch('gmmobjg',theta2, options)


	%meanval.m has the command, "disp(['# of iterations for delta convergence:  ' num2str(i)])", which shows the iterations it takes to converge in computing the mean utility for each iteration of the main optimization routine. 

	%gmmobjg.m has the command "disp(['GMM objective:  ' num2str(f1)])", which shows the value of the moment expression at each meanval.m iteration

	%Either (1) or (2) shows various things at each optimization iteration

comp_t = toc/60;

	%toc is  the elapsed optimization time that tic started


disp('                                                                    ')
disp('**********************************************************************')
disp('                                                                   ')

	% computing the s.e.

vcov = var_cov(theta2);
se = sqrt(diag(vcov));

	%var_cov is a separate matlab file

 

disp(['Object vcov dimensions: ' num2str(size(vcov)) ])
disp(['Object se dimensions: ' num2str(size(se)) ])



theta2w = full(sparse(theti,thetj,theta2));
t = size(se,1) - size(theta2,1);
se2w = full(sparse(theti,thetj,se(t+1:size(se,1))));

disp(['Object theta2w dimensions:     ' num2str(size(theta2w)) ])
disp(['Object t dimensions:     ' num2str(size(t)) ])
disp(['Object se2w dimensions:     ' num2str(size(se2w)) ])


disp('                                                                    ')
disp('**********************************************************************')
disp('                                                                   ')

	% the MD estimates

omega = inv(vcov(2:25,2:25));
xmd = [x2(1:24,1) x2(1:24,3:4)];
ymd = theta1(2:25);

disp(['Object omega dimensions: ' num2str(size(omega)) ])
disp(['Object xmd dimensions: ' num2str(size(xmd)) ])
disp(['Object ymd dimensions: ' num2str(size(ymd)) ])


beta = inv(xmd'*omega*xmd)*xmd'*omega*ymd;
resmd = ymd - xmd*beta;
semd = sqrt(diag(inv(xmd'*omega*xmd)));
mcoef = [beta(1); theta1(1); beta(2:3)];
semcoef = [semd(1); se(1); semd];

disp(['Object omega dimensions: ' num2str(size(omega)) ])
disp(['Object beta dimensions: ' num2str(size(beta)) ])
disp(['Object resmd dimensions: ' num2str(size(resmd)) ])
disp(['Object semd dimensions: ' num2str(size(semd)) ])
disp(['Object mcoef dimensions: ' num2str(size(mcoef)) ])
disp(['Objectsemcoef dimensions: ' num2str(size(semcoef)) ])


Rsq = 1-((resmd-mean(resmd))'*(resmd-mean(resmd)))/((ymd-mean(ymd))'*(ymd-mean(ymd)));
Rsq_G = 1-(resmd'*omega*resmd)/((ymd-mean(ymd))'*omega*(ymd-mean(ymd)));
Chisq = size(id,1)*resmd'*omega*resmd;


disp('                                                                    ')
disp('**********************************************************************')
disp('                                                                   ')

diary results.txt
disp(horz)
disp('  ')
for i=1:size(theta2w,1)
     disp(vert(i,:))
     disp([mcoef(i) theta2w(i,:)])
     disp([semcoef(i) se2w(i,:)])
end

disp('                                                                    ')
disp('**********************************************************************')
disp('                                                                   ')

	%disp(['GMM objective:  ' num2str(options(8))]) old matlab version, obsolete

disp(['GMM objective:  ' num2str(fval)])
disp(['MD R-squared:  ' num2str(Rsq)])
disp(['MD weighted R-squared:  ' num2str(Rsq_G)])

	%disp(['# of objective function evaluations:  ' num2str(options(10))]) old matlab, obsolete

disp(['run time (minutes):  ' num2str(comp_t)])

 

diary off

disp('                                                                    ')
disp('**********************************************************************')
disp('                                                                   ')

delete gmmresid.mat 
delete mvalold.mat











