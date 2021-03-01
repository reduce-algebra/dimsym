module dimsym23;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% package:                DIMSYM.                            %
%                                                            %
%          Symmetries of Differential Equations              %
%                          and                               %
%          Linear Differential Equation Solver.              %
%                                                            %
% Author:  James Sherring                                    %
%                                                            %
% modified by  Geoff Prince and Michael Jerie.               %
%                                                            %
% dim2ode module (interface for the odesolve package)        %
% Author:  Michael Jerie, Oct 1999.                          %
% Last modified: February 6th, 2004                                                           %
%                                                            %
%  Note: Dimsym redfines some Reduce functions.              %
%                                                            %
%      You can expect to see the following in the            %
%      compile log.                                          %
%                                                            %
%      *** Function `ordop' has been redefined               %
%      *** Function `ordpp' has been redefined               %
%                                                            %
%      *** @ redefined                                       %
%      *** Function `intersection' has been redefined        %
%                                                            %
%                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

create!-package('(dimsym23 dim2ode),'nil);


%% global decs for dim2ode module %%

% flags for already loaded package and package not found.
fluid '(!*loadodesolve !*noodepackage);

fluid '(!*odehardints);

%% end of global decs for dim2ode %%


%% start of the original (pre '99) dimsym code (including subsequent changes).
write "Dimsym 2.3, 12-October-1999. ";terpri();
write "Symmetry determination and linear D.E. package";terpri();
write "(c) 1992, 1993 James Sherring; 1999 James Sherring, Geoff Prince and Michael Jerie";terpri();
write "Any publication resulting from these calculations must reference this program.";
terpri();
write "Users are free to modify this program but it is not to be redistributed";terpri();
write "in modified form.";terpri();

algebraic(for all x let cos(x)**2=1-sin(x)**2);

%algebraic operator x,u,xi,phi,ideq,deq,deteqn,c,gen,sign,liebacklund,dothis;
algebraic operator x,u,xi,phi,ideq,deq,deteqn,c,gen,liebacklund,dothis;
share prosymvec,symvec,!*p,!*q,!*r;

fluid '(!*funord ldepsfound fdepsfound funknsfound alldepsfound subs
        subdepth bvars neqnflag madesub restvars !*dunkns genfns eqnsused
        ratio solveexit solvesuccess op!*sub2sftries restvar npvars);

global '(!*u !*x !*c !*xi !*phi);
!*x:='x; !*u:='u; !*c:='c; !*xi:='xi; !*phi:='phi;
global '(!*detunknowns); !*detunknowns := '(phi xi c);

% these are dynamic and keep track of things throughout a program run

global '(!*dets !*ccount !*liesubcount !*stopatnum !*dothisatnum !*traceatnum
 !*deteqns !*vars !*deqlis !*predeqlis !*basefns !*dunknsfound !*gentypes
 !*eghistory !*sqsublis !*algsublis !*usubvars !*firstusubvars
 !*freerestrictions !*linindepconditions !*specialfns !*divides !*hidenbefore
 !*hideneqns !*longhideneqns !*freeunknownatoms !*freeunknownops !*recnumdeps
 !*canonlist !*hardints !*startdets !*startccount !*op!*xdpshr !*opusage
 !*intfac!*opusage !*op!*xdp!*opusage !*splitvars !*proeqnsmade);

symbolic( !*dets := !*ccount := !*startdets := !*startccount :=0 );
symbolic( symvec := nil         ); % so that we dont do op!*rmred1/2 
                                   % if the eqns we are solving aren't
                                   % made by mkdets
symbolic( prosymvec := nil      ); % so that mkprgens(); doesn't return an error
                                   % if no eqns solved yet
global '(!*mkdetstypes !*mkdetsargtypes);
!*mkdetstypes:='(point custom1 custom2 liebacklund);
!*mkdetsargtypes:='(liebacklund);




%%%%%%%%%%% Some algorithm flags

global '( !*splitonsimp !*op!*sub2sfby1 !*op!*findfrby1 !*op!*findlgby1 
  !*allics  !*factordivides !*intonlyslvbl  !*sqliesub  !*autosplitdets
  !*subeqnsinprovec !*subeqnsinprocoefs !*mkproatstart  !*simpsymvecwithsub
  !*logspdiv  !*intcon1subnew !*subonintcon1 !*intcon2subnew !*subonintcon2
  !*op!*exp1tmexpand !*gcdreorder !*op!*intfacfirstorderonly !*holdvalues
  !*ignoreall !*keepeqngrps !*myezgcd !*subeqnsinreverse !*op!*get1tmonless
  !*intcon2onlowdunkn !*intfdeps !*keepints !*keep1tms !*canonbyshortest
  !*slvforhigh !*op!*xdpon0 !*keephiding !*keepsubing  !*intfaconfree
  !*intcon2bylowsdfdunkn !*subeqnsbydunkn !*i2highlong !*intcon2onmin
  !*intcon1by1 !*intcon2by1 !*splitonpro !*slvpro);

fluid '(!*slvdfs !*forcelowerprovalues);
symbolic (!*slvdfs:=nil);           % Permit solving for values
                                    % with df of dunkns in val?
                                    % This _should_ be set/cleared by ops...
symbolic (!*forcelowerprovalues:=t);% Unless using op!*proeqn
                                    % or pro phis have dependence set 

%% ok for user to change these, but could alter algorithm quite a lot...

symbolic (!*keepints:=nil);         % Do we keep explicit integrals
                                    % (not including dunkns)
symbolic (!*keepsubing:=t);         % Do we keep subbing in op!*sub2sf ?
symbolic (!*intfaconfree:=t);
symbolic (!*intcon2bylowsdfdunkn:=nil); % Set to true after testing...
symbolic (!*i2highlong:=nil);
symbolic (!*intcon2onmin:=t);       % We only want to calculate minimal
                                    % int cons, ala Reid.
symbolic (!*intcon1by1:=nil);       % Do we want to calculate more than
                                    % 1 intcon1 at a time?
symbolic (!*intcon2by1:=nil);       % Do we want to calculate more than
                                    % 1 intcon2 at a time?
symbolic (!*intcon2onlowdunkn:=nil);
symbolic (!*op!*sub2sfby1:=t);  
symbolic (!*op!*findfrby1:=nil);
symbolic (!*subeqnsbydunkn:=t);
symbolic (!*op!*findlgby1:=t);
symbolic (!*factordivides:=t);      % to pull out only important information
%symbolic (!*logspdiv:=t);          % log all divides with specialfns


%% I don't recomend changing these, could destabilise algorithm
%% Some of these are now redundant...

symbolic (!*splitonpro:=nil);       % Do we split eqns with unexp'd pro phis?
                                    % Only if their dependance is set...
symbolic (!*slvpro:=nil);           % Do we solve eqns with unexp'd pro phis?
symbolic (!*op!*get1tmonless:=t);   % Only get 1tm eqns with sep'n of vars?
symbolic (!*allics:=t);             % Add all integrability cond'ns,
                                    % not just one level
symbolic (!*sqliesub:=t);           % Doesn't seem to make much difference...
symbolic (!*autosplitdets:=nil);    % Whether splitdets is run after making
                                    % determining equations if type='point.
symbolic (!*subeqnsinprovec:=nil);  % Dont need to if we sub in procoeffs
symbolic (!*subeqnsinprocoefs:=t);  % Maybe faster???
symbolic (!*mkproatstart:=nil);       %
symbolic (if null !*mkproatstart then !*subeqnsinprocoefs:=t);
                                    % Need to be more sure than this ...
symbolic (!*simpsymvecwithsub:=t);  % Maybe faster???
symbolic (!*subonintcon1:=nil);     % Give new eqns a chance to go around,
                                    % but could be redund.
symbolic (!*subonintcon2:=nil);     % Give new eqns a chance to go around,
                                    % but could be redund.
symbolic (!*op!*exp1tmexpand:=nil);  % This doesnt help to solve,
                                    % and leads to redundancies 
symbolic (!*gcdreorder:=nil);
symbolic (!*op!*intfacfirstorderonly:=nil);
symbolic (!*holdvalues:=nil);       % In mycleark, only used by redund.
                                    % check ops
symbolic (!*ignoreall:=nil);        % Flag to override logging of divides
symbolic (!*keepeqngrps:=nil);      % Used in debugging only
symbolic (!*keep1tms:=nil);         % Too Messy to keep 1 term equations
                                    % (after they have subed)
symbolic (!*myezgcd:=t);            % Determines whether ezgcd is used in
                                    % gcd and factoring calculations. Much
                                    % faster for long gcds if space avail.
symbolic (!*subeqnsinreverse:=nil); % For putting in standard form.
symbolic (!*intfdeps:=t);           % Can do this if we have intpatch loaded;
symbolic (!*canonbyshortest:=nil);  % 
symbolic (!*slvforhigh:=t);         % So we dont need to simp !*canonlist.

%% NEVER change these flags!!!

symbolic (!*splitonsimp:=nil);      % Dont permit splitting before simping yet,
                                    % it stuffs up.
symbolic (!*op!*xdpon0:=nil);       % Invalid to sub 0/1 for explicit vars
                                    % we want to restrict by.
symbolic (!*intonlyslvbl:=t);       % Only integrate eqns if then solveable
symbolic (!*intcon1subnew:=nil);    % Don't sub this into others, its not
                                    % in stand. form yet
symbolic (!*intcon2subnew:=nil);    % Don't sub this into others, its not
                                    % in stand. form yet

%%%%%%%%%%% Some algorithm parameters, ok for user to change

global '(!*eqnlengthlimit !*factorizelimit !*liesublimit !*intnnumlimit
 !*intndenlimit !*asmtnumlimit !*asmtdenlimit !*exp1tmordlimit
 !*intfacnumlimit !*intfacdenlimit);

symbolic (!*eqnlengthlimit:=100000);% Put aside anything longer than this...
symbolic (!*factorizelimit:=10000); % Factorize equations with fewer terms
                                    % than this.
symbolic (!*liesublimit:=100);       % As high as neccessary, not infinite...
                                    % Recursion limit for levels of subbing
                                    % eqns into themselves when making
                                    % determining eqns.
symbolic (!*exp1tmordlimit:=20);    % Lim on ord of polys made by op!*exp1tm
symbolic (!*intnnumlimit:=100);     % Limit on the num. of what we integrate
symbolic (!*intndenlimit:=4);       % Limit on the den. of what we integrate
symbolic (!*asmtnumlimit:=100);     % Limit on the num. of assignments made
                                    % by op!*slvall
symbolic (!*asmtdenlimit:=5);       % Limit on the den. of assignments made
                                    % by op!*slvall
symbolic (!*intfacnumlimit:=4); 
symbolic (!*intfacdenlimit:=2);


%%%%%%%%% declarations for solve operations...

% op!*shr1tm   solve single term eqn- set to zero or remove dependence
% op!*simpeq   simplify eqn if flaged to do so
% op!*splitEc  split eqn- only by lvars not causing restrictions on funkns
% op!*exp1tm   solve single term eqn-df by one var only, set dunkn to poly
% slv          solve for a dunkn and set that value
% op!*intslv   solve for a dfdunkn and integrate to get value for dunkn
% op!*trgexp   2nd order const coeff ode with pure trig or pure exp sol'ns
% op!*intfac   first order integrating factor method
% op!*splitEd  split eqn- posibly by fvars causing restrictions on funkns
% op!*sub2sf   put !*deteqns into standard form
% addallintcons  add all integrability conditions to !*deteqns
% op!*inttbta  term by term integration

global '(!*solvedetsalgorithms !*op!*xdpslvops !*op!*intfacslvops);

!*op!*xdpslvops:='(op!*shr1tm  op!*slvtwo  op!*exp1tm
                   op!*trgexp  op!*intfac  op!*intslv  op!*slvall);

!*op!*intfacslvops:='(op!*shr1tm  op!*slvtwo  op!*exp1tm  op!*trgexp
                      op!*intfac  op!*intslv  op!*slvall  ); 

% wow, recursive calls to op!*intfac! :-)


!*solvedetsalgorithms:='(simp   proxpnd  stdsplit  split    alansplit  std   
                         std1   easy     nosub     stdform  stdform1   nofxp 
                         hidef  findf    hidel     findl    inttbt     rmvred
                         new    pro2sf   proeqns   standardform   stdformx);


put('simp,     'solveops, '(op!*simpeq  ));
put('proxpnd,  'solveops, '(op!*simpeq  op!*proexp  ));
put('stdsplit, 'solveops, '(op!*simpeq  op!*splitEb op!*proexp  ));
put('split,    'solveops, '(op!*shr1tm  op!*simpeq  op!*splitEc op!*proexp)); 
put('alansplit,'solveops, '(op!*proexp  op!*shr1tm   op!*simpeq   op!*splitEb));

put('new,'solveops,'(op!*shr1tm   op!*splitEc  op!*simpeq     op!*proeqn
                     op!*slvtwo   op!*get1tm   op!*hidefr     op!*hidelg
                     op!*exp1tm   op!*slvspl   op!*trgexp
                     op!*intfac   op!*intslv   op!*xdpxpd     op!*slvall
                     op!*splitEd  %op!*slvalldfsearly
                     op!*sub2sf   addintcons
                     op!*findlg   op!*findfr   op!*slvalldfs  ));

put('std,'solveops,'(op!*shr1tm   op!*splitEc  op!*simpeq     op!*proexp
                     op!*slvtwo   op!*get1tm   op!*hidefr     op!*hidelg
                     op!*exp1tm   op!*slvspl   op!*trgexp
                     op!*intfac   op!*intslv   op!*xdpxpd     op!*slvall
                     op!*splitEd  %addintcon1   op!*slvalldfsearly
                     op!*sub2sf   addintcons
                     op!*findlg   op!*findfr   op!*slvalldfs  ));

put('std1,'solveops,'(op!*shr1tm   op!*splitEc  op!*simpeq     op!*proexp
                      op!*slvtwo   op!*get1tm   op!*hidefr     op!*hidelg
                      op!*exp1tm   op!*slvspl   op!*trgexp
                      op!*intfac   op!*intslv   op!*get1tm9    op!*xdpshr
                      op!*xdpxpd   op!*slvall   op!*slvalldfs  op!*splitEd  
                      addintcon1   op!*sub2sf   addintcon2     op!*findlg   
                      op!*findfr   ));

put('easy, 'solveops,'(op!*shr1tm   op!*splitEc  op!*simpeq     op!*proexp
                       op!*slvtwo   op6          ));  

put('nosub, 'solveops,'(op!*shr1tm     op!*splitEc  op!*simpeq   op!*proexp
                        op!*slvtwo     op!*get1tm   op!*hidefr   op!*hidelg
                        op!*exp1tm     op!*slvspl   op!*trgexp
                        op!*intfac     op!*intslv   op!*xdpxpd   op!*slvall
                        op!*splitEd    addintcon1   op!*findlg   op!*findfr
                        op!*slvalldfs  ));

put('standardform,  'solveops,'(op!*sub2sf  addintcons  ) );

put('stdform,  'solveops,'(op!*simpeq  op!*proexp  op!*sub2sf  addintcons  ) );

% New solvedets algorithm, stdformx, placed on the !*solvedetsalgorithms list.
% It does splitting before op!*sub2sf.

put('stdformx,  'solveops,'(op!*simpeq  op!*splitEb op!*proexp  op!*sub2sf  
                            addintcons  ) );

put('stdform1, 'solveops,'(op!*shr1tm  op!*simpeq  op!*splitEc  op!*splitEd
                           op!*proexp  op!*slvtwo  op!*get1tm  %op!*exp1tm
                           op!*slvspl  op!*hidefr  op!*hidelg  op!*sub2sf 
                           addintcons  op!*findlg  op!*findfr   ));

put('nofxp, 'solveops,'(op!*shr1tm   op!*splitEc  op!*simpeq     op!*proexp
                        op!*slvtwo   op!*get1tm   op!*exp1tm     op!*slvspl   
                        op!*trgexp   op!*intfac   op!*intslv   
                        op!*xdpxpd   op!*slvall   %op!*splitEd    addintcon1  
                        %op!*sub2sf   addintcon2   op!*slvalldfs
                        ));

put('hidef,  'solveops,'(op!*hidefr) );
put('findf,  'solveops,'(op!*findfr) );
put('hidel,  'solveops,'(op!*hidelg) );
put('findl,  'solveops,'(op!*findlg) );
put('inttbt, 'solveops,'(op!*inttbta));
%put('pro2sf, 'solveops,'(op!*proeqn op!*splitEc op!*splitEd op!*sub2sf addintcons));
put('pro2sf, 'solveops,'(op!*shr1tm  op!*simpeq  op!*splitEc  op!*splitEd
                           op!*proeqn  op!*slvtwo  op!*get1tm  %op!*exp1tm
                           op!*slvspl  op!*hidefr  op!*hidelg  op!*sub2sf 
                           addintcons  op!*findlg  op!*findfr   ));
put('proeqns,'solveops,'(op!*proeqn));


flag('(op!*proexp op!*linchk op!*shr1tm op!*simpeq op!*splitE op!*splitEa
  op!*splitEb op!*splitEc op!*splitEd op!*exp1tm op6 op!*xdpxpd op!*xdpshr 
  op!*xdp op!*get1tm op!*get1tm9 slv op!*slvshr op!*slvall op!*slvalldfs
  op!*intslv op!*inttbt op!*inttbta op!*intfac op!*trgexp op!*hidelg
  op!*slvtwo op!*drpdep op!*slvspl op!*odetst op!*slvalldfsearly op!*proeqn),
 'onebyoneop);

flag('(op!*proexp op!*linchk op!*shr1tm op!*simpeq op!*splitE op!*splitEa
 op!*splitEb op!*splitEc op!*splitEd op6 op!*xdpxpd op!*xdpshr op!*xdp
 op!*get1tm op!*get1tm9  op!*exp1tm slv op!*slvshr op!*slvall op!*slvalldfs
 op!*intslv op!*inttbt op!*inttbta op!*intfac op!*slvalldfsearly op!*trgexp
 op!*findfr op!*findlg op!*slvtwo op!*drpdep op!*slvspl op!*rmred1 op!*rmred2
 op!*sub2sf addintcons addallintcons addintcon2 addintcon1 op!*proeqn),
 'restart);

flag('(op!*proexp op6 op!*xdp op!*get1tm op!*get1tm9 op!*xdpxpd op!*xdpshr
  op!*drpdep op!*proeqn),'keep);
  % keep eqns even when solved by these ops...

flag('(op!*simpeq),'dolots); % these ops are repeated (lots)
%put('op?,'selectfn,'lowordlowlen);


%%%%%%%%%%%%%%%%%%%  Trace stuff

global '(!*traceliesub !*tracesymvec !*tracesoft !*tracehard !*tracecute
  !*tracetmstmp  !*traceeqngrp !*traceeqngrpsub  !*traceeqncnt !*traceeqns
  !*traceops !*traceopson !*traceselectfn !*tracesetk !*tracecleark
  !*tracenewarb !*tracedeath !*traceslow !*traceop!*splitE !*traceop!*intslv
  !*traceop!*trgexp !*tracecanon !*tracelength !*traceop!*sub2sf
  !*tracesubeqns !*tracemysubk  !*traceintcon1 !*traceintcon2 !*traceweight
  !*traceop!*rmred1 !*traceop!*rmred2 !*tracemkgens !*tracedunknorder
  !*tracekord!* !*traceintcalls !*traceslv !*tracehidden !*traceop6
  !*traceop!*get1tm !*traceop!*xdp !*tracetrysub
  !*tracedeplev !*tracetagsimp);

!*tracesoft:='t;         % so that showeqngrp at least shows something!

symbolic procedure trace;
  !*tracesoft:=!*tracesetk:=!*traceeqngrp:=!*tracedeath:=!*tracecute
  :=!*tracetmstmp:=!*traceops:=!*tracenewarb:=!*traceeqncnt
  :=!*tracedunknorder%:=!*tracesymvec
  :=!*tracecanon:=!*tracehidden:=!*traceop!*xdp:=!*tracedeplev:=t;

symbolic operator trace;

symbolic procedure traceall;
 !*traceliesub:=!*tracesymvec:=!*tracesoft:=!*tracehard:=!*tracecute
 :=!*tracetmstmp:= !*traceeqngrp:=!*traceeqngrpsub:= !*traceeqncnt
 :=!*traceeqns:=!*traceops:=!*traceopson:=!*traceselectfn:=!*tracesetk
 :=!*tracecleark:=!*tracenewarb:=!*tracedeath:=!*traceslow:=!*traceop!*splitE :=!*traceop!*intslv:=!*traceop!*trgexp:=!*tracecanon:=!*tracelength
 :=!*traceop!*sub2sf:=!*tracesubeqns:=!*tracemysubk:=!*traceintcon1
 :=!*traceintcon2:=!*traceweight:=!*traceop!*rmred1:=!*traceop!*rmred2
 :=!*tracemkgens:=!*tracedunknorder:=!*tracekord!*:=!*traceintcalls
 :=!*traceslv:=!*traceop6:=!*traceop!*get1tm:=!*traceop!*xdp:=!*tracetagsimp
 :=!*tracedeplev:=t;

symbolic operator traceall;

symbolic procedure tracecute;
  !*tracecute:=!*tracecanon:=!*tracehidden:=!*tracedeplev:=!*tracetmstmp:=t;

symbolic operator tracecute;

symbolic procedure notrace;
  !*tracesetk:=!*traceeqngrp:=!*tracedeath:=!*tracecute
  :=!*traceops:=!*tracenewarb:=!*traceeqncnt:=!*tracedunknorder
  :=!*tracesymvec:=!*tracehidden:=nil;
 
symbolic operator notrace;

symbolic procedure traceat n;
% start tracing at equation n
<<if not fixp n then rederr("Bad integer"); 
  if n=0 then !*traceatnum:=nil else !*traceatnum:=n-1;
  n>>;

symbolic operator traceat;

symbolic procedure stopat n;
% stop solving equations when up to eqn n
<<if not fixp n then rederr("Bad integer"); 
  if n=0 then !*stopatnum:=nil else !*stopatnum:=n-1;
  n>>;

symbolic operator stopat;

symbolic procedure dothisat n;
% evaluate dothis() when up to eqn n
<<if not fixp n then rederr("Bad integer"); 
  if n=0 then !*dothisatnum:=nil else !*dothisatnum:=n;
  n>>;

symbolic operator dothisat;
 
%%%%%%%%%%% Verify stuff

global '(!*verifyop!*splitE !*verifyint !*verifygens !*verifyop!*xdp
 !*verifymysetk !*verifylin);

symbolic procedure verify;
  !*verifyop!*splitE := !*verifyint := !*verifygens := !*logspdiv 
  := !*verifymysetk := t ;

symbolic operator verify;

symbolic procedure noverify;
  !*verifyop!*splitE := !*verifyint := !*verifygens := !*verifymysetk
  := !*verifylin := nil ;

symbolic operator noverify;

lisp(!*verifylin := t);                  % Just to be quite sure...

%%%%%%%%%

symbolic procedure help;
begin
  terpri();terpri();
  write "Welcome to Dimsym 2.0 By James Sherring (J.Sherring@latrobe.edu.au)"; terpri();
  terpri();
  write "This Help is under development";terpri();
  terpri();
  write "Some info about this version of dimsym...";terpri();
  write "mkdets takes the single arguments ",!*mkdetstypes;terpri();
  write " or specify the order desired with ",!*mkdetsargtypes;terpri();
  terpri();
  write "solvedets uses one the following algorithms as argument: ",!*solvedetsalgorithms;
  terpri();terpri();
  end;

symbolic operator help;

%*****end of init

symbolic procedure freeunknown u;
% delares free unknown functions that we are not trying to determine
begin scalar v;
  for each v in u do if pairp v 
    then !*freeunknownops:=union(!*freeunknownops,list !*a2k car v)
    else !*freeunknownatoms:=union(!*freeunknownatoms,list !*a2k v);
  end;

symbolic operator freeunknown;

rlistat '(freeunknown);

symbolic procedure totder(j,u); 
mk!*sq mysimpq totderq(j,simp!* u);
symbolic operator totder;

symbolic procedure totderq(j,u); 
begin scalar vars,var,sum;
  sum:=nil ./ 1;
  vars:=union(alldepsf numr u,alldepsf denr u);
  for each var in vars do
    sum:=addsq(sum,multsq(totderk(var,j),diffsq(u,var)));
  return sum;
  end;

symbolic procedure totderf(j,u); totderq(j, u ./ 1);

symbolic procedure totderk(u,j);
if eqcar(u,!*x) then 
    if eqn(cadr u,j) 
        then (1 ./ 1) 
        else (nil ./ 1)
else if eqcar(u,!*u) then mysimpq !*k2q !*a2k (!*u.cadr u.addndx(j,cddr u))
else (nil ./ 1);

symbolic procedure addndx(j,ndx);
%
% add j to the u-index ndx.
% It is important that this is done in a consistent manner, so that there is
% a unique representation for u-indicies... ie u(1,1,2)=u(1,2,1) etc
% I use the convention of having higest to lowest, ie I use u(1,2,1).
% This must be consistent with the procedure cleanndx which cleans the input
% differential equations to use this convention.
%
if null ndx then list j
%else if atom ndx then rederr list(ndx,"invalid as index")
else if j < car ndx then (car ndx) . addndx(j,cdr ndx)
%         ^ use ordp here instead?
else j . ndx;

symbolic procedure vecder(v,u);
mk!*sq mysimpq vecderq(simp!* v, simp!* u);
symbolic operator vecder;

symbolic procedure vecderq(v,u);
% v is a standard quotient vector. u is a standard quotient.
begin scalar sum;
  sum:= nil ./ 1;
  for each term in numr v do 
    sum:=addsq(sum,multsq(diffsq(u,cadaar term),(cdr term)./1));
  return multsq(sum, 1 ./ denr v);
  end;

symbolic procedure vecderf(v,u); vecderq(v,u./1);

symbolic procedure adddfvar(u,v);
% u is a df var
% v is a df varlist.
% insert u into v
if null v then list u
else if numberp car v then  (car v).adddfvar(u,cdr v)
else if (car v)=u
  then if cdr v and numberp cadr v
    then u.(1+cadr v).(cddr v)
    else u. 2 . cdr v
  else if ordop(u,car v) 
    then u.v
    else (car v).adddfvar(u,cdr v);

symbolic procedure mydiffk(kern,var);
% *** assuming that diff of this ker will still be a kern ***
%!*q2k diffp(!*k2p kern,var);
if not pairp kern then rederr list("not dunkn in mydiffk of",kern)
else if not(var member depsk kern)
 then rederr list("bad mydiffk: ",kern," not depending on ",var)
else if (car kern) neq 'df then list('df,kern,var)
else (car kern).(cadr kern).(adddfvar(var,cddr kern));

algebraic procedure comm(a,b);
% for commutators of vectorfields
% - not the most efficient way of doing it, but easy :+)
vecder(a,b) - vecder(b,a);

symbolic operator comm;

%****** end of differentiation stuff


%****** General support for vectorfields

% must have 'partdf ordered highest in !*funorder 

symbolic procedure p2alistvector V;
% external prefix representation of vectors to internal alist representation
%
% V is a prefix form, with its vector structure given by linear combinations
% of the basis vectorfields given explicitly as "partdf(variable)"
% meaning dey by dey variable. 
%
% Convert to an alist representation of 
% ( basis vectorfield variable  . its coefficient)
% basis vectorfield variable is a kernel (!*x j) or (!*u j i1 i2 i3 ... i(n)) 
% where the j represents enumeration of the variables, and the i1 ... i(n)
% represent derivatives of the u j wrt the x i1 up to the x i(n).
% "its coefficient" is a ('!*sq expression ~)
%
% The conversion relies on vectorfield identifier 'partdf being highest in
% the ordering of the kernels in the vector V.
%
begin scalar vecnum,vecden,newvec;
  V:=simp!* V;
  vecnum:=mynumr V;
  vecden:=denr V;
  while vecnum do begin
    if  car mvar vecnum neq 'partdf then rederr list("Bad vector");
    newvec:=
      list(cadr mvar vecnum,
           mk!*sq multsq(lc vecnum ./ 1, 1 ./ vecden)).newvec;
    vecnum:=cdr vecnum;
    end;
  return newvec;
  end;

symbolic procedure alist2pvector aV;
% internal alist representation of vectors to external prefix representation
begin scalar pr,vec;
  vec:=(nil ./ 1);
  for each pr in aV do vec:=addsq(vec,multsq(!*k2q list('partdf,car pr),
                                              simp!* cadr pr
                                              ));
  return mk!*sq vec;
  end;

symbolic procedure prolong(n,V);
% V is a prefix form vectorfield
% n is the order of the prolongation being performed. 
% add all neccessary terms to V to make it the nth prolongation of V.
begin scalar new,ndx;
  V := p2alistvector V;
% for every possible index up to length n and possible values from 1 to !*p 
  ndx:='(1);
  while (count ndx)<(n+1) do begin % calculate up to the nth prolongation
    for i:=1:!*q do if null assoc(new:=cleanndx1 !*u.i.sortndx ndx,V) then
      V := list(new,mkprolongcoeff(new,V)) . V;
    ndx:=incndx ndx;
    end;        % ie if the coeff isn't there then make it
  return alist2pvector V;
  end;

symbolic operator prolong;

symbolic procedure incndx ndx;
% increment the index ndx by one. order is unimportant, so don't repeat.
if null ndx then  1 . ndx
else if (car ndx)=!*p then 1 . incndx cdr ndx
  % which will put 1's right through
else if null cdr ndx then ((car ndx) +1).cdr ndx
else if (cadr ndx)=!*p then mkndx(ndx,car ndx +1) % ie carry through
else (car ndx) . incndx cdr ndx;

symbolic procedure mkndx(u,i);
% make an indexlist of same length as u, made up of i's
if null u then nil
else i . mkndx(cdr u,i);

symbolic procedure mkprolongcoeff(uvar,V);
% uvar is u(i,ndx), with the multi-index (list) (k,J)
% V is internal alist vector.
% Calculate the coefficient of  partdf( u(i,ndx)) in the prolongation of V
% if a basis vectorfield doesn't appear in V then its coeff needs
% to be calculated as a prolongation, or it is assumed to be zero if it is
% a base-manifold basis vectorfield.
%
% just use the recursive definition of prolongation, ie
% phi(i)^(k,J) = D_k phi(i)^(J) - (sum m) u(i)_(m,J) * D_k xi(m)
%
% If the flag !*subeqnsinprocoefs then substitute the original differential
% equations into the result. This saves a lot of time later.
%
begin scalar val,old,k,i,J;
  if not (car uvar=!*u and pairp cdr uvar and pairp cddr uvar)
    then rederr list("bad vector?");
  i:=cadr uvar;
  k:=caddr uvar;
  J:=cdddr uvar;  % should I clean this first? shouldn't have to...
  old:=assoc(!*u . i . J, V);    % ie phi(i)^(J)
  if old then old:=simp!* cadr old;
  if old then val:=totderq(k,old) else val:= (nil ./ 1);
  for m:=1:!*p do begin         % !*p is the number of x vars
    old:=assoc(list(!*x,m),V);      % ie xi(m)
    if old then old:=simp!* cadr old else old:= (nil ./ 1);
    val := addsq(val,
                 multsq(negf(!*k2f !*a2k (!*u.i.addndx(m,J))) ./ 1,
                        totderq(k,old)
                        )
                 );
    end;
  if !*subeqnsinprocoefs then val:=liesubq val;
  return mk!*sq val;
  end;

%****** end vector support

%****** Some general support for eqngrp data structure

symbolic procedure egnums eglis;
% takes a list of eqngrps and returns the list of their egnums
if null eglis then nil
else (eg!*no car eglis).(egnums cdr eglis);

symbolic procedure dcs;
% show the death certificate of each eqngrp made.
% ie how and why the eqnrp has been solved or dropped.
% For debugging purposes.
begin scalar eqngrp;
  for each eqngrp in !*eghistory do 
    <<terpri();
      write(tmstmp(),"Equation ",eg!*no eqngrp," : ",eg!*dc eqngrp);
      terpri()>>;
   end;

symbolic operator dcs;

symbolic procedure nonzeroeqngrps;
% show any eqngrp made which doesnt siplify to zero.
% For debugging purposes: all eqns should be zero in the end (almost)
begin scalar eqngrp,x;
  for each eqngrp in !*eghistory do
    if (x:=mynumr mysimpf eg!*eqn eqngrp) 
      then <<write("Eqngrp(",eg!*no eqngrp,") : ");terpri();showf x >>;
  end;

symbolic operator nonzeroeqngrps;

symbolic procedure unsolvedeqns eqngrps;
% the eqngrps w/o a death certificate
% returns a list of equation numbers.
% For debugging purposes.
if null eqngrps then nil
else if (eg!*dc car eqngrps)='() 
  then (eg!*no car eqngrps).unsolvedeqns cdr eqngrps
else unsolvedeqns cdr eqngrps;

symbolic operator unsolvedeqns;

% *** just different ways to show an eqngrp

symbolic procedure showf sf;
% just shows a standard form as an algebraic structure
<<if !*tracehard then <<write sf; terpri() >>;
  if !*tracesoft then begin scalar eqn;
    eqn:=mk!*sq(sf./ 1);
    algebraic write eqn;
    end>>;

symbolic procedure showsq sq;
% just shows a standard quotient as an algebraic structure
<<if !*tracehard then <<write sq; terpri() >>;
  if !*tracesoft then begin scalar eqn;
  eqn:=mk!*sq sq;
  algebraic write eqn;
  end>>;

symbolic procedure showalg sq;
% just shows an algebraic structure
<<if !*tracehard then <<write sq; terpri() >>;
  if !*tracesoft then algebraic write sq>>;
 
symbolic procedure showeqngrp eqngrp;
begin 
  terpri();
  write(tmstmp(),"Equation ",eg!*no eqngrp," made by ",eg!*srcop eqngrp,
          " from ",eg!*bc eqngrp," has length ",eg!*len eqngrp);
  terpri();
  write ("and has leading derivative ",eg!*highdfdunkn eqngrp);
  terpri();
  if !*tracesoft then <<showf eg!*eqn eqngrp; terpri() >>;
  if !*traceeqngrpsub then 
    <<write(eg!*highdfdunkn eqngrp," = "); 
     showsq eg!*subval eqngrp; terpri() >>;
  if !*tracehard then <<write(tmstmp(),eqngrp); terpri(); terpri() >>;
  showdepsineqn eg!*eqn eqngrp;
  end;

symbolic procedure showdepsineqn eqn;
% show the dependencies of each dunkn in eqn
begin scalar dunkn,dep;
  for each dunkn in dunknsE eqn do if dep:=assoc(dunkn,depl!*)
    then <<write(dunkn," depends on ",cdr dep); terpri() >>;
  end;

%******** operations on single standard form equations ******

symbolic procedure subvalf(kern,eqn);
% the standard quotient value for kern if eqn is solved for kern
% *** it is assumed that kern appears as a toplevel mvar.**1 in list eqn
mysimpq multsq( addf(negf eqn, !*t2f((kern .** 1) .* x))./1, 1 ./ x) 
  where x=coeffE(kern,eqn);

%*** some coefficient and selector type functions

symbolic procedure coeffq(kern,u);
% find the coefficient of kern in standard quotient u
multsq(coefff(kern,numr u) ./ 1, 1 ./ denr u);

symbolic procedure coefff(kern,eqn);
% find the coefficient of kern in standard form eqn.
% It is assumed that eqn is linear in kern
if null eqn or domainp eqn or atom eqn then nil
else if mvar eqn=kern then lc eqn
% else if ordop(kern,mvar eqn) then nil % speed??
else addf( multf( !*p2f lpow eqn,coefff(kern,lc eqn) ),
           coefff(kern,red eqn));

symbolic procedure coeffE(kern,eqn);
% find the coefficient of kern in standard form eqn.
% *** it is assumed that kern appears as a toplevel mvar.**1 in list eqn
% faster then coefff
if null eqn then nil else if domainp eqn
  then rederr list("bad eqn in coeffE",kern,eqn)
else if (mvar eqn)=kern then lc eqn
else coeffE(kern,red eqn);
 
symbolic procedure dcoeffq(kern,u);
% find the coefficient of kern in standard quotient u
% **** include df of kern too
multsq(dcoefff(kern,numr u) ./ 1, 1 ./ denr u);
 
symbolic procedure dcoefff(kern,eqn);
% find the coefficient of kern in standard form eqn.
% **** include df of kern too
% It is assumed that eqn is linear in kern
if null eqn or domainp eqn then nil
else if dunknk mvar eqn and (dunknin mvar eqn)=kern then lc eqn
else addf( multf( !*p2f lpow eqn,dcoefff(kern,lc eqn) ),
           dcoefff(kern,red eqn));
 
symbolic procedure dtermq(kern,u);
% find the terms involving kern in standard quotient u
multsq(dtermf(kern,numr u) ./ 1, 1 ./ denr u);
 
symbolic procedure dtermf(kern,eqn);
% find the terms involving kern in standard quotient u
% It is assumed that eqn is linear in kern
if null eqn or domainp eqn then nil
else if dunknk mvar eqn and (dunknin mvar eqn)=kern 
  then (lt eqn) .+ dtermf(kern,red eqn)
else addf( multf( !*p2f lpow eqn,dtermf(kern,lc eqn) ),
           dtermf(kern,red eqn));
 
%**** some funcions for deciding variable dependence

symbolic procedure dunknk kern;
% is kernel kern an dunknown or derivative of one?
not atom kern and (edunknk kern  or ((car kern)='df and edunknk cadr kern));

symbolic procedure edunknk kern;
% is kernel kern explicitly a dunknown?
pairp kern and car kern member !*detunknowns;

symbolic procedure funknk kern;
% is kernel kern an funknown or derivative of one?
efunknk kern or (pairp kern and (car kern)='df and efunknk cadr kern);

symbolic procedure efunknk kern;
% is kernel kern explicitly an funkn?
% funkns can be atoms or pairs...
if atom kern
  then member(kern,!*freeunknownatoms)
  else member(car kern,!*freeunknownops);

symbolic procedure funknsp u;
% does u have any freeunkns in it anywhere?
if null u then nil
else if efunknk u then t
else if pairp u then (funknsp car u) or (funknsp cdr u);

symbolic procedure funkndepsp u;
% does u have any freeunkns with dependence in it anywhere?
if null u then nil
else if efunknk u then (depsk u or (pairp u and alldepsp u))
else if pairp u then (funkndepsp car u) or (funkndepsp cdr u);

symbolic procedure boundvarsE eqn;
% Find all variables in depl!* for the detunknowns
% appearing in standardform eqn
% Remember that the detunknowns appear linearly ie either as 
% themselves or in df statements, and at the top level of ordering.
if eqn then union(depsk mvar eqn, boundvarsE red eqn);

symbolic procedure funknsf eqn;
% find all funkns that eqn involves
if null eqn or domainp eqn
  or (null !*freeunknownatoms and null !*freeunknownops) then nil
else begin scalar funknsfound;
  funknsf1 eqn;
  return funknsfound;
  end;

symbolic procedure funknsf1 v;
if v and not domainp v
  then <<funknsk mvar v; funknsf1 lc v; funknsf1 red v>>;

symbolic procedure funknsk v;
if not dunknk v then
  if funknk v then <<funknsfound:=union(list v,funknsfound);
                     if pairp v then <<funknsa car v;funknsa cdr v>> >>
else if pairp v and ((car v)=!*x or (car v)=!*u) then nil
else funknsa v;

symbolic procedure funknsa v;
if v then
  if funknk v then <<funknsfound:=union(list v,funknsfound);
                     if pairp v then <<funknsa car v; funknsa cdr v>> >>
  else if pairp v then <<funknsa car v; funknsa cdr v>>;

symbolic procedure fdepsf eqn;
% find all variables that eqn involves from funkn terms
if null eqn or domainp eqn
  or (null !*freeunknownatoms and null !*freeunknownops) then nil
else begin scalar fdepsfound;
  fdepsf1 eqn;
  return fdepsfound;
  end;

symbolic procedure fdepsf1 v;
if v and not domainp v then <<fdepsk mvar v; fdepsf1 lc v; fdepsf1 red v>>;

symbolic procedure fdepsk v;
if not dunknk v then
if funknk v
  then <<fdepsfound:=union(depsk v,fdepsfound);
         if pairp v then fdepsa cdr v>>
else if pairp v and ((car v)=!*x or (car v)=!*u) then nil
else fdepsa v;

symbolic procedure fdepsa v;
if v then
  if funknk v then fdepsfound:=union(depsk v,fdepsfound)
  else if pairp v then <<fdepsa car v; fdepsa cdr v>>;

symbolic procedure funknatomswithdepsf eqn;
% find any funkns which are atoms (not operators)
% and which have implicit dependence.
begin scalar keep;
  for each funkn in funknsf eqn do
    if atom funkn and depsk funkn then keep:=funkn.keep;
  return keep;
  end;

symbolic procedure ldepsf u;
%find all x and u variables that standardform u involves explicitly.
begin scalar ldepsfound;
  ldepsf1 u;
  return ldepsfound;
  end;

symbolic procedure ldepsf1 v;
if v and not domainp v then <<if not dunknk mvar v then ldepsa mvar v;
                              ldepsf1 lc v; ldepsf1 red v>>;

symbolic procedure ldepsa v;
if not pairp v then nil
else if ((car v)=!*x or (car v)=!*u)
  then ldepsfound:=union(list v,ldepsfound)
else <<ldepsa car v; ldepsa cdr v>>;

symbolic procedure alldepsp u;
% find all x and u variables that prefix form u involves,
% explicitly or implicitly
begin scalar alldepsfound;
  alldepsa u;
  return alldepsfound;
  end;

symbolic procedure alldepsf u;
% find all x and u variables that standardform u involves,
% explicitly or implicitly
begin scalar alldepsfound;
  alldepsf1 u;
  return alldepsfound;
  end;
 
symbolic procedure alldepsf1 v;
if v and not domainp v 
  then <<alldepsa mvar v; alldepsf1 lc v; alldepsf1 red v>>;

symbolic procedure alldepsa v;
if v % and not domainp v 
     % ?why not this ? domainp(mk!*sq !*k2q '(x 1))=t !
     % only a test for sf's!
  then begin scalar x;
  if pairp v
    then if ((car v)=!*x or (car v)=!*u) 
      then return (alldepsfound:=union(list v,alldepsfound))
    else if (car v)='df 
      then return alldepsa cadr v;
  if (x:=assoc(v,depl!*))
    then return (alldepsfound:=union(cdr x,alldepsfound));
  if pairp v then <<alldepsa car v; alldepsa cdr v>>;
  end$

%*** some 'whats in' selector functions

symbolic procedure dunknin dfdunkn;
% dfdunkn is either a dunkn or df of one. Return the dunkn.
if (car dfdunkn)='df then cadr dfdunkn else dfdunkn;

symbolic procedure depsk kern;
% the variables the kern (or the kern in the df) depends on
(if xx then cdr xx) 
  where xx=if pairp kern and (car kern)='df
             then assoc(cadr kern,depl!*)
             else assoc(kern,depl!*);

symbolic procedure diffd(eqn,dunkn);
% Is there a df of dunkn term in eqn?
eqn and ( (car mvar eqn='df and cadr mvar eqn=dunkn)
  or diffd(red eqn,dunkn) );

symbolic procedure dunknsE eqn;
% just the actual dunkns
begin scalar !*dunkns;
  dunknsE1 eqn;
  return !*dunkns;
  end;

symbolic procedure dunknsE1 eqn;
if eqn then 
  <<!*dunkns:=union(list dunknin mvar eqn,!*dunkns); dunknsE1 red eqn>>;

symbolic procedure dfdunknsE eqn;
% the df form of the dunkns
% why bother with a routine for this? consistancy, my friend :-)
for each u on eqn collect mvar u;

symbolic procedure highphisE eqn;
% the high order phis in the eqn (ie prolongations)
begin scalar !*dunkns;
  highphisE1 eqn;
  return !*dunkns;
  end;

symbolic procedure highphisE1 eqn;
if eqn then 
  <<if car dunknin mvar eqn=!*phi and cddr dunknin mvar eqn 
      then !*dunkns:=union(list dunknin mvar eqn,!*dunkns);
    highphisE1 red eqn>>;

symbolic procedure dforderE eqn;
% the highest order derivative in eqn
begin integer ord;
  while eqn do begin
    ord:=max2(ord,dford mvar eqn);
    eqn:=red eqn;
    end;
  return ord;
  end;    

symbolic procedure highphiorderE eqn;
% the highest order phi prolongation in eqn
begin integer ord;
  while eqn do begin
    if car dunknin mvar eqn=!*phi
      then ord:=max2(ord,count cddr dunknin mvar eqn);
    eqn:=red eqn;
    end;
  return ord;
  end;    

symbolic procedure specialfnsE u;
% does standard form equation u involve any special functions?
% special functions are any kernels which are pairs,
% and which arent on the list !*ignorefns
u and ((specialfnsf lc u) or (specialfnsE red u));

symbolic procedure specialfnsf u;
% are there any special functions in standard form u?
u and pairp u and 
  ((specialfnsk mvar u) or (specialfnsf lc u) or (specialfnsf red u));

symbolic procedure specialfnsk u;
% is kernel u a special function?
if pairp u and not 
  ((car u)=!*u or (car u)=!*x or dunknk u or funknk u
  or (u member !*specialfns))
  then begin
    write(tmstmp(),"While making eqngrp ",!*dets,
          ", found special function:"); terpri();
    if !*tracehard then <<write(u); terpri() >>;
    if !*tracesoft then <<algebraic write(u); terpri() >>;
    !*specialfns:=u.!*specialfns;
    end;

symbolic procedure dropeqnfacE eqn;
% drop any factors from standard form eqn
% ie divide by GCD of coeffs
if eqn then begin scalar fac,ans,m,mm;
  mm:=count eqn;
  m:=mytermsf eqn;
  if m>!*factorizelimit 
    then <<if !*tracecute then write "$",mm,"%",m,":-(";
           return eqn>>;
  if !*tracecute then write("$",mm,"%",m);
  dividelog (fac:=eqngcdE eqn);
  if null red eqn
    then ans:=!*k2f mvar eqn
    else ans:=myquotf1E(eqn,fac);
  if !*tracecute then if fac=1
    then write "(0)"
    else write("(",mytermsf fac,")",mytermsf ans);
  if null ans then rederr("Division failed in dropeqnfacE");
  return ans;
  end;

symbolic procedure dropfacsqE u; 
% remove any factors in standard quotient u;
if null numr u then u
else begin scalar fac;
  fac:=eqngcdE numr u;
  fac:=mygcdf!*(fac,denr u);
  dividelog fac;
  if !*tracecute and (fac neq 1) then write "(\",mytermsf fac,"\)";
%    then <<write "When forming sq, got factor of length ",
%                  mytermsf fac,": ";terpri();
%           showf fac>>;
  return myquotf1E(numr u,fac)./quotf1(denr u,fac);
  end;

symbolic procedure shorterp(u,v);
%are there less terms in u than v?
mytermsf(u) < mytermsf(v);

symbolic procedure eqngcdE eqn;
% the greatest common divisor of the coefficients of eqn.
% Use this instead of gcdf(lc eqn,red eqn)
% because dunkns are funordered, which stuffs gcd up.
if null eqn then rederr("Taking eqngcdE of null eqn???")
else if domainp eqn then rederr list("Taking eqngcdE of bad eqn ",eqn)
else if not dunknk mvar eqn 
  then rederr list("Taking eqngcdE of eqn with bad mvar ",mvar eqn,
                   " in ",eqn)
else if null red eqn then lc eqn
%else mygcdf!*(lc eqn,eqngcdE red eqn);
else begin scalar lclist,rest,gcd;
  lclist:=for each rest on eqn collect lc rest;
  lclist:=sort(lclist,'shorterp);
  gcd:=car lclist;
  lclist:=cdr lclist;
  while lclist do begin
%    if !*tracecute and !*traceslow
%      then write ".",mytermsf gcd,"%",mytermsf car lclist;
    gcd:=mygcdf!*(gcd,car lclist);
    lclist:=cdr lclist;
%    if !*tracecute and !*traceslow then write ".";
    end;
  return gcd;
  end;

% Bugfix for mygcdf!*.
% Bug found with ezgcd returning different kord!* after
% a call.  This can change the ordering in an eqn.
% we must use !*funord ordering.  The new mygcdf!* (below)
% just takes care of the kord!* problem.  MJ.

%symbolic procedure mygcdf!*(u,v);
%begin scalar ans,funord;
%  if !*gcdreorder then <<funord:=!*funord; !*funord:=nil>>; %could be faster
%  if !*traceslow then write "{",count u,",",count v,"%";
%  ans:=if !*myezgcd then ezgcdf(u,v) where !*ezgcd=t else gcdf!*(u,v);
%  if !*traceslow then write "(",count ans,")}";
%  if null ans then rederr("gcd failed in mygcdf!*");
%  if !*gcdreorder then !*funord:=funord; 
%   % we dont need to reorder anything do we?
%   % not if args are parsed properly :-)
%  return ans;
%  end;


symbolic procedure mygcdf!*(u,v);
 begin scalar kord,ans;
  kord:=kord!*;
  ans:= xmygcdf!*(u,v);
  setkorder kord;
  return ans;
 end;

symbolic procedure xmygcdf!*(u,v);
begin scalar ans,funord;
  if !*gcdreorder then <<funord:=!*funord; !*funord:=nil>>; %could be faster
  if !*traceslow then write "{",count u,",",count v,"%";
  ans:=if !*myezgcd then ezgcdf(u,v) where !*ezgcd=t else gcdf!*(u,v);
  if !*traceslow then write "(",count ans,")}";
  if null ans then rederr("gcd failed in mygcdf!*");
  if !*gcdreorder then !*funord:=funord; 
   % we dont need to reorder anything do we?
   % not if args are parsed properly :-)
  return ans;
  end;
% end of bugfix for mygcdf!*.


symbolic procedure myquotf1E(eqn,u);
% ie eqn/u
% distribute the division through each lc.
if null eqn then eqn
else if null u then rederr "null divisor in myquotf1E"
else if not dunknk mvar eqn then rederr "bad eqn in myquotf1E"
else ((lpow eqn) .* quotf1(lc eqn,u))
     .+ myquotf1E(red eqn,u);

symbolic procedure dependlevels;
begin scalar deplevele,deplevels,thisdepcount,dunkns,j;
  dunkns:=dunknsineqns();
  j:=0;
  while dunkns do begin
    thisdepcount:=0;
    for each dunkn in dunkns do 
      if count(depsk dunkn)=j 
        then <<dunkns:=delete(dunkn,dunkns);
               thisdepcount:=thisdepcount+1>>;
    deplevele:=list(j,thisdepcount).deplevele;
    j:=j+1;
    if j>50
      then <<terpri();
             write" are there really dunkns with > 50 dependencies???";
             terpri()>>;
    end;
  if !*splitonpro
    then dunkns:=dunknsinprosymvec()
    else dunkns:=dunknsinsymvec();
  j:=0;
  while dunkns do begin
    thisdepcount:=0;
    for each dunkn in dunkns do 
      if count(depsk dunkn)=j 
        then <<dunkns:=delete(dunkn,dunkns);
               thisdepcount:=thisdepcount+1>>;
    deplevels:=list(j,thisdepcount).deplevels;
    j:=j+1;
    if j>50
      then <<terpri();
             write" are there really dunkns with > 50 dependencies???";
             terpri()>>;
    end;
  return list('symvec.deplevels,'eqns.deplevele);
  end;

symbolic procedure showdeplevs;
if !*tracedeplev then <<
  terpri();terpri();
  write "::",dependlevels(),":: ";terpri();
  terpri()>>;

%****** the datatype eqngrp is born here!

% an eqngrp is an eqn and a collection of information about that eqn.
% it is a list of:
%
%1  eqnno               The unique number of this eqn,
%                       which is a count of eqns.
%    (eg!*no)
%2  eqn                 The actual standardform equation
%                       which is implicitly zero
%    (eg!*eqn)          It is a differential equation for the determinable
%                       unknown functions (dunkns) it contains and it is
%                       structured so that each toplevel mvar is a dunkn or
%                       a derivative of a dunkn (dfdunkn).
%3  opstried            Operations which have been tried on this eqn
%    (eg!*ops)
%4  cause-of-death      Method this eqn was solved by, 
%                       or reason for dropping it
%    (eg!*dc)
%5  highdfdunkn         The highest ordered dunkn or dfdunkn in the eqn,
%    (eg!*highdfdunkn)  according to the special ordering in dfdunknorder
%6  subval              The value of highdfdunkn if the eqn is solved for it
%    (eg!*subval)
%7  simpflag            True if the eqn needs to be simplified
%    (eg!*simpstatus)
%8  length              The number of terms? in the eqn
%    (eg!*len)
%9  srceqns             The eqn(s) that this eqn was derived from
%    (eg!*bc)
%10  srcop              The operation that this eqn was derived from
%    (eg!*scrop)
%11  order              The order of the highest derivative in eqn
%    (eg!*ord)
%12  op!*sub2sftries         The eqn #s of any equations that have been tried
%    (eg!*op!*sub2sftries)   to be substituted into this equation
%                            (by op!*sub2sf or ?)
%13  intcon2tries
%    (eg!*i2tries)
%14  hiphiord           The highest order phi prolongation in the equation
%                       (normaly 0)
%    (eg!*hiphiord)
%15  slvterm            dunkn.coeff if eqn can be solved for dunkn
%    (eg!*slvterm)
%16  numdeps            number of dependencies of dunkns in eqn
%    (eg!*deps)
%17  nexthigh           next highest dfdunknf (after highdfdunkn)
%    (eg!*nexthigh)

symbolic procedure mkeqngrp(eqn,srceqns,srcop);
%
% takes a standardform eqn and returns the initial eqngroup
begin scalar eqngrp;
  if null eqn 
    then <<write tmstmp(),"adding zero equation???";terpri();return nil>>;
  !*dets:=!*dets+1;
  if !*tracecute then write tmstmp(),"B"; 
  if !*stopatnum and !*dets>!*stopatnum then solveexit:=t;
  if !*dothisatnum and !*dets=!*dothisatnum then mysimpq !*k2q '(dothis);
  if !*traceatnum and !*dets>!*traceatnum then trace();
  eqn:=dropeqnfacE eqn;
%  if null red eqn then <<dividelog lc eqn; eqn:=!*k2f mvar eqn>>;
  if (lnc eqn)<0 then eqn:=negf eqn;
  eqngrp:=list(!*dets,          eqn, list(),       list(),            nil,
                nil,            nil, mytermsf eqn, srceqns,           srcop,
                list(),         nil, nil,          highphiorderE eqn, list(),
                count boundvarsE eqn, nil );
  if !*traceeqngrp then showeqngrp eqngrp;
  specialfnsE eqn;
  if !*tracedunknorder then logdunknsE eqn;
  if !*tracelength then write(tmstmp(),"(",eg!*len eqngrp,")");
  if !*keepeqngrps then !*eghistory:=eqngrp.!*eghistory;
%  if !*tracecanon then <<write "!*recnumdeps=",!*recnumdeps;terpri()>>;
  if !*verifylin then op!*linchk eqngrp;
  return eqngrp;
  end;

symbolic procedure logdunknsE eqn;
% keep a list of dfdunkns in eqns.
if eqn then <<
  if not ((mvar eqn) member !*dunknsfound) 
    then !*dunknsfound:=(mvar eqn).!*dunknsfound;
  logdunknsE red eqn>>;

%******** selectors for standard equation groups ******

% Much of the information in an eqngrp might never be needed (if it is solved
% by an early op), so some information is only calculated and stored when 
% first called for. This should save time...

symbolic procedure eg!*no eqngrp;
% the equation number of the equation group
car eqngrp;
 
symbolic procedure eg!*eqn eqngrp;
% the actual equation in the equation group
cadr eqngrp;

symbolic procedure eg!*ops eqngrp;
% the operations that have been tried to solve this eqn
caddr eqngrp;

symbolic procedure eg!*dc eqngrp;
% Death certificate- the cause of death of the equation.
% Ie, the method it was solved by, if it has been.
cadddr eqngrp;

symbolic procedure eg!*highdfdunkn eqngrp;
% the highest ordered dfdunkn in the eqn
(if car x then car x else car rplaca(x,highdfdunknf eg!*eqn eqngrp))
  where x=cddddr eqngrp;

symbolic procedure eg!*subval eqngrp;
% the value of highdfdunkn if eqn is solved for it
(if car x 
  then car x 
  else car rplaca(x, subvalf(eg!*highdfdunkn eqngrp,eg!*eqn eqngrp)))
where x=cdr cddddr eqngrp;

symbolic procedure eg!*simpstatus eqngrp;
% does the equation need to be simped?
caddr cddddr eqngrp;

symbolic procedure eg!*len eqngrp;
% number of terms in the equation 
cadddr cddddr eqngrp;
 
symbolic procedure eg!*bc eqngrp;
% birth certificate- where did the eqn come from?
car cddddr cddddr eqngrp;

symbolic procedure eg!*srcop eqngrp;
% source operation: what operation(s) generated this eqn
cadr cddddr cddddr eqngrp;

symbolic procedure eg!*ord eqngrp;
% order of the highest derivative in eqn
(if car x 
  then car x 
  else car rplaca(x,dforderE eg!*eqn eqngrp))
where x=cddr cddddr cddddr eqngrp;

symbolic procedure eg!*op!*sub2sftries eqngrp;
% eqn #s of eqns that have been tried for subing into this eqn
cadddr cddddr cddddr eqngrp;

symbolic procedure eg!*i2tries eqngrp;
% eqn #s of eqns that have been tried to make incon2 eqns with
car cddddr cddddr cddddr eqngrp;

symbolic procedure eg!*hiphiord eqngrp;
% the highest order phi prolongation in the equation (normaly 0)
cadr cddddr cddddr cddddr eqngrp;

symbolic procedure eg!*slvterm eqngrp;
% dunkn.coeff if eqn can be solved for dunkn
(if car x 
  then car x 
  else car rplaca(x,sterm eqngrp))
where x=cddr cddddr cddddr cddddr eqngrp;

symbolic procedure eg!*deps eqngrp;
% the number of dependencies of dunkns in eqn
cadddr cddddr cddddr cddddr eqngrp;

symbolic procedure eg!*nexthigh eqngrp;
% the next highest dfdunknf (after highdfdunkn)
(if car x 
  then caar x 
  else (caar rplaca(x,if y then list highdfdunknf y else list nil)
       ) where y=numr eg!*subval eqngrp
)
where x=cddddr cddddr cddddr cddddr eqngrp;

% ****** functions that can change an eqngrp

symbolic procedure addopstried(op,eqngrp);
% we have tried the operation op to solve eqngrp
rplaca(cddr eqngrp, op.eg!*ops eqngrp);
  
symbolic procedure deathcert(eqngrp,causeofdeath);
begin
 if !*tracecute then write tmstmp(),"D";
 if !*tracecute then write "[",count alleqngrps(),"]";
 if !*tracedeath
   then <<terpri();
          write tmstmp(),("Equation ",eg!*no eqngrp,": ",causeofdeath);
          terpri() >>;
 rplaca(cdddr eqngrp,list causeofdeath);
 end;

symbolic procedure eg!*tagtosimp eqngrp;
  rplaca(cddr cddddr eqngrp,t);

symbolic procedure eg!*clearops eqngrp;
rplaca(cddr eqngrp, nil);

symbolic procedure clearops;
for each eqngrp in alleqngrps() do eg!*clearops eqngrp;

symbolic operator clearops;

symbolic procedure addop!*sub2sftries(trynums,eqngrp);
% op!*sub2sf has tried to sub eg#s trynums into eqngrp
rplaca(cdddr cddddr cddddr eqngrp,
       union(trynums,eg!*op!*sub2sftries eqngrp));

symbolic procedure addi2tries(trynum,eqngrp);
rplaca(cddddr cddddr cddddr eqngrp, trynum.eg!*i2tries eqngrp);
  
symbolic procedure simpdets;
for each eqngrp in alleqngrps() do eg!*tagtosimp eqngrp;
 
symbolic operator simpdets;

%*****end of eqngrp support

symbolic procedure addeqngrp eqngrp;
% add the eqngrp to !*deteqns, provided the eqn is different to the others.
begin scalar eqn,srceqns,repeated,checkeqngrps,dc;
  if (dc:=eg!*dc eqngrp) 
    then <<write tmstmp(),"**** Hey, looks like this fellow is dead!!!",
                          " His story is that";terpri();
           write dc;terpri();
           showeqngrp eqngrp;
           rederr "I don't like this, I'm confused!">>;
  eqn:=eg!*eqn eqngrp;
  srceqns:=eg!*bc eqngrp;
  if not pairp srceqns then srceqns:=list srceqns;
  checkeqngrps:=!*deteqns;
  while checkeqngrps and not repeated do
    if (eqn=eg!*eqn car checkeqngrps) 
      and not (eg!*no car checkeqngrps member srceqns)
        % above line prevents dropping an eqn because it mirrors its parent,
        % which can happen if an eqn simplifies to itself,
        % and the new eqn is added before the parent is dropped...
        then <<repeated:=t;
               deathcert(eqngrp,list("Dropped like a hotcake by addeqngrp ",
                                     "because it is identical to equation ",
                                      eg!*no car checkeqngrps))>>
        else checkeqngrps:=cdr checkeqngrps;
  if not repeated then !*deteqns:=inseqngrp(eqngrp,!*deteqns);
  end;

symbolic procedure inseqngrp(eqngrp,eqngrplis);
% sorted-insert eqngrp into the sorted eqngrplis
%if null eqngrplis or ordeqngrpbydepsdfordlen(eqngrp,car eqngrplis)
%if null eqngrplis or ordeqngrpbydfordlen(eqngrp,car eqngrplis)
if null eqngrplis or ordeqngrpbylendford(eqngrp,car eqngrplis)
  then eqngrp.eqngrplis
  else (car eqngrplis).inseqngrp(eqngrp,cdr eqngrplis);

symbolic procedure ordeqngrpbydfordlen(eqngrp1,eqngrp2);
% is eqngrp1 'easier' than eqngrp2?
% ie lower dford or equal dford and lower length
if eg!*ord(eqngrp1) neq eg!*ord(eqngrp2)
  then eg!*ord(eqngrp1)<eg!*ord(eqngrp2) 
else if eg!*len(eqngrp1) neq eg!*len(eqngrp2)
       then eg!*len(eqngrp1)<eg!*len(eqngrp2)
else if eg!*deps(eqngrp1) neq eg!*deps(eqngrp2)
       then eg!*deps(eqngrp1)<eg!*deps(eqngrp2)
else ordop(mvar eg!*eqn eqngrp1,mvar eg!*eqn eqngrp2);

symbolic procedure ordeqngrpbylendford(eqngrp1,eqngrp2);
% is eqngrp1 'easier' than eqngrp2?
% ie lower dford or equal dford and lower length
(eg!*len(eqngrp1)+2*eg!*ord(eqngrp1)) < (eg!*len(eqngrp2)+2*eg!*ord(eqngrp2))
 or ( (eg!*len(eqngrp1)+2*eg!*ord(eqngrp1)) 
        = (eg!*len(eqngrp2)+2*eg!*ord(eqngrp2)) 
      and eg!*ord(eqngrp1)<eg!*ord(eqngrp2) );

symbolic procedure ordeqngrpbydepsdfordlen(eqngrp1,eqngrp2);
% is eqngrp1 'more important' or 'easier' than eqngrp2?
% ie more deps, or same deps and lower dford or equal dford and lower length
eg!*deps(eqngrp1) > eg!*deps(eqngrp2) 
or ( eg!*deps(eqngrp1)=eg!*deps(eqngrp2) 
     and ( eg!*ord(eqngrp1) < eg!*ord(eqngrp2) 
           or ( eg!*ord(eqngrp1)=eg!*ord(eqngrp2)
                and eg!*len(eqngrp1)<eg!*len(eqngrp2) )
         )
   );

symbolic procedure alleqngrps();
union(!*deteqns,union(!*hideneqns,!*longhideneqns));

%***** funorder routines:

% permit kernel ordering by function _name_ , not function value

% KORDER and associated sorting routines permit sorting by set variables or
% kernels (set funtions with set arguments ) only.
%
% We want to include functions with unspecified args, the ordering for equal
% functions then being determined by the ordering of the arguments.
%
% This means having to change the reduce ORDOP procedure.
%
% this is the original version:
%
%symbolic procedure ordop(u,v);
% Is u ordered ahead of v or does u=v?
%   begin scalar x;
%        x := kord!*;
%    a:  if null x then return ordp(u,v)
%         else if u eq car x then return t
%         else if v eq car x then return;
%        x := cdr x;
%        go to a
%   end;

% this is to let us redefine this function.
remflag('(ordop),'lose);

% this is the new version to do extra sorting:
symbolic procedure ordop(u,v);
   begin scalar x;
        x := kord!*;
    a:  if null x then if !*funord then return ordop1(u,v)
                                   else return ordp(u,v)
         else if u eq car x then return t
         else if v eq car x then return;
        x := cdr x;
        go to a
   end;
% Added 4/5/99 at david hartley's suggestion.
% A patches.red changed ordpp to include a copy
% of ordop and then further extended it.  Since dimsym
% uses an edited version of ordop the copy inside 
% the new ordpp needed to be altered.  MJ.
% This is to let us redefine this function.
remflag('(ordpp),'lose);
symbolic procedure ordpp(u,v);
   begin scalar x;
        if car u eq car v then return cdr u>cdr v;
        x := kord!*;
        u := car u;
        v := car v;
%    a:  if null x then return ordpa(u,v)
    a:  if null x then if !*funord then return ordop1(u,v)
                                   else return ordpa(u,v)
         else if u eq car x then return t
         else if v eq car x then return nil;
        x := cdr x;
        go to a
   end;
% Is u ordered before v?
% this compares the functions w/o args, if they are on the list
% if they are and they are equal then compare their argments (at great cost?)
%
% !*funord is a list of functions from highest ordering to lowest.
%
symbolic procedure ordop1(u,v);
if atom u then if atom v then ordp(u,v) else nil
else if atom v then t         
%else if null !*funord then ordp(u,v)
  else begin scalar x;
        x := !*funord;
        if (car u) eq (car v) then return ordop2(cdr u,cdr v);
    a:  if null x then return ordp(u,v)
        else if (car u) eq car x then return t
        else if (car v) eq car x then return;
        x := cdr x;
        go to a
   end;

symbolic procedure ordop2(u,v);
% only used for comparing inside kernels
if null u then null v 
else if null v then t
else if atom u then if atom v then if numberp u then numberp v and  u geq v
                                   else if numberp v then t 
                                   else orderp(u,v)
                    else nil
else if atom v then t
else if (car u)=car v then ordop2(cdr u,cdr v)
else begin scalar x;
        x := !*funord;
    a:  if null x then return ordop2(car u,car v)
          else if (car u) eq car x then return t
          else if (car v) eq car x then return;
        x := cdr x;
        go to a
   end;

% We must have dunkns ordered before df, and df ordered before
% anything else in equations. Do we need partdf ordered before dunkns?
% As it stands, we cant just add new dunkns into !*funord, so we wont
% let the user add dunkns
% To modify !*funord, any new dunkns must go between partdf and df.

!*funord:='(partdf c phi xi df);

%*****end of funorder

% kludge rules for integrating a derivative...

% 1 var
%for all x,y 
%  let int(df(x,y),y)=x;

% % 1 var with num
% for all x,y 
%   let int(df(x,y,2),y)=df(x,y);
% for all x,y,z such that fixp z and z > 2 
%   let int(df(x,y,z),y)=df(x,y,z-1);

% % 2 vars
% for all x,y,z 
%   let int(df(x,y,z),z)=df(x,y);
% for all x,y,z such that not fixp z
%   let int(df(x,y,z),y)=df(x,z);

% % 2 vars with nums
% for all x,y,z
%   let int(df(x,y,z,2),z)=df(x,y,z),
%       int(df(x,y,2,z),y)=df(x,y,z);
% for all w,x,y,z such that fixp z and z>2
%   let int(df(w,x,y,z),y)=df(w,x,y,z-1),
%       int(df(w,x,z,y),x)=df(w,x,z-1,y);

% % 3 vars
% for all w,x,y,z 
%   let int(df(w,x,y,z),z)=df(w,x,y);
% for all w,x,y,z such that not fixp z
%   let int(df(w,x,y,z),y)=df(w,x,z),
%       int(df(w,x,z,y),x)=df(w,z,y);

% % 3 vars with nums
% for all w,x,y,z
%   let int(df(w,x,y,z,2),z)=df(w,x,y,z),
%       int(df(w,x,y,2,z),y)=df(w,x,y,z),
%       int(df(w,x,2,y,z),x)=df(w,x,y,z);
% for all v,w,x,y,z such that fixp z and z>2
%   let int(df(v,w,x,y,z),y)=df(v,w,x,y,z-1),
%       int(df(v,w,x,z,y),x)=df(v,w,x,z-1,y),
%       int(df(v,w,z,x,y),w)=df(v,w,z-1,x,y);

% % 4 vars
% for all v,w,x,y,z
%   let int(df(v,w,x,y,z),z)=df(v,w,x,y);
% for all v,w,x,y,z such that not fixp z 
%   let int(df(v,w,x,y,z),y)=df(v,w,x,z),
%       int(df(v,w,x,z,y),x)=df(v,w,z,y),
%       int(df(v,w,z,x,y),w)=df(v,z,x,y);

% % 4 vars with nums
% for all v,w,x,y,z
%   let int(df(v,w,x,y,z,2),z)=df(v,w,x,y,z),
%       int(df(v,w,x,y,2,z),y)=df(v,w,x,y,z),
%       int(df(v,w,x,2,y,z),x)=df(v,w,x,y,z),
%       int(df(v,w,2,x,y,z),w)=df(v,w,x,y,z);
% for all u,v,w,x,y,z such that numberp u and u>2
%   let int(df(v,w,x,y,z,u),z)=df(v,w,x,y,z,u-1),
%       int(df(v,w,x,y,u,z),y)=df(v,w,x,y,u-1,z),
%       int(df(v,w,x,u,y,z),x)=df(v,w,x,u-1,y,z),
%       int(df(v,w,u,x,y,z),w)=df(v,w,u-1,x,y,z);

% % 5 vars
% for all u,v,w,x,y,z
%   let int(df(u,v,w,x,y,z),z)=df(u,v,w,x,y);
% for all u,v,w,x,y,z such that not fixp z 
%   let int(df(u,v,w,x,y,z),y)=df(u,v,w,x,z),
%       int(df(u,v,w,x,z,y),x)=df(u,v,w,z,y),
%       int(df(u,v,w,z,x,y),w)=df(u,v,z,x,y),
%       int(df(u,v,z,w,x,y),v)=df(u,z,w,x,y);

% % 5 vars with nums
% for all u,v,w,x,y,z
%   let int(df(u,v,w,x,y,z,2),z)=df(u,v,w,x,y,z),
%       int(df(u,v,w,x,y,2,z),y)=df(u,v,w,x,y,z),
%       int(df(u,v,w,x,2,y,z),x)=df(u,v,w,x,y,z),
%       int(df(u,v,w,2,x,y,z),w)=df(u,v,w,x,y,z),
%       int(df(u,v,2,w,x,y,z),v)=df(u,v,w,x,y,z);
% for all s,u,v,w,x,y,z such that numberp s and s>2
%   let int(df(u,v,w,x,y,z,s),z)=df(u,v,w,x,y,z,s-1),
%       int(df(u,v,w,x,y,s,z),y)=df(u,v,w,x,y,s-1,z),
%       int(df(u,v,w,x,s,y,z),x)=df(u,v,w,x,s-1,y,z),
%       int(df(u,v,w,s,x,y,z),w)=df(u,v,w,s-1,x,y,z),
%       int(df(u,v,s,w,x,y,z),v)=df(u,v,s-1,w,x,y,z);

%***** end of kludge

symbolic procedure checktype type;
if not (if pairp type then member(car type,!*mkdetsargtypes)
                      else member(type,!*mkdetstypes))
  then rederr list(type," not supported by mkdets, try one of ",
                   !*mkdetsargtypes," with argument or one of ",
                   !*mkdetstypes," without argument.");


symbolic procedure mkdets type;
% read the diff eqns, do all initialisation, calculate determining equations.
% type is the type of symmetries to look for, eg point, handled by setxiphi.
begin scalar eqn;
  checktype type;
  !*deteqns:=nil;   % drop any equations hanging around from a previus run
  !*hideneqns:=nil; % shouldn't be any of these, but you never know...
  !*opusage:=nil;   % we haven't used any ops to solve yet...
  !*intfac!*opusage:=!*op!*xdp!*opusage:=nil; % ditto
  !*splitvars:=!*proeqnsmade:=nil;
%  if not (type='custom2 or type='custom1) then !*p:=!*q:=!*r:=0;
  !*p:=!*q:=!*r:=0;
  if !*predeqlis then readequations(); 
               % might already have been called if we checked a vec first.
               % also sets !*deqlis, !*algsublis and !*sqsublis
  setxiphi(type);
  mkprosymvec(type);
  if type='liebacklund 
      or (pairp type and car type='liebacklund and numberp cadr type)
    then lbrmdeps dunknsinsymvec();
  if !*mkproatstart then mkprolongations();
    % calculates the needed prolongations of phi
  if !*subeqnsinprovec 
    then prosymvec:=mk!*sq simp!* liesubq simp!* prosymvec;
     % if we substitute the equations in here, maybe faster?
  for each deqn in !*deqlis do begin
    eqn := mynumr mysimpq vecderf(simp!* prosymvec,deqn);
    eqn := mynumr mysimpq liesubq (eqn ./ 1);
      % substitute in the original diff eqns
    if null eqn then rederr("Base determining equation is zero???");
    addeqngrp mkeqngrp(eqn,'deq,'lie);
    if type='point and !*autosplitdets and !*mkproatstart

      then solvedets 'stdsplit;
        % splits !*deteqns by high u's and solves single term eqns.
        % Do this now because it could make big simplifications for following
        % base determining equations if we can remove dependencies...
    end;
  end;

symbolic operator mkdets;

symbolic procedure lbrmdeps dunkns;
% Remove the dependence of dunkns on any vars
% that would be subbed out by liesub,
% ie any (totders of) vars on !*usubvars
for each dunkn in dunkns do
  for each depu in depsk dunkn do
    for each subu in !*usubvars do
      if subu=depu or ((cadr subu)=(cadr depu)
                       and (minusudfs(cddr depu,cddr subu) neq -1))
        then begin
          write "We have an equation for ",subu,
                " so any dependence on ",depu,
                "would be substituted out. ",
                "So I am removing the dependence of ",dunkn,
                " on ",depu,".";terpri();
          depend1(dunkn,depu,nil);
          end;

symbolic procedure mkprosymvec type;
begin scalar ndx,cndx,i,locsymvec,locprosymvec;
%   if not (type='custom2 or type='custom1) 
%     then begin % clear all references to xi's and phi's
%       put(!*xi,'kvalue,nil);
%       put(!*phi,'kvalue,nil);
%       end;
  if not (type='custom2) 
    then begin % create symvec and prosymvec  
      locsymvec:=nil;
      for i:=1:!*p do
        locsymvec:=addf(locsymvec,
                        multf(!*k2f !*a2k list('partdf,list(!*x,i)),
                              !*k2f  list(!*xi,i)
                              )
                        );
      for i:=1:!*q do
        locsymvec:=addf(locsymvec,
                        multf(!*k2f !*a2k list('partdf,list(!*u,i)),
                              !*k2f  list(!*phi,i)
                              )
                           );
      locprosymvec:=locsymvec;
      symvec:=mk!*sq(locsymvec ./ 1);
      ndx:='(1); 
      while (count ndx) < (!*r+1) do begin
        cndx:=sortndx ndx; 
        for i:=1:!*q do
          locprosymvec:=addf(locprosymvec,
                             multf(!*k2f !*a2k list('partdf,!*u.i.cndx),
                                        !*k2f  (!*phi.i.cndx)
                                        )
                                  );
        ndx:=incndx ndx;
        end;
      prosymvec:=mk!*sq(locprosymvec ./ 1);
      end;
  end;

symbolic procedure loaddeq eqn;
% load the differential equation onto !*predeqlis.
% this is a kludge to avoid bug? in 3.4
<<!*predeqlis:=eqn.!*predeqlis;count !*predeqlis>>;

symbolic operator loaddeq;

symbolic procedure readequations;
%
% read the differential equations stored as values of the operator deq(j)
% these must have the form deq(j):=u(?,...,?)=whatever
% with no u derivative appearing as the left hand side of one of these
% equations allowed to be used in any of the other equations.
%
% another role of this procedure is to set up
% the alists !*algsublis and !*sqsublis
% which are lists of substitution rules to be used by a sublis stsatement
% to reinsert the differential equations into the base determining equations.
% !*algsublis has algebraic values while !*sqsublis has the sq equivalents.
%
begin scalar neweqn,val,kval,thiseqn,uvar,madesub,uu,subpr;integer j;
  depl!*:=cleanndx depl!*;
  !*deqlis:=!*usubvars:=!*algsublis:=!*sqsublis:=!*vars:=nil;
  if null !*predeqlis then 
     !*predeqlis:=for each kval in get('deq,'kvalue) collect cadr kval;
%  for each kval in get('deq,'kvalue) do begin
%    thiseqn:=cadr kval;
  !*predeqlis := cleanndx !*predeqlis;
     % make sure all indicies in u's are ordered
  for each thiseqn in !*predeqlis do begin         % set up the substitutions
    if (car thiseqn) neq 'equal 
      then rederr 
               list("equations not set up as substitution for u derivative");
    uvar:=simp!* cadr thiseqn;
    if not((denr uvar=1) and (null red numr uvar) 
           and (lc numr uvar=1) and (ldeg numr uvar=1)
           and pairp(uvar:=mvar numr uvar) and (car uvar)=!*u)
      then rederr 
             list("equations not set up as substitution for u derivative");
      % *** also have to check that eqns don't resub eachother...
    val:=simp!* caddr thiseqn;
    !*usubvars:=uvar.!*usubvars;
    !*algsublis:=(uvar. prepsq val) . !*algsublis;
    !*sqsublis:=(uvar. mk!*sq val) . !*sqsublis;
      % for substituting back in to determining eqns
  end;
  !*firstusubvars:=!*usubvars:=reverse sort(!*usubvars,'ordop);
  checksemistandardform(!*usubvars,!*sqsublis);
                          % check that the equations are in standard form
  !*predeqlis:=nil; % so that any new equations loaded start fresh
  !*liesubcount:=0; % count the number of times we resub eqns.
  madesub:=t; % just to push into loop first time.
  while madesub do begin
    madesub:=nil;
    for each subpr in !*sqsublis do if not madesub then begin 
                                           % make substitutions into this eqn
      val:=liesubq1 simp!* cdr subpr;
      if madesub then << % replace eqnsubs for new eqn
        write(tmstmp(),"Made substitution into given equation for ",
              car subpr);terpri();
%        !*sqsublis:=rplacd(subpr,mk!*sq val);
%        !*algsublis:=rplacd(assoc(car subpr,!*algsublis),prepsq val)>>;
        !*sqsublis :=((car subpr). mk!*sq val).delete(subpr,!*sqsublis );
        !*algsublis:=((car subpr). prepsq val).delete(subpr,!*algsublis)>>;
      end;
    end;
  !*sqsublis:=reverse sort(!*sqsublis,'compsublisprs); 
    % so that !*deqlis is sorted by (reverse of) the same
    % so that we do vecder on easiest eqns first- *MUCH* faster
    % (eg x5 for del4x)
  for each uu in !*firstusubvars do begin 
                      % just use original eqns, not any prolongations made
    subpr:=assoc(uu,!*sqsublis);
    neweqn:=mynumr mysimpq addsq(negsq !*k2q uu,simp!* cdr subpr);
    !*deqlis:=neweqn . !*deqlis; % list of equations to find syms for
    mysetk(list('ideq,j:=j+1),mk!*sq (neweqn ./ 1) );
    end;
  for each eqn in !*deqlis do !*vars:=union(!*vars,alldepsf eqn);
  findpqr !*vars;                % sets !*p,!*q and !*r
end;

symbolic procedure compsublisprs(u,v);
% u,v are prs on !*sqsublis, ie uval . !*sq val
% is u 'simpler' based on diff'n order then length?
if (count car u)<(count car v) then t % compare df order
else if (count car u)=(count car v) then (mywsqlen cdr u)<(mywsqlen cdr v) 
else nil;                                        % compare length of value


symbolic procedure mywsqlen u;
% a weighted length function for !*sq's
% based on length of eqn after subbing this expr in
% ( approx = [length of orig eqn + length of numr * #of subs]*length of denr )
begin scalar n,d;
  u:=simp!* u;
  d:=tmsf denr u;
  if d=0 then d:=1;
  n:=tmsf numr u;
  if n=0 then n:=1;
  return n*d*d;
  end;
  
symbolic procedure cleanndx v;
%
% looks for an occurence of operator u at any level of v,
% if it finds !*u then it sorts its derivative list,
% which will be the cddr of the list with car !*u.
% returns v with the appropriate sorting done on u's
%
% this procedure must be consistent with totder
%
% *** it would be good to match this with u's used here...
%
if null v or atom v then v
else if (car v)=!*u then !*u.(cadr v).sortndx cddr v
else (cleanndx car v).(cleanndx cdr v);

symbolic procedure cleanndx1 v;
% v is an indexed 'u variable
% this is to match it with other u's used.
!*q2k mysimpq !*k2q v;

symbolic procedure sortndx v;
% takes a list of integers and sorts them from highest to lowest
% returns the sorted list
if null v or null cdr v then v
else begin scalar j,newndx;
  for each j in v do newndx:=addndx(j,newndx);
  return newndx;
  end;

symbolic procedure findpqr v;
%
% find the highest order derivative of !*u 
% which occurs in the list v of x's and u's.
% Also finds the number of independent and dependent variables used and sets them as !*p,!*q,!*r
%
% v is a list of x,u variables
%
begin scalar var;
  !*p:=!*q:=!*r:=0;
  for each var in v do 
    if car var=!*x then !*p:=max2(!*p,cadr var)
    else if car var=!*u then begin
      !*q:=max2(!*q,cadr var);
      !*r:=max2(!*r,count cddr var);
%      !*p:=eval expand(!*p.cddr var,'max2);% the greatest integer in the list p.cddr var
      !*p:=mymax(!*p . cddr var);
      end;
  end;

symbolic procedure mymax u;
% u is a list of integers
% return the greatest of them
if null u then rederr("Bad list of integers.")
else if null cdr u then car u
else max2(car u,mymax cdr u);

symbolic procedure checksemistandardform(usubvars,sqsublis);
%
% usubvars are the leading derivatives of the original equations.
% check that the equations are in semi-standard form,
% ie that no leading derivatives is a total derivative of another leading derivative,
% and that the leading derivative of an eqn is ordered (by what?) higher than other 
% derivatives in the equation.
%
% just check the first condition for now ... 
%
begin scalar highu,ueqn;
  for each u1 in usubvars do for each u2 in usubvars do
    if not (u1=u2) and (cadr u1)=(cadr u2) then
      if (minusudfs(cddr u1,cddr u2) neq -1) then begin
        terpri();
        write "*** Equation for ",u1," is not in semi-standard form, ",
              "because leading derivative ",u1," is a total derivative of ",u2,
              " which we have an equation for ";terpri();
        terpri();
        end;
  for each u in usubvars do begin
    ueqn  := simp!* cdr assoc(u,sqsublis);
    if (ueqn=0)
      then highu:=nil
      else highu:=highu union(alldepsf numr ueqn, alldepsf denr ueqn);
    if highu then
      if highu=u
        then <<rederr "Equation for ",u," is not in semi-standard form, ",
                      "because leading derivative ",highu,
                      " appears on RHS of same given eqn.";terpri();terpri()>>
      else if (count highu) > (count u)
        then <<write "Equation for ",u," is not in semi-standard form, ",
                     "because it should be solved for ",highu,
                     " as leading derivative, ",
                     "which is higher order than ",u,".";terpri();terpri()>>
      else if (count highu)=(count u) and uorder(highu,u)
        then <<write "Equation for ",u,
                     " is not strictly in semi-standard form, ",
                     "because it should be solved for ",highu,
                     " as leading derivative, ",
                     "which is nominally ordered higher than ",u,",",
                     " although this should be ok.";terpri();terpri()>>;
    end;
  end;

symbolic procedure highu varlist;
% what is the highest ordered u in varlist?
begin scalar hu;
  for each var in varlist do
    if (car var)=!*u and (null hu or (count var) > (count hu)
                          or ( (count var)=(count hu) and uorder(var,hu) ))
      then hu:=var;
  return hu;
  end;

symbolic procedure uorder(u,v); 
% u and v are !*u derivatives, of equal order.
% is u ordered before v?
% (we could have just used uderivorder for this fn)
(cadr u) > (cadr v)
or ((cadr u)=(cadr v) and uderivorder(cddr u,cddr v));

symbolic procedure uderivorder(u,v);
if null u or null v or atom u or atom v
  then rederr ("Bad uderivorder comparison of ",u," and ",v,".")
  % this would mean u and v were equal or of different diff'n order.
else if (car u)>(car v) then t
else if (car u)=(car v) then uderivorder(cdr u,cdr v);

symbolic procedure proeqnifnd varlist;
% varlist is a list of x and u variables appearing in an expression
% requiring the original equations to be substituted.
% Look for any vars u(j,...) which are used,
% and which are total derivatives of lhs u(i,...) of original equations,
% and prolong that equation.
begin scalar chksubs, lastsub; %1
  for each var in varlist do
    if (car var)=!*u and not (var member !*usubvars) then begin %2
      lastsub:=nil;                              % last
      chksubs:=reverse !*usubvars;               
      while chksubs do begin %3
        if (cadr var)=(cadr car chksubs) 
              %^ j in u(j,...)  ^ i in u(i,...))
            and (minusudfs(cddr var,cddar chksubs) neq -1)
                %  is u(j,...)  a df of u(i,...
          then if not lastsub
            then <<chk!#subs();
                   addproeqn(var,car chksubs);
                   lastsub:=car chksubs>>
            else begin %4
              if (minusudfs(cddr lastsub,cddar chksubs) neq -1)
                  or (minusudfs(cddar chksubs,cddr lastsub) neq -1)
                then %write "I dont think the equations ",
                     %      "are in standard form..."
                     % this case is probably ok, we are just diff'ing
                     % an equation (we) made
                else <<
                  write "*** I got equation for ",var,
                        " by differentiating equation for ",lastsub,",";
                  write "but I could also differentiate equation for ",
                        car chksubs," !!!";terpri();
                  write "I think this means integrability ",
                        "conditions exist,",
                        " which _should_ be included by _you_ ",
                        "if they are nontrivial!";terpri();
                      >>
              end; %4
        chksubs:=cdr chksubs;
        end; %3
      end; %2
  end; %1

symbolic procedure addproeqn(udfvar,usbvar);
% udfvar appears in an expression requiring the original equations to be substituted.
% And udfvar is a totder of usbvar, which appears as the lhs of a given equation (or pro done)
begin scalar val;
  write(tmstmp(),"Forming prolongation of equation for ",usbvar,
                 " to get equation for ",udfvar);terpri(); 
  val:=cadr cdr assoc(usbvar,!*sqsublis);
  for each j in minusudfs(cddr udfvar,cddr usbvar) do val:=totderq(j,val);
  val:=mysimpq val;
  val:=liesubq1 val; % dont use liesub because it resets !*liesubcount
% resub the equations in here if needed!
  !*usubvars:=sort(udfvar.!*usubvars,'rordop);
  !*algsublis:=(udfvar. prepsq val).!*algsublis;
  !*sqsublis:=sort((udfvar. mk!*sq val).!*sqsublis,'rcompsublisprs);
  end;

symbolic procedure rordop(u,v);
ordop(v,u);

symbolic procedure rcompsublisprs(u,v);
compsublisprs(v,u);

symbolic procedure liesub u;
mk!*sq mysimpq liesubq simp!* u;

symbolic operator liesub;

symbolic procedure liesubq v; 
begin scalar madesub;
  if !*predeqlis then readequations(); % needed if mkdets has not been called
  !*liesubcount:=0;
  if !*traceliesub then <<
    write(tmstmp(),"doing liesubq on "); terpri(); showsq v;
    write("!*algsublis is: ",!*algsublis); terpri();
    write("!*sqsublis is: ",!*sqsublis); terpri() >>;
  return liesubq1 v;
  end;

symbolic procedure liesubq1 v;
begin scalar vv,varsinv,madesub1;
  madesub1:=t; % to push into loop first time
  while madesub1 do begin;
    madesub1:=nil;
    varsinv:=sort(union(alldepsf numr v,alldepsf denr v),'ordop);
    proeqnifnd varsinv;
    if intersect(varsinv,!*usubvars) then begin 
                                       % actually have to do a substitution
      chk!#subs();
      if !*tracecute then write "?";
      if !*traceliesub
         then <<write("doing liesubq loop on "); terpri(); showsq v>>;
      if !*sqliesub
        then vv:=mysimpq quotsq(sqliesubf numr v, sqliesubf denr v)
        else vv:=simp!* liesuba prepsq v;
      if !*traceliesub then <<
        write("liesubqloop on "); terpri(); showsq v;
        write("gave:"); terpri(); showsq vv>>;
      v:=vv;
      madesub1:=madesub;
      end;
    end;
  if !*tracecute then if madesub 
    then write ("+",!*liesubcount)
    else write "|";
  return v;
  end;

symbolic procedure sqliesubf v;
% v is a standard form
% return a sf
if null v or domainp v then (v ./ 1)
else addsq(multsq(exptsq(sqliesubk mvar v,ldeg v),
                  sqliesubf lc v),
           sqliesubf red v);

symbolic procedure sqliesubk v;
% v is kernel
% return sq 
%simp!* liesuba v;
(if pairp x and (car x)='!*sq then cadr x else !*k2q x) where x=liesuba v;

symbolic procedure liesuba v;
% v is an algebraic expresion
if null v or atom v then v
else if (car v)=!*u then tryliesuba v
%else if (car v)='df then (car v).(ftliesuba cadr v).(cddr v) 
     % dont sub df vars! 
     % dont sub dfs at all for now... what if unsimplified df(u 1,u 1) etc...
else if (car v)='df then v
else (liesuba car v).(liesuba cdr v);

symbolic procedure tryliesuba v;
% v is a !*u variable
% Return an algebraic expresion: (!*sq ~ ~) if !*sqliesub:=t,
% prefix form otherwise
if v member !*usubvars then begin scalar val;
  madesub:=t;
%  chk!#subs();
  val:=cdr assoc(v,if !*sqliesub then !*sqsublis else !*algsublis);
  if !*traceliesub
    then <<write("In tryliesuba, substituting ",v," as ");terpri();
           showalg val>>;
  return val;
  end
else v;
 
symbolic procedure chk!#subs;
<<!*liesubcount:=!*liesubcount+1;
  if !*liesubcount>!*liesublimit
     then rederr("Substitutions too deep when substituting given differential equation.")>>;

symbolic procedure minusudfs(y,z);
% y and z are the (ordered!) dfvarlists from u variables
% eg  y=(4,4,4,2,1,1) from yy=u(3,4,4,4,2,1,1)
% and z=(2,2,1) from zz=u(3,2,2,1)
% if yy is a totder of zz then return the difference y-z
% else return -1 
if null z then y
else if null y then -1
else if (car y)=(car z)
  then ((lambda x; if x=-1 then -1 else x) minusudfs(cdr y,cdr z))
else if (car y>car z) % from ordering in sortndx
  then ((lambda x; if x=-1 then -1 else (car y).x) minusudfs(cdr y,z))
else -1;

symbolic procedure setxiphi type;
%
% type is the type of symmetries to look for: supported types are:
%  point        Sets phi(i), xi(i) to depend on x(j)'s and u(j)'s.
%  liebacklund  Sets phi(i) to depend on x(j)'s, u(j)'s and u(j,k)'s, xi(i)=0
%  liebacklund(n)  Sets phi(i) to depend on x(j)'s, u(j)'s and
%                derivatives of u(j)'s up to order n, xi(i)=0
%  custom1      Uses user defined values or dependencies for phi(i), xi(i).
%               Warning- type 'custom1' does not reset any previously defined
%               values or dependencies for xi(i), phi(i) or higher order phis
%               which might have
%               been set by previous calls to mkdets or solvedets.
%  custom2      uses user defined values for symvec and prosymvec.
%               note the importance of seting prosymvec as well as symvec
%               eg, just do PROSYMVEC:=PROLONG(3,SYMVEC);
%               This call of prolong is less efficient than the usual one...
%
begin scalar x, basevars, lbvars, n, m,ordnvars, ordn!-1vars;
  checktype type;
  if type='custom2 then return;
  if not (type='custom2 or type='custom1) 
    then begin % clear all references to xi's and phi's
      put(!*xi,'kvalue,nil); 
      put(!*phi,'kvalue,nil);
      for j:=1:!*p do if (x:=assoc(list(!*xi , j),depl!*))
        then depl!* := delete(x, depl!*);
      for j:=1:!*q do if (x:=assoc(list(!*phi , j),depl!*))
        then depl!* := delete(x, depl!*);      
      end;
  basevars := nil;
  for j:=1:!*p do 
    <<basevars := list(!*x,j) . basevars; !*basefns:=list(!*xi,j) . !*basefns>>;
  for j:=1:!*q do 
    <<basevars := list(!*u,j) . basevars; !*basefns:=list(!*phi,j) . !*basefns>>;
  if type='point then begin
    for j:=1:!*p do depl!* := ( list(!*xi , j). basevars) . depl!*;
    for j:=1:!*q do depl!* := ( list(!*phi, j). basevars) . depl!*;
    end;
  if type='liebacklund then begin
    lbvars:=basevars;
    for j:=1:!*q do for k:=1:!*p do lbvars := list(!*u,j,k) . lbvars;
%    for j:=1:!*p do depl!* := ( list(!*xi , j). lbvars) . depl!*;
    for j:=1:!*p do mysetk(list(!*xi , j), 0);
    for j:=1:!*q do depl!* := ( list(!*phi, j). lbvars) . depl!*;
    end;                             % of type='liebacklund
  if pairp type and car type='liebacklund and numberp(m:=cadr type) then begin
    for j:=1:!*q do ordnvars:= list(!*u,j) . ordnvars;
    lbvars:=basevars;
    for n:=1:m do begin % add u vars of order n to lbvars
      ordn!-1vars:=ordnvars;
      ordnvars:=nil;
      for j:=1:!*p do for each var in ordn!-1vars do 
        ordnvars:=union(ordnvars,list((car var).(cadr var).addndx(j,cddr var)));
        % ie add df(var,x j) to ordnvars if its not already there.
      lbvars:=union(lbvars,ordnvars);
      end;
%    for j:=1:!*p do depl!* := ( list(!*xi , j). lbvars) . depl!*;
    for j:=1:!*p do mysetk(list(!*xi , j), 0);
    for j:=1:!*q do depl!* := ( list(!*phi, j). lbvars) . depl!*;
    end;                            % of type='liebacklund n
  end;

symbolic procedure mkprolongations;
% calculate the necessary prolongations
% we need one corresponding to each u in !*vars
% for q:=1:!*q do forallndx(!*p,!*r,function mkproifused,q); 
begin scalar var;
  for each var in !*vars do if car var=!*u and cddr var 
%    then mkpro(cadr var,cddr var);
    then mysetk(!*a2k 'phi.(cdr var),mk!*sq mkpro(cadr var,cddr var));
  end;

symbolic procedure mkpro(i,ndx);
% calculate the characteristic phi(j,ndx) if neccessary
%
% just use the recursive definition of prolongation, ie
% phi(i)^(k,J) = D_k phi(i)^(J) - (sum m) u(i)_(m,J) * D_k xi(m)
% where ndx=(k,J)
%
begin scalar k,J,m,kval,phival,xival,newval,xideriv;
  if null i then rederr "cant do mkpro on null phi!";
  if null ndx then return !*k2q list(!*phi,i);
  k:=car ndx;
  J:=cdr ndx;
  kval:=get(!*phi,'kvalue);
  if assoc(!*phi. i . ndx,kval) 
    then return simp!* cadr assoc(!*phi.i.ndx,kval);  % already made!
  phival:=assoc(!*phi. i . J, kval);                  % ie pair for phi(i)^(J)
  phival:=if phival then 
    simp!* cadr phival                               % ie phi(i)^(J)
    else if null J
      then !*k2q list(!*phi,i) 
%      else mkpro(i,J);
      else if !*forcelowerprovalues
        then <<if !*tracecute then write "forced setting value for ",'phi.i.J;
               phival:=mkpro(i,J);
               mysetk(!*a2k ('phi.i.J),mk!*sq phival);
               phival>>
        else !*k2q !*a2k ('phi.i.J);
  newval:=totderq(k,phival);
  for m:=1:!*p do begin                         % !*p is the number of x vars
%    xideriv:=simp!* list(!*xi,m,k);             % needs to be set somewhere, all at once?
    xival:=simp!* list(!*xi,m);
    xideriv:=totderq(k,xival);
    newval := addsq(newval,
                    multsq(negf(!*k2f !*a2k (!*u.i.addndx(m,J))) ./ 1,
                           xideriv
                           )
                    );
    end;
  if !*subeqnsinprocoefs then newval:=liesubq newval;
%  mysetk(!*phi . i . ndx, mk!*sq newval);
  return newval
  end;


%********end of mkdets

symbolic procedure solvedets alg;
% alg is the algoritm to use to solve the equations
begin
  if not(alg member !*solvedetsalgorithms)
    then <<write "*** ",alg," is not a recognised algorithm. Must be one of"; terpri();
           write !*solvedetsalgorithms; terpri()>>
    else <<write "Solving equations using ",alg," algorithm.";terpri();terpri();
           solvedetsbyops get(alg,'solveops)>>;
  end;

symbolic operator solvedets;

symbolic procedure addintcon1;
% add intcon1 integrability conditions to !*deteqns
begin integer dets,eqngrp;
  dets:=!*dets;
  for each eqngrp in intcon1() do addeqngrp eqngrp;
  if dets neq !*dets then solvesuccess:=t;
  end;

symbolic procedure addintcon2;
% add intcon2 integrability conditions to !*deteqns
begin integer dets,eqngrp;
  dets:=!*dets;
  for each eqngrp in intcon2() do addeqngrp eqngrp;
  if dets neq !*dets then solvesuccess:=t;
  end;

symbolic procedure addintcons;
% add integrability conditions to !*deteqns
begin integer dets,eqngrp;
  dets:=!*dets;
  for each eqngrp in intcon1() do addeqngrp eqngrp;
  for each eqngrp in intcon2() do addeqngrp eqngrp;
  if dets neq !*dets then solvesuccess:=t;
  end;

symbolic procedure addallintcons;
% add _all_ integrability conditions to !*deteqns
begin integer dets,eqngrp;
  while not (dets=!*dets) do begin 
    dets:=!*dets;
    for each eqngrp in intcon1() do addeqngrp eqngrp;
    op!*sub2sf();
    for each eqngrp in intcon2() do addeqngrp eqngrp;
    op!*sub2sf();
    if dets neq !*dets then solvesuccess:=t;
    if null !*allics then dets:=!*dets; 
       % so that we dont loop: ie only calculate one level of ICs
    end;
  end;
 
symbolic procedure solvedetsbyops ops;
% apply the ops to !*deteqns
% if op is flagged as onebyoneop 
%   then onebyoneop does op on each eqngrp in !*deteqns
%   else op works on the whole list
% any op can stop the solve loop by setting the flag solveexit
% Typically, this is set if !*dets, the number of deteqns, exceeds a limit
begin scalar op,opstotry,solveexit,solvesuccess,opstried;
  if !*tracecute then <<write("Solving equations ",egnums !*deteqns," by ",ops); terpri()>>;
  solveexit:=nil;
  opstotry:=ops;
  !*hidenbefore:=nil; % we havent hidden any equations yet.
  while (!*deteqns or !*hideneqns) and opstotry and not solveexit do begin
    solvesuccess:=nil;
    op:=car opstotry;
    opstotry:=cdr opstotry;
    if kord!* then
       <<write(tmstmp(),"*****!!!We have kord!*=",kord!*,
               " in solvedetsbyops, hope this is ok...");terpri()>>;
    if !*tracecute and not(op member opstried) then <<write op," "; opstried:=op.opstried>>;
    if !*traceops then <<write("trying ",op); terpri() >>;
    if flagp(op,'onebyoneop) 
      then onebyoneop op       % does op on each eqngrp in !*deteqns
      else apply(op,nil);      % does op on the whole of !*deteqns
    if solvesuccess then recordop(op,'!*opusage);
    if solvesuccess and flagp(op,'restart) then begin
      opstotry:=ops; % flag solvesuccess is set if op succeeds on any eqngrp
      if !*traceeqncnt then << terpri();
        write("There are ",count !*deteqns," equations remaining."); terpri() >>; 
      if !*traceeqns then << terpri();
        write("The equations remaining are ",egnums !*deteqns); terpri() >>; 
      end;
    end;
    terpri();
    if solveexit then <<
      write("Halting solve algorithm with early exit condition."); terpri()>>; 
    write("There are ",count !*deteqns," equations remaining."); terpri(); 
    if !*hideneqns then <<
      write("There are ",count !*hideneqns," hidden equations remaining."); terpri()>>; 
    if !*longhideneqns then <<
      write("There are ",count !*longhideneqns," long hidden equations remaining."); terpri()>>; 
    if !*deteqns
      then write("The equationgroup numbers of equations remaining are ", egnums !*deteqns);
    terpri();
    rmsubs();
  end;

symbolic procedure solve1byops(eqngrp,ops2try,opusagelist);
% solve a single eqngrp by trying lots of ops
begin scalar res,ops,op;
  ops:=ops2try;
  while ops and null res do begin
    op:=car ops;
    ops:=cdr ops;
    res:=apply(op,list eqngrp);
    end;
  if res then recordop(op,opusagelist);
  return res;
  end;

symbolic procedure recordop(op,oplist);
% oplist is an alist of list( ops , number of succesfull applications )
(if x
  then set(oplist, list(op,(cadr x)+1).delete(x,eval oplist))
  else set(oplist, list(op,1).eval oplist))
 where x=assoc(op,eval oplist);

symbolic procedure onebyoneop op;
% apply operation op to each eqngrp in !*deteqns that hasn't already been 
% done
% if succsessfull then delete eqn, & if the op is flagged as restart, then 
% exit
% It is left to the op to add any new eqns to !*deteqns
begin scalar eqngrpstotry,eqngrp,substosolve,selectfn;
  eqngrpstotry := !*deteqns;
  selectfn:=get(op,'selectfn); % a function to choose best eqn to try first
  if !*traceselectfn then <<write(tmstmp(),op," has selectfn ",selectfn); terpri()>>;
  while eqngrpstotry and not solvesuccess and not solveexit do begin
    if selectfn
      then eqngrp := apply(selectfn,list eqngrpstotry)
      else eqngrp := car eqngrpstotry;
      if !*traceselectfn then <<write(tmstmp(),"selected eqn#",eg!*no eqngrp,
                                      " from ",egnums eqngrpstotry);
                                terpri() >>;
    eqngrpstotry := delete(eqngrp,eqngrpstotry);
    if not (op member eg!*ops eqngrp) then begin
      if !*traceopson 
        then <<write(tmstmp(),"trying ",op," on eqn ",eg!*no eqngrp); terpri() >>;
      substosolve:=apply(op,list eqngrp); %returns t.truesubstosolve
      if not flagp(op,'dolots) then addopstried(op,eqngrp);
      if substosolve then begin
        tagsimpeqns(alleqngrps(), cdr substosolve);
        if not flagp(op,'keep) then !*deteqns:=delete(eqngrp,!*deteqns);
        solvesuccess:=t;
        end;
      end;
    end;
  end;

symbolic procedure lowordlowlen eqngrps;
% takes a list of equationgroups
% returns the one with lowest order derivatives, and fewest terms
begin scalar ez,loword,lowlen,eg2;
  ez:=car eqngrps;
  loword:=eg!*ord ez;
  lowlen:=eg!*len ez;
  for each eg2 in cdr eqngrps do
    if eg!*ord eg2<loword or (eg!*ord eg2=loword and eg!*len eg2<lowlen)
      then <<ez:=eg2; loword:=eg!*ord ez; lowlen:=eg!*len ez>>;
  return ez;
  end;

symbolic procedure lowlen eqngrps;
% takes a list of equationgroups
% returns the one with fewest terms
begin scalar ez,lowlen,eg2;
  ez:=car eqngrps;
  lowlen:=eg!*len ez;
  for each eg2 in cdr eqngrps do if eg!*len eg2 < lowlen
    then <<ez:=eg2; lowlen:=eg!*len ez>>;
  return ez;
  end;
 
symbolic procedure tagsimpeqns(eqngrplis, subs2solve);
begin scalar eqngrp;
  if subs2solve then for each eqngrp in eqngrplis do
%    if dunknsineqn(eg!*eqn eqngrp,subs2solve) 
    if intersect(dunknsE eg!*eqn eqngrp,subs2solve)
      then <<eg!*tagtosimp eqngrp;
             if !*tracetagsimp then write " tagtosimp eqn ",eg!*no eqngrp>>;
  end;

symbolic procedure simpeqns(eqngrplis, subs2solve);
% take the eqnlis and simplify any eqns involving the unknowns in subs2solve
% return the simplified list
begin scalar holdlis,seqngrp,eqngrp;
    if null subs2solve then return eqngrplis;
    rmsubs();
    holdlis :=nil;
    for each eqngrp in eqngrplis do 
      if dunknsineqn(eg!*eqn eqngrp,subs2solve) % only simp if needed
    then begin
        seqngrp := simpeqngrp eqngrp;
        if seqngrp then holdlis := seqngrp . holdlis;
        end
    else holdlis:=eqngrp . holdlis;
    return holdlis;
    end;

symbolic procedure dunknsineqn(eqn, dunkns);
% Does eqn involve any of the dunkns?
eqn and ( member(mvar eqn,dunkns)
          or (car mvar eqn='df and member(cadr mvar eqn,dunkns))
          or dunknsineqn(red eqn, dunkns));


%*******end of solve


% op!*proexp
%
% this operation expands phi prolongations in the equation.
% normaly(?) we want to have no unexpanded phi prolongations, 
% but it can be computationaly expensive to calculate them,
% and we might get some simplifications (single term equations)
% by doing this slowly rather than all at once.
% for now, we just calculate all pros in this eqn, but later, just do one level.
symbolic procedure op!*proexp eqngrp;
if highphisE eg!*eqn eqngrp then begin scalar pro,pros,sp;
  pros:=highphisE eg!*eqn eqngrp;
  for each pro in pros do begin
    sp:=numr simp!* pro;
    if sp and (mvar sp)=pro
      then mysetk(pro,mk!*sq mkpro(cadr pro, cddr pro));
    end;
  return t.pros;
  end;

%******end of op!*proexp


% op!*proeqn
%
% this operation expands phi prolongations as new equations.
%
symbolic procedure op!*proeqn eqngrp;
if complement(highphisE eg!*eqn eqngrp,!*proeqnsmade) then begin
  scalar pro,pros,proval,neqn,neqngrp,!*forcelowerprovalues;
  !*forcelowerprovalues:=nil;
  pros:=highphisE eg!*eqn eqngrp;
  for each pro in pros do
    if not(pro member !*proeqnsmade) then begin
      proval:=mkpro(cadr pro,cddr pro);
      neqn:=addf(multf(!*k2f pro,denr proval),
                 negf numr proval);
      neqngrp:=mkeqngrp(neqn,eg!*no eqngrp,'op!*proeqn);
      addeqngrp neqngrp;
      !*proeqnsmade:=pro.!*proeqnsmade;
      end;
    return list t %return t; %this should be a list
  end;

%******end of op!*proeqn

% op!*linchk
%
% this operation checks the validity of the determining equations,
% ie that the given equation is linear in the dunkns (with these appearing
% at the top level of ordering)
%
symbolic procedure op!*linchk eqngrp;
if not linchk1 eg!*eqn eqngrp then
  <<showf eg!*eqn eqngrp; 
     rederr list("not linear in determinable unknowns detected in equation ",
                  eg!*no eqngrp," ",eg!*eqn eqngrp)>>;


symbolic procedure linchk1 eqn;
% True if the equation is linear in dunkns.
% We only check mvars, not too likely to be more obscure than that.
null eqn or
  dunknk mvar eqn and nodunknsf lc eqn and linchk1 red eqn;

symbolic procedure nodunknsf u;
% u is a standard form: check that it doesn't have any dunkn functions in it.
% just check that mvars arent dunkns, not likely to get in such a mess
% that dunkns appear inside kernels so don't waste a lot of time looking :-)
null u or domainp u or
  not dunknk mvar u and nodunknsf lc u and nodunknsf red u;

%******end of op!*linchk

symbolic procedure op!*shr1tm eqngrp;
% solve a single term=0 type eqn. We only solve the ones here which will
% cause the other equations to contract: ie dunkn=0 or first order type.
% The equation is some factor times a derivative
% of a detunknown function. Strip the factor and put it on the divide list
%
if (not eg!*simpstatus eqngrp) and null red eg!*eqn eqngrp 
      and (!*slvpro or eg!*hiphiord eqngrp=0)
then begin scalar dfdunkn;
  if car (dfdunkn:= mvar eg!*eqn eqngrp)='df then
    if null cdddr dfdunkn % ie a first order derivative
    then <<depend1(cadr dfdunkn,caddr dfdunkn,nil);
           dividelog lc eg!*eqn eqngrp;
           if !*tracecute then write "!";
           showdeplevs();
           deathcert(eqngrp,list("Solved by op!*shr1tm: removed the ",
                     caddr dfdunkn," dependence in ",cadr dfdunkn));
           return list(t,dunknin dfdunkn)>>
    else return nil
  else <<mysetk(dfdunkn,0); % an unknown=0 type eqn
         dividelog lc eg!*eqn eqngrp;
         deathcert(eqngrp,list("Solved by op!*shr1tm: ",dfdunkn," set to zero."));
         return list(t,dfdunkn)>>;
  end;

symbolic procedure newmultpolyn(varords,deps);
% varords is a list of var1,ord1,var2,ord2, ...,var_i,deg_i
% deps is a list of dep1,dep2,...
% return a sf multivariate polynomial of degree deg_i in var_i
% and with coefficients which are arbitrary functions with dependence list deps.
begin
  if null varords then return !*k2f newarb1 deps;
  if not pairp varords or null cdr varords or not pairp cdr varords
    or not fixp cadr varords or (cadr varords)<0
      then rederr "bad power list";
  if (cadr varords)=0 then return newmultpolyn(cddr varords,deps);
  return addf(multf(!*p2f ((car varords).(cadr varords)),                          %-Highest ord 
                    newmultpolyn(cddr varords,deps)),                              %/  in 1st
              newmultpolyn((car varords).((cadr varords)-1).(cddr varords),deps)); %-Rest
  end;

symbolic procedure newmultpoly u;
begin scalar powlist,deps,x;
  while u and pairp u and cdr u and pairp cdr u and fixp cadr u do begin
    if not fixp cadr u and (cadr u)>0 then rederr "bad power list";
    powlist:=(!*a2k car u).(cadr u).powlist;
    u:=cddr u;
    end;
  deps:=for each x in u collect !*a2k x;
  return mk!*sq !*f2q newmultpolyn(powlist,deps);
  end;

symbolic operator newmultpoly;
rlistat '(newmultpoly);

symbolic procedure newpolyn(var,ord,deps);
% make a standard form polynomial of order ord in variable var
% with coefficients which are arbitrary functions with dependence list deps.
begin scalar sum;
    sum := !*k2f newarb1 deps;
    for j:=1:ord do 
      sum:=addf(sum, multf( !*k2f newarb1 deps,!*p2f (var .** j) ));  
    return sum;
    end;

symbolic procedure newpoly u;
begin scalar x;
  return mk!*sq !*f2q newpolyn(!*a2k car u,cadr u,
                               for each x in cddr u collect !*a2k x);
  end;

symbolic operator newpoly;
rlistat '(newpoly);

symbolic procedure newarb1 deps;
% return a new arbitrary function with dependence list deps as a kernel.
% the identifier C is reserved for this purpose, so that arbitrary
% functions can be labeled c(1) ... c(!*ccount) where !*ccount is a global variable
% whose value is the current number of arbitrary functions
begin scalar var;
    !*ccount:=!*ccount+1;
    if deps then for each var in deps do depend1(list(!*c,!*ccount),var,t);
    rmsubs();
    if !*tracenewarb 
      then <<terpri(); 
             write("Made new arbitrary function ",list(!*c,!*ccount),
                   " depending on ",deps); terpri() >>;
    return list(!*c,!*ccount);
    end;

symbolic procedure newarb u;
begin scalar x;
  return mk!*sq !*k2q newarb1 (for each x in u collect !*a2k x);
  end;

symbolic operator newarb;
algebraic (let newarb = newarb(nil));
rlistat '(newarb);

symbolic procedure newconst;
newarb1 nil;

symbolic operator newconst;
algebraic (let newconst=newconst());

symbolic procedure dets2poly(polyvar,deg);
begin scalar val,dunkn,deps,dunknsineqns;
  dunknsineqns:=dunknsineqns();
  for each dunkn in dunknsineqns do
    if val:=simp!* dunkn neq dunkn
      then dunknsineqns:=union(delete(dunkn,dunknsineqns),dunknsinE numr val);
  for each dunkn in dunknsineqns do
    if polyvar member (deps:=depsk dunkn)
      then mysetk(dunkn, mk!*sq !*f2q newpolyn(polyvar,deg,delete(polyvar,deps)));
  end;

symbolic operator dets2poly;

%******end of op!*shr1tm

symbolic procedure op!*simpeq eqngrp;
% if the eqngrp has been flagged for simplification then do it
if eg!*simpstatus eqngrp then begin scalar seqn;
  if !*tracecute then write "S";
  if (seqn:=simpeqngrp eqngrp) then addeqngrp seqn;
  return list t;
  end;

symbolic procedure simpeqngrp eqngrp;
% simplify the equation group.
% if the new equation is zero then return nothing
% otherwise make a new eqngrp
begin scalar seqn,bvars;
  rmsubs();
  bvars:=boundvarsE eg!*eqn eqngrp;
  seqn:=mynumr mysimpf eg!*eqn eqngrp;
  if intersect(bvars,!*usubvars) 
    then <<write "Found I can substitute original equations again ",
                 "when making eqngrp#",!*dets+1;terpri();
           seqn:=mynumr mysimpq sqliesubf seqn>>;
  if seqn then <<deathcert(eqngrp,list("Simplified to eqngrp#",!*dets+1));
                 return mkeqngrp(seqn,eg!*no eqngrp,'op!*simpeq)>>
          else deathcert(eqngrp,
                         list("Simplified trivialy -  at ",!*dets,
                              " equations."));
  end;

symbolic procedure simpeqn eqn;
begin scalar seqn;
  rmsubs();
  seqn := simp!* list('!*sq, eqn ./ 1, nil);
  return mynumr seqn;
  end;
  
%*******end of op!*simpeq

global '(!*op!*splitEonhighusonly !*op!*splitEonlvars !*op!*splitEonfunkns !*splitbyone);


%   Here we try to split an equation by the explicit apearance of an x or u 
% variable 'z', which doesn't appear in the dependence list of any of
% the determinable unknowns in the equation. Ie, z is a 'free' variable
% because it isn't a 'bound' variable. 
%
%   This means that we can write the equation as the sum
%
%   E = E(i)*alpha(i)
%
% where the alpha(i)'s don't contain any dunkns and are linearly independent
% in z, and the E(i)'s are independent of z. We then 'split' the equation E
% into the equations E(i).
%
%   To acheive this, we simply reorder the equation E so that any terms
% involving z are leading terms.
%
%   If z appears in any freeunknown f in E, then we have to make certain
% assumptions: the assumption we make is that each of the leading terms
% (or products of terms) are linearly independent in z- this effectivly gives
% an equation which we assume to be never zero. This is the most general
% assumption that could be made. If you are considering a group classification
% problem, then you will want to know what assumptions have been made, and so
% they should be reported...

symbolic procedure op!*splitEa eqngrp;
begin
  !*op!*splitEonhighusonly:=nil;
  !*op!*splitEonlvars:=t;
  !*op!*splitEonfunkns:=t;
  !*splitbyone:=t;
  return op!*splitE eqngrp;
  end;

symbolic procedure op!*splitEb eqngrp;
begin
  !*op!*splitEonhighusonly:=t;
  !*op!*splitEonlvars:=t;
  !*op!*splitEonfunkns:=nil;
  !*splitbyone:=nil;
  return op!*splitE eqngrp;
  end;

symbolic procedure op!*splitEc eqngrp;
begin
  !*op!*splitEonhighusonly:=nil;
  !*op!*splitEonlvars:=t;
  !*op!*splitEonfunkns:=nil;
  !*splitbyone:=nil;
  return op!*splitE eqngrp;
  end;

symbolic procedure op!*splitEd eqngrp;
begin
  !*op!*splitEonhighusonly:=nil;
  !*op!*splitEonlvars:=nil;
  !*op!*splitEonfunkns:=t;
  !*splitbyone:=t;
  return op!*splitE eqngrp;
  end;

symbolic procedure op!*splitE eqngrp;
%
% diferent types of op!*splitE possible, depending on values of
% !*op!*splitEonhighusonly, !*op!*splitEonlvars !*op!*splitEonfunkns !*splitbyone
%
if (not eg!*simpstatus eqngrp or !*splitonsimp)
    and (eg!*hiphiord eqngrp=0 or !*splitonpro)then
begin scalar eqn,bvars,fvar,feqnlis,linindeplis,term,neqngrp;
  eqn:=eg!*eqn eqngrp;
  bvars:=boundvarsE eqn; % the vars bound to dunkns. pass this as fluid to lkfvar
  if !*op!*splitEonlvars then begin
    bvars:=union(bvars,fdepsf eqn); % dont cause any restrictions here...
    lkfvarf eqn;         % look for a freevar in eqn. Sets fvar if it finds one
    end;
  if !*op!*splitEonfunkns and not fvar then begin
    fvar:=complement(fdepsf eqn,bvars);  % list of vars in funkns but not dunkns
    if fvar then fvar:=car fvar;         % just take one of them
    end;
  if not fvar then return;
  feqnlis:=f2jns eqn;
  linindeplis:=for each term in feqnlis collect car term;
  if !*traceop!*splitE then begin
    write(tmstmp(),"Split by ",fvar," in op!*splitE"); terpri();
    if !*tracehard then <<write(tmstmp(),"found feqnlis "); terpri();
                          write(feqnlis); terpri() >>;
    end;
  showindepifneed linindeplis;
  for each freepair in feqnlis do begin
    neqngrp:=mkeqngrp(cadr freepair,eg!*no eqngrp,
                       list("op!*splitE by taking the coefficient of ",car freepair));
    if eg!*simpstatus eqngrp then eg!*tagtosimp neqngrp; % if old one needed simping, this does 2
    addeqngrp neqngrp;
    end;
  deathcert(eqngrp,list("Split by ",linindeplis," in op!*splitE to give new equations ",
            !*dets+1 - count feqnlis," to ",!*dets));
%  solvesuccess:=t;
  if (!*verifyop!*splitE and op!*splitEcheck(eqn,feqnlis)) or null cdr feqnlis 
    then <<write "op!*splitE on eqngrp ",eg!*no eqngrp;terpri();
           if !*traceop!*splitE then showeqngrp eqngrp;
           write "Got fvar ",fvar," and split to get ",count feqnlis,"new eqns";terpri();
           if !*traceop!*splitE then showeqngrp eqngrp;
           showindep(fvar,linindeplis);
           if !*tracehard then write "got feqnlis",feqnlis;
           rederr("Bad op!*splitE")>>;
  if not (fvar member !*splitvars) then !*splitvars:=fvar.!*splitvars;
  return list t;
  end;

symbolic procedure op!*splitEcheck(eqn,feqnlis);
begin scalar term;
  for each term in feqnlis do eqn:=addf(eqn,negf multf(car term,cadr term));
  return eqn;
  end;

symbolic procedure lkfvarf eqn;
% look for a freevar in standardform eqn. 
% Set fvar if we find one
if null eqn or domainp eqn or fvar then nil % stop looking if one already found
else <<lkfvar mvar eqn; lkfvarf lc eqn; lkfvarf red eqn>>;

symbolic procedure lkfvar u;
% could look in funkns as well, but best to do that if all else fails
if atom u or fvar then nil
else if car u=!*x or car u=!*u 
  then if !*op!*splitEonhighusonly 
    then if car u=!*u and cddr u and not (u member bvars) 
      then fvar:=u 
      else nil
    else if not (u member bvars) 
      then fvar:=u 
      else nil
  else <<lkfvar car u; lkfvar cdr u>>;
  
symbolic procedure f2jns u;
% u is a standardform
% normalise u to a jnseqn
% ie, bring the fvar factors to the front
% This may split the term...
%
%   An seqn is a list of list(freeterm,coeff)
% where freeterm is a product of factors each involving fvar (passed as fluid)
% and coeff is the coefficient of that product, possibly involving fvar.
% Both freeterm and coeff are represented as standardforms.
% 
%   A jnseqn is an seqn where the coeffs dont involve fvar;
% there are no repeated freeterms; each freeterm is either a power of fvar,
% a function of fvar, or a product of such with each product ordered by factor;
% and the seqn has no repeated freeterms and is orded by freeterm.
%
if null u then nil
else if domainp u then list list (1,u)
else addjns(multjns(p2tjns lpow u,f2jns lc u),f2jns red u);
 
symbolic procedure p2tjns p;
% p is a standardpower, convert it to a single term of an jnseqn
if !*op!*splitEonhighusonly 
  then if highvar car p
    then list(!*p2f p,1)
    else list(1,!*p2f p)
  else if (car p=fvar or (not dunknk car p and fvarin car p)) 
%  else if fvarin car p    % why not just use this instead?
    then list(!*p2f p,1)
    else list(1,!*p2f p);
 
symbolic procedure multjns(u,v);
% u is a single term of an jnseqn
% v is a full jnseqn
% multiply them
if not u or not v then nil
else addtjns(list(multf(car u,caar v),multf(cadr u,cadar v)),
            multjns(u,cdr v));

symbolic procedure addjns(u,v);
% u and v are both jnseqjns. add them.
if null v then if null u then nil else u
else if null u then v
else <<for each x in u do v:=addtjns(x,v); v>>;
 
symbolic procedure addtjns(u,v);
% u is a single term of an jnseqn.
% v is a full jnseqn.
% add them
if null v then if null u then nil else list u
else if null u then v
else if car u=caar v then list(car u,addf(cadr u,cadar v)).cdr v
else if ordp(car u,caar v) then u.v
else (car v).addtjns(u,cdr v);
 
symbolic procedure fvarin u;
% does fvar appear in u?
u=fvar or (funknk u and fvar member depsk u)
or (pairp u and (fvarin car u or fvarin cdr u));


% Bugfix for highvar.  
% If eqn involves freeunknowns then it is possible 
% that a kernel is just an idenitfier and not a list. MJ.

%symbolic procedure highvar v;
%%is kernel v a df of u variable?
%car v=!*u and cddr v;

symbolic procedure highvar v;
%is kernel v a df of u variable?
% ie. is v = (u,a,b) for some integers a,b. eg (u,1,1).
if pairp v then car v=!*u and cddr v;

%end of bugfix for highvar.


symbolic procedure boundkern u;
% is u a function of bvars? 
% ie does it explicity involve any x or u vars
% which are in bvars?
% we do not consider funkns here, because this is used for
% splitting by manyvars at once, and funkns should be split 1 at a time.
pairp u and (
  (if car u=!*u or car u=!*x then u member bvars else nil)
  or boundkern car u or boundkern cdr u );

symbolic procedure showindepifneed u;
if (specfnsl u and not ({fvar,u} member !*linindepconditions)) or !*traceop!*splitE
  then <<showindep(fvar,u); !*linindepconditions:={fvar,u}.!*linindepconditions>>;
  
symbolic procedure showindep(fv,u);
% u is a list of terms to be lin indep.
% show the linear independence conditions, if they are new ones.
begin
  terpri();
  write "Must have all of"; terpri();
  for each x in u do showf x;
  terpri();
  write("linearly independent in ",fv); terpri(); terpri();
  end;

symbolic procedure showindeps;
for each u in !*linindepconditions do showindep(car u, cadr u);

symbolic operator showindeps;

symbolic procedure specfnsl u;
% u is a list of terms used to split the eqn.
% are there any special functions in there that we should pay attention to?
% pay attention to anything except powers of x,u
% Ie notice functions like sin,cos,
% and especially notice freeunknow functions.
u and (specfnsf car u or specfnsl cdr u);

symbolic procedure specfnsf u;
u and not domainp u and (specfnsk mvar u or specfnsf lc u or specfnsf red u);

symbolic procedure specfnsk u;
not (pairp u and (car u=!*u or car u=!*x));

%******end of op!*splitE

symbolic procedure op!*exp1tm eqngrp;
% Solve a single term=0 type eqn. 
% Substitute an appropriate polynomial.
%
% Only solve if there is just one df variable, otherwise expressions
% become too complicated.
%
if (not eg!*simpstatus eqngrp)  and null red eg!*eqn eqngrp 
      and (!*slvpro or eg!*hiphiord eqngrp=0) and (car mvar eg!*eqn eqngrp)='df
  then begin
    scalar dunkn,dfdunkn,var,vars,ord,deps,dfvars,val,toobig;
    % the actual eqn, the unknown function, var to int wrt, order of
    dividelog lc eg!*eqn eqngrp;
    dfdunkn := mvar eg!*eqn eqngrp;
    dunkn := cadr dfdunkn;
    dfvars:= cddr dfdunkn;
    if not !*op!*exp1tmexpand and 
      cdr dfvars and ((not numberp cadr dfvars) or (cddr dfvars)) then return;
      % ie df'ed by more than one var...
    deps := depsk dunkn;
%    if deps then deps:=cdr deps; % what is this doing here? get with it james :-)
    val:= nil ./ 1;
    while dfvars and not toobig do begin
      var := car dfvars;
      vars:= var.vars;
      if cdr dfvars and numberp cadr dfvars 
        then <<ord:=cadr dfvars; dfvars:=cddr dfvars>>
        else <<ord:= 1;          dfvars:=cdr  dfvars>>;
      if !*tracecute then write "^",ord-1;
      if ord > !*exp1tmordlimit 
        then toobig:=t
        else val := addsq(val,!*f2q newpolyn(var,ord-1,delete(var,deps)));
      end;
    if toobig then << 
      write "in op!*exp1tm found polynomial of order ",ord-1," which is too big ";terpri();
      write "(!*exp1tmordlimit=",!*exp1tmordlimit,")";terpri();
      return>>;
    deathcert(eqngrp,list("Solved by op!*exp1tm: ",dunkn," polynomial in ",vars));
    mysetk( dunkn, mk!*sq val);
    return list (t,dunkn);
    end;

%*******end of op!*exp1tm

%**** slv operations
%
% Here we solve an equation by making an assignment for a dunkn in the equation.
%
% There are several ways of doing this:
%
% 1. We have sdunkn1=dunkn2 as the eqn.
%    Then make the assignment, reducing the dependence of dunkn2 if needed.
%    This is done as op op!*slvtwo
%
% 2. We have sdunkn=val, where val does not depend on any vars more than dunkn.
%    We also have val depends on less than sdunkn
%    Then remove the dependence of sdunkn on these vars. (we don't really solve by assignment)
%    this is done as op op!*drpdep
%
% 3. We have sdunkn=val, where val does not depend on any vars more than dunkn.
%    We also have sdunkn depends on xdeps, which only appear in val explicitly,
%    ie xdeps=complement(depsk sdunkn,union(fdeps val, boundvarsE val))
%    so val=sum_i Ai * fi(xdeps) where Ai does not depend on xdeps,
%    and the fi are linearly independent explicict functions of xdeps w/o dunkns.
%    Then make the assignment sdunkn=sum_i ci * fi(xdeps) and include the equations ci=Ai
%    where the ci are new dunkns with appropriate dependence.
%    this is done as op op!*slvspl
%
% 4. We have sdunkn=val, where val does not depend on any vars more than dunkn.
%    Just make the assignment.
%    this is done as op op!*slvall

symbolic procedure op!*slvtwo eqngrp;   
if red eg!*eqn eqngrp and null red red eg!*eqn eqngrp then
begin scalar eqn,sdunkn,sval,xdeps;
  eqn:=eg!*eqn eqngrp;
  if % (null red eqn) or (red red eqn) or
         not (domainp lc eqn)     or (car mvar eqn='df)
      or not (domainp lc red eqn) or (car mvar red eqn='df) 
    then return;
  sdunkn:=eg!*highdfdunkn eqngrp;
  sval:=eg!*subval eqngrp;
  if !*traceslv then begin 
    terpri();
%    write "In op!*slvtwo on eqngrp ",eg!*no eqngrp," with xdeps ",
%           complement(depsk sval,depsk sdunkn); terpri();
    showeqngrp eqngrp;
    write "op!*slvtwo: setting ",sdunkn," to ";terpri(); showsq sval;
    write "removing the dependence of ",mvar numr sval," on ",
          complement(depsk mvar numr sval,depsk sdunkn); terpri();
    end;
  for each var in complement(depsk mvar numr sval,depsk sdunkn) do
    <<depend1(mvar numr sval,var,nil);
      if !*tracecute then write "!";
      showdeplevs();
      xdeps:=t;
     >>;
  if xdeps
    then deathcert(eqngrp,list("Solved by slvtwo: solved for ",sdunkn,
                               " and removed dependence of ",mvar numr sval,
                               " on ",complement(depsk mvar numr sval,depsk sdunkn)))
    else deathcert(eqngrp,list("Solved by slvtwo: solved for ",sdunkn));
  mysetk(sdunkn, mk!*sq sval );
  if xdeps
    then return list (t,sdunkn,mvar numr sval)
    else return list (t,sdunkn);
  end;  

symbolic procedure op!*drpdep eqngrp;   
if eg!*slvterm eqngrp then begin scalar sterm,valn,xdeps;
  if !*traceslv then <<write tmstmp(),"drpdep on eqngrp ",eg!*no eqngrp;terpri()>>;
  sterm:=eg!*slvterm eqngrp;
  valn:= addf(negf eg!*eqn eqngrp, !*t2f sterm);
  xdeps:=complement(depsk caar sterm,union(alldepsf cdr sterm,alldepsf valn));
  if !*traceslv then begin
    showeqngrp eqngrp;
    write "removing the dependence of ",caar sterm," on ",xdeps; terpri();
    end;
  if null xdeps then return;
  for each var in xdeps do <<
    depend1(caar sterm,var,nil);
    if !*tracecute then write "!";
    showdeplevs();
    >>;
  return list (t,car sterm);
  end;  
  
symbolic procedure op!*slvspl eqngrp;
%
% The eqn can be solved for sdunkn with xdeps dependence in sdunkn
% only appearing explicitly.
%
% Use new dunkns for the coefficients of sdunkn if they are long enough,
% (more than 1 term, or non-numerical coeff)
% and then include equations for those dunkns.
%
begin scalar eqn, sterm, sdunkn,sdunkncoeff,subvaln, xdeps,jnsval,
             newvaln,thiscoeff,thisxtm,modcoeff,neweqnlis;
%  if !*traceslv then <<write tmstmp(),"slvspl on eqngrp ",eg!*no eqngrp;terpri()>>;
  eqn:=eg!*eqn eqngrp;
  if not (sterm:=eg!*slvterm eqngrp) then return;
  sdunkn:=caar sterm;
  sdunkncoeff:=cdr sterm;
  subvaln := addf(negf eqn, !*t2f sterm);
%  if !*traceslv then begin
%    write "got sdunkn ",sdunkn;terpri();
%    write "got sdunkncoeff ";terpri();showf sdunkncoeff;
%    write "got subvaln ";terpri();showf subvaln;
%    end;
  xdeps:=complement(depsk sdunkn, union(union(fdepsf subvaln,fdepsf sdunkncoeff),
                                        boundvarsE subvaln)                     );
  if not xdeps then return;
  if !*traceslv
     then <<terpri();
            write "In op!*slvspl on eqngrp ",eg!*no eqngrp," with xdeps ",xdeps; terpri();
            showeqngrp eqngrp>>;
  % now split valn/cdr sterm...??? how
  fvar:=car xdeps;            % We can only split by one at the moment... fvar is passed as fluid
  jnsval:=f2jns subvaln;      % This is the list of { free terms , coeffs of free terms }
%  if !*traceslv then <<write "got jnsval "; terpri(); write jnsval;terpri()>>;
  for each pr in jnsval do begin
    thisxtm:=car pr;
    thiscoeff:=cadr pr;
    if null red thiscoeff and (null lc thiscoeff or domainp lc thiscoeff) 
        and (!*slvdfs or (car mvar thiscoeff) neq 'df)
      then   modcoeff:=thiscoeff
      else <<modcoeff:=!*k2f newarb1 alldepsf thiscoeff;
             if !*tracecute then write ";";
             neweqnlis:=addf(negf modcoeff,thiscoeff).neweqnlis>>;
    newvaln:=addf(newvaln, multf(modcoeff,thisxtm));
    end;
  %
  % it might be good to factorise the denr (cdr sterm) and place as many factors as possible
  % onto the neweqns, rather than the newval.
  %
  if !*tracecute then write "{",mytermsf newvaln,"/",mytermsf sdunkncoeff,"}";
  if !*traceslv then begin
    write "solved eqn for ",sdunkn," to get value "; terpri();
    showsq (newvaln ./ sdunkncoeff);
    if neweqnlis
      then <<
        write "subject to the new equations";terpri();
        for each neweqn in neweqnlis do <<showf neweqn; showdepsineqn neweqn>>;
        >>
      else <<write "subject to no new equations";terpri()>>;
    terpri();
    end;
  for each neweqn in neweqnlis do addeqngrp mkeqngrp(neweqn,list eg!*no eqngrp,'slvspl);
  deathcert(eqngrp,list("Solved by slvspl: solved for ",sdunkn));
  mysetk(sdunkn, mk!*sq (newvaln ./ sdunkncoeff));
  return list (t,sdunkn);
  end;


global'(!*slvonless);

symbolic procedure op!*slvshr eqngrp;   
begin %scalar !*slvonless;
  !*slvonless:=t;
  return slv eqngrp;
  end;

symbolic procedure slv eqngrp;
begin scalar eqn, sterm, subval;
  if !*traceslv then <<write tmstmp(),"slv on eqngrp ",eg!*no eqngrp;terpri()>>;
  eqn:=eg!*eqn eqngrp;
  if not (sterm:=eg!*slvterm eqngrp) then return;
  subval := addf( negf eqn, !*t2f sterm);
  subval := multsq(subval./ 1,1 ./ cdr sterm);
  if !*traceslv then begin
    write tmstmp(),"slv on eqngrp ",eg!*no eqngrp;terpri();
    write "sterm is ",sterm; terpri();
    write "depsk caar sterm is ",depsk caar sterm;terpri();
    write "boundvarsE numr subval is ", boundvarsE numr subval;terpri();
    write "complement(depsk caar sterm, boundvarsE numr subval) is ",
           complement(depsk caar sterm, boundvarsE numr subval);terpri();
    end;
  if !*slvonless and 
     (funknsp eqn or not complement(depsk caar sterm, boundvarsE numr subval))
    then return;
               % ie only do the substitution if it gives less implicit dependence
               % and involves no freeunknowns
  deathcert(eqngrp,list("Solved by slv: solved for ",caar sterm));
  mysetk(caar sterm,mk!*sq subval);
  return list (t,caar sterm);
  end;

symbolic procedure op!*slvalldfs eqngrp;
begin scalar !*slvdfs;
  !*slvdfs:='t;
  return op!*slvall eqngrp;
  end;

symbolic procedure op!*slvalldfsearly eqngrp;
if null !*canonlist then begin scalar !*slvdfs;
  !*slvdfs:='t;
  return op!*slvall eqngrp;
  end;

symbolic procedure op!*slvall eqngrp;
begin scalar eqn, sterm, subval;
%  if !*traceslv then <<write tmstmp(),"slvall on eqngrp ",eg!*no eqngrp;terpri()>>;
  eqn:=eg!*eqn eqngrp;
  if not (sterm:=eg!*slvterm eqngrp) then return;
  subval := addf( negf eqn, !*t2f sterm);
  subval := multsq(subval./ 1,1 ./ cdr sterm);
  if (null !*slvdfs) and dfsinEp numr subval
    then <<if !*traceslv then <<
             write "In op!*slvall found value with numerator of length ",mytermsf numr subval,
                   " and denominator of length ",mytermsf denr subval,
                   " which had DFs of dunkns, which we aren't using in slv values just now.";
             terpri();
             return>>;
           if !*tracecute then write tmstmp(),"P",mytermsf numr subval,"/",mytermsf denr subval;
         >>;
  if (mytermsf numr subval) > !*asmtnumlimit or (mytermsf denr subval) > !*asmtdenlimit
    then <<write tmstmp(),
                 "In op!*slvall found value with numerator of length ",mytermsf numr subval,
                 " and denominator of length ",mytermsf denr subval," which is too long.";
           terpri();
           return>>;
  deathcert(eqngrp,list("Solved by slvall: solved for ",caar sterm));
  mysetk(caar sterm,mk!*sq subval);
  return list (t,caar sterm);
  end;

symbolic procedure sterm eqngrp;
% returns term (dunkn.**1).*coeff if eqn can be solved for dunkn
% move through looking for the one with shortest lc 
% or if it is the highdfdunkn if flag !*slvforhigh is true
begin scalar eqn, highdfdunkn, reqn,stm, stmlctms, alldeps;
  eqn:=eg!*eqn eqngrp;
  highdfdunkn:=eg!*highdfdunkn eqngrp;
  alldeps:=alldepsf eqn;
  reqn:=eqn;
  while reqn do
    if not (car mvar reqn='df) and null complement(alldeps,depsk mvar reqn)
      and not diffd(reqn,mvar reqn)
        then if !*slvforhigh and intersect(!*canonlist,alleqngrps()) and (mvar reqn)=highdfdunkn
          then <<stm:=lt reqn;reqn:=nil>>
          else if (null stm) or (mytermsf lc reqn) < stmlctms
            then <<stm:=lt reqn;                % ie the term dunkn.*coeff
                   stmlctms:=mytermsf lc reqn;
                   reqn:=red reqn>>
            else reqn:=red reqn
        else reqn:=red reqn;
  return stm;
  end;

symbolic procedure dfsinEp eqn;
% are there any df of dunkn terms in E?
eqn and ((car mvar eqn)='df or dfsinEp red eqn);

%*******end of slv

symbolic procedure otherboundvarsE(dunkn,eqn)$
% Find all variables in depl!* for the detunknowns appearing in standardform eqn
% except for terms involving dunkn.
% Remember that the detunknowns appear linearly ie either as 
% themselves or in df statements, and at the top level of ordering.
if null eqn then nil
else if (dunknin mvar eqn)=dunkn then otherboundvarsE(dunkn,red eqn)
else union(depsk mvar eqn, otherboundvarsE(dunkn,red eqn) )$

symbolic procedure op6 eqngrp$
%
% Reduce the dependence of an equation:
% If a dunkn (NOT a dfdunkn though) depends on less than lfdeps and otherbvars
% then restrict by the complement-
% ie for any dunkn which does depend on one or more of the restricted vars
% then replace it with a new dunkn w/o that dependence.
%
% Ie, a 1-term eqn(s) would be found by diffing once by (each) xdep
% So why not just leave this to op!*get1tm?
%
begin scalar eqn, resteqn, dunkn, lfdeps, xdeps, v$
  eqn:=eg!*eqn eqngrp$
  resteqn:=eqn$
  lfdeps:=union(ldepsf eqn,fdepsf eqn)$
  while resteqn do  % ie for each dfdunkn     
    if  (car (dunkn:=mvar resteqn)) neq 'df and not diffd(eqn,dunkn) 
        and (xdeps:=complement(depsk dunkn,union(otherboundvarsE(dunkn,eqn),lfdeps)))
      then resteqn:=nil          % so dunkn does not depend only on xdeps
      else resteqn:=cdr resteqn$ 
  if not xdeps then return$
  if !*traceop6 then <<write "op6 on eqngrp";terpri();showeqngrp eqngrp>>;
  for each v in xdeps do depend1(dunkn,v,nil)$
  if !*tracecute then write "!";
  showdeplevs();
  if !*tracecute then
    <<write ("Equation ",eg!*no eqngrp," op6: removing the dependence of ", dunkn," on ",xdeps)$
     terpri()>>;
  return list(t,dunkn)$
  end$

%******end of op6

symbolic procedure op!*get1tm eqngrp$
%
% Can we get a one-term equation df(R,x,n+m+1)=0
% by differentiating eqn m+1 times?
%
% Must have 
%
% (1) term df(R,x,n) is only dunkn term in eqn involving x
% (2) x does not appear explicitly non-polynomialy
% (3) x appears explicitly polynomialy of degree at most m       { (2) => (3) for some m }
% (4) x does not appear in coeff C of term df(R,x,n)
% (5) there is another dfdunkn in eqn depending on something R doesnt (if !*op!*get1tmonless)
%
% This is an important shortcut to just differentiating eqn by x m times.
%
% We would then get eqn df( C*df(R,x,n) ,x,m+1)       (even if (4) is not met).
%
% We only look for 1tm eqns with one df-var because that is all we like to solve. 
%
% If n=0 then we have the choice of each var in depsk R for x
%
% If condition (4) is not met, then we could put eqn C*df(R,x,n) + poly(x,m,otherdeps) =0 
% Dont do this yet.
%
% car='df
% cadr=dunkn
% caddr=x
% cadddr=n
%
begin scalar eqn, dfdunkns, dfdunkn, var, vars, R, x, n, m, neqn, neqngrp$
  eqn:=eg!*eqn eqngrp$
  if null red eqn then return;                      % already a one-term eqn?
  dfdunkns:=onceonlydfdunkns eqn;                   % the dfdunkns whose dunkns only appear once.
                                                    % dfdunkns now satisfy (1.i)
  for each dfdunkn in dfdunkns do if mixeddf dfdunkn
    then dfdunkns:=delete(dfdunkn,dfdunkns);        % dfdunkns now satisfy (1.ii)
  while dfdunkns and null R do begin
    dfdunkn :=car dfdunkns;
    dfdunkns:=cdr dfdunkns;
    var:=vars:=nil;
    if (car dfdunkn)='df then var:=caddr dfdunkn else vars:=depsk dfdunkn;
    if !*op!*get1tmonless and null complement(alldepsf eqn,depsk dfdunkn)
      then var:=vars:=nil;                          % now satisfies (5)
    if var and not (var member otherboundvarsE(dunknin dfdunkn,eqn))
                                                    % now satisfies (1.iii)
           and not (var member nonpolyvarsE eqn)    % not likely to call nonpolyvars
                                                    % now satisfies (2)
           and not (var member alldepsf coeffE(dfdunkn,eqn))
                                                    % now satisfies (4)
      then <<R:=dunknin dfdunkn; x:=var;
             if cdddr dfdunkn then n:=cadddr dfdunkn else n:=1>>;
    if vars and (vars:=complement(vars,otherboundvarsE(dfdunkn,eqn)))
                                                    % now satisfies (1.iii)
            and (vars:=complement(vars,nonpolyvarsE eqn))
                                                    % now satisfies (2)
            and (vars:=complement(vars,alldepsf coeffE(dfdunkn,eqn)))
                                                    % now satisfies (4)
      then <<R:=dfdunkn; x:=car vars; n:=0>>;
    end;
  if not R then return;
  m:=highpolydegf(eqn, x);
  neqn:=!*k2f  list('df, R, x, n+m+1);
  dividelog coeffE(dfdunkn,eqn);
  addeqngrp (neqngrp:=mkeqngrp(neqn,eg!*no eqngrp,'op!*get1tm));    
  if !*traceop!*get1tm then <<showeqngrp eqngrp; showeqngrp neqngrp>>;
  return list(t);
  end$


symbolic procedure op!*get1tm9 eqngrp$
%
% Can we get a one-term equation df( C*df(R,x,n) ,x,m+1)=0
% by differentiating eqn m+1 times?
%
% Must have 
%
% (1) term df(R,x,n) is only dunkn term in eqn involving x
% (2) x does not appear explicitly non-polynomialy
% (3) x appears explicitly polynomialy of degree at most m       { (2) => (3) for some m }
% (5) there is another dfdunkn in eqn depending on something R doesnt
%
% This is an important shortcut to just differentiating eqn by x m times.
%
% We would then get eqn df( C*df(R,x,n) ,x,m+1)       (even if (4) is not met).
%
% We only look for 1tm eqns with one df-var because that is all we like to solve. 
%
% If n=0 then we have the choice of each var in depsk R for x
%
% We put eqn C*df(R,x,n) + poly(x,m,otherdeps) =0 
%
% car='df
% cadr=dunkn
% caddr=x
% cadddr=n
%
begin scalar eqn, dfdunkns, dfdunkn, var, vars, R, x, n, m, neqn, neqngrp $
  eqn:=eg!*eqn eqngrp$
  dfdunkns:=onceonlydfdunkns eqn;                   % the dfdunkns whose dunkns only appear once.
                                                    % dfdunkns now satisfy (1.i)
  for each dfdunkn in dfdunkns do if mixeddf dfdunkn
    then dfdunkns:=delete(dfdunkn,dfdunkns);        % dfdunkns now satisfy (1.ii)
  while dfdunkns and null R do begin
    dfdunkn :=car dfdunkns;
    dfdunkns:=cdr dfdunkns;
    var:=vars:=nil;
    if (car dfdunkn)='df then var:=caddr dfdunkn else vars:=depsk dfdunkn;
    if !*op!*get1tmonless and null complement(alldepsf eqn,depsk dfdunkn)
      then var:=vars:=nil;                          % now satisfies (5)
    if var and not (var member otherboundvarsE(dunknin dfdunkn,eqn))
                                                    % now satisfies (1.iii)
           and not (var member nonpolyvarsE eqn)    % not likely to call nonpolyvars
                                                    % now satisfies (2)
           and not (var member alldepsf coeffE(dfdunkn,eqn))
                                                    % now satisfies (4)
      then <<R:=dunknin dfdunkn; x:=var;
             if cdddr dfdunkn then n:=cadddr dfdunkn else n:=1>>;
    if vars and (vars:=complement(vars,otherboundvarsE(dfdunkn,eqn)))
                                                    % now satisfies (1.iii)
            and (vars:=complement(vars,nonpolyvarsE eqn))
                                                    % now satisfies (2)
      then <<R:=dfdunkn; x:=car vars; n:=0>>;
    end;
  if not R then return;
  m:=highpolydegf(eqn,x);
  neqn:=addf( list( (dfdunkn .** 1) .* coeffE(dfdunkn,eqn) ),
              newpolyn(x,m,delete(x,depsk R)) );
  addeqngrp (neqngrp:=mkeqngrp(neqn,eg!*no eqngrp,'op!*get1tm9));    
  if !*traceop!*get1tm then <<showeqngrp eqngrp; showeqngrp neqngrp>>;
  return list(t);
  end$


symbolic procedure onceonlydfdunkns eqn;
% the dfdunkns whose dunkns only appear once.
begin scalar dfdunkns, dfdunkn, m1dunkns, ans;
  dfdunkns:=dfdunknsE eqn;
  while dfdunkns do begin
    dfdunkn :=car dfdunkns;
    dfdunkns:=cdr dfdunkns;
    if not (dunknin dfdunkn member m1dunkns) then
      for each d in dfdunkns do
        if (dunknin d)=(dunknin dfdunkn) then m1dunkns:=(dunknin d).m1dunkns;
    end;
  for each d in dfdunknsE eqn do if not (dunknin d member m1dunkns) then ans:=d.ans;
  return reverse ans;
  end;
  
symbolic procedure highpolydegf(eqn,var);
% we rely heavily here on the representation in standard forms:
% if var appears as mvar (a kernel) of u then var _wont_ appear as a mvar of lc u,
% and only of a lesser degree in red u.
if null eqn or domainp eqn then 0
else if (mvar eqn)=var then ldeg eqn
else max2(highpolydegf(lc eqn,var), highpolydegf(red eqn,var));

%******end of op!*get1tm

symbolic procedure op!*xdpshr eqngrp;
<<!*op!*xdpshr:=t; op!*xdp eqngrp>>;

symbolic procedure op!*xdpxpd eqngrp;
<<!*op!*xdpshr:=nil; op!*xdp eqngrp>>;

symbolic procedure op!*xdp eqngrp$
%
% op!*xdp (new name soon)
%
% restrict by xdeps
%
% Condition:
%
% dunkn D does not depend on z's.
% Every dunkn in eqn which _does_ depends on z's doesnt depend on something D does,
% and so it is ok to solve the eqn for D, after we make apropriate restrictions
% Also, we expect D to depend on less than bvars, else the restriction would not be needed.
%
% (at the moment, we have D depends on zb which _no_ other dunkn depends on)
%
% procedure xdeps returns list of possibilities: { {(z's), D}, ... }
%
begin scalar xdepslis, mostx, res;
  if !*op!*xdpshr
    then xdepslis:=xdepsshr eg!*eqn eqngrp
    else xdepslis:=xdepsxpd eg!*eqn eqngrp;
  while xdepslis and null res do begin
    mostx:=findmostx xdepslis;                      % the xdunkn combo with most xdeps      
    xdepslis:=delete(mostx, xdepslis);
    res:=op!*xdp1(eqngrp,mostx);
    end;                
  return res;
  end;  

symbolic procedure findmostx xdepslis;
if null xdepslis then nil
else if null cdr xdepslis then car xdepslis
else begin scalar mostx;
  mostx:=car xdepslis;
  xdepslis:=cdr xdepslis;
  while xdepslis do begin
    if count(car mostx) > count(caar xdepslis) then mostx:=car xdepslis;
    xdepslis:=cdr xdepslis;
    end;
  return mostx;
  end;

symbolic procedure op!*xdp1(eqngrp,xdeps)$
%
% Given one possibility xdeps= {(z's), D}, make restricions and solve if possible
%
begin scalar eqn, xdunkn, allxdeps, lxdeps, neqn, mixed$
  eqn:=eg!*eqn eqngrp$
  xdunkn:=dunknin cadr xdeps;                            
  allxdeps:=car xdeps;
  allxdeps:=complement(allxdeps,nonpolyvarsE eqn);      % we cant restrict by these vars (yet)
  lxdeps:=intersection(allxdeps,ldepsf eqn);            % we will need to set values for these
  xdeps:=complement(allxdeps,fdepsf eqn);               % we dont restrict by fdeps yet...
  if null !*op!*xdpon0 then xdeps:=complement(xdeps,lxdeps); 
  if null xdeps then return$
  if !*traceop!*xdp 
    then <<write "op!*xdp got xdunkn ",xdunkn," with xdeps ",xdeps," on eqngrp";terpri();
           showeqngrp eqngrp>>;
  neqn:=eqn;
  if !*op!*xdpon0 and lxdeps 
    then neqn:=sublxdeps01(lxdeps,neqn,xdunkn);          % choose values for explicit xdeps...
  if null neqn then return;                              % (needed to sub other than 0 or 1)
  for each d in dunknsE neqn do 
    if (dunknin d)=xdunkn and mixeddf d then mixed:='t;
  if mixed then return;                                  % we cant solve these eqns
  neqn:=restricteqn(eqn,neqn,xdunkn,xdeps);
  if null neqn 
   then <<write "op!*xdp got 0"; terpri(); return>>;        % must be redundancies left
  return op!*xdpsolvebyops(neqn,eqngrp,xdunkn,xdeps);       % solve this eqn if we can, else drop it
  end;


symbolic procedure mixeddf(dunkn);
% is dunkn a mixed derivative?
% car = df
% cadr = dunkn
% caddr = var1
% cdddr = nil/ cadddr = number & cddddr = nil then good
% cadddr = var2/ cadddr = number & caddddr = var2 then bad
(car dunkn)='df and (cdddr dunkn) and (not (numberp cadddr dunkn) or cddddr dunkn);
% car d='df and (null cdddr d or (numberp cadddr d and null cddddr d)); % ???

% symbolic procedure xdeps(eqn);
% %
% % Condition:
% %
% % dunkn D does not depend on z's.
% % Every dunkn in eqn which _does_ depends on z's doesnt depend on something D does,
% % and so it is ok to solve the eqn for D, after we make apropriate restrictions
% %
% % (at the moment, we have D depends on zb which _no_ other dunkn depends on)
% %
% % return (z's).D
% %
% begin scalar xdeps, sxdeps, dfdunkn, resteqn;
%   bvars:=boundvarsE eqn$
%   resteqn:=eqn$
%   while resteqn do begin
%     xdeps:=complement(bvars,depsk mvar resteqn)$      % 
% %          ^ what this dunkn doesnt depend on but something else does
%     if xdeps and complement(depsk mvar resteqn, otherboundvarsE(dunknin mvar resteqn,eqn))% 
% %                ^ what this dunkn depends on but nothing else does
%        and (null sxdeps) or count(sxdeps) > (count xdeps)
%          then <<sxdeps:=xdeps; dfdunkn:=mvar resteqn>>;
%     resteqn:=red resteqn$
%     end$
%   return sxdeps.dfdunkn;
%   end;


symbolic procedure xdepsshr(eqn);
%
% dunkn D does not depend on z's.
% D depends on something that nothing else does.
% and so it is ok to solve the eqn for D, after we make appropriate restrictions
%
% return list of possibilities: { {(z's), D}, ... }
%
begin scalar xdeps, sxdepslis, dunkn, dunkns, ddeps, resteqn;
  bvars:=boundvarsE eqn$
  resteqn:=eqn$
  dunkns:=dunknsE eqn;
  while resteqn do begin
    dunkn:=dunknin mvar resteqn;
    ddeps:=depsk dunkn;
    xdeps:=complement(bvars,ddeps);
    if xdeps and complement(ddeps, otherboundvarsE(dunkn,eqn))
      then sxdepslis:=list(xdeps, mvar resteqn).sxdepslis;
    resteqn:=red resteqn$
    end$
  return sxdepslis;
  end;

symbolic procedure xdepsxpd(eqn);
%
% Condition:
%
% dunkn D does not depend on z's.
% Every dunkn in eqn which _does_ depends on z's doesnt depend on something D does,
% and so it is ok to solve the eqn for D, after we make appropriate restrictions
%
% return list of possibilities: { {(z's), D}, ... }
%
begin scalar xdeps, sxdepslis, dunkn, dunkns, ddeps, resteqn;
  bvars:=boundvarsE eqn$
  resteqn:=eqn$
  dunkns:=dunknsE eqn;
  while resteqn do begin
    dunkn:=dunknin mvar resteqn;
    ddeps:=depsk dunkn;
    xdeps:=complement(bvars,ddeps);
    for each d1 in dunkns do if xdeps and not (d1=dunkn) and null complement(ddeps,depsk d1)
      then xdeps:=complement(xdeps,depsk d1); 
    if xdeps then sxdepslis:=list(xdeps, mvar resteqn).sxdepslis;
    resteqn:=red resteqn$
    end$
  return sxdepslis;
  end;

symbolic procedure nonpolyvarsE eqn;
% find all vars !*x and !*u which appear nonpolynomialy in eqn.
% Pass npvars as fluid
begin scalar npvars;
  for each red on eqn do nonpolyvarsf lc red;
  return npvars;
  end;

symbolic procedure nonpolyvarsf u;
if null u or domainp u then nil
else <<nonpolyvarsk mvar u,nonpolyvarsf lc u;nonpolyvarsf red u>>;

symbolic procedure nonpolyvarsk u;
% find all vars !*x and !*u which appear nonpolynomialy in kern,
% ie not at top level 
if pairp u and not((car u)=!*u) and not ((car u)=!*x)
then nonpolyvarsa u;

symbolic procedure nonpolyvarsa u;
if not pairp u then nil
else if (car u=!*u) or (car u=!*x) then npvars:=union(list u, npvars)
else <<nonpolyvarsk car u;nonpolyvarsk cdr u>>;

symbolic procedure sublxdeps01(lxdeps,eqn,xdunkn);
begin scalar neqn,neqn1;
  if !*traceop!*xdp then <<write "op!*xdp, subbing ",lxdeps," =0 in eqn.";terpri()>>;
  neqn:=eqn;
  for each var in lxdeps do begin
    neqn1:=mysub0f(neqn,var);
    if neqn and not (xdunkn member dunknsE neqn1) then <<
      if !*traceop!*xdp then <<write "Substituting ",var," =1 instead.";terpri()>>;
      neqn1:=mysub1f(neqn,var)>>;
    if neqn and not (xdunkn member dunknsE neqn1) then <<
      if !*traceop!*xdp then <<write "Failed to choose good value for ",var;terpri()>>;
      neqn1:=nil>>;
    neqn:=neqn1;
    end;
  neqn:=mynumr mysimpf neqn;
  if !*traceop!*xdp then <<write "op!*xdpl: after subbing, got ";terpri(); showf neqn>>;
  return neqn;
  end;

symbolic procedure mysub0f(eqn,u);
% sub u=0 in mvars of eqn
if null eqn or domainp eqn then eqn
else if (mvar eqn)=u then mysub0f(red eqn,u)
else addf(multf(!*p2f lpow eqn, mysub0f(lc eqn,u)),
          mysub0f(red eqn,u));

symbolic procedure mysub1f(eqn,u);
% sub u=1 in mvars
if null eqn or domainp eqn then eqn
else if (mvar eqn)=u then addf(mysub1f(red eqn,u), mysub1f(lc eqn,u))
else addf(multf(!*p2f lpow eqn, mysub1f(lc eqn,u)),
          mysub1f(red eqn,u));


symbolic procedure restricteqn(eqn,neqn,xdunkn,xdeps);
% restrict by neqn by xdeps,
% absorb any redundancies
begin scalar restvars, redlist, newdunkns, olddfdunkns, badnd; 
  restvars:=xdeps;               % restrict by these vars , pass as fluid
  neqn:=restricteqnf(neqn);      % restrict by restvars
  if !*traceop!*xdp then <<write "After restriction, got new equation";terpri();
                       showf neqn;
                       showdepsineqn neqn>>;
  %
  % wow, this can introduce redundancies here...
  %
  % first look for redundancies between old and new dunkns,
  % and absorb them into the new dunkns.
  % ( _not_ into the old dunkns, and _not_ by setting either to zero.)
  %
  newdunkns:=complement(dunknsE neqn,dunknsE eqn);        % newdunkns will not be dfdunkns (yet!)
  olddfdunkns:=complement(dfdunknsE neqn, newdunkns);
  redlist:=nil; % the list of redunkndant dfdunkns in neqn (between ods and nds)
  for each nd in newdunkns do for each od in olddfdunkns do % check for redundancies between them
    if not (od member redlist) and not complement(depsk od, depsk nd)
      and not complement(op!*xdpratiodeps(od,nd,neqn), depsk nd)
        then <<redlist:=od.redlist; if (dunknin od)=xdunkn then badnd:=nd>>;
  if !*traceop!*xdp and redlist then <<write "op!*xdp removing redundancies ",redlist,
                                         " with old dunkns from neqn.";terpri()>>;
  % if we can absorb xdunkn then we will have big problems, so we quit on this case
  % this could only happen if there was a dfdunkn in the original eqn which depended on _more_
  % than xdep, which is not meant to happen. So we ring the alarm bells as well...
  if badnd then begin
    terpri();
    write "op!*xdp with xdunkn ",xdunkn," and xdeps ",xdeps," on eqn ";terpri();showf eqn;
    write "We have restricted to neqn ";terpri();showf neqn;
    write "But we got xdunkn redundant to ",badnd," which should not have happened!";terpri();
    neqn:=nil; % we dont want to play with this
    end;
  for each od in redlist do neqn:=mysub0f(neqn,od);
  %
  % now look for redundancies between new dunkns
  %
  redlist:=nil; % the list of redunkndant dfdunkns in neqn (between new dunkns)
  for each nd in newdunkns do for each nd1 in newdunkns do
    if  nd neq nd1 and not(nd1 member redlist) and not complement(depsk nd1, depsk nd)
       and not complement(op!*xdpratiodeps(nd,nd1,neqn), depsk nd1)
        then redlist:=nd.redlist;
  if !*traceop!*xdp and redlist then <<write "op!*xdp removing redunancies ",redlist,
                                         " between new dunkns from neqn.";terpri()>>;
  for each nd in redlist do neqn:=mysub0f(neqn,nd); % how about clearing dependencies...
  newdunkns:=complement(newdunkns,redlist);
  %
  % Now its time to check that there are no possible redundancies of the type
  % where coeffs of terms in xdunkns involve more than 1 dunkn.
  % We can't remove these redundancies, so we quit on this.
  % (we could extend this op to relate new dunkns,
  % ie, by restricting c(1)(x,y,z) to c(2)(x,y) means restricting df(c(1),x) to df(c(2),x)
  % _and_ restricting df(c(1),z)(x,y,z) to c(3)(x,y) means restricting df(c(1),x,z) to df(c(2),x)
  % (this last bit is too complex to do for now)
  %
  redlist:=nil; % the list of unremovable redunkndant dfdunkns in neqn (between new dunkns)
  for each nd in newdunkns do for each nd1 in newdunkns do
    if  nd neq nd1 and not complement(op!*xdpratiodeps(nd,nd1,neqn), union(depsk nd,depsk nd1))
        then redlist:=nd.redlist;
  if !*traceop!*xdp and redlist
    then <<write "op!*xdp: cant remove redunancies ",redlist,
                 " between new dunkns from neqn, so quitting operation.";terpri()>>;
  if redlist then neqn:=nil;
  %
  % Hey, we still might be able to absorb some factors into nd's ... :-) not a real big issue
  % but we might also be able to absorb some nd's into eachother... 
  % and there still might be some redundancies left...
  %
  return neqn;
  end;

symbolic procedure op!*xdpratiodeps(dunkn1,dunkn2,eqn);
begin scalar x;
 x:=op!*xdpratio(dunkn1,dunkn2,eqn);
 return union(alldepsf numr x,alldepsf denr x);
 end;

symbolic procedure op!*xdpratio(od,nd,eqn);
% ratio of the coefficients of od and nd in eqn
begin scalar oc,nc,gcd;
  oc:=coeffE(od,eqn);
  nc:=coeffE(nd,eqn);
  gcd:=mygcdf!*(oc,nc);
  return quotf1(oc,gcd) ./ quotf1(nc,gcd);
  end;

symbolic procedure restricteqnf(eqn)$
% restricts eqn so that dunkns dont depend on restvars,
% by making new dunkns if neccessary
% **pass restvar as fluid
% return sf result (was sq)
if null eqn then nil
else addf( list(restricteqnk(mvar eqn) .** 1 .* lc eqn),restricteqnf(red eqn))$


symbolic procedure restricteqnk(dfdunkn)$
% **takes restvar as fluid
if intersect(restvars, depsk dfdunkn)
  then newarb1 complement(depsk dfdunkn,restvars)
  else dfdunkn$

symbolic procedure op!*xdpsolvebyops(neqn,eqngrp,xdunkn,xdeps);
%  We have got a  new eqn, now solve it for xdunkn if we can
begin scalar neqngrp, res, olddunkns$
  neqngrp:=mkeqngrp(neqn,eg!*no eqngrp,"op!*xdp");
  olddunkns:=dunknsE eg!*eqn eqngrp;                      % the dunkns in the original eqn
  res:=nil;
  if !*traceop!*xdp then <<write "op!*xdp got eqngrp"; terpri(); showeqngrp neqngrp>>;
  res:=solve1byops(neqngrp, !*op!*xdpslvops, '!*op!*xdp!*opusage);
  if null res then deathcert(neqngrp,list("created by op!*xdp then dropped and",
                             " original eqn kept as I cant solve this any further"));
  if null res and !*traceop!*xdp then <<
    write "op!*xdp failed to solve the eqn once we found it, so dropping it :-("; terpri()>>;
  if null res then return;
  if not (xdunkn member cdr res)
    and (not member(cadr res, olddunkns) or not equalset(depsk xdunkn,depsk cadr res))
      then begin                                 % solved eqn for some other dunkn ! bad move...
        write "***op!*xdp got xdunkn ",xdunkn," with xdeps ",xdeps," on eqngrp";terpri();
        showeqngrp eqngrp;
        write "restricted by ",xdeps," to get neqngrp";terpri();
        showeqngrp neqngrp;
        write "solved this eqngrp with reason ",eg!*dc neqngrp;terpri();
        write "but xdunkn ",xdunkn," is not in sdunkns ",cdr res," from solved equation neqn";
        terpri();
        write "I am worried that I might have solved for some other dunkn which is not in symvec";
        terpri();
        write "That is a potentially very dangerous thing to do...";terpri();terpri();
        end;
  if !*verifyop!*xdp then
    if simpeqn neqn % ie has a substition been made which will simplify neqn (wont work on eqn!)
      then <<write "op!*xdp failed"; terpri()>>
      else <<write "op!*xdp succeeded"; terpri()>>;
  if !*traceop!*xdp 
    then <<write "op!*xdp, equation made by op!*xdp solved with reason ",eg!*dc neqngrp;terpri()>>;
  return res;
%      addeqngrp neqngrp;
%  if neqn then return list(t)$
  end$

%******end of op!*xdp



% Here we solve eqn for a dfdunkn, take the expression for it and integrate it
% until we have the expression for the actual dunkn. Then make an assignment.
%
% We must have:
% 1. dfdunkn must be diff'd by only one var, which will be the intvar.
% 2. dunkn must depend on all vars in eqn.
% 3. dunkn must appear only once in eqn.
% 4. dfdunkn must be the (equal) lowest order derivative by intvar.
% 5. No (sq) coeff of a dfdunkn depending on intvar can involve intvar.

symbolic procedure op!*intslv eqngrp;
begin scalar eqn, reqn, neqn, alldeps, dfdunkn, intvar, intord, intcount,
             badint, funknatomswithdeps;
  reqn:=eqn:=eg!*eqn eqngrp;
  if null red eqn then return; % Dont solve single term equations here.
                               % diff by 2 vars would fail test anyway.
  alldeps:=alldepsf eqn;
  while reqn and not intvar do
    if car (dfdunkn:=mvar reqn)='df
      and (null cdddr dfdunkn or (numberp cadddr dfdunkn and null cddddr dfdunkn))
        % ie df'd by only one var
      and not complement(alldeps, depsk dfdunkn)
        % ie depends all vars in eqn
      and onlyonce(dfdunkn,eqn)
      and (dford dfdunkn)=minord(eqn, caddr dfdunkn) % do we really need this condition?
%      and not (minord(eqn, caddr dfdunkn)=0)        % or cxan we use this instead
        then intvar:=caddr dfdunkn
        else reqn:=red reqn;
  if not intvar then return;
  if null cdddr dfdunkn then intord:=1 else intord:=cadddr dfdunkn;
  neqn:=subvalf(dfdunkn,eqn); % this is the value for dfdunkn
%  if mytermsf denr neqn>!*intndenlimit or mytermsf numr neqn>!*intnnumlimit then return;
  % now we have to test to see if we can integrate this expression.
  % we need to test each coeff of dfdunkns depending on intvar - these
  % coefficients can not involve intvar. If the den doesnt involve intvar
  % then do test just on numr... this can be worked on later, for now
  % just use existing test...
  if not sqinttest(neqn,intvar) then return;
  intcount:=0; % we have integrated 0 times so far.
  funknatomswithdeps:=union(funknatomswithdepsf denr neqn,funknatomswithdepsf numr neqn);
  if !*traceop!*intslv 
    then <<write(tmstmp(),"op!*intslv: found funkns ",union(funknsf denr neqn,funknsf numr neqn),
                 " and funknatomswithdeps ",funknatomswithdeps);terpri();
           showeqngrp eqngrp>>;
  if funknatomswithdeps and not !*intfdeps then 
     <<write("op!*intslv on eqngrp ",eg!*no eqngrp,". Found freeunkns ",
             funknatomswithdeps," with implicit dependence.");terpri();
       write("Integrator will fail on this, see documentation for fix");terpri();terpri();
       return>>;
  if !*traceop!*intslv then <<write("op!*intslv: integration by ",intvar);terpri()>>;
  if (mytermsf numr neqn)>!*intnnumlimit
    then <<write "op!*intslv on eqngrp ",eg!*no eqngrp,": ",
                 "integration seems possible but I decided the numr is too long as it has ",
                 "length ", mytermsf numr neqn,"."; terpri();
           badint:=t>>;
  if mytermsf denr neqn>!*intndenlimit
    then <<write "op!*intslv on eqngrp ",eg!*no eqngrp,": ",
                 "integration seems possible but I decided the denr is too long as it has ",
                  "length ",mytermsf denr neqn,"."; terpri();
           badint:=t>>;
  while intcount neq intord and not badint do begin
    if !*traceop!*intslv then <<write("op!*intslv: integration ",intcount," got:");
                                showsq neqn;
                                write("denominator has length ",mytermsf denr neqn);terpri()>>;
    if !*traceintcalls then <<write "op!*intslv: making integration call ";terpri()>>;
    if !*tracecute then write "I",mytermsf numr neqn,"/",mytermsf denr neqn;
    if (mytermsf numr neqn)=0 then <<showeqngrp eqngrp;
                                 write "in op!*intslv: integrating ";terpri();
                                 showsq neqn>>;
    neqn:=myintsqE(neqn,intvar);
    if !*traceintcalls then <<write "out of integration call ";terpri()>>;
    intcount:=intcount+1;
    if (funkndepsp neqn) and not !*intfdeps then 
      <<write ("op!*intslv: After ",intcount," integrations, found freeunknowns in ");terpri();
       showsq neqn>>;
    if (badint:=intleft numr neqn) and not ((numr neqn) member !*hardints) then 
      <<write ("op!*intslv: After ",intcount," integrations, integrals left in ");terpri();
       showsq neqn;
       !*hardints:=(numr neqn).!*hardints>>;
    end;
  if badint then return;
  if !*traceop!*intslv then <<write(tmstmp(),"op!*intslv: integration ",intcount," got:"); showsq neqn>>;
  neqn:=addsq(neqn,!*f2q newpolyn(intvar,intord-1,delete(intvar,depsk dfdunkn)) );
  if !*traceop!*intslv then << write(tmstmp(),"op!*intslv: Expression for ",dunknin dfdunkn," is "); showsq neqn>>;
  deathcert(eqngrp,list("Solved by op!*intslv: integrated and solved for ",dunknin dfdunkn));
  mysetk(dunknin dfdunkn,mk!*sq neqn);
  if !*verifyint then
    if simpeqn eqn
      then begin
        write(tmstmp(),"**** op!*intslv: Unable to verify integration"); terpri();
        write("Integrating equation number ",eg!*no eqngrp," by ",intvar,
              " ",intcount," times to solve for ",dfdunkn); terpri();
        write("Was:"); terpri(); showf eqn;
        write("Got value for ",dunknin dfdunkn," :"); terpri(); showsq neqn;
        showdepsineqn numr neqn;
        write("Equation evaluates to"); terpri(); showf simpeqn eqn;
        end
      else <<write("op!*intslv integration succeeded..."); terpri() >>;
  return list (t,dunknin dfdunkn);
  end;

symbolic procedure onlyonce(dfdunkn,eqn);
% does dunknin dfdunkn appear only once in eqn?
null eqn 
or ( (dunknin dfdunkn) neq (dunknin mvar eqn) or mvar eqn=dfdunkn)
and onlyonce(dfdunkn,red eqn);

symbolic procedure minord(eqn, var);
% what is the lowest order derivative by var of the dunkns in eqn which
% depend on var ?
begin scalar x;integer k;
  k:=10000; % infinity- we would never tackle equations of this order!
  while eqn and not (k=0) do begin
    if (car mvar eqn='df) and (x:=member(var,mvar eqn))
      then k:=(if cdr x and numberp cadr x then min2(k,cadr x) else 1)
      else if var member depsk mvar eqn then k:=0;
    eqn:=red eqn;
    end;
  return k;
  end;

% ****** end of op!*intslv

symbolic procedure intcheck(neweqnq,oldeqn,intvar);
% Check that integration worked. Return t if it did, nil if it didn't
%
% df of new eqn should be some multiple of old eqn.
% ( but new eqn has been multiplied by some factor, this has
%   to be dropped otherwise const of int'n wont dissapear...)
% So drop factors, compare for equality.
% This is not foolproof: it is possible that two eqns can be equivalent
% with factors not removable by dropfac, eg a factor of cos... so simp!*
%
begin scalar ne;
 neweqnq:=multsq(neweqnq, (denr neweqnq) ./ (lc numr neweqnq) );
 ne:=dropeqnfacE numr diffsq(neweqnq,intvar);
 if mvar ne neq mvar oldeqn then return nil;
 return null numr mysimpf addf(multf(coefff(mvar ne,ne),oldeqn),
                               negf multf(coefff(mvar oldeqn,oldeqn),ne));
 end;

% Term by term integration.
% Each dunkn term must be integrable in its own right,
% Ie the dfdunkn is diffed wrt var and the coeff is indep of var,
% or the dfdunkn is indep of var 

symbolic procedure op!*inttbt eqngrp;
begin scalar eqn, dunkn, dfvars, ivars, val, intvar, badint, rm, den, funknatomswithdeps;
  eqn := eg!*eqn eqngrp;
  if !*intonlyslvbl
    then for each dunkn in isdunkns eqn do 
      if (dfvars:=getdfvars dunkn) and not cdr dfvars % ie only one dfvar
        then ivars:=union(ivars,dfvars) % consider integrating by each of these
    else for each rm on eqn do
      ivars:=union(ivars,getdfvars mvar rm);
  while ivars do if inttestf(eqn, car ivars) 
    then <<intvar:=car ivars; ivars:=nil>>
    else ivars:=cdr ivars;
  if null intvar then return;
  if (funknatomswithdeps:=funknatomswithdepsf eqn) and not !*intfdeps then 
     <<write(tmstmp(),"op!*inttbt: Pre integration, found freeunkns with implicit dependence ",
              funknatomswithdeps);terpri();
       write("Integrator will fail on this, see documentation for fix");terpri();
       return>>;
  if !*traceintcalls then <<write "op!*inttbt: making integration call ";terpri()>>;
  if !*tracecute then write "I",mytermsf eqn;
  val:=myintsqE(eqn ./ 1,intvar);
  if !*traceintcalls then <<write "out of integration call ";terpri()>>;
  rmsubs();
  val := simp!* val;
  den := denr val;   % keep for multiplying constant
  val := mynumr val;
  if (badint:=intleft val) and not (val member !*hardints) then begin
    terpri();
    write ("op!*inttbt: integrals left in "); terpri();
    showf val;
    !*hardints:=val.!*hardints;
    end;
  if badint then return;
  val:=addf( val,multf(!*k2f newarb1 delete(intvar, alldepsf eqn), den) );
  if !*verifyint then
    if not intcheck(val ./ 1,eqn,intvar) 
      then begin
        write(tmstmp(),"Integrating equation number ",eg!*no eqngrp," by ",intvar);
        terpri();
        write("Was:"); terpri();
        showf eqn;
        write("Got:"); terpri();

        showf val;
        showdepsineqn val;
        write("**** op!*inttbt: Unable to verify integration");
        end
      else write(tmstmp(),"op!*inttbt integration succeeded...");
    terpri();
  addeqngrp mkeqngrp(val,eg!*no eqngrp,'op!*inttbt);
  deathcert(eqngrp,list("Integrated wrt ",intvar,
            " to give new equation ",!*dets));
  return list t;
  end;

% term by term integration with division
symbolic procedure op!*inttbta eqngrp;
begin scalar eqn, neqn, dfdunkn, a, dfvars, 
             val, intvar, badint, rm, den, funknatomswithdeps;
  eqn := eg!*eqn eqngrp;
%  isdunkns:=isdunkns eqn;
  rm:=eqn;
  while rm and not intvar do begin
    dfdunkn:=mvar rm;
    a:=coefff(dfdunkn,eqn);
    if% dunknin dfdunkn member isdunkns
      %and
      (dfvars:=getdfvars dfdunkn) %and not cdr dfvars % ie only one dfvar
      and sqinttest(neqn:=multsq(eqn ./ 1, 1 ./ A),car dfvars)
        then intvar:=car dfvars
        else rm:=red rm;
    end;
  if null intvar then return;
  if (funknatomswithdeps:=union(funknatomswithdepsf numr neqn,funknatomswithdepsf denr neqn))
     and not !*intfdeps then
     <<write(tmstmp(),"op!*inttbta: Pre integration, found freeunkns with implicit dependence ",
              funknatomswithdeps);terpri();
       write("Integrator will fail on this, see documentation for fix");terpri();
       return>>;
  if !*traceintcalls then <<write "op!*inttbta: making integration call ";terpri()>>;
  if !*tracecute then write "I",mytermsf numr neqn,"/",mytermsf denr neqn;
  val:=myintsqE(neqn, intvar);
  if !*traceintcalls then <<write "out of integration call ";terpri()>>;
  den := denr val;
  val := mynumr val;
  if (badint:=intleft val) and not (val member !*hardints) then 
    <<terpri(); write ("op!*inttbta: integrals left in "); terpri(); showf val;
      !*hardints:=val.!*hardints>>;
  if badint then return;
  val:=addf( val,multf(!*k2f newarb1 delete(intvar, alldepsf eqn), den) ); 
  if !*verifyint then
    if not intcheck(val ./ 1,eqn,intvar) 
      then begin
        write(tmstmp(),"op!*inttbta: Integrating equation number ",eg!*no eqngrp,
              " by ",intvar); terpri();
        write("Was:"); terpri(); showf eqn;
        write("Got:"); terpri(); showf val;
        showdepsineqn val;
        write("**** op!*inttbta: Unable to verify integration");
        end
      else write("op!*inttbta integration succeeded...");
      terpri();
  addeqngrp mkeqngrp(val,eg!*no eqngrp,'op!*inttbta);
  deathcert(eqngrp,list("op!*inttbta: Integrated wrt ",intvar," to give new equation ",!*dets));
  return list t;
  end;

symbolic procedure isdunkns eqn;
% Find all dunkns which are differentiated, appear only once,
% and which depend on (at least) alldepsf eqn and boundvarsE eqn.
% Ie, find all dfdunkns which if the eqn is integrated by their derivatives,
% then the equation could be solved for.
begin scalar resteqn, dunkn, dunknsonce, dunknsmore, isdunkns, alldeps;
  alldeps:=alldepsf eqn;
  resteqn:=eqn;
  while resteqn do begin % find all dunkns appearing only once
    dunkn:=dunknin mvar resteqn;
    if dunkn member dunknsonce then <<
      dunknsmore:=dunkn.dunknsmore;
      dunknsonce:=delete(dunkn,dunknsonce)>>
    else if not (dunkn member dunknsmore) then dunknsonce:=dunkn.dunknsonce;
    resteqn:=cdr resteqn;
    end;
  for each dunkn in dunknsonce do
    if not complement(alldeps,depsk dunkn) then isdunkns:=dunkn.isdunkns;
  return isdunkns;
  end;

symbolic procedure getdfvars kern;
% kern is a df of a dunkn kernel.
% find all vars that the dunkn is diffd by.
begin scalar vars;
    kern:=cddr kern; % ie the derivative list
    while kern do begin
        vars:=union(list car kern,vars);
        if cdr kern and fixp cadr kern then kern:=cddr kern else kern:=cdr kern;
        end;
    return vars;
    end;
  
%****** just some general integration stuff

fluid '(difnum failed);

symbolic procedure dfvarlist u;
% u is a dfdunkn kernel in the case of an edunkn (explicit determinable unknown function)
% or 'df . edunkn . dfvarlist 
% return the dfvarlist
if (car u)='df then cddr u;

symbolic procedure minusdfdunkns(u,v);
% u and v are dfdunkns,
% return the diference of their dfvarlist if u is a derivative of v
if dunknin(u)=dunknin v then minusdfvarlist(dfvarlist u,dfvarlist v) else -1;

symbolic procedure minusdfvarlist(u,v);
% u and v are the (ordered!) dfvarlists from dfdunkns variables
% eg u=(x 2,2,x 3,x 4,u(1,1)) from ,say, uu=df(phi 3,x 2,2,x 3,x 4,u(1,1))
% if uu is a deriv of vv then return the difference y-z
% else return -1 
% use reduce ordering of dfvarlists: higher args are higher orderd,
% ie (simplified) df(anything,a,b) satisfies ordop(a,b)='t
if null v then u
else if null u then -1
else if (car u)=(car v)
  then if pairp cdr u and numberp cadr u
    then if pairp cdr v and numberp cadr v 
      then if (cadr u)=(cadr v) then minusdfvarlist(cddr u,cddr v)
      else if (cadr u)-(cadr v)=1 then 
        ((lambda z; if z=-1 then z else (car u).z)
          minusdfvarlist(cddr u,cddr v))
      else if (cadr u)-(cadr v)>0 then
        ((lambda z; if z=-1 then z else (car u).((cadr u)-(cadr v)).z)
         minusdfvarlist(cddr u,cddr v))
      else -1 % ie (cadr u)-(cadr v)<0
      % end of both nums
    else if (cadr u)=2 then ((lambda z; if z=-1 then z else (car u).z)
                              minusdfvarlist(cddr u,cdr v))
    else ((lambda z; if z=-1 then z else (car u).((cadr u)-1).z)
           minusdfvarlist(cddr u,cdr v))
    % end of u num
  else if pairp cdr v and numberp cadr v then -1    % v with num, u w/o num
  else minusdfvarlist(cdr u,cdr v)                % u&v same, w/o num
else if not ordop(car u,car v) then -1 % ok:-)
else if pairp cdr u and numberp cadr u then 
  ((lambda z; if z=-1 then z else (car u).(cadr u).z)
    minusdfvarlist(cddr u,v) )
else ((lambda z; if z=-1 then z else (car u).z)
       minusdfvarlist(cdr u,v) );


symbolic procedure minusdfvarlisttwo(u,v);
% u and v are the (ordered!) dfvarlists from dfdunkns variables
% eg u=(x 2,2,x 3,x 4,u(1,1)) from ,say, uu=df(phi 3,x 2,2,x 3,x 4,u(1,1))
%
% Subtract any common derivatives from the two dfvarlists
% and return them as a list of two dfvarlists.
%
% Use reduce ordering of dfvarlists: higher args are higher orderd,
% ie (simplified) df(anything,a,b) satisfies ordop(a,b)=t
%
if null v or null u then list(u,v)
else begin scalar uu,vv,u2do,v2do,var;
  u2do:=uu:=u;
  v2do:=vv:=v;
  while u2do and v2do do
    if (car u2do)=(car v2do)
      then <<var:=car u2do;
             u2do:=deldfvar1(var,u2do);
             v2do:=deldfvar1(var,v2do);
             uu:=deldfvar1(var,uu);
             vv:=deldfvar1(var,vv)>>
    else if ordop(car u2do,car v2do)
      then u2do:=deldfvar1(car u2do,u2do)
    else v2do:=deldfvar1(car v2do,v2do);
  return list(uu,vv);
  end; 

symbolic procedure intleft u;
%% are there any explicit integrals left?
%if pairp u then if car u='int then t
%                else (intleft car u) or (intleft cdr u);
if pairp u 
  then if car u='int 
    then if not atom cadr u and (not atom nm and dunknk mvar nm where nm=numr simp!* cadr u)
      then t                  % integral of expresion with dunkns in it
      else null !*keepints    % integral of expresion w/o dunkns in it
    else (intleft car u) or (intleft cdr u);

symbolic procedure inttestf(u,intvar);
% is standard form u integrable wrt intvar?
null u or inttestt(lt u,intvar) and inttestf(red u,intvar);

symbolic procedure inttestt(ltm,intvar);
not (intvar member depsk caar ltm)
or (intvar member getdfvars caar ltm) and not (intvar member alldepsf cdr ltm);

symbolic procedure myintf(u,intvar);
% integrate standard form u by variable intvar
% return standard quotient.
if null u then (nil ./ 1)
%else simp!* list('int, mk!*sq , intvar); 
else simp!* list('int, prepsq (u ./ 1), intvar); 

symbolic procedure myintsq(u,intvar);
% integrate sq u by variable intvar
% return standard quotient.
if null numr u then (nil ./ 1)
%else simp!* list('int, mk!*sq u, intvar); 
else simp!* list('int, prepsq u, intvar); 

symbolic procedure myintsqE(sqeqn,intvar);
% integrate sq-eqn by variable intvar
% return standard quotient.
myintperv(sq2perv sqeqn,intvar);
%myintsq(sqeqn,intvar);

symbolic procedure myintperv(perveqn,intvar);
% integrate sq-eqn by variable intvar
% return standard quotient.
% perveqn is list of list(dfdunkn, sq coeff)
%
% we can optimise code by stoping integration when it fails on one part...
%
% car perveqn is the first term
% caar perveqn is the first dfdunkn
% cadar perveqn is the first dfdunkn's sq coefficient
% cdr is the rest
%
       % We are relying heavily on the ordering of dfdunkns being higher than any
       % mvar in the coefficient, as has been set by funorder. Saves expensive multsq
%
if null perveqn then (nil ./ 1)
else if (intvar member depsk caar perveqn)
  and (intvar member union(alldepsf numr cadar perveqn, alldepsf denr cadar perveqn))
%  then rederr "badint" % set up something to handle this
  then << %if !*trace? then write "bad integral in ",caar perveqn," term.";
          !*k2q list('int, mk!*sq perv2sq perveqn, intvar)>>
else if intvar member depsk caar perveqn
       % Integrate the dfdunkn part
  then addsq( list( (deldfvar(intvar,caar perveqn) .** 1)
                    .* numr cadar perveqn
                   )
               ./ denr cadar perveqn ,
              myintperv(cdr perveqn,intvar)
             )
       % integrate the coeff part
  else addsq( ( list( ((caar perveqn) .** 1)
                       .* numr intres
                     ) 
                 ./ denr intres 
                where intres=myintsq(cadar perveqn, intvar)
               ),
              myintperv(cdr perveqn,intvar)
             );


symbolic procedure deldfvar(var,dfdunkn);
% dfdunkn has at least one occurance of var in its deriv list.
% remove one occurance of var
% car dfdunkn = 'df
% cadr dfdunkn = the unknown
% cddr dfdunkn = the deriv list.
if not(var member cddr dfdunkn) 
  then rederr list("cant remove var ",var," from derivative list ",cddr dfdunkn)
  else if null cdddr dfdunkn then cadr dfdunkn
  else (car dfdunkn) . (cadr dfdunkn) . deldfvar1(var,cddr dfdunkn);

symbolic procedure deldfvar1(var,u);
if car u=var 
  then if cdr u 
    then if numberp cadr u 
      then if cadr u = 2 
        then (car u).cddr u
        else (car u).(cadr u-1).cddr u
      else cdr u
    else nil
  else if numberp cadr u % implies u MUST have a cdr 
    then (car u).(cadr u).deldfvar1(var,cddr u) % implies u MUST have a cddr
    else (car u).deldfvar1(var,cdr u);

%******end of integration stuff

%******end of op!*inttbt

% first order integrating factor method

symbolic procedure op!*intfac eqngrp;
%
% if we have an eqn
%
% A(x)*y' + B(x)*y + C(x) = 0
%
% then we divide by A(x) to put in standard form,
%
% ie   y' + A(x)/B(x)*y + C(x)/A(x) = 0
%
% This than has Integrating factor I=exp int(B/A,x)
%
% So we have (Iy)' + I C(x)/A(x) = 0
%
% This integrates as Iy + int(I C(x)/A(x),x) + const = 0
%
% We don't do this if A(x) is too long, as determined by !*intndenlimit
%
% We have to find a dunkn phi in eqn, appearing twice:
% once with dfvarlist Z, once with dfvarlist (Z,z), 
% so that df(phi,(Z,z)) serves as y' and df(phi,Z) serves as y.
% (Z is a list of variables, z is a single variable).
% We then insist that the equation is solvable for phi after this and possibly
% some further term-by-term integrations: So Z must only involve z's,
% and phi must depend on at least as much as the rest of the eqn.
begin scalar eqn, teqn, dfdunkns, dfdunkn, ord1dunkn, ord0dunkn, intvar, 
             A, B, C,intfac, intfac1, dfvarlist,
             cdeps, alldeps, neweqn, badint, funknatomswithdeps,
             intfaceqngrp, res;
  teqn:=eqn:=eg!*eqn eqngrp;
  dfdunkns:=dfdunknsE eqn;
  alldeps:=alldepsf eqn;
  while teqn and null intvar do begin 
    % is (mvar teqn) the ord1dunkn? 
    dfdunkn:=mvar teqn;
    dfvarlist:=dfvarlist dfdunkn;
    if dfvarlist
      and (null cdr dfvarlist or 
           (numberp cadr dfvarlist and null cddr dfvarlist))
        % ie dfdunkn is just df'd wrt one var
      and ( (null cdr dfvarlist and member(dunknin dfdunkn,dfdunkns))
           or (red teqn and (mvar red teqn=(deldfvar(caddr dfdunkn,dfdunkn)))))
        % ie does the ord0dunkn appear in eqn?
        % uses fact that if ord0dunkn is a df, then it has to be next dfdunkn...
      and not complement(alldeps,depsk dfdunkn)
      and ord1op!*intfactest1(eqn,dfdunkn,deldfvar(caddr dfdunkn,dfdunkn)
                          ,caddr dfdunkn)
        then <<ord1dunkn:=dfdunkn;
              ord0dunkn:=deldfvar(caddr dfdunkn,dfdunkn);
              intvar:=caddr dfdunkn>>
        else teqn:=red teqn;
    end;
  if not intvar then return;
  if !*op!*intfacfirstorderonly and car ord0dunkn='df then return;
  A:=coeffE(ord1dunkn,eqn);                                          % A(x)
  if (mytermsf A)>!*intndenlimit then return;
  B:=coeffE(ord0dunkn,eqn);                                          % B(x)
  if null A or null B then rederr("Bad op!*intfac");
  C:=addf(addf(eqn,negf !*t2f ((ord0dunkn .** 1) .* B)),
          negf !*t2f ((ord1dunkn .** 1) .* A));
  if (funknatomswithdeps:=union(funknatomswithdepsf A,funknatomswithdepsf B))
     and not !*intfdeps then 
     <<write(tmstmp(),"op!*intfac: Pre integration, found freeunkns with implicit dependence ",
              funknatomswithdeps);terpri();
       write("Integrator will fail on this, see documentation for fix");terpri();
       return>>;
  if (mytermsf B)>!*intnnumlimit then return;  
  if !*traceintcalls
    then <<write "op!*intfac: making integration call to make intfac";terpri()>>;
  if !*tracecute then write "I",mytermsf B,"/",mytermsf A;
  intfac:=myintsq(multsq(B ./ 1,1 ./ A), intvar);                      % I
  if !*traceintcalls then <<write "out of integration call ";terpri()>>;
  if (badint:=intleft intfac) and not ((numr intfac) member !*hardints) then
    <<terpri();
      write ("op!*intfac part a: integrals left in "); terpri();
      showsq intfac;
      !*hardints:=(numr intfac).!*hardints>>;
  if badint then return;
  if (mytermsf numr intfac)>!*intfacnumlimit then <<
    write "In op!*intfac, got integrating factor with ",mytermsf numr intfac,
          " terms in the numerator, which is too long:";terpri();
    showsq intfac;
    return>>;
  if (mytermsf denr intfac)>!*intfacdenlimit then <<
    write "In op!*intfac, got integrating factor with ",mytermsf denr intfac,
          " terms in the denominator, which is too long:";terpri();
    showsq intfac;
    return>>;
  if null !*intfaconfree and funknsp intfac then <<
    write "In op!*intfac, got integrating factor with freeunknowns,",
          " which we aren't using just now:";terpri();
    showsq intfac;
    return>>;
%  intfac:=list('exp,mk!*sq intfac);
  intfac:=list('exp,prepsq intfac);
  intfac:=simp!* intfac;
  intfac1:=multsq( 1 ./ A,intfac);
  neweqn:=multsq(intfac1,C ./ 1);                               % I/A(x) * C(x)
  if not sqinttest(neweqn,intvar) and (!*traceintcalls or !*tracecute) then
    <<terpri();
      write(tmstmp(),"Halfway through integrating factor method on eqngrp ",
            eg!*no eqngrp," but unable to integrate by ",intvar); terpri();
      return>>; 
  if (funknatomswithdeps:=funknatomswithdepsf C) and not !*intfdeps then 
     <<write("op!*intfac: Pre integration, found freeunkns with implicit dependence ",
              funknatomswithdeps);terpri();
       write("Integrator will fail on this, see documentation for fix");terpri();
       return>>;
  if (mytermsf numr neweqn)>!*intnnumlimit or (mytermsf denr neweqn)>!*intndenlimit then return;  
  if !*traceintcalls then <<write "op!*intfac: making integration call for val";terpri()>>;
  if !*tracecute then write "I",mytermsf numr neweqn,"/",mytermsf denr neweqn;
  neweqn:=mysimpq myintsqE(neweqn, intvar);                      % int( I*C(x)/A(x) , x)
  if !*traceintcalls then <<write "out of integration call ";terpri()>>;
  if (badint:=intleft neweqn) and not (neweqn member !*hardints) then
    <<terpri();
      write ("op!*intfac part b: integrals left in "); terpri();
      showsq neweqn;
      !*hardints:=neweqn.!*hardints>>;
  if badint then return;
  cdeps:=delete(intvar,alldepsf eqn);
  neweqn:=addsq(multsq(!*k2q ord0dunkn,intfac),       % Iy + int(I C(x)/A(x),x)
                addsq(!*k2q newarb1 cdeps,neweqn));   % + const
  neweqn:=mysimpq neweqn;
  if !*verifyint then
    if not intcheck(neweqn,eqn,intvar) 
      then begin
        write("Integrating factor method on equation number ",
              eg!*no eqngrp," by ",intvar); terpri();
        write("Was:"); terpri(); showf eqn;  terpri();
        write("A(x)="); showf A; terpri();
        write("B(x)="); showf B; terpri();
        write("C(x)="); showf C; terpri();
        write("Integrating factor is"); showsq intfac; terpri();
        write("Got neweqn"); terpri(); showsq neweqn;
        write("df of neqn is"); terpri();
        showf dropeqnfacE numr diffsq(neweqn,intvar);
        showdepsineqn numr neweqn;
        write("**** op!*intfac: Unable to verify integration"); terpri();
        end
      else <<write("op!*intfac (primary) integration succeeded ..."); terpri()>>;
  intfaceqngrp := mkeqngrp(mynumr neweqn,eg!*no eqngrp,'op!*intfac);
%  if (car ord0dunkn)='df 
%    then res:=op!*intslv intfaceqngrp
%    else res:=op!*slvall intfaceqngrp;
  res:=nil;
  % now we will try lots of ops on this, importantly we solve single term eqns now
  res:=solve1byops(intfaceqngrp, !*op!*intfacslvops, '!*intfac!*opusage);
  if null res then deathcert(intfaceqngrp,list("created by op!*intfac then dropped and",
                             " original eqn kept as I cant solve this any further"));
  if null res then return; % ie it couldnt be integrated any further
%  if !*tracecute then <<
%    write "eqn made by op!*intfac was solved with reason ",eg!*dc intfaceqngrp;terpri()>>;
  if !*traceintcalls then <<
    write "We solved this eqn with reason ",eg!*dc intfaceqngrp;terpri()>>;
  if not (dunknin ord0dunkn member cdr res) 
      and (cdr res and not equalset(depsk cadr res,depsk ord0dunkn) ) then begin
    terpri();terpri();
    write("*** Integrating factor method on equation number ",
          eg!*no eqngrp," by ",intvar); terpri();
    write("Was:"); terpri(); showf eqn;  terpri();
    showdepsineqn eqn;
    write("A(x)="); showf A; terpri();
    write("B(x)="); showf B; terpri();
    write("C(x)="); showf C; terpri();
    write("Integrating factor is"); showsq intfac; terpri();
    write("Got neweqn"); terpri(); showsq neweqn;
    write "But when we solved this eqn, it was solved with reason ",eg!*dc intfaceqngrp;terpri();
    write "I expected it to be solved for ",dunknin ord0dunkn;terpri();
    write "This is potentially dangerous...";terpri();terpri();
    end;
  deathcert(eqngrp,list("Integrated wrt ",intvar,
              " by op!*intfac using the first order integrating factor method on the ",
               ord0dunkn," term to give new equation ",!*dets));
%  return list(t);
  if !*verifyint then
    if simpeqn eqn 
      then begin
        write("Integrating factor method on equation number ",
              eg!*no eqngrp," by ",intvar); terpri();
        write("Was:"); terpri(); showf eqn;  terpri();
        write("A(x)="); showf A; terpri();
        write("B(x)="); showf B; terpri();
        write("C(x)="); showf C; terpri();
        write("Integrating factor is"); showsq intfac; terpri();
        write("Got neweqn"); terpri(); showsq neweqn;
        write("df of neqn is"); terpri();
        showf dropeqnfacE numr diffsq(neweqn,intvar);
        showdepsineqn numr neweqn;
        write("**** op!*intfac: Unable to verify integration after solving associated equation."
              ); terpri();
        end
      else <<write("op!*intfac integration (secondary assignment) succeeded ..."); terpri()>>;
  return res; % ie the result of solving intfaceqngrp
  end;

symbolic procedure ord1op!*intfactest1(eqn,ord1dunkn,ord0dunkn,intvar);
%
% Can we solve eqn using the first order integrating factor method
% with ord0dunkn as the function and intvar as the differentiation variable?
% Ie, we want ord0dunkn to appear differentiated by intvar once,
% and to appear undifferentiated by intvar,
% and every other term to integrable by intvar.
%
% For generality, we may allow ord0dunkn to be a dunkn df'd by other vars.
%
null eqn
or (mvar eqn=ord1dunkn or mvar eqn=ord0dunkn 
      or (intvar member getdfvars mvar eqn) 
      or not (intvar member depsk mvar eqn)
   ) and ord1op!*intfactest1(red eqn,ord1dunkn,ord0dunkn,intvar);


symbolic procedure dfdunkns eqn;
if null eqn then nil
else if car mvar eqn='df then union(list mvar eqn,dfdunkns red eqn)
else dfdunkns red eqn;

symbolic procedure ord1varsk dfdunkn;
% find all variables that the dunkn in the kernel dfdunkn is
% differentiated by exactly once. 
if car dfdunkn neq 'df then nil
else begin scalar varlist,ord1vars;
  varlist:=cddr dfdunkn;
  ord1vars:=nil;
  while varlist do 
    if null cdr varlist then <<ord1vars:=(car varlist).ord1vars; varlist:=nil>> 
    else if not numberp cadr varlist 
      then <<ord1vars:=(car varlist).ord1vars; varlist:=cdr varlist>>
    else varlist:=cddr varlist;
  return ord1vars;
  end;
  
symbolic procedure delvar(u,v);
% u is a df var, v is a dfdunkn containing u.
% remove u from dflist of v.
% If v is ONLY dfd by u then it is nolonger a df dunkn...
if null cdddr v
  then if caddr v=u 
    then dunknin v 
    else rederr list("derivative ",v," doesnt contain ",u)
  %else delete(u,v);
  else (car v).(cadr v).deldfvarlis(u,cddr v);

symbolic procedure deldfvarlis(u,v);
if null v then rederr list("derivative ",v," doesnt contain ",u)
else if car v=u
  then if (cdr v and numberp cadr v)
    then if cadr v=1 then cddr v else (cadr v).(cadr v - 1).cddr v
    else cdr v
  else (car v).deldfvarlis(u,cdr v);

symbolic procedure ord1op!*intfactest(eqn,ord0dunkn,ord1var);
%
% Can we solve eqn using the first order integrating factor method
% with ord0dunkn as the function and ord1var as the differentiation variable?
% Ie, we want ord0dunkn to appear differentiated by ord1var once,
% and to appear undifferentiated by ord1var,
% and every other term to integrable by ord1var.
% We dont actually have to look for the zero order term because w/o it, the
% equation would be term by term integrable by ord1var.
%
% For generality, we may allow ord0dunkn to be a dunkn df'd by other vars.
%
null eqn
or (mvar eqn=ord0dunkn or (ord1var member getdfvars mvar eqn) or 
      not (ord1var member depsk mvar eqn)
   ) and ord1op!*intfactest(red eqn,ord0dunkn,ord1var);

symbolic procedure sqinttest(u,intvar);
% is standard quotient u integrable wrt intvar?
null numr u
or sqtinttest(lt numr u,denr u,intvar)
   and sqinttest((red numr u)./denr u,intvar);


symbolic procedure sqtinttest(ltm,d,intvar);
not (intvar member depsk caar ltm)
or (intvar member getdfvars caar ltm)
   and not (intvar member alldepsf multsq((cdr ltm)./1,1 ./d));

%******end of op!*intfac

symbolic procedure op!*trgexp eqngrp;
% solve second order const coeff type eqns- ie giving trig or exp solns
%
% A*f''(x) + B*f(x) + C = 0
% 
% f:=arbtrig( sqrt(b/a)*x, restdeps) - C/B
%
begin scalar eqn, trigdunkn2, trigdunkn, ord2co, phase, sign,
      val, trigdeps, trigvar, trigarg, simpsign;
  eqn:=eg!*eqn eqngrp;
  if !*traceop!*trgexp then << write tmstmp(),"in op!*trgexp on "; terpri(); showeqngrp eqngrp>>;
  trigdunkn2:=trigdunkn eqn; % find a dunkn that is diff'd twice, and which
                         % depends on more than rest of eqn.
  if not trigdunkn2 then return;
  if !*traceop!*trgexp then <<write(tmstmp(),"got trigdunkn ",trigdunkn2); terpri() >>;
  trigdunkn:=dunknin trigdunkn2;
  trigvar:=caddr trigdunkn2;
  trigdeps:=delete(trigvar,depsk trigdunkn);
  ord2co:=coeffE(trigdunkn2,eqn);                    % = A
  phase:=coeffE(trigdunkn,eqn);                      % = B
  val:=addf(addf(negf eqn,                            
                 !*t2f(trigdunkn2 .** 1 .* ord2co)),  % = -C
            !*t2f(trigdunkn .** 1 .* phase));
  val:=multsq( val ./ 1, 1 ./ phase);                 % = -C/B
  phase:=multsq(phase ./ 1, 1 ./ ord2co);             % = B/A
    %
    % Now it's time to take the sqrt... decisions about sign???
    %
    %
  if numberp numr phase and numberp denr phase
    then if (mynumr phase>0) = (denr phase>0) 
      then sign:=1 else sign:=-1
    else if (simpsign:=simp!* list('sign,prepsq phase))=(1 ./ 1) or simpsign=(-1 ./ 1)
      then sign:=numr simpsign;
%  if sign= -1 then return; % should use exp functions here...
  if !*traceop!*trgexp and (simpsign=(1 ./ 1) or simpsign=(-1 ./ 1))
    then <<write "Using SIGN of ";showsq phase; write " = ",sign; terpri()>>;
  if sign= -1 then phase:=negsq phase;
  phase:=mk!*sq phase;
  if null sign then <<terpri();
    algebraic write("***In op!*trgexp, can't determine sign of ",phase);
    write "Try using a LET rule for SIGN of that expression.";terpri();
    terpri();
    return   >>;
  phase:=simp!* list('sqrt,phase);                    % = sqrt B/A
  trigarg:=mk!*sq multsq(phase,!*k2q trigvar);
  val:=mk!*sq mysimpq addsq(val,
      if sign=1 then sqarbtrig (trigarg,trigdeps)
      else if sign=-1 then sqarbexp (trigarg,trigdeps) 
    );
  if !*traceop!*trgexp then begin
    showeqngrp eqngrp;
    algebraic write("op!*trgexp: setting ",trigdunkn," to arbtrig/exp in ",trigarg);
    terpri();
    write(trigdunkn," = "); terpri();
    algebraic write val;
    end;
  deathcert(eqngrp,list("Solved by op!*trgexp: ",trigdunkn,
                        " set to arbtrig/exp in ",trigarg));
  mysetk(trigdunkn,val);
  if !*verifyint 
    then if numr mysimpf eqn 
      then rederr("Bad op!*trgexp")
      else write ("op!*trgexp integration succeded");
%  solvesuccess:=t;
  return list(t,trigdunkn);
  end;

  
symbolic procedure trigdunkn eqn;
% find dunkn in eqn which appears only diff'd twice by some var and undiff'd
% by that var,
% with alldeps and otherboundvars less than depsk trigdunkn w/o var.
% return the dfdunkn
% caddr mvar resteqn=trigvar
begin scalar resteqn,tdunkn,alldeps;
  alldeps:=alldepsf eqn;
  resteqn:=eqn;
  while resteqn and not tdunkn do begin
    if (car mvar resteqn)='df and diff2dunkn mvar resteqn
      and not complement(alldeps,depsk mvar resteqn)
      and indf0(dunknin mvar resteqn,eqn)
      and jindf02(dunknin mvar resteqn,caddr mvar resteqn,eqn)
      and restnointvar(eqn,dunknin mvar resteqn,caddr mvar resteqn)
      and intvarnocoeff(eqn,caddr mvar resteqn)
      and not ((caddr mvar resteqn) member otherboundvarsE(dunknin mvar resteqn,eqn))
        then tdunkn:=mvar resteqn
        else resteqn:=cdr resteqn;
      end;
  return tdunkn;
  end;

symbolic procedure diff2dunkn u;
% u is a df of a dunkn
% is u a dunkn differentiated exactly twice by some var?
% car u='df; cadr u=dunkn; caddr u=var; cadddr u=2
cdddr u and numberp cadddr u and (cadddr u)=2 and null cddddr u;

symbolic procedure indf0(dunkn,eqn);
% does dunkn appear in eqn undifferentiated?
eqn and ((mvar eqn)=dunkn or indf0(dunkn,red eqn));

symbolic procedure jindf02(dunkn,dfvar,eqn);
% does dunkn appear in eqn 
% only as itself undifferentiated
% or differentiated by dfvar twice??
null eqn 
or ( (
     (mvar eqn)=dunkn                                      % as itsef
     or (dunknin mvar eqn) neq dunkn                       % not it
     or (diff2dunkn mvar eqn and (caddr mvar eqn)=dfvar)   % df'd by dfvar twice
     )
     and jindf02(dunkn,dfvar,red eqn) 
   );

symbolic procedure restnointvar(eqn,dunkn,intvar);
% is dunkn the only dunkn in eqn involving intvar?
% should we check coefficients here too? done in intvarnocoeff.
null eqn or ( (not member(intvar,mvar eqn) or dunknin mvar eqn=dunkn)
              and restnointvar(red eqn,dunkn,intvar) );

symbolic procedure intvarnocoeff(eqn,intvar);
% ensure that intvar does not appear in any coefficients.
null eqn or ( not member(intvar,alldepsf lc eqn) 
              and intvarnocoeff(red eqn,intvar) );

symbolic procedure sqarbtrig(u,deps);
% make a standard quotient linear combination of sin and cos of u
% with coefficient which are arbitrary functions with dependence list deps.
mysimpq (
  addf( multf(!*k2f newarb1 deps,!*k2f list('sin,u)),
        multf(!*k2f newarb1 deps,!*k2f list('cos,u))
      ) ./ 1);

symbolic procedure arbtrig u;
begin scalar x;
  return mk!*sq sqarbtrig(car u,for each x in cdr u collect !*a2k x);
  end;
 
symbolic operator arbtrig;

symbolic procedure arbtrig1 u;
begin scalar x;
  return mk!*sq  addsq(sqarbtrig(car u,for each x in cdr u collect !*a2k x),
                       !*k2q newarb1 (for each x in cdr u collect !*a2k x)
                       );
  end;

symbolic operator arbtrig1;

rlistat '(arbtrig arbtrig1);

symbolic procedure sqarbexp(u,deps);
% make a standard quotient linear combination of exp u and exp -u
% with coefficient which are arbitrary functions with dependence list deps.
mysimpq (
  addf( multf(!*k2f newarb1 deps,!*k2f list('exp,u)),
        multf(!*k2f newarb1 deps,!*k2f list('exp,list('minus,u)))
      ) ./ 1);

symbolic procedure arbexp u;
begin scalar x;
  return mk!*sq sqarbexp(car u,for each x in cdr u collect !*a2k x);
  end;
 
symbolic operator arbexp;

symbolic procedure arbexp1 u;
begin scalar x;
  return mk!*sq  addsq(sqarbexp(car u,for each x in cdr u collect !*a2k x),
                       !*k2q newarb1 (for each x in cdr u collect !*a2k x)
                       );
  end;

symbolic operator arbexp1;

rlistat '(arbexp arbexp1);

%******end of op!*trgexp

symbolic procedure op!*rmred1;
% look for redundant dunkns and eqns.
%
% Two dunkns are related if they appear in the same eqn,
% or if they are related to another dunkn.
% Form a list of sets of related dunkns.
%
% Any dunkn not related to one in symvec must be redundant, so set it to 0
%
% dunkns appearing in eqns, but not related to one in symvec,
% are type a redundant. (computationally expensive to have around)
% dunkns with known dependencies (ie on depl!*) but not related to one in symvec
% are type b redundant. (only annoying to have around in showing dependencies)
%
% These arise (only?) as constants of integration in eqns which are 
% substituted out when including integrability conditions.
%
% Ie, we integrate an eqn and then differentiate it... Silly, hey?
%
% Type c redundancies are explained, and found, in op!*rmred2.
%
if symvec then begin 
  scalar dunknsineqns, dunknsinsymvec, dunknlinks, link, biglink,var,simpdunkns;
  dunknsineqns:=dunknsineqns();
  dunknsinsymvec:=dunknsinsymvec();
    % we are only interested in dunkns linked to those in symvec
    % so set all others to 0
  for each eqngrp in alleqngrps() do
    dunknlinks:=addlinks(dunknsE eg!*eqn eqngrp,dunknlinks);
  if !*traceop!*rmred1 then begin
    write(tmstmp(),"op!*rmred1,found dunkns used in symvec are ",dunknsinsymvec); terpri();
    write("op!*rmred1,found dunknlinks are ",dunknlinks); terpri();
    end;
  for each link in dunknlinks do if not intersect(dunknsinsymvec,link)
    then begin
      write("op!*rmred1: found redundant functions type a ",link,
            " and setting them to 0"); terpri();
      for each dunkn in link do <<
        simpdunkns:=dunkn.simpdunkns;
        mysetk(dunkn,0)>>;
      tagsimpeqns(alleqngrps(), link);
      solvesuccess:=t;
      end
    else biglink:=link.biglink;
  for each fn in depl!* do if pairp fn 
    and edunknk fn and not member(fn,dunknsinsymvec)
    and not member(fn,biglink) then 
      if null assoc(fn,get(car fn,'kvalue)) 
        then begin 
            % not set by user...???
          terpri();
          write("found redundant function type b ",car fn,
                " and clearing it.");terpri();
          mycleark fn;
          end
        else for each var in cdr assoc(fn,depl!*) do begin
          terpri();
          write("found redundant function type b ",car fn,
                " with value, and deleting its dependence on depl!* (???).");terpri();
          depend1(fn,var,nil);
          end;
       terpri();
  end;

symbolic procedure addlinks(u,v);
% V is a list of lists of dunkns linked by eqns,
% each list in V is a set of dunkns linked by eqns.
% u is a list of dunkns which should be linked. Add them into V
begin scalar link,newlink,newv;
  newv:=v;
  for each link in v do if intersect(u,link) then begin
    newlink:=union(newlink,link);
    newv:=delete(link,newv);
    end;
  newv:=union(newlink,u).newv;
  return newv;
  end;

%******end of op!*rmred1
 
symbolic procedure sorteddepsk kern;
sort(depsk kern,function ordop);

symbolic procedure op!*rmred2;
%
% This procedure attempts to remove type c redundancies
%
% This is when there are redundant solutions, ie there appear 2 arbitrary
% functions in the solution with the same dependencies and the same
% coefficients.
%
% for each pair of dunkns in symvec, pull out the symmetry vector corresponding
% to those dunkns, ie sum of dfdunkn*coeff*vector
% Compare each term for a common ratio.
%
if symvec then begin % #1
  scalar nsvec, dunkn1, dunkn2, dunknsinsymvec, deps1, vec1, vec2, reds;
  dunknsinsymvec:=dunknsinsymvec();
  nsvec:=numr simp!* symvec;
  for each dunkn1 in dunknsinsymvec do if not (dunkn1 member reds) then begin % #2
    deps1:=sorteddepsk dunkn1;
    vec1:=nil;
    for each dunkn2 in dunknsinsymvec do if dunkn1 neq dunkn2 
      and not (dunkn2 member reds) and deps1=sorteddepsk dunkn2 then begin %#3
        if null vec1 then vec1:=dtermf(dunkn1,nsvec);
          % only calculate if needed, but only calculate once.
        vec2:=dtermf(dunkn2,nsvec);
        if vectorequivf(vec1,vec2) then begin % #4
          terpri();
          if !*traceop!*rmred2 then algebraic write symvec;
          write("op!*rmred2: Found redundant functions ",dunkn1," and ",dunkn2,
                " in symvec."); terpri();
          reds:=dunkn2 . reds;
          end; % #4
        end; % #3
    end; % #2
    if reds then begin % #6
      write("found redundant functions type c ",reds,
            " and setting them to 0"); terpri();
      for each dunkn1 in reds do mysetk(dunkn1,0);
      tagsimpeqns(alleqngrps(), reds);
      solvesuccess:=t;
      end; % #6
  end; % #1

symbolic procedure vectorequivf(vec1,vec2);
begin scalar ratio;
  if !*traceop!*rmred2 then begin
    write(tmstmp(),"testing for equivalence between"); terpri();
    showf vec1;
    write("and"); terpri();
    showf vec2;
    end;
  return vectorequivf1(vec1,vec2);
  end;

symbolic procedure vectorequivf1(vec1,vec2);
% vec1 and vec2 are standardform vectors both linear in different dunkns (or
% their dfdunkns). Are they equivalent except for the name of their dunkn,
% up to some constant multiple?
% basisvectors are toplevel mvars here.
if null vec1 then null vec2
else if null vec2 then nil
else ((mvar vec1)=(mvar vec2)) and equivf(lc vec1,lc vec2) % ie equiv lt
      and vectorequivf1(red vec1,red vec2);

symbolic procedure equivf(u,v);
% u and v are standardforms both linear in different dunkns (or their dfdunkns).
% Are they equivalent except for the name of their dunkn, 
% up to some constant multiple?
% dunkns are toplevel mvars here.
if null u then null v
else if null v then nil
else equivk(u,v) and ratiop(lc u,lc v) and equivf(red u,red v);

symbolic procedure equivk(u,v);
% u and v are different dfdunkns.
% Are they equivalent except for the name of their dunkn?
% ie do they have the same derivative list?
if not ((car u)='df) then not ((car v)='df)
else if not ((car v)='df) then nil
else (cddr u)=(cddr v);

symbolic procedure ratiop(u,v);
% u and v are standard forms
% does u = ratio * v?
% *** can be done much better than this...
begin scalar r;
  r:=multsq(u ./ 1, 1 ./ v);
  if not (domainp numr r and domainp denr r) 
    or (ratio and r neq ratio)
      then return nil;
  if not ratio then ratio:=r;
  return t;
  end;

%******end of op!*rmred2

symbolic procedure op!*sub2sf;
% Takes !*deteqns
% returns them in canonical form, so that each equation has been substituted
% into all the others.
%
% Now we want op!*sub2sf to be able to work by returning after it has done one substitution.
% We need to keep track of where all the equations are:
% they will be found on !*canonlist and holdlis.
%
if !*deteqns and null cdr !*deteqns then !*canonlist:=!*deteqns
else if !*deteqns then begin 
  scalar holdlis,eqn,oldeqngrp,neweqngrp,madesub,
         neweqnoncanon,neqnflag,eqnsused,op!*sub2sftries,oldcan,br;
  integer subdepth;
  neweqnoncanon:=nil; % subdepth:=0;
  if !*traceop!*sub2sf then <<write(tmstmp(),"op!*sub2sf on ",egnums !*deteqns);terpri() >>;
  !*canonlist:=intersection(!*canonlist,!*deteqns);
  holdlis:=complement(!*deteqns,!*canonlist);
  oldcan:=!*canonlist;
  while holdlis and not (!*op!*sub2sfby1 and neweqnoncanon) do begin
    if !*canonbyshortest                                     % start again
      then oldeqngrp:=shortesteg holdlis
      else oldeqngrp:=lowdunkneqngrp holdlis;
    holdlis:=delete(oldeqngrp,holdlis);
    eqn:=eg!*eqn oldeqngrp;
    op!*sub2sftries:=eg!*op!*sub2sftries oldeqngrp; % dont need to try to sub these in...
    eqnsused:=nil; subdepth:=nil;
    if !*tracesubeqns then <<
      write "About to sub eqngrps ",egnums !*canonlist," into ";terpri();
      showeqngrp oldeqngrp;
      write "We don't use op!*sub2sftries= ",op!*sub2sftries;terpri()>>;
    if !*subeqnsbydunkn
      then eqn:=newsubeqns(eqn,!*canonlist)
      else eqn:=subeqns(eqn,!*canonlist);
    addop!*sub2sftries(egnums !*canonlist,oldeqngrp);
    if neqnflag then begin % subs made
      solvesuccess:=t;
      if eqn 
        then begin
          if (red eg!*eqn oldeqngrp or null !*keep1tms) then << % ie keep 1-term eqns
            deathcert(oldeqngrp,list("subbed in op!*sub2sf using equations ",eqnsused,
                                     " and became eqngrp ",!*dets+1));
            !*deteqns:=delete(oldeqngrp,!*deteqns)>>;
          neweqngrp := mkeqngrp(eqn,eg!*no oldeqngrp,"op!*sub2sf from subbing eqns".eqnsused);
          addeqngrp neweqngrp;
          if !*keepsubing then addop!*sub2sftries(egnums !*canonlist,neweqngrp);
          neweqnoncanon:=t;
          end
        else <<deathcert(oldeqngrp,list("subbed in op!*sub2sf using equations ",eqnsused,
                                   " and simplified to zero."));
               !*deteqns:=delete(oldeqngrp,!*deteqns)>>;
               % delete it now or else we might move on by...
      end % of if eqn then << >> else 
      else neweqngrp:=oldeqngrp; % no subs made
    if eqn % then begin % insert it into !*canonlist, return the top to holdlis
      then <<br:=ins2canon(neweqngrp,!*canonlist,holdlis);
             !*canonlist:=car br;
             holdlis:=cadr br>>;
    end; % of while holdlis and not (!*op!*sub2sfby1 and neweqnoncanon)
  if !*tracecanon then <<
    write(tmstmp(),"!*canonlist= ",egnums !*canonlist);
    if holdlis then <<write(" holdlis= ",egnums holdlis); terpri()>>;
    write "!*canonlist is now "; terpri();
    write "   ",for each x in !*canonlist collect eg!*highdfdunkn x;terpri();
%    write "holdlist is now "; terpri();
%    write "   ",for each x in holdlis collect eg!*highdfdunkn x;terpri()
    >>;
  if !*traceslow then <<write "<";terpri()>>;
  if not(oldcan=!*canonlist) or neweqngrp then solvesuccess:=t;
  end;

symbolic procedure ins2canon(neqngrp,clist,hlist);
% insert neqngrp into clist (canonlist) and return rest on hlist
% ie neqngrp will be last on clist (or first on hlist if equal highdunkns)
%
begin scalar newc,insflag;
  if !*tracecute then write "&";
  while clist and not insflag do % find where this eqn slots in
    if dunknorder(eg!*highdfdunkn car clist,eg!*highdfdunkn neqngrp)
        or (eg!*highdfdunkn car clist)=(eg!*highdfdunkn neqngrp)
      then insflag:=t
      else <<newc:=(car clist).newc; clist:=cdr clist>>;
   %
   % Now we have neqngrp comes after all those on newc,
   % and before those on clist (unless hdunkns eq with car clist)
   % (we could have highdfdunkns equal here, 'cause they wont sub out w/o !*keepsubing)
   %
   if !*tracecute and clist then write ("_",count clist);
   if !*tracecanon and clist then write "Got equation for ",eg!*highdfdunkn neqngrp;
%   if !*tracecanon and clist then <<
%     write "Got equation for ",eg!*highdfdunkn neqngrp," which comes before"; terpri();
%     write "   ",for each x in clist collect eg!*highdfdunkn x;terpri()
%      >>;
   if null !*keepsubing or null clist
       or (eg!*highdfdunkn car clist) neq (eg!*highdfdunkn neqngrp)
     then <<newc:=reverse (neqngrp.newc);
            hlist:=append (clist,hlist)>>
     else <<newc:=reverse ((car clist).newc);
            hlist:=append (neqngrp.(cdr clist),hlist)>>;
%   if !*tracecanon then <<
%     write "!*canonlist is now "; terpri();
%     write "   ",for each x in newc collect eg!*highdfdunkn x;terpri();
%     write "holdlist is now "; terpri();
%     write "   ",for each x in hlist collect eg!*highdfdunkn x;terpri()>>;
   if !*traceop!*sub2sf then <<
     terpri();
     write(tmstmp(),"after sub, !*canonlist:= ",egnums !*canonlist); terpri()>>;
   verifycanonorder newc;
   return list(newc,hlist);
   end;

symbolic procedure verifycanonorder clist;
%
% Make sure that each egngrp in clist (!*canonlist)
% is ordered before each one after it
%
% It would also be good to test that leading derivs
% dont appear in lhs's
%
for each eglist on clist
  do for each eg in (cdr eglist) do
    if not dunknorder(eg!*highdfdunkn eg,eg!*highdfdunkn car eglist)
      then <<terpri();terpri();
        write "***Failed to verify ordering of clist(=!*canonlist?):";terpri();
        write egnums clist;terpri();
        write "   ",for each x in clist collect eg!*highdfdunkn x;terpri();
        write "Got bad order with ",eg!*highdfdunkn car eglist,
              " and ",eg!*highdfdunkn eg;terpri();
        rederr "Bad ordering!">>;

symbolic procedure shortesteg eqngrplis;
% which is the shortest eqngrp?
% Use length instead of termsf because this is
% a more accurate estimate of subbing difficulty.
if null eqngrplis then nil
else if null cdr eqngrplis then car eqngrplis
else begin scalar short,sl;
  short:=car eqngrplis;
  sl:=slength short;
  for each eqngrp in cdr eqngrplis do if slength(eqngrp)<sl 
    then <<short:=eqngrp; sl:=count eg!*eqn short>>;
  return short;
  end;
  
symbolic procedure slength eqngrp;
begin scalar f,fd,l;
  f:=funknsf eg!*eqn eqngrp;
  if f then for each funkn in f do fd:=union(fd,depsk f);
  if f then l:=3+count fd else l:=0;
  return l+count eg!*eqn eqngrp;
  end;
   
symbolic procedure subeqns(eqn,subeqngrplis);
%
% take standard form eqn.
% substitute all eqns in subeqngrplis.
% return standardform eqn.
%
begin scalar sublis,madesub;
  % eqnsused, neqnflag and op!*sub2sftries are fluid
  % subdepth passed as fluid, just to keep score of how hard it was to sub this eqn in.
  eqnsused:=nil; % we havent yet subbed any eqns into eqn
  neqnflag:=nil;
  if !*subeqnsinreverse then subeqngrplis:=reverse subeqngrplis; 
  sublis:=subeqngrplis;
  while sublis and eqn do begin
    madesub:=nil; subdepth:=0;
    if (eg!*no car sublis member op!*sub2sftries)
      then <<if !*tracecute then write "#">>
      else begin % actually try this sub
      if !*tracecute then write "?";
      if !*tracesubeqns 
        then <<terpri();
               write "about to sub "; terpri();
               showeqngrp car sublis; 
               write "into"; terpri();
               showf eqn; terpri()
             >>;
      eqn:=mynumr mysubf(eqn,(eg!*highdfdunkn car sublis),
                         (eg!*subval car sublis) );
      if !*tracesubeqns then <<write "Got :"; terpri(); showf eqn>>;
      if madesub then begin
          if !*tracecute then if eqn then write("+(",subdepth,")") else write "0";
          eqnsused:=(eg!*no car sublis).eqnsused;
          op!*sub2sftries:=nil; % this is a new eqn, so we can sub them in again :-(
          eqn:=mynumr mysimpf eqn;
          if !*keepsubing then eqn:=dropeqnfacE eqn; % ie, if it wont be done by mkeqngrp
          if !*tracesubeqns then begin
            write(tmstmp(),"subbed ",eg!*highdfdunkn car sublis,
                  " from equation ",eg!*no car sublis," and got:"); terpri();
            showf eqn;
            end;
          if !*keepsubing
            then sublis:=subeqngrplis
            else sublis:=nil;               % so that we dont do any more subs...
          neqnflag:=t;
          end;
      end; % of trying to actually sub this eqn in
      if not madesub then sublis:=cdr sublis;
    end;
  if !*tracecute and eqn and not neqnflag then write "|";
  return eqn;
  end;

symbolic procedure newsubeqns(eqn,subeqngrplis);
%
% take standard form eqn.
% Substitute all eqns in subeqngrplis,
% except those on op!*sub2sftries for the first attempt,
% because these have already been fully subbed before.
% (actually, if subdepth=0 then we might be able to keep avoiding them???)
% return standardform eqn.
%
begin scalar dfdunkns2sub,dfdunkn2sub,shortsublis,trsub,subeqngrp,madesub;            % #1
  % subeqngrplis, eqnsused, neqnflag and op!*sub2sftries are fluid
  % subdepth passed as fluid, just to keep score of how hard it was to sub this eqn in.
  eqnsused:=nil; % we havent yet subbed any eqns into eqn
  neqnflag:=nil;
  if !*subeqnsinreverse then subeqngrplis:=reverse subeqngrplis; 
  for each eqngrp in reverse subeqngrplis do
    if (eg!*no eqngrp member op!*sub2sftries)
      then <<if !*tracecute then write "#">>
      else <<if !*tracecute then write "?";
             shortsublis:=eqngrp.shortsublis>>;
  dfdunkns2sub:=sort(dfdunknsE eqn,'dunknorder);
  while dfdunkns2sub do begin                            % #2 try to sub for each dfdunkn in eqn
    madesub:=nil; subdepth:=0;
    dfdunkn2sub:=car dfdunkns2sub;
    dfdunkns2sub:=cdr dfdunkns2sub;
    trsub:=trysub(dfdunkn2sub,eqn,shortsublis);
       % this just subs subval for highest order dfdunkn possible
       % returns list(madesubflag,eqn,subeqngrp)
    if car trsub then begin
        eqn:=cadr trsub;
        subeqngrp:=caddr trsub;
        if !*tracecute then if eqn then write("+(",subdepth,")") else write "0";
        eqnsused:=(eg!*no subeqngrp).eqnsused;
        op!*sub2sftries:=nil; % this is a new eqn, so we can sub them all in again :-(
        shortsublis:=subeqngrplis; % ditto
      if !*keepsubing
        then dfdunkns2sub:=sort(dfdunknsE eqn,'dunknorder)
        else dfdunkns2sub:=nil;                     % so that we dont do any more subs...
      neqnflag:=t;
      end;
    end;                                              % #2
  if !*tracecute and eqn and not neqnflag then write "|";
  return eqn;
  end;

fluid '(onlysubthisdfdunkn);

symbolic procedure trysub(dfdunkn2sub,eqn,sublis);
% try to sub each eqngrp in sublis for dfdunkn2sub in eqn.
% just do one substitution.
%
% return list(madesubflag,eqn,subeqngrp)
begin scalar madesubflag,subeqngrp;
  if !*tracetrysub then <<write "trying to sub for ",dfdunkn2sub;terpri()>>;
  while sublis and not madesubflag do begin            % #1 try to sub each eqn for dfdunkn2sub
    subeqngrp:=car sublis;
    sublis:=cdr sublis;
%    if !*tracecute then if not (eg!*no subeqngrp member op!*sub2sftries) 
%      then write "#"
%      else write "?";
    if %not (eg!*no subeqngrp member op!*sub2sftries) and
       minusdfdunkns(dfdunkn2sub,eg!*highdfdunkn subeqngrp) neq -1
      then begin                                       % #2 actually make the sub
        madesubflag:=t;
        if !*tracecute then write "]?[";
        if !*tracesubeqns 
          then <<terpri();
                 write "about to sub "; terpri();
                 showeqngrp subeqngrp; 
                 write "into"; terpri();
                 showf eqn; terpri()
               >>;
        onlysubthisdfdunkn:=dfdunkn2sub;
        eqn:=dropeqnfacE mynumr mysimpq mysubf(eqn,(eg!*highdfdunkn subeqngrp),
                                      (eg!*subval subeqngrp) );
        if !*tracesubeqns then <<write "Got :"; terpri(); showf eqn>>;
        end;                                            % #2
    end;
  if madesubflag
    then return list(t,eqn,subeqngrp)
    else return list(nil);
  end;

symbolic procedure mysubf(eqn,u,v);
% makes substitutions of kern u for sq v on the linear mvars in eqn.
% returns sq
if null eqn then (nil ./ 1)
else if !*subeqnsbydunkn then mynewsubf(eqn,u,v)
else addsq(dropfacsqE multsq(mysubk(mvar eqn,u,v),(lc eqn) ./ 1),
           mysubf(red eqn,u,v) );
  
symbolic procedure mynewsubf(eqn,u,v);
% makes substitutions of kern u for sq v on the linear mvars in eqn.
% returns sq
% Here we just sub for onlysubthisdfdunkn, which is passed as fluid
if null eqn then (nil ./ 1)
else begin scalar coeff,rest,sval,nsval,dsval,neqn,fac;
  coeff:=coeffE(onlysubthisdfdunkn,eqn); 
  rest:=addf(negf multf(!*k2f onlysubthisdfdunkn,coeff),eqn);
  sval:=mysubk(onlysubthisdfdunkn,u,v);
  nsval:=mynumr sval;
  dsval:=denr sval;
  fac:=mygcdf!*(coeff,dsval);
  dsval:=quotf1(dsval,fac);
  coeff:=quotf1(coeff,fac);
  if nsval and rest then <<
    fac:=mygcdf!*(eqngcdE nsval,eqngcdE rest);
    % important not to have dunkns in fac!!
    nsval:=myquotf1E(nsval,fac);
    rest :=myquotf1E(rest ,fac);
    >>;
  neqn:=addf(multf(nsval,coeff),
             multf(dsval,rest));
  % We dont need to (dropeqnfacE) here-- already dropped all posible factors
  % We are relying on there being no factors between nsval & dsval or coeff & rest,
  % which is of if they have been (dropeqnfacE)'d and (dropeqnfacsq)'d respectivly.
  return neqn ./ 1 ;                      % sq for compatibility
  end; 
  
symbolic procedure mysubk(x,y,val);
%
% makes substitutions of kern y for sq val on kernel x
% returns the sq value after substitution
%
% subdepth passed as fluid, just to keep score of how hard it was
% to sub this eqn in
%
if  (dunknin x) neq (dunknin y) then !*k2q x
else if x=y then <<madesub:=t; val>>
else if  (car x) neq 'df then !*k2q x
else if !*subeqnsbydunkn and x neq onlysubthisdfdunkn then !*k2q x
else begin scalar z;
  if !*tracemysubk then <<write(tmstmp(),"df subtracting ",y," from ",x); terpri() >>;
  z:=minusdfvarlist(dfvarlist x,dfvarlist y);
  if !*tracemysubk then <<write "got z= ",z; terpri() >>;
  if !*tracemysubk and z=-1 then <<write(tmstmp(),"can't sub here...");terpri() >>;
  if z=-1 then return !*k2q x;
  madesub:=t;
  if !*tracemysubk then <<write "subval is "; terpri(); showsq val>>;
  while z do <<
    subdepth:=subdepth+1;
    val:=dropfacsqE mysimpq diffsq(val,car z);
    if cdr z and numberp cadr z 
      then if cadr z=1 then z:=cddr z else z:=(car z).(cadr z-1).cddr z
      else z:=cdr z>>;
  if !*tracemysubk then <<write "subval is now "; terpri(); showsq val>>;
  return val;
  end;

symbolic procedure lowdunkneqngrp eqngrplis;
if eqngrplis then begin scalar low,rest;
  low := car eqngrplis;
  rest:= cdr eqngrplis;
  while rest do if dunknorder(eg!*highdfdunkn low,eg!*highdfdunkn car rest)
    then <<low := car rest; rest:= cdr rest>>
    else rest:= cdr rest;
  return low
  end;

symbolic procedure highdfdunknf eqn;
% find the highest ordered dfdunkn
% **must be nonzero eqn!
if eqn then begin scalar highdfdunkn, resteqn;
  highdfdunkn:=mvar eqn;
  resteqn:=eqn;
  while red resteqn do if not dunknorder(highdfdunkn,mvar (resteqn:=red resteqn))
     then highdfdunkn := mvar resteqn;
  return highdfdunkn;
  end;

% Equations are considered to be solved for their dfdunkn which is highest
% according to the ordering >_s, which is implemented by the procedure
% dunknorder.
%
% This is used to substitute one equation into another, and tends to result in
% some equations just involving dunkns which are _lower_ in the ordering >_s.


symbolic procedure dunknorder(dunkn1,dunkn2);
% Is dunkn1 ordered ahead of dunkn2?  (dfdunkns, actualy...)
% ie, the predicate dunkn1 >_s dunkn2
%
% Eg, want df(a,x,2) >_s df(a,x)
%
begin scalar x,y;
  if dunkn1=dunkn2 then return nil;
  if null dunkn1 then return nil;
  if null dunkn2 then return t;
  if (x:=count dunknin dunkn1) neq (y:=count dunknin dunkn2)
    then return x>y;
  if  (x:=recnumdeps dunknin dunkn1) neq (y:=recnumdeps dunknin dunkn2)
    then return x>y;
  if (x:=dford dunkn1) neq (y:=dford dunkn2)
    then return x>y;
  if  (x:=dfvarlist dunkn1) neq (y:=dfvarlist dunkn2)
    then return dfvarorder(x,y);
  if (x:=dunknin dunkn1) neq (y:=dunknin dunkn2)
    then return not ordop(x,y); % oldest first...
  rederr list("dunknorder failed to distinguish between ",
               dunkn1," and ",dunkn2);
  end;

symbolic procedure recnumdeps dunkn;
% how many variables (x's and u's) does dunkn depend on?
% (at the time this function was first called with this arg for consistancy)
begin scalar x;
  if null (x:=assoc(dunkn,!*recnumdeps))
    then <<!*recnumdeps:=list(dunkn,count depsk dunkn).!*recnumdeps;
           return count depsk dunkn>>
    else return cadr x;
  end;

symbolic procedure dford u;
% u is a dfdunkn kernel.
% what order derivative is u?
if car u neq 'df then 0
else begin integer k;
  u:=cddr u;
  while u do
    if cdr u and numberp cadr u then <<k:=k+cadr u; u:=cddr u>>
    else <<k:=k+1; u:=cdr u>>;
  return k
  end;

symbolic procedure dfvarorder(dfvars1,dfvars2);
% Is the derivative list dfvars1 ordered ahead of dfvars2?
% The highest ordered
% one will be the one differentiated more by a higher ordered var
% Note that this ordering is very dependent on the default df ordering...
% so go with ordop...
%
%df(y,var1,var2,var3) will satisfy ordop(var1,var2)=ordop(var1,var3)=ordop(var2,var3)=t
%
% Eg, want df(a,x,2) >_s df(a,x),
%          df(a,x)   >_s df(a,y)   if ordop(x,y)  (so we can compare first derivatives first)
%
if null dfvars1 then nil
else if null dfvars2 then t % neither of these cases should ever occur...
%else if not car dfvars1=car dfvars2 then not ordop(car dfvars1,car dfvars2)
else if not ((car dfvars1)=(car dfvars2)) then ordop(car dfvars1,car dfvars2)
else if numberp cadr dfvars1 
  then if numberp cadr dfvars2 
    then if not (cadr dfvars1=cadr dfvars2)
          then (cadr dfvars1 > cadr dfvars2)
          else dfvarorder(cddr dfvars1,cddr dfvars2)
    else t
else if numberp cadr dfvars2 then nil
else dfvarorder(cdr dfvars1,cdr dfvars2);

%*****end of op!*sub2sf

symbolic procedure intcon1;
% For each equation on canonlist, consider the integrability conditions relating to the
% dependence of the highdfdunkn term (lack of dependence is equiv to a first order eqn),
% and calculate only the one with lowest new sdfdunkn
% Return a list of all new equations. 
begin scalar eqngrp, intlis, pr2do, new, eqngrps2do;
  if !*traceops then <<write "in intcon1"; terpri() >>;
  for each eqngrp in intersection(!*deteqns,!*canonlist) do 
    if not member('intcon1,eg!*ops eqngrp)
        and not (!*stopatnum and !*dets > !*stopatnum)
        and (new:=intcon1newsdunkns eqngrp)                 % so that there are intcons to get
      then if !*intcon1by1
        then <<if (null pr2do or dunknorder(cadr pr2do,car new))
                then pr2do:=list(eqngrp,car new)
              >>
        else eqngrps2do:=eqngrp.eqngrps2do;
  if pr2do then intlis:=calcintcon1 car pr2do;
  for each eqngrp in eqngrps2do do intlis:=append(calcintcon1 eqngrp,intlis);
  if intlis then solvesuccess:=t;
  return intlis;
  end;

symbolic procedure intcon1newsdunkns eqngrp;
% what are the (likely) sdfdunkns in the equations that will result from intcon1 eqngrp?
% return them sorted by dunknorder
begin scalar deps,sdunkn,sdunkns,new,eqn;
  deps:=alldepsf eg!*eqn eqngrp;
  eqn:=eg!*eqn eqngrp;
  for each var in complement(deps,depsk eg!*highdfdunkn eqngrp) do begin
    sdunkn:=nil;
    while eqn do begin
      if var member alldepsf lc eqn then new:=mvar eqn
        else if var member depsk mvar eqn then new:=mydiffk(mvar eqn,var)
        else new:=nil;
      if new and (null sdunkn or dunknorder(sdunkn,new))
        then sdunkn:=new;
      eqn:=red eqn;
      end;
    if sdunkn then sdunkns:=sdunkn.sdunkns;
    end;
  return reverse sort(sdunkns,'dunknorder); % ordered from lowest to highest 
  end;

symbolic procedure calcintcon1 eqngrp;
begin scalar neqngrp, var, deps, neqn, madesub, intlis, sublis;
  addopstried('intcon1,eqngrp);
  deps:=alldepsf eg!*eqn eqngrp;
  if complement(deps,boundvarsE eg!*eqn eqngrp)
    then <<write "Intcon1 on eqngrp ",eg!*no eqngrp," with free vars ",
                  complement(deps,boundvarsE eg!*eqn eqngrp),".";terpri();
           write "Much faster to use op!*splitE.";terpri();
           !*splitvars:=union(!*splitvars,complement(deps,boundvarsE eg!*eqn eqngrp))>>;
           % this is a slower form of spliting equations...
  for each var in complement(deps,depsk eg!*highdfdunkn eqngrp) do begin
    % ie for each var in the eqn, except ones which would just sub out again
    neqn:=difff(eg!*eqn eqngrp,var);
    if neqn then neqn:=mynumr mysimpf mynumr neqn;
    if !*traceintcon1 then begin
      write(tmstmp(),"intcon1: Differentiated equation ",eg!*no eqngrp,
            " wrt ",var," to get"); terpri();
      showf neqn;
      end;
    if neqn then begin
      neqngrp := mkeqngrp(neqn,eg!*no eqngrp,'intcon1);
      write("made eqn ",eg!*no neqngrp," by differentiating eqn ",eg!*no eqngrp,
            " by ",var); terpri();
      intlis:=neqngrp.intlis;
      if !*intcon1subnew then sublis:=neqngrp.sublis;
      end;
%    solvesuccess:=t;
%        % ?? Even if we haven't added a new condition, we have eliminated one possibility
%        % ?? and we need to try another possibility if there is one.
    end; % of this var
  return intlis;
  end;

fluid '(gotnew sublis intlis intcontwo);

symbolic procedure addintcontwo2;
begin scalar intcontwo;
  intcontwo:=t;
  addintcon2();
  end;

symbolic procedure intcon2;
% Takes !*deteqns in canonical form on !*canonlist.
% For each pair of equations, calculate the integrability conditions,
% and then sub all the eqns on the !*canonlist into this eqn.
% Return a list of all new equations. These have to be added into the canonical
% list, and any new integrability conditions then calculated.
begin scalar eqngrp, eqngrplis, eqngrp1, eqngrp2, eqn,
      madesub, intlis, sublis, cons2try, minpr, pr1, pr2, pr, notminlist,
      dunkns, lowdunkncons2try;
  if !*traceops then <<write tmstmp(),"in intcon2"; terpri() >>;
  if !*traceslow then <<write ">intcon2"; terpri()>>;
  sublis:=intersection(!*deteqns,!*canonlist);
  for each eqngrplis on intersection(!*deteqns,!*canonlist) do
    for each eqngrp2 in cdr eqngrplis do
      if not member((eg!*no (eqngrp1:=car eqngrplis)),eg!*i2tries eqngrp2)
        and dunknin(eg!*highdfdunkn eqngrp1)=dunknin(eg!*highdfdunkn eqngrp2)
        and not (!*stopatnum and !*dets > !*stopatnum) 
          then cons2try:=list(intcon2newsdunkn(eqngrp1,eqngrp2),
                              intconweight(eqngrp1,eqngrp2),
                              eqngrp1,
                              eqngrp2,
                              intcon2cdunkn(eqngrp1,eqngrp2))
                         .cons2try;
  if null cons2try and !*traceslow then <<write "<"; terpri()>>;
  if null cons2try then return;
  if !*intcon2onmin then begin
    for each pr1 in cons2try do for each pr2 in cons2try do
      if not (pr1=pr2 or minusdfvarlist(car cddddr pr1,car cddddr pr2)=-1)
        then notminlist:=pr1.notminlist;       % ie not a minimal int con
    cons2try:=complement(cons2try,notminlist); % we only want to calculate minimal int cons
    if !*traceweight or !*traceintcon2 then 
      for each pr in notminlist do 
        <<write "weight for ",car cddddr pr," between ",eg!*highdfdunkn caddr pr," and ",
                eg!*highdfdunkn cadddr pr," is ",cadr pr; terpri();
          write "expected sdunkn from intcon between ",eg!*highdfdunkn caddr pr," and ",
                eg!*highdfdunkn cadddr pr," is  ",car pr; terpri();
          write "But this isn't a minimal intcon so I won't consider it (for now)."; terpri()>>;
    end;
  if null cons2try then <<
     cons2try:=notminlist;
     terpri();write "All minimal intcons evaluated, now we do the rest";terpri();terpri();
     >>;
  if null cons2try then rederr "failed to find max intcon2s";
  if !*intcon2onlowdunkn then begin
    for each pr in cons2try do
      if not member((dunknin car cddddr pr),dunkns)
        then dunkns := (dunknin car cddddr pr).dunkns;
    dunkns:=reverse sort(dunkns,'dunknorder); % sort them from lowest ordered to highest
    while null lowdunkncons2try and dunkns do begin
      for each pr in cons2try do
        if (dunknin car cddddr pr) = car dunkns then lowdunkncons2try :=pr.lowdunkncons2try;
      dunkns:=cdr dunkns; % if none of these, then look for intcons in next dunkn
      end;
    cons2try:=lowdunkncons2try;
    end;
  if null cons2try then rederr "failed to find intcon2s in any min dunkn";
  minpr:=car cons2try; 
  if !*intcon2bylowsdfdunkn

    then for each pr in cdr cons2try do
     if car minpr and (dunknorder(car minpr,car pr)

                       or ((car pr)=(car minpr) and ((cadr pr) < cadr minpr)))
       then minpr:=pr
    else for each pr in cdr cons2try do
      if (cadr pr) neq 0 and (cadr pr) < cadr minpr then minpr:=pr;
  if !*traceweight or !*traceintcon2 then begin
    for each pr in cons2try do 
      <<write "weight for ",car cddddr pr," between ",eg!*highdfdunkn caddr pr," and ",
              eg!*highdfdunkn cadddr pr," is ",cadr pr; terpri();
        write "expected sdunkn from intcon between ",eg!*highdfdunkn caddr pr," and ",
              eg!*highdfdunkn cadddr pr," is  ",car pr; terpri()>>;
    if minpr then <<
      write "weight of min intcon is ",cadr minpr; terpri();
      write "expected sdunkn is ",car minpr; terpri()>>;
    end;
  if !*intcon2by1 then cons2try:=list minpr;
  for each pr in cons2try do begin
    eqngrp1:=caddr pr;
    eqngrp2:=cadddr pr;
    eqn:=calcintcon2(eqngrp1,eqngrp2);
    if eqn then eqn:=mynumr mysimpf eqn;
    if eqn then if !*subonintcon2 then eqn:=subeqns(eqn,sublis);
    if eqn then begin
      eqngrp:=mkeqngrp(eqn,list(eg!*no eqngrp1, eg!*no eqngrp2),'int2);
      terpri();
      write("made eqn ",eg!*no eqngrp," from integrability conditions on eqns ",
            eg!*no eqngrp1," and ",eg!*no eqngrp2," for ",car cddddr pr); terpri();
      intlis:=eqngrp.intlis;
      if !*intcon2subnew then sublis:=eqngrp.sublis;
      end;
    if !*traceintcon2 then <<write "Integrability condition after subs is";terpri();showf eqn>>;
    end; % of making this intcon
  if !*traceslow then <<write "<";terpri()>>;
  solvesuccess:=t; % even if we haven't added a new condition, we have eliminated one possibility
                   % and we need to try another possibility if there is one
  return intlis;
  end;

symbolic procedure intcon2newsdunkn(eqngrp1,eqngrp2);
% what is the (likely) sdfdunkn in the equation that will result from 
% calculating the integrability conditions between these two eqngrps?
begin scalar s1,s2,z,temp; % #3
  temp:=minusdfvarlisttwo(dfvarlist eg!*highdfdunkn eqngrp1,dfvarlist eg!*highdfdunkn eqngrp2);
%  if numr eg!*subval eqngrp1 then s1:=highdfdunknf numr eg!*subval eqngrp1;
%  if numr eg!*subval eqngrp2 then s2:=highdfdunknf numr eg!*subval eqngrp2;
  if (not onetegp eqngrp1 and null car temp)
     or (not onetegp eqngrp1 and null cadr temp)
        then rederr "in intcon2newsdunkn, seen sub2sf failed";
  s1:=eg!*nexthigh eqngrp1;
  s2:=eg!*nexthigh eqngrp2;
  if !*i2highlong
    then begin
      s2:=intcon2newsdunkn1(eg!*subval eqngrp2,car temp);
      s1:=intcon2newsdunkn1(eg!*subval eqngrp1,cadr temp);
      end
    else begin % this is just a rough estimate
      z:=car temp;
      if s2 then while z do <<
        if s2 and car z member depsk s2 then s2:=mydiffk(s2,car z);
        if cdr z and numberp cadr z 
          then if (cadr z)=1 then z:=cddr z else z:=(car z).((cadr z)-1).(cddr z)
          else z:=cdr z>>;
      z:=cadr temp;
      if s1 then while z do <<
        if s1 and car z member depsk s1 then s1:=mydiffk(s1,car z);
        if cdr z and numberp cadr z 
          then if (cadr z)=1 then z:=cddr z else z:=(car z).((cadr z)-1).(cddr z)
          else z:=cdr z>>;
      end;
  if s1 and s2 
    then if dunknorder(s1,s2) then return s1 else return s2
    else if s1 then return s1
    else if s2 then return s2;
  end;

symbolic procedure intcon2newsdunkn1(sqeqn,z);
% what are the (likely) sdfdunkns in the equations that will result from diff sqeqn by varlist z?
% find them and take the lowest
begin scalar pveqn, npveqn, var, term, z, sdunkn;
  pveqn:=sq2perv sqeqn; % perverse eqn;
  while z do begin
    var:=var z;
    npveqn:=nil;
    for each term in pveqn do
      if var member alldepsf numr cadr term or var member alldepsf denr cadr term
        then npveqn:=list(car term,
                          multsq(difff(numr cadr term,var), 1 ./ denr cadr term))
                     .npveqn
      else if var member depsk car term
        then npveqn:=list(mydiffk(car term,var), cadr term).npveqn
    if cdr z and numberp cadr z 
      then if (cadr z)=1 then z:=cddr z else z:=(car z).((cadr z)-1).(cddr z)
      else z:=cdr z;
    pveqn:=npveqn;
    end;
  if npveqn then <<sdunkn:=caar npveqn; npveqn:=cdr npveqn>>;
  while npveqn do begin
    if dunknorder(sdunkn, caar npveqn) then sdunkn:=caar npveqn;
    npveqn:=cdr npveqn;
    end;
  return sdunkn;
  end;

symbolic procedure sq2perv sqeqn;
% take a sq eqn, return it as a pveqn,
% ie a list of list(dunkn, sq coeff).
if null numr sqeqn then nil
else list(mvar numr sqeqn, multsq((lc numr sqeqn) ./ 1, 1 ./ denr sqeqn))
      . sq2perv( (red numr sqeqn) ./ denr sqeqn);

symbolic procedure perv2sq pveqn;
% take a pveqn, return it as a sq eqn
% we need to multiply out because factors might have been dropped...
if null pveqn then ( nil ./ 1)
else addsq( multsq( !*k2q caar pveqn, cadar pveqn),
            perv2sq cdr pveqn
           );
  
symbolic procedure intconweight(eqngrp1,eqngrp2);
% how costly is it to calculate the integrability conditions between these two eqngrps?
% This will depend on the number of differentiations involved,
% how long the equations are, and what order derivatives the new eqn will have.
% For now, we use w = (ord of highest deriv (of dunkn) in resulting eqn) 
%                   + (n:=number of derivs needed to get it)
%                   + (3 + #deps funkns) if there are any funkns
%
begin scalar s1,s2,ss1,ss2,o1,o2,sdunkn1,sdunkn2,f1,f2,fd1,fd2;
      integer n,m,mm,nfd1,nfd2,w;
  s1:=numr eg!*subval eqngrp1;
  s2:=numr eg!*subval eqngrp2;
  if null s1 and null s2 then return 1000000;  % integrability condition will be trivial
  f1:=funknsf eg!*eqn eqngrp1;                 % the funkns in eqngrp1
  f2:=funknsf eg!*eqn eqngrp2;                 % the funkns in eqngrp2
  for each funkn in f1 do fd1:=union(fd1,depsk funkn);  
  for each funkn in f2 do fd2:=union(fd2,depsk funkn);
  nfd1:=count fd1;                             % the number of funkn dependencies in eqngrp1  
  nfd1:=count fd1;                             % the number of funkn dependencies in eqngrp2  
%  ss1:=(if s1 then highdfdunknf s1 else nil);
%  ss2:=(if s2 then highdfdunknf s2 else nil);
  ss1:=eg!*nexthigh eqngrp1;
  ss2:=eg!*nexthigh eqngrp2;
  o1:=(if s1 then dford ss1 else 0);           % order of highest deriv after eg!*highdfdunkn
  o2:=(if s2 then dford ss2 else 0);           % order of highest deriv after eg!*highdfdunkn
  sdunkn1:=eg!*highdfdunkn eqngrp1;
  sdunkn2:=eg!*highdfdunkn eqngrp2;
  for each var in union(dfvars sdunkn1,dfvars sdunkn2) do begin
    m:=vardforder(sdunkn1,var)-vardforder(sdunkn2,var);
    if m<0 then mm:=-m else mm:=m;
    n:=n+mm;                                   % number of differentiations needed to get intcon
    % want to work out likely sdunkn in new eqn...
    if m>0 then o1:=o1+m else o2:=o2-m;        % order of highest deriv (dunkn) in resulting eqn
    end;
  %
  %
  w:=min2(o1,o2);            % (expected) ord of highest deriv (of dunkn) in resulting eqn
  w:=w+n;                    % number of derivs needed to get it
  if f1 or f2 then w:=w+3;   % if there are any funkns
  w:=w+nfd1+nfd2;            % + #deps funkns
  return w;
  end;

symbolic procedure dfvars dfdunkn;
%  what vars is dfdunkn differentiated wrt.
if (car dfdunkn) neq 'df then nil
else begin scalar v,vars;
  v:=cddr dfdunkn;
  while v do if (cdr v and fixp cadr v)
    then <<vars:=(car v).vars; v:=cddr v>>
    else <<vars:=(car v).vars; v:=cdr v>>;
  return vars;
  end;

symbolic procedure vardforder(dfdunkn,var);
% how many times is dfdunkn df'd by var?
if  (car dfdunkn) neq 'df then 0
else begin scalar v;integer o;
  v:=cddr dfdunkn;
  while v do if (car v)=var
    then if (cdr v and fixp cadr v) 
      then <<o:= cadr v;v:=nil>>
      else <<o:= 1;v:=nil>>
    else v:=cdr v;
  return o;
  end;  
      
symbolic procedure intcon2cdunkn(eqngrp1,eqngrp2);
% the common dfdunkn found by differentiating each eqn as needed,
% so that the common dfdunkn can be equated to give integrability condition
begin scalar z,cdunkn;
  if  (dunknin eg!*highdfdunkn eqngrp1) neq (dunknin eg!*highdfdunkn eqngrp2)
    then rederr "No such integrability condition!";
  z:=car minusdfvarlisttwo(dfvarlist eg!*highdfdunkn eqngrp1,dfvarlist eg!*highdfdunkn eqngrp2);
  cdunkn:=eg!*highdfdunkn eqngrp2;
  while z do <<
    cdunkn:=mydiffk(cdunkn,car z);
    if cdr z and numberp cadr z 
      then if (cadr z)=1 then z:=cddr z else z:=(car z).((cadr z)-1).(cddr z)
      else z:=cdr z>>;
  return cdunkn;
  end;

symbolic procedure onetegp eqngrp;
% is this a onetermeqngrp?
null red eg!*eqn eqngrp;

symbolic procedure calcintcon2(eqngrp1,eqngrp2);
begin scalar e1, e2, z, eqn, sdfdunkn, fac, fac1, fac2, temp; 
  addi2tries(eg!*no eqngrp1,eqngrp2);
  addi2tries(eg!*no eqngrp2,eqngrp1);
  if !*traceintcon2 then begin
    write(tmstmp(),"intcon2 between equations ",eg!*no eqngrp1," and ",
           eg!*no eqngrp2," with common dunkn ",dunknin eg!*highdfdunkn eqngrp1);
    terpri(); showeqngrp eqngrp1; showeqngrp eqngrp2;
    end;
  temp:=minusdfvarlisttwo(dfvarlist eg!*highdfdunkn eqngrp1,dfvarlist eg!*highdfdunkn eqngrp2);
  if (not onetegp eqngrp1 and null car temp)
     or (not onetegp eqngrp1 and null cadr temp)
        then rederr "in calcintcon2, seen sub2sf failed";
  e1:=eg!*eqn eqngrp1;
  e2:=eg!*eqn eqngrp2;
  z:=car temp;
  while z do <<
    e2:=mynumr difff(e2,car z);
    if cdr z and numberp cadr z 
      then if (cadr z)=1 then z:=cddr z else z:=(car z).((cadr z)-1).(cddr z)
      else z:=cdr z>>;
  z:=cadr temp;
  while z do <<
    e1:=mynumr difff(e1,car z);
    if cdr z and numberp cadr z 
      then if (cadr z)=1 then z:=cddr z else z:=(car z).((cadr z)-1).(cddr z)
      else z:=cdr z>>;
  if null e1 or null e2 then rederr "calcintcon2 failed with zero eqn";
  sdfdunkn:=highdfdunknf e1;
  if sdfdunkn neq highdfdunknf e2 then rederr "calcintcon2 failed to match dfdunkns";
  if !*traceintcon2 then <<write "out of intcon2 diff'n loop with eqns";terpri()>>;
  e1:=mynumr mysimpf e1;
  e2:=mynumr mysimpf e2;
  if !*traceintcon2 then <<showf e1;terpri();showf e2;terpri()>>;
  fac1:=coeffE(sdfdunkn,e1);
  fac2:=coeffE(sdfdunkn,e2);
  if null fac1 or null fac2
    then rederr list("Failed to take coeffs of ",sdfdunkn);
  if !*traceslow or !*traceintcon2 
    then <<write "intcon2, going into mygcdf!*";terpri()>>;
  fac:=mygcdf!*(fac1,fac2);
  if !*traceslow or !*traceintcon2 
    then <<write "intcon2, out of mygcdf!*";terpri()>>;
  fac1:=quotf1(fac1,fac); fac2:=quotf1(fac2,fac);
  if not (fac1 and fac2) then rederr("Division failed in intcon2");
  eqn:=addf(multf(fac2,e1), negf multf(fac1,e2));
  if coeffE(sdfdunkn,eqn) 
    then rederr list("Failed to drop sdfdunkns in intcon2");
  if !*traceintcon2 then begin % #4'
    write "New equations are:"; terpri(); showsq e1; showsq e2;
    write("Common sdfdunkn is now ",sdfdunkn," with reduced coeffs");terpri();
    showf fac1; showf fac2;
    write "Integrability condition before substitution is"; terpri();showf eqn;
    end;  % #4'
  return eqn;
  end;


%***********

symbolic procedure op!*hidelg eqngrp;
%
% put any very equations onto !*longhideneqns
%
% We only do this once normally, after we have split the equations and solved any simple stuff.
%
% this op isnt a 'restart op...
%
if eg!*len eqngrp > !*eqnlengthlimit then begin
  !*longhideneqns:=inseqngrp(eqngrp,!*longhideneqns); % so it is in its sorted order.
  !*deteqns:=delete(eqngrp,!*deteqns);
  if !*tracehidden then <<
    write(tmstmp(),"Putting aside eqngrp ",eg!*no eqngrp," because it has ",
        eg!*len eqngrp," terms, which is too long.");terpri()>>;
  end;


symbolic procedure op!*findlg;
%
% pull (the easiest) equation off !*longhidenlis
% for now, just pull off the first.
% Well, maybe just take them all
if !*longhideneqns then
  if !*op!*findlgby1 then begin scalar eqngrp;
    eqngrp:=car !*longhideneqns;!*longhideneqns:=cdr !*longhideneqns;
    if !*tracehidden then <<
      write(tmstmp(),"Found equation group ",eg!*no eqngrp," which was hidden.");terpri()>>;
  %  if null !*longhideneqns then <<write "And now I've found all the lost sheep :-)";terpri()>>;
    addeqngrp eqngrp;
    solvesuccess:=t;
    if null !*longhideneqns then !*eqnlengthlimit:=2*!*eqnlengthlimit;
    end
  else begin scalar eqngrp;
    if !*tracehidden then <<
      write(tmstmp(),"finding eqngrps ",egnums !*longhideneqns);terpri()>>;
  %  write "And now I've found all the lost sheep :-)";terpri();
    for each eqngrp in !*longhideneqns do addeqngrp eqngrp;
    !*longhideneqns:=nil;
    !*eqnlengthlimit:=2*!*eqnlengthlimit;
    solvesuccess:=t;
    end;   


%**********************


symbolic procedure op!*hidefr;
%
% put any equations involving freeunknowns onto !*hideneqns
%
% We only do this once normally, after we have split the equtions and solved any simple stuff.
%
% We should also make sure there have been no assignments made involving freeunknowns...
% (or actually make sure it doesn't happen in the first place)
%
% this op isnt a 'restart op...
%if null !*hidenbefore and null !*hideneqns then begin scalar eqngrp;
if (!*keephiding and (nothiden !*deteqns)) or null !*hidenbefore and null !*hideneqns
then begin scalar eqngrp;
  !*hidenbefore:=t;
  for each eqngrp in !*deteqns do addopstried('op!*hidefr,eqngrp);
  for each eqngrp in !*deteqns do
    if funknsp eg!*eqn eqngrp then 
      <<!*hideneqns:=inseqngrp(eqngrp,!*hideneqns);
        !*deteqns:=delete(eqngrp,!*deteqns)>>;
  if null !*hideneqns then return;
%  if null !*deteqns then <<!*deteqns:=!*hideneqns; !*hideneqns:=nil; return>>;
%  ^ this line would cause new eqngrps to come back onto eqngrplis w/o restarting,
%   causing unsimped eqns to be solved! Just let the eqns come back with op!*findfr in good time
  if !*tracehidden then <<
    write(tmstmp(),"Hiding eqngrps ",egnums !*hideneqns);terpri();
    write("because they have freeunknowns in them.");terpri();
    write("But don't worry, I'll find them again later :-)");terpri()>>;  
%  !*hideneqns:=reverse !*hideneqns; % so it is in its sorted order.
  end;


symbolic procedure nothiden eqngrplis;
% is there an eqn with an funkn which hasn't been hidden before?
if null eqngrplis then nil
else if not ('op!*hidefr member eg!*ops car eqngrplis)
        and funknsp eg!*eqn car eqngrplis
  then t
else nothiden cdr eqngrplis;

symbolic procedure op!*findfr;
% pull (the easiest?) equations off !*hidenlis
if !*hideneqns then if !*op!*findfrby1
  then begin scalar eqngrp;
    eqngrp:=car !*hideneqns;!*hideneqns:=cdr !*hideneqns;
    if !*tracehidden then <<
      write(tmstmp(),"Found equation group ",eg!*no eqngrp," which was hidden.");terpri();
      if null !*hideneqns then <<write "And now I've found all the lost sheep :-)";terpri()>>;>>;
    addeqngrp eqngrp;
    solvesuccess:=t;
    end
  else begin scalar eqngrp;
    if !*tracehidden then <<
      write(tmstmp(),"finding eqngrps ",egnums !*hideneqns);terpri();
      write "And now I've found all the lost sheep :-)";terpri()>>;
    for each eqngrp in reverse !*hideneqns do addeqngrp eqngrp;
    !*hideneqns:=nil;
    solvesuccess:=t;
    end;   

symbolic procedure op!*rcvrfr;
%
% look for a free-hidden eqn which needs simping and then becomes not-free
% - should we then hide all free eqns in !*deteqns?
%
% actualy, we should test that any simp!* will shorten the eqn...
%
if !*hideneqns then begin scalar eqngrp, hideneqnstotry, neqngrp;
  hideneqnstotry:=!*hideneqns;
  while hideneqnstotry and not solvesuccess do begin
    eqngrp:=car !*hideneqns;
    hideneqnstotry:=cdr hideneqnstotry;
    if eg!*simpstatus eqngrp then begin
      neqngrp:=simpeqngrp eqngrp;
      if funknsp eg!*eqn neqngrp
        then !*hideneqns:=inseqngrp(neqngrp,!*hideneqns)
        else <<addeqngrp neqngrp; solvesuccess:=t>>;
      !*hideneqns:=delete(eqngrp,!*deteqns);
      end;
    end;
  end;

%******** end of op!*findlg

% is eqn an ode?
%
% level 1 test = GODE in dfdunkn, var: 
%   we need each dfdunkn to be df of dfdunkn by var,
%   or not depend on var.
%
% level 2 test = ODE in dfdunkn, var:
%   also need no coeffs depend on other vars
%
% level 3 test = CCODE in dfdunkn, var: 
%   also need coeffs all constant.
%
% Do we also want that dfdunkn depends on more than rest of eqn so we can solve for it ???
%

symbolic procedure op!*odetst eqngrp;
begin scalar odepr;
  if null (odepr:=godetest eg!*eqn eqngrp) then return;
  terpri();
  write "Found eqngrp ",eg!*no eqngrp," is a generalised ode for ",car odepr," in ",cdr odepr;
  terpri();
%  showeqngrp eqngrp;
  end;

symbolic procedure godetest eqn;
% is eqn a GODE in dfdunkn, var?
% If so, then return dfdunkn . var
begin scalar dfdunkns,alldfdunkns,otherdfdunkns,dfvars,GODEdfdunkn,thisdfdunkn,thisvar;
  dfdunkns:=alldfdunkns:=dfdunknsE eqn;
  while dfdunkns and not GODEdfdunkn do begin             % look at each dfdunkn,
    thisdfdunkn:=car dfdunkns;
    otherdfdunkns:=delete(thisdfdunkn,alldfdunkns);
    dfvars:=depsk thisdfdunkn;
    while dfvars and not GODEdfdunkn do begin             % and look at each var in each dfdunkn
      thisvar:=car dfvars;
      if godetest1(thisdfdunkn,thisvar,otherdfdunkns) then GODEdfdunkn:=thisdfdunkn;
      dfvars:=cdr dfvars;
      end;
    dfdunkns:=cdr dfdunkns;
    end;
  if GODEdfdunkn then return (GODEdfdunkn.thisvar);
  end;


symbolic procedure godetest1(odedfdunkn,odevar,otherdfdunkns);
% want to check that each dfdunkn in otherdfdunkns is df of odedfdunkn by odevar,
% or else that it doesnt depend on odevar.
begin scalar ans;
  ans:=t;                                     % assume true until we find out otherwise
  for each dfdunkn in otherdfdunkns do
    if ans and odevar member depsk dfdunkn and not godetest2(dfdunkn,odedfdunkn,odevar)
      then ans:=nil;
  return ans;
  end;

symbolic procedure godetest2(dfdunkn,odedfdunkn,odevar);
% is dfdunkn a df of odedfdunkn by odevar?
dunknin(dfdunkn)=dunknin(odedfdunkn) and
  ((lambda diff; (diff neq -1) and (car diff)=odevar and 
                (null cdr diff or (numberp cadr diff and null cddr diff)))
    minusdfdunkns(dfdunkn,odedfdunkn));
 

%************************** end of solve ops

%****** some output routines

symbolic procedure showdets;
begin scalar j,lead,para;
  j:=setdeteqns(); % sets the kvalues of deteqn for access in algebraic mode
                   % and returns the number of determining equations=deteqn 0
  if not (j=0)
    then write("There are ",j,
                 " determining equations remaining, which are...")
    else write "There are no determining equations remaining.";
  terpri();
  for k:=1:j do <<algebraic write "deteqn(",k,")=",deteqn(k);>>;
  if depl!* then <<terpri();write("The remaining dependencies are ..."); terpri()>>;
  showdepsused();
  terpri();
  if not (j=0) then begin
    for each eqngrp in alleqngrps() do lead:=(eg!*highdfdunkn eqngrp).lead;
    for each eqngrp in alleqngrps() do para:=union(dfdunknsE eg!*eqn eqngrp,para);
    for each var in lead do para:=delete(var,para);
    write "The dunkns in the remaining equations are: ",dunknsineqns();terpri();terpri();
    write "The leading derivatives are: ",lead; terpri();terpri();
    write "The parametric derviatives in the remaining equations are:";terpri();
    write para;terpri();terpri();
    end;
  end;

symbolic operator showdets;

symbolic procedure showeqngrps;
begin scalar eqngrp;
  for each eqngrp in !*deteqns do showeqngrp eqngrp;
  if !*hideneqns then <<write "And the hidden equationgroups are:";terpri()>>;
  for each eqngrp in !*hideneqns do showeqngrp eqngrp;
  showdeps();
  if !*longhideneqns then write "And there are also longhidden equationgroups.";terpri();
  end;

symbolic operator showeqngrps;

symbolic procedure setdeteqns;
% sets the kvalues of deteqn so they can be accessed in algebraic mode
begin integer j; scalar eqngrp;
  remprop('deteqn,'kvalue); % clears the old values?
  for each eqngrp in alleqngrps() do begin
    mysetk(list('deteqn,j:=j+1),mk!*sq mysimpf eg!*eqn eqngrp);
%    mysetk(list('deteqn,eg!*no eqngrp),mk!*sq mysimpf eg!*eqn eqngrp);
    mysetk(list('deteqn,minus j),eg!*no eqngrp); % the internal equation number
    end;
  mysetk('(deteqn 0),j);
  return j;
  end;

symbolic procedure mkgens;
begin 
  scalar svec, genfns, dunknsineqns, dunkn, gen, link, coeff,
         fingens, infgens, eqngens, dunknlinks;
  integer n,m,l;
  genfns:=dunknsinsymvec();
  dunknsineqns:=dunknsineqns();
  if !*tracecute then <<write("dunkns in eqns are ",dunknsineqns); terpri()>>;
  rmsubs(); % to force simplification
  svec:=simp!* symvec;
  for each eqngrp in alleqngrps() do
    dunknlinks:=addlinks(dunknsE eg!*eqn eqngrp,dunknlinks);
  if !*tracemkgens then <<write("found genfns ",genfns); terpri() >>;
  %
  % the generators from dunkns still in eqns (possibly with dependencies)
  %
  if not (svec=0) then for each link in dunknlinks do begin
    gen:=nil ./ 1 ;
    for each dunkn in genfns do 
      if dunknin dunkn member link and (coeff:=dtermq(dunkn,svec))
        then gen := addsq(gen,coeff); % the coeff here includes the dunkn etc
    if domainp denr gen then gen := (mynumr gen) ./ 1;
    if numr gen then eqngens := gen . eqngens;
    end;
  %
  % the generators from dunkns not in eqns but with dependencies
  %
  if svec neq 0 then for each dunkn in genfns do 
    if assoc(dunkn,depl!*) and not (dunkn member dunknsineqns) then begin
      gen := dtermq(dunkn,svec);
      if domainp denr gen then gen := (mynumr gen) ./ 1;
      if numr gen then infgens := gen . infgens;
      end;
  %
  % the generators from dunkns without dependencies,
  % except for dunkns still in eqns
  %
  if svec neq 0 then for each dunkn in genfns do 
    if not assoc(dunkn,depl!*) and  not member(dunkn,dunknsineqns)
      and numr (gen:=coeffq(dunkn,svec)) then begin
        if domainp denr gen then gen := (mynumr gen) ./ 1;
        fingens := gen . fingens;
        end;
  remprop('gen,'kvalue);
  while fingens do begin
    gen := car fingens;
    fingens := cdr fingens;
    mysetk( list('gen,(n:=n+1)),mk!*sq gen);
    end;
  while infgens do begin
    gen := car infgens;
    infgens := cdr infgens;
    mysetk( list('gen,(m:=m+1)+n),mk!*sq gen);
    end;
  while eqngens do begin
    gen := car eqngens;
    eqngens := cdr eqngens;
    mysetk( list('gen,(l:=l+1)+n+m),mk!*sq gen);
    end;
  !*gentypes:=list(n,m,l);
  showgens();
  mysetk(list('gen,0),(l+n+m));
  end;

symbolic operator mkgens;

symbolic procedure mkprgens;
begin
scalar dummy;
dummy:=symvec;
symvec:=prosymvec;
mkgens();
symvec:=dummy;
end;
symbolic operator mkprgens;

symbolic procedure check u;
% check that symmetry generator (!*sq) u is valid...
begin scalar upr,x,eqnlis,tsteqn;
  if !*predeqlis then readequations(); % needed if mkdets has not been called
  if !*p=0 or !*q=0 or not (fixp !*p and fixp !*q and fixp !*r)
    then <<write "*** Can't find original equations to verify symmetry";terpri();
           return>>;
  upr:=mysimpq simp!* prolong(!*r,u);  % this may make more prolongation coeffs than needed...
  for each tsteqn in !*deqlis do
    (if numr x then eqnlis:=(mk!*sq x).eqnlis)
       where x=mysimpq liesubq ((mynumr vecderf(upr, tsteqn)) ./ 1);
  if eqnlis then return 'list.eqnlis;
  end;

symbolic operator check;

algebraic procedure showgens;
% n is number of finite symms, ie const coeffs, not involved in any eqns
% m is number of other syms (with inf dim algebra)
begin scalar l,m,n,cond;
  if not !*gentypes and cdr !*gentypes and cddr !*gentypes then mkgens();
  n:=lisp car !*gentypes; m:=lisp cadr !*gentypes; l:=lisp caddr !*gentypes;
  write "There are ",n+m+l," symmetries found.";
  if n>0 then write "The generators of the finite algebra are:";
  for j:=1:n do write "Gen(",j,") = ",gen(j);
  if m>0 then write "The generators of the infinite algebra are:";
  for j:=(n+1):(n+m) do write "Gen (",j,") = ",gen(j);
  if l>0 then <<
    write ("The generators for the remaining equations are:");
    write ("(The unknowns in these generators satisfy the remaining determining equations.)")>>;
  for j:=(n+m+1):(n+m+l) do write "Gen (",j,") = ",gen(j);
  if (symbolic !*verifygens) 
    then for j:=1:n+m+l do if (cond:=check gen j) then begin 
      write "gen(",j,
            ") has the condition(s) that the following expressions be zero:";
      symbolic terpri(); write cond;
    end;
  showdepsused();
  end;
 
symbolic procedure showdeps();
begin scalar dep;
  terpri();
  for each dep in depl!* do if dunknk car dep
    then <<write(car dep," depends on ",cdr dep); terpri() >>;
  end;  

symbolic operator showdeps;

symbolic procedure showdepsused();
begin scalar dep,dunknsineqns,dunknsinsymvec;
  dunknsineqns:=dunknsineqns();
  dunknsinsymvec:=dunknsinsymvec();
  terpri();
  for each dep in depl!* do 
    if dunknk car dep 
      and ((car dep member dunknsinsymvec) or (car dep member dunknsineqns))
        then <<write(car dep," depends on ",cdr dep); terpri() >>;
  end;
 
symbolic operator showdepsused;
 
symbolic procedure dfdunknsineqns;
begin scalar eg, rm, dfdunknsineqns;
  for each eg in alleqngrps() do
    for each rm on eg!*eqn eg do
      dfdunknsineqns:=union(dfdunknsineqns,mvar rm);
  return dfdunknsineqns;
  end;
 
symbolic procedure dunknsineqns;
begin scalar eg, dunknsineqns;
  for each eg in alleqngrps() do
    dunknsineqns:=union(dunknsineqns,dunknsE eg!*eqn eg);
  return dunknsineqns;
  end;
  
symbolic procedure dunknsinsymvec; 
begin scalar svec,term,dunknsinsymvec;
  rmsubs(); % to force simplification
  svec:=simp!* symvec;
  symvec:=mk!*sq svec;
  for each term in numr svec do 
    dunknsinsymvec:=union(dunknsinsymvec,dunknsE cdr term);
  return dunknsinsymvec;
  end;

symbolic procedure dunknsinprosymvec; 
begin scalar svec,term,dunknsinsymvec;
  rmsubs(); % to force simplification
  svec:=simp!* prosymvec;
  prosymvec:=mk!*sq svec;
  for each term in numr svec do 
    dunknsinsymvec:=union(dunknsinsymvec,dunknsE cdr term);
  return dunknsinsymvec;
  end;

symbolic procedure dfdunknsinsymvec;
begin scalar svec,term,dfdunknsinsymvec;
  svec:=simp!* symvec;
  symvec:=mk!*sq svec;
  for each term in numr svec do
    dfdunknsinsymvec:=union(dfdunknsinsymvec,dfdunknsE cdr term);
  return dfdunknsinsymvec;
  end;


symbolic procedure stats;
begin
  showtime; terpri();
  write ("Total of ",!*dets-!*startdets," equations used, with ",!*ccount-!*startccount,
         " new arbitrary functions made."); terpri();
  terpri();
  write("successful operations were :",!*opusage);terpri();
  if !*intfac!*opusage then
    <<write("successful operations on op!*intfac eqns were :",!*intfac!*opusage); terpri()>>;
  if !*op!*xdp!*opusage then
    <<write("successful operations on op!*xdp eqns were :",!*op!*xdp!*opusage);terpri()>>;
  terpri();
  write "Variables used to split determining equations were ",!*splitvars;terpri();
  if !*tracedunknorder then showdunkns();
  !*startccount:=!*ccount;   % so that we show the stats relevant since the last call
  !*startdets:=!*dets;
  !*opusage:=!*intfac!*opusage:=!*op!*xdp!*opusage:=nil;
  end;

symbolic operator stats;

symbolic procedure showdunkns;
begin scalar dunkn;
  !*dunknsfound:=sort(!*dunknsfound,'dunknorder);
  write "ordered list dunkns found is"; terpri();
  for each dunkn in !*dunknsfound do write dunkn; terpri();
  end;

algebraic procedure showcomms;
begin scalar n;
  if not !*gentypes and cdr !*gentypes and cddr !*gentypes then mkgens();
  n:=(lisp car !*gentypes)+(lisp cadr !*gentypes)+(lisp caddr !*gentypes);
  for j:=1:n-1 do for k:=j+1:n do
    write("comm( gen(",j,"), gen(",k,") )= ",comm(gen j,gen k));
  end;

symbolic operator showcomms;

%******end of output

%********* The following code is borrowed from ExCalc (with permission).

algebraic operator partdf;
% changed from "operator partdf" 4/5/99
newtok '((!@) partdf);

symbolic procedure flatindxl u;
   for each j in u collect if atom j then j else cadr j;

symbolic procedure partdfprn u;
    if null !*nat then <<prin2!* '!@;
                         prin2!* "(";
                         if cddr u then inprint('!*comma!*,0,cdr u)
                          else maprin cadr u;
                         prin2!* ")" >>
     else begin scalar y; integer l;
            l := flatsizec flatindxl cdr u+1;
            if l>(linelength nil-spare!*)-posn!* then terpri!* t;
            % avoids breaking of the operator over a line;
            y := ycoord!*;
            prin2!* '!@;
            ycoord!* :=  y - if (null cddr u and indexvp cadr u) or
                                (cddr u and indexvp caddr u) then 2
                              else 1;
                if ycoord!*<ymin!* then ymin!* := ycoord!*;
                if null cddr u then <<maprin cadr u;
                                     ycoord!* := y>>
                 else <<for each j on cddr u do
                          <<maprin car j;
                            if cdr j then prin2!* " ">>;
                        ycoord!* := y;
                        if atom cadr u then prin2!* cadr u
                         else <<prin2!* "(";
                                maprin cadr u;
                                prin2!* ")">>>>
          end;

put('partdf,'prifn,'partdfprn);

symbolic procedure indexvp u;
   null atom u and flagp(car u,'indexvar);

%******end of excalc stuff


algebraic procedure readdets;
begin
  lisp(!*deteqns:=nil);
  adddets();
  end;

symbolic procedure adddets;
%
% Read the determining equations stored as values of the operator deteqn(j)
% for j>0. j<0 stores eqngrp numbers for those eqns...
%
% Place the list of newlyformed eqngrps onto the list !*deteqns
%
begin scalar kvals,kval,thiseqn,neweqn,oldnum;
  rmsubs();
  kvals:=get('deteqn,'kvalue);
  for each kval in kvals do
    if caar kval='deteqn and fixp cadar kval and cadar kval>0 then begin
      oldnum:=assoc(list('deteqn,-cadar kval),kvals);
      if oldnum then oldnum:=cadr oldnum;
      thiseqn:=cadr kval;
      neweqn:=mynumr simp!* thiseqn;
      if neweqn then 
        addeqngrp  mkeqngrp(neweqn,car kval,
                            list("readdets from eqngrp ",oldnum) 
                            );
      end
    else if not (caar kval='deteqn and fixp cadar kval)
      then <<write("Unable to read ",car kval); terpri() >>;
  end;

symbolic operator adddets;

%*******end of readdets

%***** just some utility functions


symbolic procedure mynumr u;
if not (u=0) then <<dividelog denr u; numr u>>;

symbolic procedure dividelog u;
% keep a record of interesting divides
% and return the list of those found this time
if u and (not domainp u) and not !*ignoreall
  and ( ((!*freeunknownatoms or !*freeunknownops) and not(u member !*divides) 
         and (funknsinf u) and noticedivf u)
       or (!*logspdiv and noticedivf u)            )
  then begin
    scalar factors,fac,!*ezgcd,thesedivides;
    if !*myezgcd then !*ezgcd:=t;
    if !*factordivides
      then <<if !*tracecute or !*traceslow then write "~",count u;%,"%",mytermsf u;
             factors:=fctrf u;
             if !*tracecute or !*traceslow then write "~";>>
      else factors:= 1 . list(u.1);
    %
    %   **********************************
    %   fctrf is a reduce 3.4 function
    %   not supported in reduce 3.3 ???
    %   It returns a factored form, not supported in reduce 3.3 ???
    %   **********************************
    %
    for each fac in cdr factors do 
      if not (car fac member !*divides) 
        and ((!*logspdiv and noticedivf car fac) or funknsinf car fac) then begin 
          write(tmstmp(),"*** free or special functions found when dividing by ");terpri();
          showf car fac; terpri();
          !*divides:=(car fac).!*divides;
          end;
    for each fac in cdr factors do 
      if not (car fac member thesedivides) 
        and ((!*logspdiv and noticedivf u) or funknsinf car fac)
          then thesedivides:=(car fac).thesedivides;
    return thesedivides;   % for keeping track of divides related to eqngrps, when implemented...
    end;

symbolic procedure showdivides;
for each u in reverse !*divides do begin
  write("Free or special functions found when dividing by ");
  terpri();
  showf u; terpri();
  end;

symbolic operator showdivides;

symbolic procedure showsimpdivides;
for each u in reverse !*divides do begin
  write("Free or special functions found when dividing by ");
  terpri();
  showf u; 
  write "which simplifies to";terpri();
  showsq mysimpf u;
  end;

symbolic operator showsimpdivides;

symbolic procedure noticedivf u;
% should we notice standardform u as interesting?
u and not domainp u and 
  (noticedivk mvar u or noticedivf lc u or noticedivf red u);
 
symbolic procedure noticedivk u;
% should we notice kernel u as interesting?
% yes if it is not just an x or u variable, 
not (pairp u and (car u=!*u or car u=!*x));
%and (funknsink u 
%  or (pairp u and (cddr u or pairp cadr u)));

symbolic procedure funknsinf u;
% are there any funkns in sf u? 
u and not domainp u and
(funknsink mvar u or funknsinf lc u or funknsinf red u);

symbolic procedure funknsink u;
% are there any funkns in u?
funknk u or (pairp u and (funknsink car u or funknsink cdr u));

symbolic procedure complement(u,v);
% returns the complement of v in u, ie elements of u not in v
begin scalar comp,w;
  for each w in u do if not member(w,v) then comp:=w.comp;
  return comp
  end;

symbolic procedure equalset(u,v);
% are sets u and v equal (but maybe ordered differently)?
null complement(u,v) and null complement(v,u);  

symbolic procedure intersection(u,v);
% returns the intersection of sets u and v
% in the order of u!
if u then begin scalar w,x;
  for each x in u do if x member v then w:=x.w;
  return reverse w;
  end;

symbolic procedure intersect(u,v);
% predicate, does  u intersect v?
% marginaly faster than intersection.
begin scalar w;
  while u do
    if (car u) member v
      then <<w:=t; u:=nil>>
      else u:=cdr u;
  return w;
  end;

symbolic procedure mysetk(u,v);
% set the value of u to v,
% car u is the name of an operator
% v is a prefix value
begin scalar su,propu,sqval;
  rmsubs();
  sqval:=simp!* v;
  if edunknk u then begin
    su:=numr simp!* u;
    if null su or  (propu:=mvar su) neq u
      then <<terpri();terpri(); write "*** value for ",u," already set ???";terpri();terpri()>>;
    % This next check is an important check, yet is costly in simping symvec.
    % However, it should never occur, so only test if we are being fussy. :-)
    if !*verifymysetk and not (u member dunknsinsymvec()) and null cddr u then begin
      terpri();terpri();
      write "*** Setting value for ",u," which does not appear in symvec ???";terpri();
      write " This is could be a very dangerous thing to do...";terpri();
      terpri();
      end;
    if assoc(u,depl!*) and v=0
      then for each var in cdr assoc(u,depl!*) do depend1(u,var,nil);
    if assoc(u,!*recnumdeps) % and v=0
      then !*recnumdeps:=delete(assoc(u,!*recnumdeps),!*recnumdeps);
    if !*tracesetk then begin
      if !*tracehard then <<write ("setting ",u," to"); terpri();
                            write v; terpri() >>;
      if !*tracesoft then algebraic (write ("setting ",u," to ",v));
      end;
    if !*tracecute then write ("=",mytermsf numr sqval, "/",mytermsf denr sqval," ");
                              % "(",mytermsf numr symvec,"/",mytermsf denr symvec,")");
    end;
  setk1(u,v,t);
  if edunknk u then showdeplevs();
  if !*simpsymvecwithsub then symvec:=mk!*sq simp!* symvec;
  if !*tracesymvec then algebraic( write "symvec := ",symvec );  
  dividelog denr simp!* v;
  end;

symbolic procedure mycleark u;
% clear the value of u,
% car u is the name of an operator
% u should not be used...
if not !*holdvalues then begin scalar var;
  if assoc(u,depl!*)
    then for each var in cdr assoc(u,depl!*) do depend1(u,var,nil);
    if assoc(u,!*recnumdeps) 
      then !*recnumdeps:=delete(assoc(u,!*recnumdeps),!*recnumdeps);
  if !*tracecleark then <<
    if !*tracehard then <<write ("clearing ",u); terpri() >>;
    if !*tracesoft then algebraic write("clearing ",u)
    >>;
  setk1(u,nil,nil);
  end;

symbolic procedure mysimpq u;
<<rmsubs(); simp!* list('!*sq,u,nil)>>;

symbolic procedure mysimpf u;
<<rmsubs(); simp!* list('!*sq, u ./ 1, nil)>>;

symbolic procedure count v;
% count the number of elements in list v
%if pairp v then 1+count cdr v else 0;
if null v then 0
else begin integer i;
  while pairp v do <<i:=i+1;v:=cdr v>>;
  return i;
  end;


%******end of util

symbolic procedure tmstmp;
% used as an argument to write statements
% Show !*dets, as a time measure
% Show count !*deteqns for interest
if !*tracetmstmp
  then compress append(append('!".'!<.explode !*dets,
                              '! .explode (count !*deteqns)),
                              list('!>,'!"))
  else " ";

symbolic procedure cute u;
if !*tracecute then 
  if !*tracetmstmp 
    then write("<",!*dets," ",count !*deteqns,">",u)
    else write u;

symbolic procedure cute1 u;
if !*tracecute then write u;

symbolic procedure mywrite1 u;
<<write u;terpri()>>;

rlistat '(cute cute1 mywrite1);

% The following code is included with permission from the author

module intpatch;  % Integrate dependent variables & rational powers

% Author: Francis J. Wright <F.J.Wright@QMW.ac.uk>
% Date: 19 June 1992  (Very minor update 16/11/92)
% REDUCE version: 3.4 or 3.4.1; PSL

% This version (19 June 1992) fixes a bug that integrals that remained
% symbolic were not returned as unique kernels - I now use !*kk2q.
% This bug caused rather obscure symptoms, such as failures with on
% factor.

% DOCUMENTATION
% =============
% This patch has two separate functions:

% 1: It allow integrals containing IMPLICITLY dependent variables,
% created using the DEPEND command, to remain symbolic rather than
% cause an error, whilst preserving other error handling as normal.
% ON FAILHARD turns this facility off.
% This facility was developed from a patch by Herbert Melenk,
% which this patch is intended to replace.

% 2: It integrates simple rational powers of the integration
% variable that the integrator currently fails to integrate.


% This file is intended to be compiled using faslout and then loaded
% using load_package, but alternatively the source code can be input
% when required and either compiled or interpreted.  It could also be
% edited into the standard source file INT.RED and built permanently
% into REDUCE.

% Because of the autoload definitions for int (see ENTRY.RED) it is
% safe to load this module before the integrator has been loaded, and
% it does not cause the integrator to load.  The integrator can safely
% be AUTOLOADED by calling INT afterwards, either from algebraic mode
% or from symbolic mode using any of the methods described below, but
% MUST NOT BE EXPLICITLY LOADED AFTER this patch, otherwise this patch
% will not work because the simpfn property of int will get
% overwritten.  I cannot see any solution to this, other than avoiding
% doing it!

% WARNING: To call the integrator from symbolic mode when this patch
% is loaded DO NOT CALL SIMPINT DIRECTLY.  Instead, there are several
% possibilities that should work, including any of the following:

% 1.  Use the form ????{'int, integrand, int_var}, where ???? is
% reval, aeval or simp, etc, depending on the form of result you want,
% and the arguments are in prefix (or pseudo-prefix) form.  This is
% equivalent to calling the integrator from algebraic mode, and I
% currently think this is probably the best approach.

% 2.  Call simpintpatch directly, but this is not general and so not
% recommended.

% 3.  Replace calls of "simpint" by (say) "integrator", and do
% essentially the following:
% global'(integration_function);
% symbolic(integration_function := get('int, 'simpfn));
%
% symbolic procedure integrator u;
%    apply1(integration_function, u);
% Note: apply1 is defined in RLISP.RED.


fluid '(!*failhard);

%% fluid '(soft!-rerror!-number);
% Hope not necessary - it seems not to be.

put('int, 'simpfn, 'SimpIntPatch);

symbolic procedure SimpIntPatch u;
   % Driver for various patches:
   % 1: Catch errors from SimpInt, trap error number 7 only,
   % and pass on all other errors as normal hard REDUCE errors.
   % 2: Post-process unintegrated rational powers.
   begin scalar r, !*redefmsg, !*uncached; !*uncached := t;
      % !*redefmsg rebound to avoid PSL messages
      % about redefinition of rerror.
      %% integer soft!-rerror!-number;    % defaults to 0, not nil
      put('int, 'simpfn, 'SimpInt);       % assumed & reset by SimpInt
      copyd('rerror, 'intpatch!-rerror);  % redefine rerror
      r := errorset!*({'SimpInt, mkquote u}, nil);
      copyd('rerror, 'rerror!%);          % restore rerror
      put('int, 'simpfn, 'SimpIntPatch);  % reset INT interface
      if pairp r then <<
         % First call of SimpInt succeeded -
         % try to reprocess any integrals left:
         put('int, 'simpfn, 'SimpIntRatPow);
         u := resimp car r;  % this works ONLY with !*uncached := t;
         put('int, 'simpfn, 'SimpIntPatch);
         return u
      >>
      else if !*failhard or not(r eq 7) then
         rederr EMSG!*              % Error failure
      else return !*kk2q('int . u)  % Remain symbolic
   end;


% Integrator error trap patch to allow controlled error handling
% ==============================================================

% The error numbers generated by SIMPINT and the corresponding
% error message that would be output by INT are the following,
% collected from the INT source code:

%  1  =  "Improper number of arguments to INT"
%  2  =  "Improper number of arguments to INT"
%  3  =  "Too many arguments to INT"
%  4  =  "FAILHARD switch set"
%  5  =  "Invalid polynomial in int-quadterm"
%  6  =  "Empty list to mapply"
%  7  =  "Can't integrate in the presence of side-relations" (TRAPPED)
%  8  =  "Invalid exponent"
%  9  =  "FAILHARD switch set"

% If any other error number, such as 0, should occur then it
% corresponds to some other non-specific error.


symbolic procedure rerror!%(message);
   % This is precisely the definition of rerror in RLISP.RED,
   % but redefining it here makes sure it is loaded,
   % and also avoids the need to save it.
   % Precisely this procedure is also defined in SOFTSOLV.
   rederr message;

symbolic procedure intpatch!-rerror(number,message);
   %%   << soft!-rerror!-number := number; error1() >>;
   % The following will suffice provided errorsets
   % are not nested in the integrator.
   % It makes error message text available in EMSG!*.
   error(number, message);


% Integrator postprocessor patch to integrate simple rational
% powers that the integrator currently fails to integrate.
% =======================================================

symbolic procedure SimpIntRatPow u;  % u = (integ var)
   % Integrate integrands of the form var**(m/n),
   % which the integrator leaves in a rather bizarre form -
   % hence the precise form of the following code.
   % Returns original integral if it has the wrong form.
   begin scalar integ, var, power;
      integ := car u;  var := cadr u;
         % assumes true prefix forms, already evaluated by SimpInt.
%     power := errorset!*(
%        {'FindRatPow, mkquote integ, mkquote var}, nil);
%     errorset!*(u,v) == errorset(u,v,!*backtrace)
%     Backtrace from this is unlikely to be interesting, so ...
      power := errorset(
         {'FindRatPow, mkquote integ, mkquote var}, nil, nil);
      if errorp power then return !*kk2q('int . u);
      power := car power;  % correct form of integrand found.
      % integrand = var**power, so return integral:
      power := reval {'plus, power, 1};
      return simp!* {'quotient, {'expt, var, power}, power}
   end;

symbolic procedure FindRatPow(monom, var);
   % Return power of a monomial in var, as a
   % rational number in UNSIMPLIFIED prefix form
   % or cause error return to enclosing errorset.
   if eqcar(monom, 'quotient) then
      {'plus, FindRatPow(cadr monom, var), 
         {'minus, FindRatPow(caddr monom, var)}}
   else if eqcar(monom, 'times) then
      'plus . for each el in cdr monom collect FindRatPow1(el, var)
   else FindRatPow1(monom, var);

symbolic procedure FindRatPow1(monom, var);
   if monom eq 1 then 0  % only possible constant by linearity
   else if monom = var then 1
   else if eqcar(monom, 'expt) and
      cadr monom = var then caddr monom
   else error1();  % wrong form

endmodule;

% end of intpatch
%******************************************

symbolic procedure mytermsf u;
if null u or domainp u then 0
else 1+(mytermsf lc u)+(mytermsf red u); 

endmodule;

end;

