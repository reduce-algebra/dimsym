module dim2ode;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% module:                    DIM2ODE                         %
%                                                            %
%           Support for dimsym <-> odesolve interface        %
%                                                            %
% Author:   Michael Jerie, October 1999.                     %
% Last modified: February 6th, 2004                                                           %
%                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
comment on this file

1.
??? appears where there is an important issue to be addressed

2.
mask dependencies in the case where deteqn is like
   Af''' + ... + Bf' + C where A,B,C are fns of x and
f is a fn of x,y but only diffed with respect to x in this det. 

3.
add in some tracing functions

4.
must put op!*odesolve on the solveops list or whatever it is.

5.
use a global var !*odehardint as the ode counterpart of !*hardint
and put in a showhardint() function for displaying unevaluated
integrals at the end of a solvedets call.;
% end comment

symbolic procedure oderestnointvar(eqn,dunkn,intvar);
% is dunkn the only dunkn in eqn involving intvar?
null eqn or 
  (  (not member(intvar,mvar eqn)
      or 
      if dunknin mvar eqn=dunkn then null 
            complement(alldepsf !*k2f mvar eqn,list(intvar))
      or
      mvar eqn = dunkn
      )
   and restnointvar(red eqn,dunkn,intvar) );

symbolic procedure oderestnodunkn(eqn,dunkn,intvar);
% return t if dunkn is diffed by intvar only or undiffed.
null eqn or ( (if dunknin mvar eqn = dunkn then
                 if car mvar eqn ='df then
                     if null (caddr mvar eqn = intvar) then 'nil
                       else 't
                   else 't
                 else 't)
              and oderestnodunkn(red eqn,dunkn,intvar)
             );

symbolic procedure findodedunkn eqn;
% return the dunkn for which eqn is an ode else nil
% check the mvars of the standard terms in the sf eqn to see
% if eqn is an ode for some dunkn in the mvars
begin scalar resteqn,tdunkn,alldeps;
  alldeps:=alldepsf eqn;
  resteqn:=eqn;
  while resteqn and not tdunkn do begin
    if (car mvar resteqn)='df
      %% mvar resteqn looks like (df dunkn var deg)
      % rule out mixed deriv's here
      and (  ( length mvar resteqn = 4 and numberp cadddr mvar resteqn)
              or
              length mvar resteqn = 3   )
      %% note dunknin mvar resteqn might depend on more than one var
      % rule out fn of more than one var ??? mask this in future ???
      and not ordp(length depsk dunknin mvar resteqn,2)
      % candidate for splitting
      and null complement(alldeps,depsk mvar resteqn)
      % no other dunkn diffed wrt the intvar 
      and oderestnointvar(eqn,dunknin mvar resteqn,caddr mvar resteqn)
      % dunkn not diffed by anything else (REDUNDANT unless we do the ??? mask above)
      and oderestnodunkn(eqn,dunknin mvar resteqn,caddr mvar resteqn)
      %% we don't want only constant coeffs
      %and intvarnocoeff(eqn,caddr mvar resteqn)
      %% possible to have a freeunknown in the coeffs depending on the intvar + others?
      %% use funknsp or funkndepsp to test?
      %and nofunknceoffdeps(eqn,caddr mvar resteqn)
      and not ((caddr mvar resteqn) member otherboundvarsE(dunknin mvar resteqn,eqn))
        then tdunkn:=mvar resteqn
        else resteqn:=cdr resteqn;
      end;
  return tdunkn;
  end;

symbolic procedure mkdummyid u;
% returns a single pair. intern is
% needed so the system recognises the new id.
 u.intern gensym(); 

symbolic procedure ds2osetdepssublist sublist;
% fix dependencies, only including those deps on depl!*.  Any others??
begin scalar deps,pair;
  for each pair in sublist do
    if (deps:=assoc(car pair,depl!*)) then begin
      if ordp(length cdr deps, 2) then 
        write "*** WARNING: more than one dependency in ode in ds2osetdepssublist",terpri()
      % we assume here that cdr pair is an id         
      else depend1(cdr pair,cdr assoc(cadr deps,sublist) ,t);
    end;
  end;

symbolic procedure ds2osbtrcttail(u,v);
% v is the tail of u. Return u minus the tail.
begin scalar tempu,tempv;
  if null v then return u;
  if null u then return nil;
  tempu:=reverse u;
  tempv:=reverse v;
  for each x in tempv do if car tempu = x then tempu:=cdr tempu
    else write "*** WARNING: bad tail in ds2osbtrcttail";
  return reverse tempu;
 end;
%trst ds2osbtrcttail;

symbolic procedure ds2osubvar(u,pair);
%returns u with occurences of (car pair) replaced with (cdr pair).
if null u then nil else
  begin scalar x,y,ans;
    if (y:=member(car pair, u)) then  
      <<x:=append(ds2osbtrcttail(u,y),
           append(if (not pairp cdr pair) then list(cdr pair)
                    else cdr pair
                  ,cdr y));
        ans:=ds2osubvar(x,pair)>>
    else return u;
    return ans;
  end;
%trst ds2osubvar;

symbolic procedure ds2osubeqn(eqn,sublist);
% return eqn with each occurence of car of each pair in sublist
% replaced with the cdr of that pair.  eqn is a sf so eqn can
% be either a (nested) list/pair, a number, an id, or null !!!.
if null eqn then nil 
 else if idp eqn then eqn 
 else if numberp eqn then eqn
 else
    begin scalar neweqn,y,x,u;
      neweqn:=eqn;
      %for each u in eqn do 
      while (eqn and pairp eqn) do begin
        u:=car eqn; eqn:= cdr eqn;    
        % u is an id, number or list
        if idp u then nil 
         else if numberp u then nil
         else if pairp u then
            % y is a pair
            if (y:=assoc(u,sublist)) then neweqn:=ds2osubvar(neweqn,y)
            else <<if (x:=ds2osubeqn(u,sublist))=u then nil 
                   else neweqn:=ds2osubvar(neweqn,list(u,x))>>;
      end;
    return neweqn;
   end;
%trst ds2osubeqn;

symbolic procedure dimsym2odesolve eqn;
% eqn is a sf.  Returns list(neweqn, dunkn, indvar, sublist).
% neweqn is nil if eqn is not an ode (according to findodedunkneqn),
% otherwise neweqn is eqn with new id's (with correct deps) subbed
% for the vars.  dunkn is the dependent var, indvar is the independent var.
% sublist is a list of pairs ( (oldvar. new id) ... )
begin scalar odedunkn,dunkn,var,dunknsanddepslist,sublist,neweqn;
  % At the moment findodedunkn only returns an odedunkn
  % if there is dependence on only one variable in eqn.  ? Too restrictive ?.
  % odedunkn is the possibly diff'd dunkn for which eqn is an ode.
  odedunkn:=findodedunkn eqn;
  if null odedunkn then return list(nil, nil, nil, nil);
  % mj
  %write "dimsym2odesolve:  odedunkn is ",odedunkn;terpri();
  dunkn:=dunknin odedunkn;
  var:=caddr odedunkn;
  %% Now we need to make a new equation for odesolve call.
  %% Need to replace dunkns c(i) by identifiers otherwise
  %% odesolve gives an error.  Any other c(i) appart from 
  %% odedunkn will be constants?
  dunknsanddepslist:=append(dunknsE eqn,alldepsf eqn);
  % sublist is list of pairs (dunkn or dep . new identifier)
  % make initial sublist without setting deps
  while dunknsanddepslist do begin
     sublist:=mkdummyid car dunknsanddepslist. sublist;
     dunknsanddepslist:=cdr dunknsanddepslist;
     end;
  % now set dependencies
  ds2osetdepssublist sublist;
  % we don't need to do this do we?
  sublist:=reverse sublist;
  neweqn:=ds2osubeqn(eqn,sublist);
  neweqn:=numr simp prepf neweqn;
  return list(neweqn,cdr assoc(dunkn,sublist),
               cdr assoc(var,sublist), sublist);
 end;

%% this declaration is now made in the header module dimsym22p
%% flag for already loaded package and package not found.
%%fluid '(!*loadodesolve !*noodepackage);

symbolic procedure myfilep u;
% u=file Check if "file.b" is in any of the loaddirectories!*  
begin scalar temp,found;
  temp:=loaddirectories!*;
  while temp and null found do
    begin scalar y;
       if null explode2 car temp then <<temp:=cdr temp;return>>;
       y:=reverse cdr reverse explode car temp;
       found:=filep compress append(y, append(explode u,'(!. b !" )) );
       temp:=cdr temp;
     end;
  return found;  
 end;

symbolic procedure myloadodesolve;
% load the package odesolve if it hasn't already been loaded
 if null !*loadodesolve then
     %<<algebraic load odesolve; !*loadodesolve:='t >>;
     begin scalar versionlist;
        % versionlist is a list of releases of the odesolve package
        % currently there are:
        % odesolve: the original package, and currently released with
        %           reduce 3.7 
        % odesolv1: updated but unofficial release, although how 
        %           much improvement to the linear ode solving has
        %           been done I'm not sure.
        versionlist:='(odesolv1 odesolve);
        while null !*loadodesolve and versionlist do 
          begin scalar x;
             x:=car versionlist;
             versionlist:=cdr versionlist;
             if myfilep x % need to do 
                          % evload(list('quote,list u))
                          % for module details but Unknown sideeffects
                          % from doing this. So at the moment if the
                          % system can't find a module then there will
                          % be a fatal error.
                          %%and  begin scalar y,found;
                          %%       y:=get(x,'package);
                          %%       while y and found do
                          %%         << found:=myfilep car y;
                          %%            y:=cdr y>>;
                          %%      end
                   then << load!-package x;
                           !*loadodesolve:='t;
                           terpri();
                           write "*** op!*odesolve has successfully loaded ",x;
                           terpri();terpri()>>;
           end;
          if null !*loadodesolve then 
             <<!*noodepackage:='t;
               terpri();
               write "*** op!*odesolve: Unable to load the package odesolve";
               terpri();
               write "    I have checked the load directories ... ";terpri();
               for each x in loaddirectories!* do <<if null explode2 x then nil 
                                                      else <<write "      ",x; terpri()>>
                                                     >>;
               write "*** The solvedets algorithm will continue without op!*odesolve";
               terpri();terpri()
               >>;
      end;  
symbolic operator myloadodesolve;

%% this declaration is now made in the header module dimsym22p
%%fluid '(!*odehardints);

symbolic procedure showodehardints;
% display ints on !*odehardints
if !*odehardints then begin scalar ints,hardint;
    ints:=!*odehardints;
    while ints do <<
      hardint:=car ints;
      algebraic write hardint;
      ints:=cdr ints>>;
    end;

symbolic operator showodehardints;
%rlistat '(showodehardints);
%equiv: put('showodehardints, 'stat, 'rlist);

symbolic procedure myodesolve(ode,dunkn,indvar);
% entry/exit to the package odesolve.
% Returns solution, where solution an sq and
% is the rhs of "y=solution" returned by odesolve, otherwise nil
% if odesolve is unsuccessful.
% ??? Need to reset flags after odesolve call ???
begin scalar odesoln,u,v,rhs;
  % odesoln is a prefix form: (list (equal lhs rhs))
  odesoln:=odesolve(prepf numr simp ode,dunkn,indvar);
  
  % ??? reset flags - No lingering effect from using odesolve ???
  
  % now check if odesolve succeeded.
  % Must have:
  % 1 explicit solution "y = function of x" ie. no implicit or parametric solns
  % 2 no unevaluated integrals in the rhs.
  % 3 The odesolve package will return multiple "soln's" in a
  % list, which we wont allow (perhaps can't happen for our 
  % linear deteqns anyway?).
  if ordp(length(odesoln),3) then return; % no multiple solutions.
  % car u must be either 'list or 'equal
  if car (u:=cadr odesoln) = 'equal then 
      if cadr u = dunkn then rhs:=simp caddr u %rhs is explicit soln
        else return % implicit solution
    else return; % parametric solution
  if odedunknleft(rhs,dunkn) then return; % just to be sure not implicit soln
  % list of psuedo prefix unevaluated integrals
  if (v:=odeintleft rhs) then
      begin 
        !*odehardints:=append(v,!*odehardints);
        for each x in v do <<terpri();
           write "op!*odesolve failed : unable to integrate ";
           terpri(); mathprint prepsq simp x; terpri()>>;
      end
  else return rhs;
  end;

symbolic procedure odeintleft u;
% return ((!*sq (int ... ) t) ... ) if 'int is in sq u
  if pairp u then 
     begin scalar ints;
      if car u ='int then ints:=mk!*sq simp u.ints;
      return append(odeintleft(car u),
                   append(odeintleft(cdr u), ints));
     end;

symbolic procedure removenils u;
% remove any 'nil from the middle of a list
begin scalar temp;
   while u do <<if null (car u = 'nil) then temp:=car u. temp; u:=cdr u>>;
   return reverse temp;
  end;

symbolic procedure odedunknleft(u,dunkn);
% return t if dunkn is in sq u
% dunkn is a gensym identifier
if  pairp u then 
    if  car u = dunkn then 't
    else odedunknleft(car u,dunkn) or odedunknleft(cdr u,dunkn);

symbolic procedure op!*odesolve eqngrp;
%
% solve eqns which are odes using the package 
% "odesolve" by Francis. J. Wright.
%
begin scalar odesoln,soln,odeeqngrp,eqn,odedunkn;
  if !*noodepackage then return;
  eqn:=eg!*eqn eqngrp;
  if null eqn then return;
  % test if eqn is ode. If it is then dimsym2odesolve 
  % returns list(neweqn, dunkn, indvar, sublist)
  %terpri();write "entered op!*odesolve with eqn = ",mathprint prepf eqn;
  %terpri();
  if null car (odeeqngrp:=dimsym2odesolve eqn) then return; % not an ode
  myloadodesolve();  % is an ode (odesoln is sq or nil)
  odesoln:= myodesolve(prepf car odeeqngrp,cadr odeeqngrp,caddr odeeqngrp);
  odedunkn:=car assoc2(cadr odeeqngrp,cadddr odeeqngrp); % the dunkn we hope to solve for
  if null odesoln then return  % unable to solve
    else  % in this case we have solution and we must sub back to dimsym form.
       begin scalar odearbconstsublist; 
          % list of arbconst(i) in odesoln
          odearbconstsublist:= odearbleft odesoln;
          soln:=odesolve2dimsym(odesoln,odearbconstsublist);
          % odesoln is an sq, therefore so is soln (unsimped)
          soln:=odesolve2dimsym(soln,cadddr odeeqngrp);
         soln:=mk!*sq simp list('!*sq, soln, 'nil); % ??? DON'T FORGET THE DIVIDES.
         write "op!*odesolve has solved eqn ",eg!*no eqngrp," :";
         terpri();
         writepri(mkquote odedunkn,'first);
         writepri(" = ",'nil);
         writepri(mkquote prepsq simp soln,'last);terpri();
         if null (denr simp soln = 1) then <<
           write "*** WARNING: unrecorded divide in op!*odesolve. ";
           terpri();
           write "             Here is the denominator: "; terpri();
           mathprint prepf denr simp soln; terpri()>>;
       end; 
  deathcert(eqngrp,list("Solved by op!*odesolve: ",odedunkn,
                        " set to ",soln));
  mysetk(odedunkn,soln);
 return list(t,odedunkn); 
 end;

symbolic procedure odearbleft u;
% look for (arbconst i) in (u =)odesoln
% return a list of pairs ( ((arbconst i). newconst()) ...)
% ready for subbing.
begin scalar arblist,uarblist;
  arblist:=arbsinE u;
  while arblist do <<if member(car arblist,cdr arblist) then nil
                       else uarblist:=car arblist.uarblist;
                     arblist:=cdr arblist>>;
  while uarblist do <<arblist:=list(newconst(),car uarblist).arblist;
                      uarblist:=cdr uarblist>>;
  return arblist;
 end;

symbolic procedure arbsinE u;
% list of every use or (arbconst i) including repeats.
    if pairp u then
       if car u='arbconst then list(u)
          else append(arbsinE(car u),arbsinE(cdr u));

symbolic procedure odesolve2dimsym(odesoln,sublist);
% odesoln is an sq !!! BE CAREFULL
% returns the solution found by odesolve with the old
% dimsym dunkns and vars subbed back in for the gensym()'s
begin scalar subbacklist;
    while sublist do <<subbacklist:=reversepair car sublist. subbacklist;
                       sublist:=cdr sublist>>;
    subbacklist:=reverse subbacklist;
    return o2dssubeqn(odesoln,subbacklist);
   end;

symbolic procedure reversepair u;
% reverse (a.b) -> (b.a)
  if pairp u then 
      if atom cdr u then cdr u.list(car u)
      else cadr u.list(car u);

symbolic procedure o2dssubeqn(eqn,sublist);
% return eqn with each occurence of car of each pair in sublist
% replaced with the cdr of that pair.  eqn is a sf so eqn can
% be either a (nested) list/pair, a number, an id, or null !!!.
if null eqn then nil else
    begin scalar neweqn,y,x,u;
      neweqn:=eqn;
      %for each u in eqn do 
      while (eqn and pairp eqn) do begin
        u:=car eqn;
        eqn:= if (atom cdr eqn and null cdr eqn ='nil)
                   then list(cdr eqn) 
                else cdr eqn;
        % u is an id, number or list
        % y is a pair (or list)
        if (y:=assoc(u,sublist)) then neweqn:=ds2osubvar(neweqn,y)
        else <<if (x:=o2dssubeqn(u,sublist))=u then nil 
              else neweqn:=ds2osubvar(neweqn,list(u,x))>>;
      end;
    return neweqn;
   end;
%trst ds2osubeqn;

!*solvedetsalgorithms:='(simp   proxpnd  stdsplit  split    alansplit  std  
                         odeslv  %% Added this one in for op!*odesolve trial !!
                         std1   easy     nosub     stdform  stdform1   nofxp 
                         hidef  findf    hidel     findl    inttbt     rmvred
                         new    pro2sf   proeqns   standardform  stdformx);

put('odeslv,'solveops,'(op!*shr1tm   op!*splitEc  op!*simpeq     op!*proexp
                      op!*slvtwo   op!*get1tm   op!*hidefr     op!*hidelg
                      op!*exp1tm   op!*slvspl   
                      op!*trgexp     
                      op!*intfac 
                      op!*odesolve % HERE IT IS  
                      op!*intslv   op!*xdpxpd     op!*slvall
                      op!*splitEd  %addintcon1   op!*slvalldfsearly
                      op!*sub2sf   addintcons
                      op!*findlg   op!*findfr   op!*slvalldfs  ));

flag('(op!*odesolve),'onebyoneop);

flag('(op!*odesolve),'restart);

endmodule;

end; % of file
