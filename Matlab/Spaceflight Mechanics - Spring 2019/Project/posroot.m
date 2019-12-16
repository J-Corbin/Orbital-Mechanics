function x = posroot(Roots)
%{ 
  This subfunction extracts the positive real roots from
  those obtained in the call to MATLAB's 'roots' function.
  If there is more than one positive root, the user is
  prompted to select the one to use.
 
  x         - the determined or selected positive root
  Roots     - the vector of roots of a polynomial 
  posroots  - vector of positive roots

  User M-functions required: none
%}
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~

%...Construct the vector of positive real roots:
posroots  = Roots(find(Roots>0 & ~imag(Roots))); 
npositive = length(posroots);

%...Exit if no positive roots exist:
if npositive == 0
    fprintf('\n\n ** There are no positive roots.  \n\n')
    return
end

%...If there is more than one positive root, output the
%   roots to the command window and prompt the user to
%   select which one to use:
if npositive == 1
    x = posroots;
else
    fprintf('\n\n ** There are two or more positive roots.\n')
    for i = 1:npositive
        fprintf('\n root #%g = %g',i,posroots(i))
    end
    fprintf('\n\n Make a choice:\n')
    nchoice = 0;
    while nchoice < 1 | nchoice > npositive
        nchoice = input(' Use root #? ');
    end
    x = posroots(nchoice);
    fprintf('\n We will use %g .\n', x)
end
end