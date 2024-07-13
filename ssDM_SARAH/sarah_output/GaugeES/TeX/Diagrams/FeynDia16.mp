% Diagrams/FeynDia16.mp -- do not edit, generated automatically by SSDMv2-GaugeES.tex
input feynmp
require_RCS_revision "1.30";
beginchar(1, 75*1.0pt#, 75*1.0pt#, 0);
"feynmf: 1";
LaTeX_unitlength:=1.0pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
vleft (__l1);
vright (__r1, __r2);
vconnect ("scalar", __l1, __v1);
vconnect ("scalar", __v1, __r1);
vconnect ("wiggly", __r2, __v1);
vlabel ("$H^0$", __l1);
vlabel ("$H^-$", __r1);
vlabel ("$W_{{k \mu }}$", __r2);
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
