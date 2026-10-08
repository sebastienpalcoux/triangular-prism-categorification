# Independent database audit; no internet access and no floating point.
Print("GAP ",GAPInfo.Version,"\n");
for o in [1..128] do
 for i in [1..NrSmallGroups(o)] do
  G:=SmallGroup(o,i);; C:=CharacterTable(G);; irr:=Irr(C);; ind:=Indicator(C,2);;
  for j in [1..Length(irr)] do
   if ind[j]=-1 then
    for a in [1..Length(irr)] do
     b:=Position(irr,ComplexConjugate(irr[a]));;
     mult:=ScalarProduct(irr[a]*irr[b],irr[j]);;
     if mult>0 and IsEvenInt(mult) then Print([o,i,j,a,b,mult],"\n");fi;
    od;
   fi;
  od;
 od;
od;
QUIT;
