# Check the nonabelian finite simple groups of order at most 126000.
# Enumeration: GAP's SimpleGroupsIterator; character tables: CTblLib.
if LoadPackage("CTblLib")<>true then Error("CTblLib required");fi;
Print("GAP ",GAPInfo.Version,"; CTblLib ",PackageInfo("CTblLib")[1].Version,"\n");
Witnesses:=function(t)
 local irr,ind,w,j,a,m;
 irr:=Irr(t); ind:=Indicator(t,2); w:=[];
 for j in [1..Length(irr)] do
  if ind[j]=-1 then
   for a in [1..Length(irr)] do
    m:=ScalarProduct(irr[a]*ComplexConjugate(irr[a]),irr[j]);
    if m>0 then Add(w,[j,a,m]);fi;
   od;
  fi;
 od;
 return w;
end;;
numberChecked:=0;; positiveTables:=[];;
for g in SimpleGroupsIterator(1,126000) do
 name:=IsomorphismTypeInfoFiniteSimpleGroup(g).shortname;;
 libname:=name;
 t:=CharacterTable(libname);;
 if t=fail then Error("Missing table: ",name);fi;
 if Size(t)<>Size(g) or not IsSimple(t) then Error("Group/table mismatch");fi;
 w:=Witnesses(t);; numberChecked:=numberChecked+1;;
 if Length(w)>0 then Add(positiveTables,Identifier(t));fi;
 Print(name,"; order=",Size(t),"; table=",Identifier(t),"; witnesses=",w,"\n");
od;
if numberChecked<>34 or positiveTables<>["U3(5)"] then Error("Unexpected simple-group result");fi;
# Cyclic simple groups have no irreducible character of indicator -1.
# Ten exceptions to the uniform defect-zero argument for dual-Burnside.
NoZeroColumns:=function(t) local chars; chars:=Irr(t);
 return Filtered([2..Length(chars)],j->ForAll(chars,x->x[j]<>0));end;;
exceptionalTables:=[];;
for name in ["M12","M22","M24","J2","HS","Suz","Ru","Co1","Co3","BM"] do
 t:=CharacterTable(name);;
 nonvanishing:=NoZeroColumns(t);;
 if Length(nonvanishing)>0 then Add(exceptionalTables,name);fi;
 Print("dual-Burnside ",name,": nonidentity columns without zero=",nonvanishing,"\n");
od;
if exceptionalTables<>["M22","M24"] then Error("Unexpected sporadic result");fi;
Print("PASS: historical simple-group and ten-sporadic character checks.\n");
QUIT;
