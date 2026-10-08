# Optional independent regeneration of the four group character-ring inputs.
# gap -q census/regenerate_group_rings.g > regenerated_group_rings.json
# GAP 4.12.1 was used for the frozen output.
EmitCharacterRing := function(g, name)
    local chars, n, dims, tensor;
    chars := Irr(g);
    n := Length(chars);
    dims := List(chars, x -> x[1]);
    tensor := List([1..n], i -> List([1..n], j -> List([1..n], k ->
        ScalarProduct(chars[i]*chars[j], chars[k]))));
    Print("{\"group\":\"", name, "\",\"order\":", Size(g),
          ",\"dimensions\":", dims, ",\"fusion_matrices\":", tensor, "}");
end;
Print("[\n");
EmitCharacterRing(AlternatingGroup(5), "A5");
Print(",\n");
EmitCharacterRing(PSL(2,7), "PSL2_7");
Print(",\n");
EmitCharacterRing(AlternatingGroup(6), "A6");
Print(",\n");
EmitCharacterRing(PSL(2,11), "PSL2_11");
Print("\n]\n");
QUIT;
