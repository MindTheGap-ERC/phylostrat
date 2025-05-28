SUFFIX1=full
SUFFIX2=cont
SUFFIX3=inc
#rb_command1="morph <- readDiscreteCharacterData(\"data/char_mat_${SUFFIX1}.nex\"); taxa <- readTaxonData(\"data/taxa_${SUFFIX1}.tsv\"); outfile <- \"output/test_${SUFFIX1}.log\"; source(\"code/fbd_inference.rev\")"
rb_command2="morph <- readDiscreteCharacterData(\"data/char_mat_${SUFFIX2}.nex\"); taxa <- readTaxonData(\"data/taxa_${SUFFIX2}.tsv\"); outfile <- \"output/test_${SUFFIX2}.log\"; source(\"code/fbd_inference.rev\")"
rb_command3="morph <- readDiscreteCharacterData(\"data/char_mat_${SUFFIX3}.nex\"); taxa <- readTaxonData(\"data/taxa_${SUFFIX3}.tsv\"); outfile <- \"output/test_${SUFFIX3}.log\"; source(\"code/fbd_inference.rev\")"

#echo $rb_command1 | "rb.exe" &
echo $rb_command2 | "rb.exe" &
echo $rb_command3 | "rb.exe" &
wait