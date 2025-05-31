SUFFIX1=cont
SUFFIX2=inc
NO=3

rb_command1="morph <- readDiscreteCharacterData(\"data/char_mat_${SUFFIX1}.nex\"); taxa <- readTaxonData(\"data/taxa_${SUFFIX1}.tsv\"); outfile <- \"output/test_${SUFFIX1}_${NO}.log\"; treefile_name <- \"output/trees_${SUFFIX1}_${NO}.log\"; source(\"code/fbd_inference.rev\")"
rb_command2="morph <- readDiscreteCharacterData(\"data/char_mat_${SUFFIX2}.nex\"); taxa <- readTaxonData(\"data/taxa_${SUFFIX2}.tsv\"); outfile <- \"output/test_${SUFFIX2}_${NO}.log\"; treefile_name <- \"output/trees_${SUFFIX2}_${NO}.log\";source(\"code/fbd_inference.rev\")"
rb_command3="morph <- readDiscreteCharacterData(\"data/char_mat_${SUFFIX2}.nex\"); taxa <- readTaxonData(\"data/taxa_${SUFFIX2}.tsv\"); outfile <- \"output/test_skyline_${NO}.log\"; treefile_name <- \"output/trees_skyline_${NO}.log\"; source(\"code/skyline_inference.rev\")"

echo $rb_command1 | "rb.exe" &
echo $rb_command2 | "rb.exe" &
echo $rb_command3 | "rb.exe" &
wait