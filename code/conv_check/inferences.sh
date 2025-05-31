SUFFIX2=cont
SUFFIX1=inc
NO=4


rb_command1="morph <- readDiscreteCharacterData(\"data/char_mat_${SUFFIX2}.nex\"); taxa <- readTaxonData(\"data/taxa_${SUFFIX2}.tsv\"); outfile <- \"output/test_${SUFFIX2}_${NO}_1.log\"; treefile_name <- \"output/trees_${SUFFIX2}_${NO}_1.log\";source(\"code/fbd_inference.rev\")"
rb_command2="morph <- readDiscreteCharacterData(\"data/char_mat_${SUFFIX2}.nex\"); taxa <- readTaxonData(\"data/taxa_${SUFFIX2}.tsv\"); outfile <- \"output/test_${SUFFIX2}_${NO}_2.log\"; treefile_name <- \"output/trees_${SUFFIX2}_${NO}_2.log\";source(\"code/fbd_inference.rev\")"
rb_command3="morph <- readDiscreteCharacterData(\"data/char_mat_${SUFFIX2}.nex\"); taxa <- readTaxonData(\"data/taxa_${SUFFIX2}.tsv\"); outfile <- \"output/test_${SUFFIX2}_${NO}_3.log\"; treefile_name <- \"output/trees_${SUFFIX2}_${NO}_3.log\";source(\"code/fbd_inference.rev\")"
rb_command4="morph <- readDiscreteCharacterData(\"data/char_mat_${SUFFIX2}.nex\"); taxa <- readTaxonData(\"data/taxa_${SUFFIX2}.tsv\"); outfile <- \"output/test_${SUFFIX2}_${NO}_4.log\"; treefile_name <- \"output/trees_${SUFFIX2}_${NO}_4.log\";source(\"code/fbd_inference.rev\")"
rb_command5="morph <- readDiscreteCharacterData(\"data/char_mat_${SUFFIX2}.nex\"); taxa <- readTaxonData(\"data/taxa_${SUFFIX2}.tsv\"); outfile <- \"output/test_${SUFFIX2}_${NO}_5.log\"; treefile_name <- \"output/trees_${SUFFIX2}_${NO}_5.log\";source(\"code/fbd_inference.rev\")"
rb_command6="morph <- readDiscreteCharacterData(\"data/char_mat_${SUFFIX1}.nex\"); taxa <- readTaxonData(\"data/taxa_${SUFFIX1}.tsv\"); outfile <- \"output/test_skyline_${NO}.log\"; treefile_name <- \"output/trees_skyline_${NO}_1.log\"; source(\"code/skyline_inference.rev\")"

echo $rb_command1 | "rb.exe" &
echo $rb_command2 | "rb.exe" &
echo $rb_command3 | "rb.exe" &
echo $rb_command4 | "rb.exe" &
echo $rb_command5 | "rb.exe" &
echo $rb_command6 | "rb.exe" &
wait