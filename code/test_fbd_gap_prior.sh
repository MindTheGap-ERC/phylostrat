ID=1
NCHAR=30
CASE=sinusoid # don't change, only sinusoid data makes sense for this analysis
bash code/inference_fbd_gap_prior.sh data/sim/fbd_gap_prior/char_mat_${ID}_nchar${NCHAR}_${CASE}.nex data/sim/fbd_gap_prior/fossils_${ID}_nchar${NCHAR}_${CASE}.csv data/sim/fbd_gap_prior/tree_${ID}_nchar${NCHAR}_${CASE}.nex output/fbd_gap_prior/num_test_${ID}_${NCHAR}_${CASE} output/fbd_gap_prior/tree_test_${ID}_${NCHAR}_${CASE} data/sim/fbd_gap_prior/mol_dat_${ID}_nchar${NCHAR}_${CASE}.nex