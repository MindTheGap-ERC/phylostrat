ID=24
NCHAR=1000
CASE=sinusoid # don't change, only sinusoid data makes sense for this analysis
bash code/inference_fbd_gap_est.sh data/fbd_gap_est/sim/char_mat_${ID}_nchar${NCHAR}_${CASE}.nex data/fbd_gap_est/sim/fossils_${ID}_nchar${NCHAR}_${CASE}.csv data/fbd_gap_est/rb_output/num_test_${ID}_${NCHAR}_${CASE} data/fbd_gap_est/rb_output/tree_test_${ID}_${NCHAR}_${CASE} data/fbd_gap_est/sim/mol_dat_${ID}_nchar${NCHAR}_${CASE}.nex