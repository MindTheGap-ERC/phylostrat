ID=1
NCHAR=30
CASE=sinusoid # don't change, only sinusoid data makes sense for this analysis
bash code/inference_fbd_gap_est.sh sim_data/fbd_gap_est/char_mat_${ID}_nchar${NCHAR}_${CASE}.nex sim_data/fbd_gap_est/fossils_${ID}_nchar${NCHAR}_${CASE}.csv sim_data/fbd_gap_est/tree_${ID}_nchar${NCHAR}_${CASE}.nex rb_output/fbd_gap_est/num_test_${ID}_${NCHAR}_${CASE} rb_output/fbd_gap_est/tree_test_${ID}_${NCHAR}_${CASE} sim_data/fbd_gap_est/mol_dat_${ID}_nchar${NCHAR}_${CASE}.nex