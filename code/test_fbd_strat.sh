ID=1
NCHAR=30
CASE=sinusoid
bash code/inference_fbd_base.sh sim_data/fbd_strat/char_mat_${ID}_nchar${NCHAR}_${CASE}.nex sim_data/fbd_strat/fossils_${ID}_nchar${NCHAR}_${CASE}.csv rb_output/fbd_strat/num_test_${ID}_${NCHAR}_${CASE} rb_output/fbd_strat/tree_test_${ID}_${NCHAR}_${CASE} sim_data/fbd_strat/mol_dat_${ID}_nchar${NCHAR}_${CASE}.nex