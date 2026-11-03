from extra_scripts.run_suppa2_full_data_comparison import pseudobulk_sample, sample_cell_type


def test_pseudobulk_labels_with_double_underscores():
    assert pseudobulk_sample("ASC__HC__m1") == "m1__ASC"
    assert pseudobulk_sample("CD16__NK_cells__LTCS__CO-10") == "CO-10__CD16__NK_cells"
    assert sample_cell_type("CO-10__CD16__NK_cells") == "CD16__NK_cells"
