import pytest
import pandas as pd 

def test_parse():
    from tenextra.parse import select_likely_receptors
    clean_clones, all_clones, ct_chains, ct_pairs =\
        select_likely_receptors(f = 'tenextra/data/filtered_contig_annotations_test.csv')
    assert isinstance(clean_clones, pd.DataFrame)
    assert clean_clones.shape[0] < all_clones.shape[0]
