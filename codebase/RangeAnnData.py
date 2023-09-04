from mudata import MuData
from anndata import AnnData
import bioframe as bf


def set_coord(adata, range_df):
    bf.is_bedframe(range_df, raise_errors=True)
    adata.varm['coord'] = range_df.set_index(adata.var_names)


def subset_by_overlap(adata, granges):
    # coord = pr.PyRanges(adata.varm['coord'].reset_index())
    varm = adata.varm['coord'].copy()
    varm['idx'] = varm.index
    idx = bf.overlap(varm, granges, how='inner')['idx']
    return adata[:, idx]


def slice_granges(adata, chrom, start, end):
    idx = bf.select(adata.varm['coord'], f"{chrom}:{start}-{end}").index
    return adata[:, idx]


class RangeAnnData(AnnData):
    # if anndata object is given initalize with anndata 
    # else use AnnData to initialize
    def __init__(self, anndata=None, **kwargs):
        if anndata is not None:
            self = anndata
        else:
            super().__init__(**kwargs)

    def set_coord(self, range_df):
        # Check if range_df is a bedframe
        bf.is_bedframe(range_df, raise_errors=True)
        self.varm['coord'] = range_df.set_index(self.var_names)

    def subset_by_overlap(self, granges):
        varm = self.varm['coord'].copy()
        varm['idx'] = varm.index
        idx = bf.overlap(varm, granges, how='inner')['idx']
        return self[:, idx]

    def slice_granges(self, chrom, start, end):
        idx = bf.select(self.varm['coord'], f"{chrom}:{start}-{end}").index
        return self[:, idx]


class RangeMuData(MuData):
    def subset_by_overlap(self, granges):
        return self.__class__({k: self.mod[k].subset_by_overlap(granges) for k in self.mod.keys()})

    def slice_granges(self, chrom, start, end):
        return self.__class__({k: self.mod[k].slice_granges(chrom, start, end) for k in self.mod.keys()})


def read_h5ad(filename):
    """Read anndata object from h5ad file and change __class__ to RangeAnnData"""
    anndata = AnnData.read_h5ad(filename)
    anndata.__class__ = RangeAnnData
    return anndata

def read_h5mu(filename):
    """Read mudata object from h5mu file and change __class__ to RangeMuData"""
    mudata = MuData.read_h5mu(filename)
    mudata.__class__ = RangeMuData
    return mudata