import mudata
import anndata
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


class RangeAnnData(anndata.AnnData):

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


class RangeMuData(mudata.MuData):
    def subset_by_overlap(self, granges):
        return self.__class__({k: self.mod[k].subset_by_overlap(granges) for k in self.mod.keys()})

    def slice_granges(self, chrom, start, end):
        return self.__class__({k: self.mod[k].slice_granges(chrom, start, end) for k in self.mod.keys()})


def read_h5ad(filename):
    """Read anndata object from h5ad file and change __class__ to RangeAnnData"""
    adata = anndata.read_h5ad(filename)
    adata.__class__ = RangeAnnData
    return adata


def read_h5mu(filename):
    """Read mudata object from h5mu file and change __class__ to RangeMuData"""
    mu = mudata.read_h5mu(filename)
    mu.__class__ = RangeMuData
    return mu
