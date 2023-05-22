from mudata import MuData
from anndata import AnnData
import bioframe as bf


def set_coord(adata, range_df):
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
    def set_coord(self, range_df):
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
        return MuData({k: self.mod[k].subset_by_overlap(granges) for k in self.mod.keys()})

    def slice_granges(self, chrom, start, end):
        return MuData({k: self.mod[k].slice_granges(chrom, start, end) for k in self.mod.keys()})
