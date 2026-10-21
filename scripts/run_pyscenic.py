"""
pySCENIC pipeline (GRNBoost2 -> cisTarget pruning -> AUCell), used for Fig 5E-G / Fig S3.
Usage: python run_pyscenic.py <expr_cells_x_genes.csv.gz> <out_prefix> [n_workers]
Outputs: <out_prefix>_adj.tsv, <out_prefix>_regulons.csv, <out_prefix>_auc.csv
"""
import sys, os, pickle
import pandas as pd

DB = os.environ.get("SCENIC_DB", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "cisTarget"))


def main():
    expr_file, out = sys.argv[1], sys.argv[2]
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    from arboreto.algo import grnboost2
    from distributed import LocalCluster, Client
    from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
    from pyscenic.utils import modules_from_adjacencies
    from pyscenic.prune import prune2df, df2regulons
    from pyscenic.aucell import aucell

    ex = pd.read_csv(expr_file, index_col=0)
    print("expression", ex.shape, flush=True)
    tfs = [t.strip() for t in open(os.path.join(DB, "hs_hgnc_tfs.txt")) if t.strip() in ex.columns]
    print("TFs in data:", len(tfs), flush=True)

    adj_file = out + "_adj.tsv"
    if not os.path.exists(adj_file):
        cluster = LocalCluster(n_workers=workers, threads_per_worker=1, memory_limit="6GB")
        client = Client(cluster)
        adj = grnboost2(ex, tf_names=tfs, verbose=True, seed=1234, client_or_address=client)
        adj.to_csv(adj_file, sep="\t", index=False)
        client.close(); cluster.close()
    adj = pd.read_csv(adj_file, sep="\t")

    dbs = [RankingDatabase(fname=os.path.join(DB, f), name=f.split(".")[0])
           for f in os.listdir(DB) if f.endswith(".feather")]
    modules = list(modules_from_adjacencies(adj, ex))
    motifs = os.path.join(DB, "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl")
    df = prune2df(dbs, modules, motifs, num_workers=workers)
    regulons = df2regulons(df)
    with open(out + "_regulons.pkl", "wb") as f:
        pickle.dump(regulons, f)
    rows = [(r.name, g, w) for r in regulons for g, w in r.gene2weight.items()]
    pd.DataFrame(rows, columns=["regulon", "gene", "weight"]).to_csv(out + "_regulons.csv", index=False)
    auc = aucell(ex, regulons, num_workers=workers, seed=1234)
    auc.to_csv(out + "_auc.csv")
    print("regulons:", len(regulons), flush=True)


if __name__ == "__main__":
    main()
