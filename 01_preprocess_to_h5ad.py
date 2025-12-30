import scanpy as sc
import os
import warnings
import gc
import matplotlib
import pandas as pd

matplotlib.use('Agg')  # 后台绘图
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

# --- 配置 ---
DATA_DIR = "./data"
OUTPUT_FILE = "./results/processed_breast_cancer.h5ad"
os.makedirs("./results", exist_ok=True)

# ⚡ 内存红线保护 ⚡
# 建议：3000 是 32GB 内存的“舒适区”上限。
MAX_CELLS_PER_SAMPLE = 3000

# ⚡ 核心 Trick：使用 inner join ⚡
JOIN_METHOD = 'inner'

# 映射表
GSE161529_GSM_MAP = {
    'GSM4909253': 'Normal', 'GSM4909254': 'Normal', 'GSM4909255': 'Normal',
    'GSM4909256': 'Normal', 'GSM4909257': 'Normal', 'GSM4909258': 'Normal',
    'GSM4909259': 'Normal', 'GSM4909260': 'Normal', 'GSM4909261': 'Normal',
    'GSM4909262': 'Normal', 'GSM4909263': 'Normal', 'GSM4909264': 'Normal',
    'GSM4909265': 'Normal', 'GSM4909266': 'Normal', 'GSM4909267': 'Normal',
    'GSM4909268': 'Normal', 'GSM4909269': 'Normal', 'GSM4909270': 'Normal',
    'GSM4909271': 'Normal', 'GSM4909272': 'Normal', 'GSM4909273': 'Normal',
    'GSM4909274': 'Normal', 'GSM4909275': 'Normal', 'GSM4909276': 'Normal',
    'GSM4909277': 'Normal', 'GSM4909278': 'Normal', 'GSM4909279': 'Normal', 'GSM4909280': 'Normal',
    'GSM4909281': 'TNBC', 'GSM4909282': 'TNBC', 'GSM4909283': 'TNBC', 'GSM4909284': 'TNBC',
    'GSM4909285': 'TNBC', 'GSM4909286': 'TNBC', 'GSM4909287': 'TNBC', 'GSM4909288': 'TNBC',
    'GSM4909289': 'HER2_Positive', 'GSM4909290': 'HER2_Positive', 'GSM4909291': 'HER2_Positive',
    'GSM4909292': 'HER2_Positive', 'GSM4909293': 'HER2_Positive', 'GSM4909294': 'HER2_Positive',
    'GSM4909295': 'Luminal_A', 'GSM4909296': 'Luminal_A', 'GSM4909297': 'Luminal_A', 'GSM4909298': 'Luminal_A',
    'GSM4909299': 'Luminal_A', 'GSM4909300': 'Luminal_A', 'GSM4909301': 'Luminal_A',
    'GSM4909302': 'Luminal_A', 'GSM4909303': 'Luminal_A', 'GSM4909304': 'Luminal_A',
    'GSM4909305': 'Luminal_A', 'GSM4909306': 'Luminal_A', 'GSM4909307': 'Luminal_A',
    'GSM4909308': 'Luminal_A', 'GSM4909309': 'Luminal_A', 'GSM4909310': 'Luminal_A', 'GSM4909311': 'Luminal_A',
    'GSM4909312': 'Luminal_A', 'GSM4909313': 'Luminal_A', 'GSM4909314': 'Luminal_A',
    'GSM4909315': 'Luminal_A', 'GSM4909316': 'Luminal_A', 'GSM4909317': 'Luminal_A',
    'GSM4909318': 'Luminal_A', 'GSM4909319': 'Luminal_A', 'GSM4909320': 'Luminal_A',
    'GSM4909321': 'Luminal_A'
}

DATASET_LEVEL_MAP = {
    "GSE240112": "Luminal_A",
    "GSE262288": "Luminal_B",
    "GSE274139": "Luminal_B",
    "GSE306201": "Luminal_A",
    "GSE289825": "TNBC"
}


def get_sample_subtype(gse_id, gsm_id):
    if gse_id == "GSE161529":
        return GSE161529_GSM_MAP.get(gsm_id, "Unknown")
    return DATASET_LEVEL_MAP.get(gse_id, "Unknown")


def safe_load_10x(directory):
    try:
        return sc.read_10x_mtx(directory, cache=False)
    except ValueError as e:
        if "Length of passed value" in str(e) and "var_names" in str(e):
            print(f"    ⚠️ 自动修复 Features 表头问题...")
            mtx_path = os.path.join(directory, "matrix.mtx.gz")
            feats_path = os.path.join(directory, "features.tsv.gz")
            if not os.path.exists(feats_path): feats_path = os.path.join(directory, "genes.tsv.gz")

            adata = sc.read_mtx(mtx_path).T
            features = pd.read_csv(feats_path, header=None, sep='\t')
            if len(features) == adata.n_vars + 1:
                features = features.iloc[1:]

            if features.shape[1] > 1:
                adata.var_names = features.iloc[:, 1].values
                adata.var['gene_ids'] = features.iloc[:, 0].values
            else:
                adata.var_names = features.iloc[:, 0].values
            return adata
        else:
            raise e


def load_and_merge_incrementally():
    merged_adata = None
    count = 0

    print(f"🚀 启动增量合并模式 (单样本上限: {MAX_CELLS_PER_SAMPLE}, Join={JOIN_METHOD})...")

    for gse_id in os.listdir(DATA_DIR):
        gse_path = os.path.join(DATA_DIR, gse_id)
        if not os.path.isdir(gse_path) or gse_id == "GSE274139": continue

        for root, dirs, files in os.walk(gse_path):
            if "matrix.mtx.gz" in files and "barcodes.tsv.gz" in files:
                try:
                    sample_id = os.path.basename(root)
                    subtype = get_sample_subtype(gse_id, sample_id)

                    # 1. 加载单个样本
                    current_adata = safe_load_10x(root)
                    current_adata.var_names_make_unique()

                    # ----------------------------------------------------
                    # 🔥 [内存保护] 下采样逻辑
                    # ----------------------------------------------------
                    if MAX_CELLS_PER_SAMPLE is not None and current_adata.n_obs > MAX_CELLS_PER_SAMPLE:
                        sc.pp.subsample(current_adata, n_obs=MAX_CELLS_PER_SAMPLE, random_state=42)
                        print(f"  [{count + 1}] 加载: {gse_id}/{sample_id} -> 下采样至 {MAX_CELLS_PER_SAMPLE} cells")
                    else:
                        print(f"  [{count + 1}] 加载: {gse_id}/{sample_id} ({current_adata.n_obs} cells)")

                    # 注入元数据
                    current_adata.obs['batch'] = gse_id
                    current_adata.obs['sample_id'] = sample_id
                    current_adata.obs['subtype'] = subtype

                    # 2. 增量合并
                    if merged_adata is None:
                        merged_adata = current_adata
                    else:
                        merged_adata = sc.concat(
                            [merged_adata, current_adata],
                            join=JOIN_METHOD,
                            index_unique=None
                        )

                    # 3. 强制垃圾回收
                    del current_adata
                    gc.collect()
                    count += 1

                except Exception as e:
                    print(f"❌ Error {sample_id}: {e}")

    if merged_adata is None: return None

    # 填补 NaN
    if JOIN_METHOD == 'outer':
        print("  Join='outer', 跳过 NaN 填补以节省内存...")
        pass

    merged_adata.obs_names_make_unique()
    print(f"✅ 合并完成！最终规模: {merged_adata.shape}")
    return merged_adata


def run_preprocess(adata):
    print("\n⚡ 开始低内存预处理...")

    # QC
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
    adata = adata[adata.obs.pct_counts_mt < 25, :]
    sc.pp.filter_cells(adata, min_genes=200)

    print(f"  QC后细胞数: {adata.n_obs}")

    print("  归一化 (Normalize)...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    print("  筛选高变基因 (HVG)...")
    # ⚡ 内存保护：只选 3000 个基因 ⚡
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, subset=False)

    # 备份全基因数据到 raw (这是关键，后续打分要用)
    adata.raw = adata

    # ⚡ 核心 Trick: 丢弃非高变基因 ⚡
    print("  切片: 只保留高变基因用于计算...")
    adata = adata[:, adata.var.highly_variable]
    gc.collect()

    print("  PCA 降维...")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')

    # Harmony
    try:
        import harmony
        print("  运行 Harmony...")
        sc.external.pp.harmony_integrate(adata, 'batch')
        use_rep = 'X_pca_harmony'
    except:
        print("  使用 PCA...")
        use_rep = 'X_pca'

    # 聚类
    print("  聚类 (Leiden)...")
    sc.pp.neighbors(adata, use_rep=use_rep)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5)

    # ----------------------------------------------------
    # 🔥 [修复点] 使用 use_raw=True 确保能找到所有 Marker
    # ----------------------------------------------------
    print("  细胞类型注释...")
    marker_genes = {
        "T_cells": ["CD3D", "CD3E", "CD2"],
        "B_cells": ["CD79A", "MS4A1"],
        "Epithelial": ["EPCAM", "KRT8", "KRT18"],
        "Macrophage": ["CD68", "CD163", "AIF1"],
        "Fibroblast": ["COL1A1", "DCN"],
        "Endothelial": ["PECAM1", "VWF"]
    }

    scored_types = []
    for cell_type, genes in marker_genes.items():
        # ⚠️ 注意：这里去 adata.raw.var_names 里找基因
        valid_genes = [g for g in genes if g in adata.raw.var_names]

        if valid_genes:
            # ⚠️ 注意：use_raw=True 确保用全基因数据打分
            sc.tl.score_genes(adata, valid_genes, score_name=cell_type, use_raw=True)
            scored_types.append(cell_type)
        else:
            print(f"  ⚠️ 警告: 未找到 {cell_type} 的 Marker 基因，将跳过。")

    if scored_types:
        # 只提取成功打分的列
        scores = adata.obs[scored_types]
        adata.obs['cell_type_major'] = scores.idxmax(axis=1)
    else:
        print("  ❌ 错误: 所有细胞类型均无法注释，标记为 Unknown")
        adata.obs['cell_type_major'] = "Unknown"

    print("  生成 UMAP 图...")
    sc.pl.umap(adata, color=['cell_type_major', 'batch', 'subtype'], save="_all_cells_integrated.png", show=False)

    return adata


if __name__ == "__main__":
    try:
        adata_merged = load_and_merge_incrementally()
        if adata_merged:
            adata_processed = run_preprocess(adata_merged)

            print(f"\n💾 正在保存结果到: {OUTPUT_FILE}")
            adata_processed.write(OUTPUT_FILE, compression="gzip")
            print("\n🎉🎉🎉 极限挑战成功！所有细胞已保存！ 🎉🎉🎉")
    except Exception as e:
        with open("error_log.txt", "w") as f:
            f.write(str(e))
            import traceback

            f.write(traceback.format_exc())
        print(f"❌ 发生错误: {e}")