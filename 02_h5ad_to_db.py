import scanpy as sc
import pandas as pd
import numpy as np
import sqlite3
import os
import warnings

warnings.filterwarnings("ignore")

# --- 配置 ---
INPUT_FILE = "./results/processed_breast_cancer.h5ad"
DB_NAME = "scRNA_data.db"


def build_database():
    if not os.path.exists(INPUT_FILE):
        print(f"❌ 未找到 {INPUT_FILE}，请先运行上一步的预处理脚本！")
        return

    print(f"📂 正在加载预处理结果: {INPUT_FILE} ...")
    adata = sc.read_h5ad(INPUT_FILE)
    print(f"   数据规模: {adata.shape}")

    # 1. 检查必要字段
    required_obs = ['subtype', 'cell_type_major']
    for col in required_obs:
        if col not in adata.obs.columns:
            raise ValueError(f"❌ 数据中缺少必要注释列: {col}，请检查预处理是否成功。")

    # 2. 构建分组键 (Group Key)
    # 格式: Subtype::CellType (例如 TNBC::T_cells)
    print("🔄 正在构建分组键...")
    adata.obs['group_key'] = adata.obs['subtype'].astype(str) + "::" + adata.obs['cell_type_major'].astype(str)

    unique_groups = adata.obs['group_key'].unique()
    print(f"   共识别出 {len(unique_groups)} 个生物学分组 (Subtype x CellType)")

    # 3. 计算 Pseudobulk (内存安全版)
    # 我们遍历每个分组，单独计算均值，而不是把整个矩阵转为 DataFrame
    print("📊 正在计算 Pseudobulk 均值 (使用 adata.raw 中的全基因数据)...")

    # 优先使用 raw (包含了所有基因，不仅仅是高变基因)
    if adata.raw is not None:
        use_adata = adata.raw.to_adata()  # 这是一个虚拟视图，不占内存
    else:
        print("⚠️ 警告: 未找到 raw 数据，将使用仅含高变基因的 X 矩阵。")
        use_adata = adata

    records = []
    total_genes = use_adata.n_vars

    for i, group in enumerate(unique_groups):
        if "Unknown" in group:
            continue

        subtype, cell_type = group.split("::")
        print(f"   [{i + 1}/{len(unique_groups)}] 处理: {subtype} - {cell_type}")

        # 获取该组细胞的索引
        cells_mask = adata.obs['group_key'] == group

        # 提取这些细胞的表达矩阵 (Sparse)
        # 注意: 只提取这几千个细胞，内存非常安全
        chunk_X = use_adata[cells_mask].X

        # 计算均值 (转为 dense array)
        # axis=0 表示沿细胞方向求均值
        if hasattr(chunk_X, "toarray"):
            mean_expression = chunk_X.mean(axis=0).A1  # .A1 转为 1D array
        else:
            mean_expression = chunk_X.mean(axis=0)

        # 提取基因名
        gene_names = use_adata.var_names

        # 组装数据 (过滤掉表达量极低 < 0.01 的基因以压缩体积)
        # 这一步能让数据库体积减小 50% 以上
        for gene_idx, exp_val in enumerate(mean_expression):
            if exp_val > 0.01:
                records.append((
                    gene_names[gene_idx],  # Gene
                    subtype,  # Subtype
                    cell_type,  # CellType
                    float(f"{exp_val:.4f}")  # Avg_Expression (保留4位小数)
                ))

    # 4. 写入 SQLite
    print(f"💾 正在写入数据库 (共 {len(records)} 条记录)...")
    conn = sqlite3.connect(DB_NAME)
    cursor = conn.cursor()

    # 建表 (符合项目规范)
    cursor.execute('DROP TABLE IF EXISTS Table_Expression')
    cursor.execute('''
        CREATE TABLE Table_Expression (
            Gene TEXT,
            Subtype TEXT,
            CellType TEXT,
            Avg_Expression REAL
        )
    ''')

    # 批量插入 (Executemany 速度极快)
    cursor.executemany('INSERT INTO Table_Expression VALUES (?,?,?,?)', records)
    conn.commit()

    # 5. 创建索引 (查询速度提升 100 倍的关键)
    print("⚡ 正在创建 B-Tree 索引...")
    cursor.execute('CREATE INDEX idx_gene ON Table_Expression (Gene)')
    cursor.execute('CREATE INDEX idx_subtype_cell ON Table_Expression (Subtype, CellType)')

    conn.close()

    print(f"\n🎉🎉🎉 数据库构建成功！")
    print(f"文件位置: {os.path.abspath(DB_NAME)}")
    print("下一步：将此文件发送给组长(你自己)用于 Streamlit 整合。")


def verify_db():
    """简单的验证函数，看看数据对不对"""
    if not os.path.exists(DB_NAME): return

    print(f"\n🔍 数据库自检 (Sample Query):")
    conn = sqlite3.connect(DB_NAME)

    # 查一个经典基因 PKM2
    df = pd.read_sql_query("SELECT * FROM Table_Expression WHERE Gene='PKM2' ORDER BY Avg_Expression DESC LIMIT 5",
                           conn)
    print("查询 Gene='PKM2' 的前5条结果:")
    print(df)

    # 统计行数
    count = conn.execute("SELECT Count(*) FROM Table_Expression").fetchone()[0]
    print(f"\n数据库总行数: {count}")
    conn.close()


if __name__ == "__main__":
    build_database()
    verify_db()