import scanpy as sc
import sqlite3
import pandas as pd
import os
import warnings

# Suppress warnings
warnings.filterwarnings("ignore")

# --- Configuration ---
INPUT_FILE = "../results/processed_breast_cancer.h5ad"
DB_NAME = "../dbs/scrna_3nf.db"  # Target 3NF Database


def build_3nf_database():
    if not os.path.exists(INPUT_FILE):
        print(f"❌ Input file not found: {INPUT_FILE}")
        return

    print(f"📂 Loading: {INPUT_FILE} ...")
    adata = sc.read_h5ad(INPUT_FILE)

    # ---------------------------------------------------------
    # CRITICAL FIX: Use RAW data to include ALL genes, not just HVGs
    # ---------------------------------------------------------
    if adata.raw is not None:
        print("⚡ Utilizing adata.raw to access full transcriptome (20k+ genes)...")
        # Creating a view of the raw data for calculation
        use_adata = adata.raw.to_adata()
    else:
        print("⚠️ Warning: adata.raw is None. Only HVGs will be included.")
        use_adata = adata

    # 1. Prepare Dimension Tables

    # A. Genes Table (From raw data)
    genes = pd.DataFrame({'Gene': use_adata.var_names})
    genes['gene_id'] = range(len(genes))
    print(f"   Total Genes identified: {len(genes)}")

    # B. Cell Groups Table (Subtype + CellType)
    # Note: .obs matches between adata and adata.raw
    adata.obs['group_key'] = adata.obs['subtype'].astype(str) + "::" + adata.obs['cell_type_major'].astype(str)
    unique_groups = adata.obs[['subtype', 'cell_type_major', 'group_key']].drop_duplicates().reset_index(drop=True)
    unique_groups['group_id'] = range(len(unique_groups))
    print(f"   Biological Groups identified: {len(unique_groups)}")

    # 2. Calculate Expression (Fact Table)
    print("📊 Calculating Pseudobulk Means (this may take a moment)...")
    records = []

    # Map group_key to group_id for speed
    group_id_map = dict(zip(unique_groups['group_key'], unique_groups['group_id']))

    # Iterate through groups
    for i, group_key in enumerate(unique_groups['group_key']):
        group_id = group_id_map[group_key]

        # Create mask for current group
        mask = adata.obs['group_key'] == group_key

        # Calculate Mean Expression using RAW data
        # Note: use_adata[mask].X ensures we pull from the full gene matrix
        chunk_X = use_adata[mask].X

        # Handle sparse matrix conversion if necessary
        if hasattr(chunk_X, "toarray"):
            mean_expr = chunk_X.mean(axis=0).A1  # .A1 converts matrix to 1D array
        else:
            mean_expr = chunk_X.mean(axis=0)

        # Optimization: Filter out very low expression to save DB space
        # Threshold 0.01 is standard for pseudobulk visualization
        gene_names_indices = range(len(mean_expr))

        # Batch append to list (Python list append is fast)
        for gene_idx, val in zip(gene_names_indices, mean_expr):
            if val > 0.01:
                records.append((gene_idx, group_id, float(f"{val:.4f}")))

        if (i + 1) % 5 == 0:
            print(f"   Processed {i + 1}/{len(unique_groups)} groups...")

    # 3. Store into SQLite (3NF Structure)
    print(f"💾 Writing {len(records)} expression records to DB...")

    # Ensure directory exists
    os.makedirs(os.path.dirname(DB_NAME), exist_ok=True)

    conn = sqlite3.connect(DB_NAME)
    cursor = conn.cursor()

    # Drop old tables to ensure fresh schema
    cursor.execute("DROP VIEW IF EXISTS Table_Expression")
    cursor.execute("DROP TABLE IF EXISTS Expression")
    cursor.execute("DROP TABLE IF EXISTS Genes")
    cursor.execute("DROP TABLE IF EXISTS CellGroups")

    # Create Tables (3NF)
    cursor.execute("CREATE TABLE Genes (gene_id INTEGER PRIMARY KEY, gene_name TEXT)")
    cursor.execute("CREATE TABLE CellGroups (group_id INTEGER PRIMARY KEY, subtype TEXT, celltype TEXT)")

    # Fact Table with Foreign Keys
    cursor.execute('''
        CREATE TABLE Expression (
            gene_id INTEGER, 
            group_id INTEGER, 
            value REAL, 
            FOREIGN KEY(gene_id) REFERENCES Genes(gene_id), 
            FOREIGN KEY(group_id) REFERENCES CellGroups(group_id)
        )
    ''')

    # Bulk Insert
    print("   Inserting Dimensions...")
    cursor.executemany("INSERT INTO Genes VALUES (?,?)", zip(genes['gene_id'], genes['Gene']))
    cursor.executemany("INSERT INTO CellGroups VALUES (?,?,?)",
                       zip(unique_groups['group_id'], unique_groups['subtype'], unique_groups['cell_type_major']))

    print("   Inserting Facts...")
    cursor.executemany("INSERT INTO Expression VALUES (?,?,?)", records)

    conn.commit()

    # 4. Create Indices (Crucial for performance)
    print("⚡ Creating Indices...")
    cursor.execute("CREATE INDEX idx_gene_id ON Expression (gene_id)")
    cursor.execute("CREATE INDEX idx_group_id ON Expression (group_id)")
    cursor.execute("CREATE INDEX idx_gene_name ON Genes (gene_name)")

    # 5. Create View for App Compatibility
    print("🔗 Creating Compatibility View...")
    cursor.execute("""
        CREATE VIEW Table_Expression AS
        SELECT g.gene_name as Gene, c.subtype as Subtype, c.celltype as CellType, e.value as Avg_Expression
        FROM Expression e
        JOIN Genes g ON e.gene_id = g.gene_id
        JOIN CellGroups c ON e.group_id = c.group_id
    """)

    conn.close()
    print(f"\n🎉 3NF Database Built Successfully: {os.path.abspath(DB_NAME)}")


if __name__ == "__main__":
    build_3nf_database()