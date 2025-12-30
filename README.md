# 🧬 BC-MOD: Breast Cancer Multi-Omics Database Platform
# 乳腺癌多组学联邦数据库平台

[![Python[dbs](dbs)](https://img.shields.io/badge/Python-3.9%2B-blue)](https://www.python.org/)
[![Streamlit](https://img.shields.io/badge/Streamlit-1.30%2B-red)](https://streamlit.io/)
[![MySQL](https://img.shields.io/badge/Database-MySQL%20%26%20SQLite-orange)](https://www.mysql.com/)

## 📖 Introduction

**BC-MOD** is a high-performance federated multi-omics data integration platform designed for breast cancer research. The project adopts a **"Federated Database"** architecture, combining **MySQL** (for user authentication, RBAC, and audit logs) with **SQLite** (for distributed storage of massive omics data) on the backend, and uses **Streamlit** for the interactive frontend visualization.

The system integrates the following five core omics datasets, supporting cross-dimensional gene retrieval and heterogeneity analysis:
1.  **🔬 scRNA-seq (Single-cell Transcriptomics)**: Expression profiles covering 300,000+ cells.
2.  **🧬 ATAC-seq (Chromatin Accessibility)**: Analysis of epigenetic regulatory regions.
3.  **⚗️ Metabolomics**: Gene-metabolite association networks.
4.  **🗺️ Spatial Transcriptomics**: In situ expression analysis of tissue slices.
5.  **🖼️ Imaging (Pathology)**: AI-assisted annotation display of H&E stained slides.

---

## ✨ Features

- **Hybrid Database Architecture**:
    - **Management Layer**: MySQL handles User Authentication (Auth), Role-Based Access Control (RBAC), and Audit Logs.
    - **Data Layer**: SQLite stores research data processed via Pseudobulk dimensionality reduction, enabling millisecond-level queries.
- **Multi-dimensional Query Modes**:
    - Supports global search by **Gene Symbol** (e.g., `FOXA1`).
    - Supports targeted analysis by **Cell Type**, **Metabolite**, **Sample ID**, and **Spatial Region**.
- **Interactive Visualization**: Built-in UMAP plots, bar charts, spatial heatmaps, and pathology annotation overlays.
- **Enterprise-grade Features**:
    - 🔐 Secure Login System (SHA-256 Encryption).
    - 🌐 Internationalization (i18n) - Toggle between Chinese/English.
    - 🛠️ Admin Dashboard (CRUD User Management).

---

## 📂 Directory Structure

```text
BC_MOD_Project/
│
├── app.py                  # [Entry] Main Application (Streamlit Frontend)
├── auth_manager.py         # [Module] MySQL Connection & Auth Logic
├── init_db.py              # [Tool] Database Initialization Script
├── requirements.txt        # [Deps] Project Dependencies
│
├── dbs/                    # [Data] Core Omics Database Files (SQLite)
│   ├── scrna.db
│   ├── atac_breast_cancer.db
│   ├── metabolomics.db
│   ├── spatial.db
│   └── imaging.db
│
├── scripts/                # [Scripts] ETL Pipeline (Crawling, Cleaning, Building DB)
│   ├── 01_crawl_geo_metadata.py
│   ├── 02_download_matrices.py
│   ├── 03_organize_files.py
│   ├── 04_process_to_h5ad.py
│   └── 05_build_sqlite_db.py
│
└── README.md               # Project Documentation
```

---

## 🚀 Quick Start

### 1. Prerequisites
Ensure the following are installed locally:
* **Python 3.8+** (Anaconda recommended)
* **MySQL Server 8.0+**

### 2. Installation

```bash
# Create virtual environment (Optional)
conda create -n bc_db python=3.9
conda activate bc_db

# Install Python dependencies
pip install -r requirements.txt
# Or manually install core libraries:
# pip install streamlit mysql-connector-python pandas matplotlib scanpy
```

### 3. Configuration
1.  Start the MySQL Service.
2.  Open `auth_manager.py` and `init_db.py`, locate the configuration section, and update the `password` to your local MySQL Root password:
    ```python
    # auth_manager.py & init_db.py
    DB_CONFIG = {
        'host': 'localhost',
        'user': 'root',
        'password': 'YOUR_PASSWORD',  # <--- Change this
        'database': 'bc_mod_admin'
    }
    ```
3.  Run the initialization script to create tables and default accounts:
    ```bash
    python init_db.py
    ```
    > Success Message: `✅✅✅ 数据库初始化完美完成！`

### 4. Run the Application
Run the following command in the project root directory:
```bash
streamlit run app.py
```
The browser will automatically open `http://localhost:8501`.

---

## 📖 User Guide

### 🔐 Login
Two default test accounts are pre-configured:

| Role | Username | Password | Permissions |
| :--- | :--- | :--- | :--- |
| **Admin** | `admin` | `admin123456` | Full Query Access + Admin Panel (User Mgmt) |
| **Guest** | `guest` | `guest123456` | Data Query Access Only |

### 🔍 Data Query (Query Mode)
Default view after login:
1.  **Sidebar**: Toggle Language (CN/EN), browse Data Dictionary (Top Genes/Metabolites).
2.  **Query Modes**:
    - **Global Search**: Enter a gene (e.g., `FOXA1`, `ESR1`) to see its profile across all omics.
    - **Advanced Modes**: Select specific analysis types via the dropdown (e.g., "scRNA Cell Type Analysis", "ATAC Sample Analysis") to explore Top N features.

### 🛠️ Admin Panel
*Visible only to Admins*
1.  Select **"🛠️ Admin Panel"** in the sidebar navigation.
2.  **User Management**: Create new users, assign roles (Admin/Guest).
3.  **Audit Logs**: View system login and operation history (stored in MySQL).

---

## ⚠️ Troubleshooting

**Q: `PermissionError: [WinError 10013]` on startup?**
* **Cause**: Default port 8501 is occupied or restricted.
* **Solution**: Run on a different port:
    ```bash
    streamlit run app.py --server.port 9527
    ```

**Q: `MySQL Connection Failed`?**
* Check if the MySQL service is running (via Task Manager).
* Verify that the password in `auth_manager.py` matches your local MySQL password.

---

## 📜 License
This project is created for a Principles of Database Design assignment and is intended for academic research and educational demonstration purposes only.

---
*© 2025 BC-MOD Project Team. All Rights Reserved.*