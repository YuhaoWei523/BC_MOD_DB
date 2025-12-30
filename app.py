import asyncio
import sys

# Windows asyncio event loop fix
if sys.platform.startswith("win"):
    asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())

import streamlit as st
import sqlite3
import pandas as pd
import os
import json
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import auth_manager as auth

# ==========================================
# 1. Internationalization (i18n) Dictionary
# ==========================================
TRANS = {
    "CN": {
        "title": "BC-MOD 乳腺癌多组学数据库",
        "sidebar_header_user": "👤 用户中心",
        "sidebar_header_nav": "🧭 功能导航",
        "sidebar_header_settings": "⚙️ 设置与工具",
        "nav_query": "🔍 数据查询",
        "nav_admin": "🛠️ 后台管理",
        "login_title": "🔐 系统登录",
        "login_btn": "登录",
        "logout_btn": "退出登录",
        "intro_title": "系统介绍",
        "intro_text": "本系统整合 scRNA-seq, ATAC-seq, Spatial, Metabolomics 以及 Imaging 数据。",
        "lbl_language": "语言 / Language",

        # Query Modes
        "qmode_label": "选择查询模式",
        "qmode_gene": "🧬 全库搜索 (By Gene)",
        "qmode_scrna": "🔬 scRNA 细胞类型分析",
        "qmode_atac": "🧬 ATAC 样本开放度分析",
        "qmode_metabo": "⚗️ 代谢物关联分析",
        "qmode_spatial": "🗺️ 空间区域异质性分析",

        # Inputs & Filters
        "input_gene": "输入基因 Symbol",
        "input_gene_ph": "尝试: FOXA1, ESR1, PKM",
        "input_metabo": "输入代谢物名称",
        "input_celltype": "选择细胞类型",
        "input_region": "选择空间区域 (Region)",
        "input_sample": "选择样本 ID",
        "input_top_n": "筛选 Top N 结果",
        "lbl_subtype": "亚型过滤",

        # Tab Titles
        "tab_scrna": "🔬 scRNA (单细胞)",
        "tab_atac": "🧬 ATAC (表观)",
        "tab_metabo": "⚗️ Metabo (代谢)",
        "tab_spatial": "🗺️ Spatial (空间)",
        "tab_imaging": "🖼️ Imaging (影像)",
        "tab_user_mgmt": "👥 用户管理",
        "tab_audit_logs": "📝 审计日志",
        "tab_data_crud": "🧬 组学数据维护 (CRUD)",

        # Admin CRUD
        "crud_create_title": "➕ 新增专家注释 (Create)",
        "crud_create_desc": "模拟专家审编流程：为特定基因添加生物学描述。",
        "crud_update_title": "📝 修正注释 (Update)",
        "crud_delete_title": "🗑️ 删除注释 (Delete)",
        "btn_submit": "提交",
        "btn_update": "更新",
        "btn_delete": "删除",
        "msg_success_create": "注释添加成功！",
        "msg_success_update": "注释更新成功！",
        "msg_success_delete": "注释删除成功！",

        # Headers & Captions
        "header_top_genes": "🔥 高表达基因排行",
        "caption_plot_limit": "注：为保证图表清晰，图表仅展示前 50 项，完整数据请见下方表格。",
        "data_browser": "📚 数据字典导览",
        "top_genes_list": "🔥 高表达基因 (scRNA)",
        "top_metas_list": "🧪 高表达代谢物",
        "atac_sim_note": "⚠️ 注：当前 ATAC 数据库缺失临床亚型标注。下图展示基于模拟元数据的分组对比。",
        "atac_raw_title": "2. 原始样本分布 (未过滤)",
        "spatial_single_note": "ℹ️ 提示：当前 Spatial 模块展示标准参考样本 V1 (HER2_Positive)。",
        "imaging_note": "💡 说明：展示 AI 辅助识别的肿瘤感兴趣区域 (ROI)。本模块为非结构化数据存储演示。",
        "warn_no_data": "未找到相关数据。",
        "success_login": "验证通过！正在进入系统...",
        "err_login": "用户名或密码错误"
    },
    "EN": {
        "title": "BC-MOD Multi-Omics Database",
        "sidebar_header_user": "👤 User Center",
        "sidebar_header_nav": "🧭 Navigation",
        "sidebar_header_settings": "⚙️ Settings & Tools",
        "nav_query": "🔍 Data Query",
        "nav_admin": "🛠️ Admin Panel",
        "login_title": "🔐 System Login",
        "login_btn": "Login",
        "logout_btn": "Logout",
        "intro_title": "Introduction",
        "intro_text": "A database integrating scRNA-seq, ATAC-seq, Spatial, Metabolomics and Imaging data.",
        "lbl_language": "Language",

        # Query Modes
        "qmode_label": "Query Mode",
        "qmode_gene": "🧬 Global Search (By Gene)",
        "qmode_scrna": "🔬 scRNA Cell Type Analysis",
        "qmode_atac": "🧬 ATAC Sample Analysis",
        "qmode_metabo": "⚗️ Metabolite Analysis",
        "qmode_spatial": "🗺️ Spatial Region Analysis",

        # Inputs & Filters
        "input_gene": "Enter Gene Symbol",
        "input_gene_ph": "Try: FOXA1, ESR1, PKM",
        "input_metabo": "Enter Metabolite Name",
        "input_celltype": "Select Cell Type",
        "input_region": "Select Spatial Region",
        "input_sample": "Select Sample ID",
        "input_top_n": "Show Top N",
        "lbl_subtype": "Subtype Filter",

        # Tab Titles
        "tab_scrna": "🔬 scRNA",
        "tab_atac": "🧬 ATAC",
        "tab_metabo": "⚗️ Metabo",
        "tab_spatial": "🗺️ Spatial",
        "tab_imaging": "🖼️ Imaging",
        "tab_user_mgmt": "👥 User Mgmt",
        "tab_audit_logs": "📝 Audit Logs",
        "tab_data_crud": "🧬 Data Maintenance (CRUD)",

        # Admin CRUD
        "crud_create_title": "➕ Add Expert Annotation (Create)",
        "crud_create_desc": "Expert biocuration: Add biological descriptions for genes.",
        "crud_update_title": "📝 Update Annotation",
        "crud_delete_title": "🗑️ Delete Annotation",
        "btn_submit": "Submit",
        "btn_update": "Update",
        "btn_delete": "Delete",
        "msg_success_create": "Annotation added successfully!",
        "msg_success_update": "Annotation updated successfully!",
        "msg_success_delete": "Annotation deleted successfully!",

        # Headers & Captions
        "header_top_genes": "🔥 Top Expressed Genes",
        "caption_plot_limit": "Note: Plot limited to top 50 items for clarity. See table for full list.",
        "data_browser": "📚 Data Dictionary",
        "top_genes_list": "🔥 Top Genes (scRNA)",
        "top_metas_list": "🧪 Top Metabolites",
        "atac_sim_note": "⚠️ Note: ATAC subtypes are simulated for demonstration.",
        "atac_raw_title": "2. Raw Sample Distribution",
        "spatial_single_note": "ℹ️ Note: Showing Reference Sample V1 (HER2_Positive).",
        "imaging_note": "💡 Note: AI-identified Tumor ROI demo.",
        "warn_no_data": "No data found.",
        "success_login": "Login successful! Redirecting...",
        "err_login": "Invalid username or password"
    }
}

# ==========================================
# 2. Configuration & Helpers
# ==========================================
st.set_page_config(page_title="BC-MOD Database", page_icon="🧬", layout="wide")

# Database File Paths
DB_PATHS = {
    "scRNA": "./dbs/scrna_3nf.db",
    "ATAC": "./dbs/atac.db",
    "Metabo": "./dbs/metabolomics.db",
    "Spatial": "./dbs/spatial.db",
    "Imaging": "./dbs/imaging.db"
}

# Session State Initialization
if 'logged_in' not in st.session_state:
    st.session_state['logged_in'] = False
    st.session_state['user_role'] = None
    st.session_state['username'] = None
if 'lang' not in st.session_state:
    st.session_state['lang'] = "CN"


def t(key):
    """Retrieve translation for key based on current language."""
    return TRANS[st.session_state['lang']].get(key, key)


def run_sqlite_query(db_key, sql):
    """Execute SQL query safely."""
    db_path = DB_PATHS.get(db_key)
    if not os.path.exists(db_path): return None
    try:
        conn = sqlite3.connect(db_path)
        df = pd.read_sql_query(sql, conn)
        conn.close()
        return df
    except:
        return None


def get_atac_meta(sample_id):
    """Simulate metadata mapping for ATAC samples."""
    types = ["TNBC", "HER2_Positive", "Luminal_A", "Luminal_B", "Normal"]
    return types[hash(sample_id) % len(types)]


def get_distinct_values(db_key, table, column):
    """Fetch unique values for dropdown menus."""
    sql = f"SELECT DISTINCT {column} FROM {table} ORDER BY {column}"
    df = run_sqlite_query(db_key, sql)
    return df[column].tolist() if df is not None else []


def get_real_top_elements(mode):
    """Fetch top elements for Data Browser."""
    try:
        if mode == "gene":
            sql = "SELECT Gene FROM Table_Expression GROUP BY Gene ORDER BY SUM(Avg_Expression) DESC LIMIT 10"
            df = run_sqlite_query("scRNA", sql)
            return df['Gene'].tolist() if df is not None else []
        elif mode == "metabo":
            sql = "SELECT Metabolite FROM Metabolite_Expression GROUP BY Metabolite ORDER BY SUM(Expression_Level) DESC LIMIT 10"
            df = run_sqlite_query("Metabo", sql)
            return df['Metabolite'].tolist() if df is not None else []
    except:
        return []
    return []


# --- MySQL Helper for Annotations ---
def ensure_annotation_table():
    """Create gene_annotations table in MySQL if not exists."""
    conn = auth.get_connection()
    if conn:
        cursor = conn.cursor()
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS gene_annotations (
                id INT AUTO_INCREMENT PRIMARY KEY,
                gene VARCHAR(50),
                note TEXT,
                author VARCHAR(50),
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """)
        conn.commit()
        conn.close()


# ==========================================
# 3. UI Modules
# ==========================================

def login_ui():
    st.markdown(f"## 🧬 {t('title')}")
    st.markdown("---")

    # Language Switcher on Login Page
    lang = st.radio("Language / 语言", ["CN", "EN"], horizontal=True)
    st.session_state['lang'] = lang

    c1, c2 = st.columns([1, 1])
    with c1:
        st.info(f"""
        **{t('intro_title')}**:
        {t('intro_text')}

        **Test Accounts**:
        - Admin: `admin` / `admin123456`
        - Guest: `guest` / `guest123456`
        """)

    with c2:
        st.subheader(t('login_title'))
        user = st.text_input("Username")
        pwd = st.text_input("Password", type="password")
        if st.button(t('login_btn'), type="primary", use_container_width=True):
            u = auth.check_login(user, pwd)
            if u:
                st.session_state['logged_in'] = True
                st.session_state['user_role'] = u['role']
                st.session_state['username'] = u['username']
                auth.log_action(user, "Login")
                ensure_annotation_table()  # Initialize annotation table on login
                st.success(t('success_login'))
                st.rerun()
            else:
                st.error(t('err_login'))


def admin_ui():
    st.header(t('nav_admin'))
    st.warning(f"Admin: {st.session_state['username']}")

    tab1, tab2, tab3 = st.tabs([t('tab_user_mgmt'), t('tab_audit_logs'), t('tab_data_crud')])

    with tab1:
        st.subheader("Create New User")
        with st.form("new_u"):
            c1, c2, c3 = st.columns(3)
            nu = c1.text_input("Username")
            np = c2.text_input("Password", type="password")
            nr = c3.selectbox("Role", ["guest", "admin"])
            if st.form_submit_button("Create"):
                if auth.create_user(nu, np, nr):
                    st.success(f"User {nu} created!")
                    auth.log_action(st.session_state['username'], f"Create user {nu}")
                else:
                    st.error("Failed. Username exists?")

    with tab2:
        st.subheader("Audit Logs")
        if st.button("Refresh"): st.rerun()
        st.dataframe(auth.get_system_logs(20), use_container_width=True)

    # --- [CRUD] Expert Annotation System ---
    with tab3:
        st.subheader(t('tab_data_crud'))

        # 1. Create Annotation
        with st.expander(t('crud_create_title'), expanded=True):
            st.caption(t('crud_create_desc'))
            c1, c2 = st.columns([1, 2])
            target_gene = c1.text_input("Target Gene", "FOXA1").strip().upper()
            note_content = c2.text_input("Annotation Content")

            if st.button(t('btn_submit')):
                if target_gene and note_content:
                    conn = auth.get_connection()
                    if conn:
                        cursor = conn.cursor()
                        cursor.execute("INSERT INTO gene_annotations (gene, note, author) VALUES (%s, %s, %s)",
                                       (target_gene, note_content, st.session_state['username']))
                        conn.commit()
                        conn.close()
                        st.success(t('msg_success_create'))
                        auth.log_action(st.session_state['username'], f"Create Annotation: {target_gene}")
                else:
                    st.warning("Please fill all fields.")

        st.divider()

        # 2. Update & Delete (Manage Existing)
        st.markdown("#### Manage Existing Annotations")
        conn = auth.get_connection()
        if conn:
            # Load all annotations
            df_notes = pd.read_sql("SELECT * FROM gene_annotations ORDER BY created_at DESC", conn)
            conn.close()

            if not df_notes.empty:
                for index, row in df_notes.iterrows():
                    with st.container(border=True):
                        c1, c2, c3 = st.columns([1, 3, 1])
                        c1.markdown(f"**{row['gene']}**")

                        # Update Logic
                        new_note = c2.text_input("Edit Note", row['note'], key=f"edit_{row['id']}")
                        if c2.button(t('btn_update'), key=f"btn_up_{row['id']}"):
                            conn = auth.get_connection()
                            cursor = conn.cursor()
                            cursor.execute("UPDATE gene_annotations SET note=%s WHERE id=%s", (new_note, row['id']))
                            conn.commit()
                            conn.close()
                            st.success(t('msg_success_update'))
                            st.rerun()

                        # Delete Logic
                        if c3.button(t('btn_delete'), key=f"btn_del_{row['id']}"):
                            conn = auth.get_connection()
                            cursor = conn.cursor()
                            cursor.execute("DELETE FROM gene_annotations WHERE id=%s", (row['id'],))
                            conn.commit()
                            conn.close()
                            st.success(t('msg_success_delete'))
                            auth.log_action(st.session_state['username'], f"Delete Annotation ID: {row['id']}")
                            st.rerun()
            else:
                st.info("No annotations found.")


# ------------------------------------------
# Main Query Logic
# ------------------------------------------
def query_ui():
    # --- Sidebar Layout (Organized) ---
    with st.sidebar:
        # 1. User Info Section
        st.markdown(f"### {t('sidebar_header_user')}")
        st.success(f"**{st.session_state['username']}** ({st.session_state['user_role']})")
        if st.button(t('logout_btn'), use_container_width=True):
            auth.log_action(st.session_state['username'], "Logout")
            st.session_state['logged_in'] = False
            st.rerun()

        st.markdown("---")

        # 2. Navigation Section
        st.markdown(f"### {t('sidebar_header_nav')}")
        nav_mode = st.radio("Go to", [t('nav_query'), t('nav_admin')], label_visibility="collapsed")

        if nav_mode == t('nav_admin'):
            if st.session_state['user_role'] == 'admin':
                return admin_ui()
            else:
                st.error("Access Denied")
                return

        st.markdown("---")

        # 3. Settings & Data Dictionary (Expanders for cleaner look)
        with st.expander(t('sidebar_header_settings'), expanded=False):
            # Language Switcher
            current_idx = 0 if st.session_state['lang'] == "CN" else 1
            new_lang = st.selectbox(t('lbl_language'), ["CN", "EN"], index=current_idx)
            if new_lang != st.session_state['lang']:
                st.session_state['lang'] = new_lang
                st.rerun()

        with st.expander(t('data_browser'), expanded=False):
            st.markdown(f"**{t('top_genes_list')}**")
            genes = get_real_top_elements("gene")
            st.code("\n".join(genes) if genes else "Loading...")

            st.markdown(f"**{t('top_metas_list')}**")
            metas = get_real_top_elements("metabo")
            st.code("\n".join(metas) if metas else "Loading...")

    # --- Main Content Area ---
    st.title(f"🧬 {t('title')}")

    # Query Mode Selection
    q_modes = {
        t('qmode_gene'): "gene",
        t('qmode_scrna'): "scrna_advanced",
        t('qmode_atac'): "atac_advanced",
        t('qmode_metabo'): "metabo_advanced",
        t('qmode_spatial'): "spatial_advanced"
    }

    col_mode, col_blank = st.columns([1, 2])  # Limit width of selector
    with col_mode:
        selected_mode_label = st.selectbox(t('qmode_label'), list(q_modes.keys()))

    mode_key = q_modes[selected_mode_label]

    # --- Mode 1: Global Gene Search (Federated) ---
    if mode_key == "gene":
        c1, c2 = st.columns([3, 1])
        gene_input = c1.text_input(t('search_label'), value="FOXA1", placeholder=t('input_gene_ph')).strip().upper()
        subtype = c2.selectbox(t('filter_label'), ["All", "TNBC", "HER2_Positive", "Luminal_A", "Luminal_B", "Normal"])

        if not gene_input: return

        # Render 5 Tabs
        tabs = st.tabs([t('tab_scrna'), t('tab_atac'), t('tab_metabo'), t('tab_spatial'), t('tab_imaging')])

        # Tab 1: scRNA
        with tabs[0]:
            # [Display Expert Annotations]
            conn_mysql = auth.get_connection()
            if conn_mysql:
                cursor = conn_mysql.cursor(dictionary=True)
                cursor.execute("SELECT * FROM gene_annotations WHERE gene = %s", (gene_input,))
                notes = cursor.fetchall()
                conn_mysql.close()
                if notes:
                    st.info(f"📋 **Expert Annotation / 专家注释 ({len(notes)})**")
                    for n in notes:
                        st.markdown(f"- {n['note']} *(By: {n['author']} @ {n['created_at']})*")

            # Main Chart
            sql = f"SELECT Subtype, CellType, Avg_Expression FROM Table_Expression WHERE Gene = '{gene_input}'"
            if subtype != "All": sql += f" AND Subtype = '{subtype}'"
            df = run_sqlite_query("scRNA", sql)
            if df is not None and not df.empty:
                st.bar_chart(df, x="CellType", y="Avg_Expression", color="Subtype")
            else:
                st.warning(t('warn_no_data'))

        # Tab 2: ATAC
        with tabs[1]:
            sql = f"SELECT sample, {gene_input} FROM sample_gene_matrix"
            df = run_sqlite_query("ATAC", sql)
            if df is not None and not df.empty:
                df['Subtype'] = df['sample'].apply(get_atac_meta)
                st.info(t('atac_sim_note'))
                df_filter = df[df['Subtype'] == subtype] if subtype != "All" else df
                avg = df_filter.groupby("Subtype")[gene_input].mean().reset_index()
                st.bar_chart(avg, x="Subtype", y=gene_input, color="Subtype")

                st.markdown("---")
                st.markdown(f"#### {t('atac_raw_title')}")
                st.bar_chart(df, x="sample", y=gene_input)
            else:
                st.info(f"Gene {gene_input} not found in ATAC DB.")

        # Tab 3: Metabo
        with tabs[2]:
            df_map = run_sqlite_query("Metabo", f"SELECT * FROM Gene_Metabolite_Map WHERE Gene = '{gene_input}'")
            if df_map is not None and not df_map.empty:
                st.dataframe(df_map)
                metas = [m.replace("'", "''") for m in df_map['Metabolite'].tolist()]
                if metas:
                    m_str = "', '".join(metas)
                    sql = f"SELECT * FROM Metabolite_Expression WHERE Metabolite IN ('{m_str}')"
                    if subtype != "All": sql += f" AND Subtype = '{subtype}'"
                    df_exp = run_sqlite_query("Metabo", sql)
                    if df_exp is not None and not df_exp.empty:
                        st.line_chart(df_exp, x="Subtype", y="Expression_Level", color="Metabolite")
            else:
                st.warning(t('warn_no_data'))

        # Tab 4: Spatial
        with tabs[3]:
            st.info(t('spatial_single_note'))
            sql = f"SELECT SampleID, Region, Avg_Expression FROM Table_SpatialExpression WHERE Gene = '{gene_input}'"
            df = run_sqlite_query("Spatial", sql)
            if df is not None and not df.empty:
                st.bar_chart(df, x="Region", y="Avg_Expression", color="SampleID")
            else:
                st.warning(t('warn_no_data'))

        # Tab 5: Imaging
        with tabs[4]:
            st.info(t('imaging_note'))
            df_anno = run_sqlite_query("Imaging", "SELECT annotation FROM annotations LIMIT 1")
            if df_anno is not None:
                try:
                    data = json.loads(df_anno.iloc[0]['annotation'])
                    fig, ax = plt.subplots(figsize=(6, 6))
                    ax.set_title("Pathology ROI Annotation")
                    ax.set_xlim(5000, 22000);
                    ax.set_ylim(13000, 8000)
                    for region in data.get('positive', []):
                        poly = patches.Polygon(region['vertices'], closed=True, facecolor='#FF4B4B', alpha=0.4,
                                               edgecolor='red')
                        ax.add_patch(poly)
                    st.pyplot(fig)
                except:
                    st.error("JSON Error")
            else:
                st.info("No Imaging Data")

    # --- Mode 2: scRNA Advanced (By CellType) ---
    elif mode_key == "scrna_advanced":
        st.subheader(t('qmode_scrna'))

        # Dynamic inputs
        cell_types = get_distinct_values("scRNA", "Table_Expression", "CellType")
        c1, c2 = st.columns(2)
        ct = c1.selectbox(t('input_celltype'), cell_types)
        top_n = c2.slider(t('input_top_n'), 3, 100, 5)

        if ct:
            # Aggregate to find top expressed genes in this cell type
            sql = f"""
                SELECT Gene, AVG(Avg_Expression) as MeanExpr 
                FROM Table_Expression 
                WHERE CellType = '{ct}' 
                GROUP BY Gene 
                ORDER BY MeanExpr DESC 
                LIMIT {top_n}
            """
            df = run_sqlite_query("scRNA", sql)
            if df is not None and not df.empty:
                st.markdown(f"#### {t('header_top_genes')} : {ct}")
                st.info(t('caption_plot_limit'))
                st.bar_chart(df.head(50), x="Gene", y="MeanExpr")
                st.dataframe(df, use_container_width=True)
            else:
                st.warning(t('warn_no_data'))

    # --- Mode 3: ATAC Advanced (By Sample) ---
    elif mode_key == "atac_advanced":
        st.subheader(t('qmode_atac'))

        # Get sample list from the wide table column names logic is tricky in SQL
        # Just use the 'sample' column
        samples = get_distinct_values("ATAC", "sample_gene_matrix", "sample")
        c1, c2 = st.columns(2)
        sid = c1.selectbox(t('input_sample'), samples)
        top_n = c2.slider(t('input_top_n'), 3, 100, 5)

        if sid:
            # Query the whole row for this sample
            df_sample = run_sqlite_query("ATAC", f"SELECT * FROM sample_gene_matrix WHERE sample = '{sid}'")
            if df_sample is not None and not df_sample.empty:
                # Transpose: Genes are columns, we need them as rows
                # Drop non-numeric 'sample' column
                df_t = df_sample.drop(columns=['sample']).T
                df_t.columns = ['Openness']
                df_t.index.name = 'Gene'
                df_t = df_t.sort_values(by='Openness', ascending=False).head(top_n)

                st.markdown(f"#### Top Open Chromatin Genes: {sid}")
                st.caption(f"Subtype (Simulated): {get_atac_meta(sid)}")
                st.bar_chart(df_t.head(50))
                st.dataframe(df_t, use_container_width=True)

    # --- Mode 4: Metabo Advanced (By Metabolite) ---
    elif mode_key == "metabo_advanced":
        st.subheader(t('qmode_metabo'))

        # Get all metabolites from Map table
        metas = get_distinct_values("Metabo", "Gene_Metabolite_Map", "Metabolite")
        selected_meta = st.selectbox(t('input_metabo'), metas)

        if selected_meta:
            meta_safe = selected_meta.replace("'", "''")

            # 1. Find Associated Genes
            sql_genes = f"SELECT DISTINCT Gene, KEGG, Note FROM Gene_Metabolite_Map WHERE Metabolite = '{meta_safe}'"
            df_genes = run_sqlite_query("Metabo", sql_genes)
            st.markdown("#### 🧬 Associated Genes")
            st.dataframe(df_genes, use_container_width=True)

            # 2. Metabolite Expression Profile
            st.markdown("#### 📊 Expression Profile (by Subtype)")
            sql_exp = f"SELECT Subtype, Expression_Level FROM Metabolite_Expression WHERE Metabolite = '{meta_safe}'"
            df_exp = run_sqlite_query("Metabo", sql_exp)
            if df_exp is not None and not df_exp.empty:
                st.bar_chart(df_exp, x="Subtype", y="Expression_Level")

    # --- Mode 5: Spatial Advanced (By Region) ---
    elif mode_key == "spatial_advanced":
        st.subheader(t('qmode_spatial'))

        regions = get_distinct_values("Spatial", "Table_SpatialExpression", "Region")
        c1, c2 = st.columns(2)
        reg = c1.selectbox(t('input_region'), regions)
        top_n = c2.slider(t('input_top_n'), 3, 100, 5)

        if reg is not None:
            # Query top genes in this region
            sql = f"""
                SELECT Gene, Avg_Expression 
                FROM Table_SpatialExpression 
                WHERE Region = '{reg}' 
                ORDER BY Avg_Expression DESC 
                LIMIT {top_n}
            """
            df = run_sqlite_query("Spatial", sql)
            if df is not None and not df.empty:
                st.markdown(f"#### Top Marker Genes in Region {reg}")
                st.info(t('spatial_single_note'))
                st.bar_chart(df.head(50), x="Gene", y="Avg_Expression")
                st.dataframe(df, use_container_width=True)


if __name__ == "__main__":
    if st.session_state['logged_in']:
        query_ui()
    else:
        login_ui()