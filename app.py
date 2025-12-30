import asyncio
import sys

# ⚠️ Windows 异步事件循环补丁
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
# 1. 中英文翻译
# ==========================================
TRANS = {
    "CN": {
        "title": "BC-MOD 乳腺癌多组学数据库",
        "sidebar_nav": "功能导航",
        "nav_query": "🔍 数据查询",
        "nav_admin": "🛠️ 后台管理",
        "login_title": "🔐 系统登录",
        "login_btn": "登录",
        "logout_btn": "退出登录",
        "intro_title": "系统介绍",
        "intro_text": "本系统采用联邦架构，整合 scRNA-seq, ATAC-seq, Spatial 及 Metabolomics 数据。",
        "search_label": "输入基因 Symbol",
        "search_placeholder": "尝试: FOXA1, ESR1, PKM",
        "filter_label": "亚型过滤",
        "tab_scrna": "🔬 scRNA (单细胞)",
        "tab_atac": "🧬 ATAC (表观)",
        "tab_metabo": "⚗️ Metabo (代谢)",
        "tab_spatial": "🗺️ Spatial (空间)",
        "tab_imaging": "🖼️ Imaging (影像)",
        "data_browser": "📚 数据字典导览",
        "top_genes": "🔥 高表达基因 (Top 10)",
        "top_metas": "🧪 高丰度代谢物 (Top 10)",
        "atac_sim_note": "⚠️ 注：当前 ATAC 数据库缺失临床亚型标注。下图展示基于模拟元数据的分组对比。",
        "atac_raw_title": "2. 原始样本分布 (未过滤)",
        "spatial_single_note": "ℹ️ 提示：当前 Spatial 模块展示标准参考样本 V1 (HER2_Positive)。亚型过滤器不适用。",
        "imaging_note": "💡 说明：展示 AI 辅助识别的肿瘤感兴趣区域 (ROI)。本模块为非结构化数据存储演示。",
        "warn_no_data": "未找到相关数据。",
        "success_login": "验证通过！正在进入系统...",
        "err_login": "用户名或密码错误"
    },
    "EN": {
        "title": "BC-MOD Breast Cancer Multi-Omics DB",
        "sidebar_nav": "Navigation",
        "nav_query": "🔍 Data Query",
        "nav_admin": "🛠️ Admin Panel",
        "login_title": "🔐 System Login",
        "login_btn": "Login",
        "logout_btn": "Logout",
        "intro_title": "Introduction",
        "intro_text": "A federated database integrating scRNA-seq, ATAC-seq, Spatial, and Metabolomics data.",
        "search_label": "Enter Gene Symbol",
        "search_placeholder": "Try: FOXA1, ESR1, PKM",
        "filter_label": "Subtype Filter",
        "tab_scrna": "🔬 scRNA (Single Cell)",
        "tab_atac": "🧬 ATAC (Chromatin)",
        "tab_metabo": "⚗️ Metabo (Metabolism)",
        "tab_spatial": "🗺️ Spatial (Transcriptomics)",
        "tab_imaging": "🖼️ Imaging (Pathology)",
        "data_browser": "📚 Data Browser",
        "top_genes": "🔥 Top Expressed Genes",
        "top_metas": "🧪 Top Abundant Metabolites",
        "atac_sim_note": "⚠️ Note: ATAC subtypes are currently simulated for demonstration purposes.",
        "atac_raw_title": "2. Raw Sample Distribution (Unfiltered)",
        "spatial_single_note": "ℹ️ Note: Spatial module shows Reference Sample V1 (HER2+). Subtype filter implies single sample.",
        "imaging_note": "💡 Note: Displaying AI-identified Tumor ROI. Demo for unstructured data storage.",
        "warn_no_data": "No data found.",
        "success_login": "Login successful! Redirecting...",
        "err_login": "Invalid username or password"
    }
}

# ==========================================
# 2. 系统配置与工具函数
# ==========================================
st.set_page_config(page_title="BC-MOD Database", page_icon="🧬", layout="wide")

# 📂 数据库路径
DB_PATHS = {
    "scRNA": "./dbs/scrna.db",
    "ATAC": "./dbs/atac.db",
    "Metabo": "./dbs/metabolomics.db",
    "Spatial": "./dbs/spatial.db",
    "Imaging": "./dbs/imaging.db"
}

# 初始化 Session State
if 'logged_in' not in st.session_state:
    st.session_state['logged_in'] = False
    st.session_state['user_role'] = None
    st.session_state['username'] = None
if 'lang' not in st.session_state:
    st.session_state['lang'] = "CN"  # 默认中文


def t(key):
    """获取当前语言的文本"""
    return TRANS[st.session_state['lang']].get(key, key)


def run_sqlite_query(db_key, sql):
    """通用查询函数"""
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
    """ATAC 模拟元数据"""
    types = ["TNBC", "HER2_Positive", "Luminal_A", "Luminal_B", "Normal"]
    return types[hash(sample_id) % len(types)]


# 🔥 核心逻辑升级：获取真实的 Top 数据
def get_real_top_elements(mode):
    """获取数据库中真实表达量最高的基因/代谢物"""
    try:
        if mode == "gene":
            # scRNA: 按平均表达量降序取前10
            sql = "SELECT Gene FROM Table_Expression GROUP BY Gene ORDER BY SUM(Avg_Expression) DESC LIMIT 10"
            df = run_sqlite_query("scRNA", sql)
            return df['Gene'].tolist() if df is not None else []
        elif mode == "metabo":
            # Metabo: 按表达水平降序取前10
            sql = "SELECT Metabolite FROM Metabolite_Expression GROUP BY Metabolite ORDER BY SUM(Expression_Level) DESC LIMIT 10"
            df = run_sqlite_query("Metabo", sql)
            return df['Metabolite'].tolist() if df is not None else []
    except:
        return []
    return []


# ==========================================
# 3. 界面模块
# ==========================================
def login_ui():
    st.markdown(f"## 🧬 {t('title')}")
    st.markdown("---")

    # 语言切换 (登录页)
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
        user = st.text_input("Username / 用户名")
        pwd = st.text_input("Password / 密码", type="password")
        if st.button(t('login_btn'), type="primary", use_container_width=True):
            u = auth.check_login(user, pwd)
            if u:
                st.session_state['logged_in'] = True
                st.session_state['user_role'] = u['role']
                st.session_state['username'] = u['username']
                auth.log_action(user, "Login")
                st.success(t('success_login'))
                st.rerun()
            else:
                st.error(t('err_login'))


def admin_ui():
    st.header(t('nav_admin'))
    st.warning(f"Admin: {st.session_state['username']}")

    tab1, tab2 = st.tabs(["User Management", "System Logs"])
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


def query_ui():
    st.title(f"🧬 {t('title')}")

    # --- 侧边栏 ---
    with st.sidebar:
        # 语言切换 (侧边栏)
        st.session_state['lang'] = st.selectbox("Language", ["CN", "EN"],
                                                index=0 if st.session_state['lang'] == "CN" else 1)

        st.markdown("---")
        st.markdown(f"### {t('data_browser')}")

        # 动态获取真实的 Top 数据
        with st.expander(t('top_genes')):
            genes = get_real_top_elements("gene")
            st.write(", ".join(genes) if genes else "Loading...")

        with st.expander(t('top_metas')):
            metas = get_real_top_elements("metabo")
            st.write(", ".join(metas) if metas else "Loading...")

    # --- 搜索区 ---
    with st.container():
        c1, c2 = st.columns([3, 1])
        # 强制设置默认值包含 FOXA1，提示语包含 ESR1, PKM
        gene_input = c1.text_input(
            t('search_label'),
            value="FOXA1",
            placeholder=t('search_placeholder'),
            help="Examples: FOXA1, ESR1, PKM"
        ).strip().upper()

        subtype = c2.selectbox(t('filter_label'), ["All", "TNBC", "HER2_Positive", "Luminal_A", "Luminal_B", "Normal"])

    # --- Tabs ---
    tabs = st.tabs([t('tab_scrna'), t('tab_atac'), t('tab_metabo'), t('tab_spatial'), t('tab_imaging')])

    # 1. scRNA
    with tabs[0]:
        sql = f"SELECT Subtype, CellType, Avg_Expression FROM Table_Expression WHERE Gene = '{gene_input}'"
        if subtype != "All": sql += f" AND Subtype = '{subtype}'"
        df = run_sqlite_query("scRNA", sql)
        if df is not None and not df.empty:
            st.bar_chart(df, x="CellType", y="Avg_Expression", color="Subtype")
        else:
            st.warning(t('warn_no_data'))

    # 2. ATAC
    with tabs[1]:
        sql = f"SELECT sample, {gene_input} FROM sample_gene_matrix"
        df = run_sqlite_query("ATAC", sql)
        if df is not None and not df.empty:
            df['Subtype'] = df['sample'].apply(get_atac_meta)

            # Grouped Chart
            st.caption(t('atac_sim_note'))
            df_filter = df[df['Subtype'] == subtype] if subtype != "All" else df
            avg = df_filter.groupby("Subtype")[gene_input].mean().reset_index()
            st.bar_chart(avg, x="Subtype", y=gene_input, color="Subtype")

            # Raw Chart
            st.markdown("---")
            st.markdown(f"#### {t('atac_raw_title')}")
            st.bar_chart(df, x="sample", y=gene_input)
        else:
            st.info(f"Gene {gene_input} not in ATAC top list.")

    # 3. Metabo
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

    # 4. Spatial
    with tabs[3]:
        # 蓝色提示框 display
        st.info(t('spatial_single_note'))

        sql = f"SELECT SampleID, Region, Avg_Expression FROM Table_SpatialExpression WHERE Gene = '{gene_input}'"
        df = run_sqlite_query("Spatial", sql)
        if df is not None and not df.empty:
            st.bar_chart(df, x="Region", y="Avg_Expression", color="SampleID")
        else:
            st.warning(t('warn_no_data'))

    # 5. Imaging
    with tabs[4]:
        # 蓝色提示框 display (Request 3)
        st.info(t('imaging_note'))

        df_anno = run_sqlite_query("Imaging", "SELECT annotation FROM annotations LIMIT 1")
        if df_anno is not None:
            try:
                data = json.loads(df_anno.iloc[0]['annotation'])
                fig, ax = plt.subplots(figsize=(6, 6))
                ax.set_title("Pathology ROI Annotation")
                ax.set_xlim(5000, 22000)
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


# ==========================================
# 6. 入口
# ==========================================
if __name__ == "__main__":
    if st.session_state['logged_in']:
        with st.sidebar:
            st.success(f"👤 {st.session_state['username']}")
            nav = st.radio(t('sidebar_nav'), [t('nav_query'), t('nav_admin')])
            st.markdown("---")
            if st.button(t('logout_btn')):
                auth.log_action(st.session_state['username'], "Logout")
                st.session_state['logged_in'] = False
                st.rerun()

        if nav == t('nav_admin'):
            if st.session_state['user_role'] == 'admin':
                admin_ui()
            else:
                st.error("Access Denied")
        else:
            query_ui()
    else:
        login_ui()