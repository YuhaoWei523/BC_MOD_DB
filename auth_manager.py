import mysql.connector
import hashlib
import streamlit as st
import pandas as pd
from datetime import datetime

# ==========================================
# 1. MySQL 数据库配置
# ==========================================
DB_CONFIG = {
    'host': 'localhost',
    'user': 'root',
    'password': '123456',
    'database': 'bc_mod_admin'
}


# ==========================================
# 2. 核心功能函数
# ==========================================
def get_connection():
    """获取数据库连接"""
    try:
        return mysql.connector.connect(**DB_CONFIG)
    except Exception as e:
        st.error(f"❌ MySQL 连接失败: {e}")
        return None


def hash_password(password):
    """密码加密 (SHA-256)"""
    return hashlib.sha256(password.encode()).hexdigest()


def check_login(username, password):
    """验证登录"""
    conn = get_connection()
    if not conn: return None

    try:
        cursor = conn.cursor(dictionary=True)
        pwd_hash = hash_password(password)

        # 防止 SQL 注入查询
        query = "SELECT * FROM users WHERE username = %s AND password_hash = %s"
        cursor.execute(query, (username, pwd_hash))
        user = cursor.fetchone()
        return user
    except Exception as e:
        st.error(f"登录查询出错: {e}")
        return None
    finally:
        if conn.is_connected(): conn.close()


def create_user(username, password, role='guest'):
    """创建新用户 (Admin功能)"""
    conn = get_connection()
    if not conn: return False

    try:
        cursor = conn.cursor()
        pwd_hash = hash_password(password)
        query = "INSERT INTO users (username, password_hash, role) VALUES (%s, %s, %s)"
        cursor.execute(query, (username, pwd_hash, role))
        conn.commit()
        return True
    except Exception as e:
        # st.error(f"创建用户失败: {e}")
        return False
    finally:
        conn.close()


def log_action(username, action):
    """记录操作日志 (触发器模拟)"""
    conn = get_connection()
    if not conn: return

    try:
        cursor = conn.cursor()
        query = "INSERT INTO system_logs (username, action) VALUES (%s, %s)"
        cursor.execute(query, (username, action))
        conn.commit()
    except Exception as e:
        print(f"日志记录失败: {e}")
    finally:
        conn.close()


def get_system_logs(limit=50):
    """获取系统日志"""
    conn = get_connection()
    if not conn: return pd.DataFrame()
    try:
        df = pd.read_sql(f"SELECT * FROM system_logs ORDER BY log_time DESC LIMIT {limit}", conn)
        return df
    finally:
        conn.close()