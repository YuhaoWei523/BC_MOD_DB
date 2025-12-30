import mysql.connector
import hashlib

# --- 配置 ---
# 这是连接 MySQL 服务器的密码 (安装时设的那个)
MYSQL_ROOT_PASSWORD = '123456'

# 这是你想要设置的系统登录账号密码
APP_USERS = [
    ('admin', 'admin123456', 'admin'),  # 用户名, 密码, 权限
    ('guest', 'guest123456', 'guest')
]


def hash_password(password):
    """跟 auth_manager.py 保持一致的加密逻辑"""
    return hashlib.sha256(password.encode()).hexdigest()


def init_database():
    print("🚀 开始初始化数据库...")

    # 1. 连接 MySQL Server (不指定数据库，因为还没建)
    try:
        conn = mysql.connector.connect(
            host='localhost',
            user='root',
            password=MYSQL_ROOT_PASSWORD
        )
        cursor = conn.cursor()
    except Exception as e:
        print(f"❌ 连接 MySQL 失败: {e}")
        print("请检查你的 MySQL 服务是否启动，以及 ROOT 密码是否正确。")
        return

    # 2. 创建库和表
    try:
        # 建库
        cursor.execute("CREATE DATABASE IF NOT EXISTS bc_mod_admin;")
        cursor.execute("USE bc_mod_admin;")

        # 建用户表
        print("📂 创建用户表...")
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS users (
                id INT AUTO_INCREMENT PRIMARY KEY,
                username VARCHAR(50) NOT NULL UNIQUE,
                password_hash VARCHAR(255) NOT NULL,
                role ENUM('admin', 'guest') DEFAULT 'guest',
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            );
        """)

        # 建日志表
        print("📂 创建日志表...")
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS system_logs (
                log_id INT AUTO_INCREMENT PRIMARY KEY,
                username VARCHAR(50),
                action VARCHAR(255),
                log_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            );
        """)

        # 3. 插入用户 (带加密)
        print("🔐 正在添加用户...")

        # 先清空旧数据(防止重复运行报错)
        cursor.execute("TRUNCATE TABLE users;")

        for user, pwd, role in APP_USERS:
            pwd_hash = hash_password(pwd)
            # 插入数据
            query = "INSERT INTO users (username, password_hash, role) VALUES (%s, %s, %s)"
            cursor.execute(query, (user, pwd_hash, role))
            print(f"   - 用户 [{user}] 创建成功 (密码: {pwd})")

        conn.commit()
        print("\n✅✅✅ 数据库初始化完美完成！")
        print("现在你可以直接运行 streamlit run app.py 了")

    except Exception as e:
        print(f"❌ 初始化过程中出错: {e}")
    finally:
        conn.close()


if __name__ == "__main__":
    init_database()