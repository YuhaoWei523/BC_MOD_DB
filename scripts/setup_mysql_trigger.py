import mysql.connector

DB_CONFIG = {
    'host': 'localhost',
    'user': 'root',
    'password': '123456',
    'database': 'bc_mod_admin'
}

SQL_TRIGGER = """
CREATE TRIGGER IF NOT EXISTS log_new_user
AFTER INSERT ON users
FOR EACH ROW
BEGIN
    INSERT INTO system_logs (username, action, log_time)
    VALUES (NEW.username, CONCAT('System Audit: New user registered - Role: ', NEW.role), NOW());
END;
"""


def install_trigger():
    try:
        conn = mysql.connector.connect(**DB_CONFIG)
        cursor = conn.cursor()

        print("🔧 正在安装 MySQL 触发器...")
        # 某些旧版驱动不支持 DELIMITER，我们直接执行单条创建语句
        # 注意：这通常需要管理员权限
        try:
            cursor.execute("DROP TRIGGER IF EXISTS log_new_user")  # 先清理旧的
            cursor.execute(SQL_TRIGGER)
            print("✅ 触发器 'log_new_user' 安装成功！")
            print("   现在每当你在后台创建新用户，system_logs 表会自动增加一条审计记录。")
        except Exception as e:
            print(f"❌ 安装失败 (SQL 错误): {e}")

        conn.commit()
        conn.close()
    except Exception as e:
        print(f"❌ 连接失败: {e}")


if __name__ == "__main__":
    install_trigger()