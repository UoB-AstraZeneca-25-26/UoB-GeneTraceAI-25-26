#!/usr/bin/env python3
from pathlib import Path
import duckdb

base_dir = Path(__file__).resolve().parent
db_path = base_dir / "db" / "celllineselector.duckdb"
out_path = base_dir / "content.txt"

con = duckdb.connect(str(db_path))

print(f"Database: {db_path}")
print(f"Output: {out_path}\n")

with out_path.open("w", encoding="utf-8") as f:
    tables = con.execute("SHOW TABLES").fetchall()
    for (table_name,) in tables:
        f.write(f"TABLE: {table_name}\n")
        print(f"TABLE: {table_name}")

        columns = con.execute(f"DESCRIBE {table_name}").fetchall()
        col_names = [col[0] for col in columns]
        f.write(f"COLUMNS: {', '.join(col_names)}\n")
        print(f"COLUMNS: {', '.join(col_names)}")

        try:
            rows = con.execute(f"SELECT * FROM \"{table_name}\" LIMIT 10").fetchall()
            if not rows:
                f.write("HEAD: (no rows)\n")
                print("HEAD: (no rows)")
            else:
                for row in rows:
                    f.write(str(row) + "\n")
                    print(row)
        except Exception as exc:
            f.write(f"HEAD: ERROR - {exc}\n")
            print(f"HEAD: ERROR - {exc}")

        f.write("\n---\n\n")
        print()

print(f"\nSaved first 10 rows for each table to: {out_path}")
con.close()
