import pandas as pd
import tqdm


def get_prefix_df(prefix, con):
    """
    Retrieve AD and RPTR loss table prefixes from a DuckDB database and
    organize them into dataframes with replicate and time information.

    Parameters
    ----------
    prefix : str
        Pool prefix to search for (e.g. "yeast_pool_A"). Currently used
        only for filtering table names.
    con : duckdb.DuckDBPyConnection
        Active DuckDB connection.

    Returns
    -------
    AD_prefixes_df : pandas.DataFrame
        DataFrame containing AD loss tables with columns:
        - id   : table name
        - rep  : replicate number
        - time : timepoint

    RT_prefixes_df : pandas.DataFrame
        DataFrame containing RPTR loss tables with columns:
        - id   : table name
        - rep  : replicate number
        - time : timepoint
    """    
    AD_prefixes = []
    RT_prefixes = []

    # Get list of all tables in the database
    tables = con.sql("SHOW TABLES").fetchall()

    # Identify loss tables belonging to this pool
    for table in tables:
        if prefix in table[0] and "loss" in table[0]:
            print(table[0])
            
            # AD tables
            if "AD" in table[0]:
                print(table[0])
                AD_prefixes.append(table[0])

            # Reporter / RT tables
            elif "RPTR" in table[0]:
                print(table[0])
                RT_prefixes.append(table[0])

            else:
                print(table[0])
                
    # Build AD prefixes dataframe
    AD_prefixes_df = pd.DataFrame({"id": AD_prefixes})
    AD_prefixes_df["id"] = AD_prefixes_df["id"].astype(str)
    
    # Extract replicate and time from table name
    AD_prefixes_df["rep"] = AD_prefixes_df["id"].str.extract("AD_(\d+)")#.astype(int)
    AD_prefixes_df["time"] = AD_prefixes_df["id"].str.extract("AD_\d+_(\d+)")#.astype(int)

    # Sort chronologically
    AD_prefixes_df = AD_prefixes_df.sort_values(by=["time", "rep"])

    # Build RPTR prefixes dataframe
    RT_prefixes_df = pd.DataFrame({"id": RT_prefixes})
    RT_prefixes_df["id"] = RT_prefixes_df["id"].astype(str)

    RT_prefixes_df["rep"] = RT_prefixes_df["id"].str.extract("RPTR_(\d+)")#.astype(int)
    RT_prefixes_df["time"] = RT_prefixes_df["id"].str.extract("RPTR_\d+_(\d+)")#.astype(int)

    RT_prefixes_df = RT_prefixes_df.sort_values(by=["time", "rep"])

    return AD_prefixes_df, RT_prefixes_df


def build_loss_table_dict(prefixes, con, step1_table="step1_AD_AD_BC_RPTR_BC_designed"):
    """
    Build loss tables with an additional Step1 row summarizing how many reads
    and unique constructs appear in the Step1 design table.
    """

    loss_table_df_dict = {}

    for table in tqdm.tqdm(prefixes):

        # Load the loss table
        loss_table_df = con.execute(f"SELECT * FROM {table}").fetchdf()

        # Determine final mapping table name
        last_table = loss_table_df["map"].iloc[-1]
        file_name = table[:-13] + "_" + last_table

        # Inspect schema
        cols = con.execute(f"PRAGMA table_info('{file_name}')").fetchdf()["name"].tolist()

        # Decide how to count reads
        if "count" in cols:
            read_expr = "SUM(m.count)"
        else:
            read_expr = "COUNT(*)"

        # Decide join logic depending on table type
        if "AD" in cols:
            join_condition = "m.AD_BC = s.AD_BC AND m.AD = s.AD"
            distinct_expr = "COUNT(DISTINCT m.AD)"
        else:
            join_condition = "m.RPTR_BC = s.RPTR_BC"
            distinct_expr = "COUNT(*)"

        # Compute Step1 statistics
        step1_counts = con.execute(f"""
            SELECT
                {read_expr} AS total_reads,
                {distinct_expr} AS unique_count
            FROM {file_name} AS m
            JOIN {step1_table} AS s
              ON {join_condition}
        """).fetchone()

        step1_total_reads = step1_counts[0]
        step1_unique = step1_counts[1]

        prev_row = loss_table_df.iloc[-1]

        # Construct Step1 row
        step1_row = {
            "map": "step1",
            "description": "In Step 1 table",
            "unique_count": step1_unique,
            "unique_AD_count": step1_unique,
            "total_reads": step1_total_reads,
            "% of previous step": step1_total_reads / prev_row["total_reads"] * 100,
            "% of total reads": step1_total_reads / loss_table_df.iloc[0]["total_reads"] * 100,
        }

        # Append Step1 row
        loss_table_df = pd.concat(
            [loss_table_df, pd.DataFrame([step1_row])],
            ignore_index=True
        )

        loss_table_df_dict[table] = loss_table_df

    return loss_table_df_dict