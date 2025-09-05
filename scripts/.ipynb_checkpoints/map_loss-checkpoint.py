import duckdb

class LossTable():
    def __init__(self, db_path, cols, reads_threshold, column_pairs):
        self.db_path = db_path
        self.con = duckdb.connect(self.db_path)
        self.cols = cols
        
    def map2_loss(self):
        # What proportion of barcodes/tiles are quality checked to be okay?          
        union_queries = [
            f"""
            SELECT 
                '{bc}' AS barcode,
                SUM({bc}_qual) AS count_ones,
                COUNT(*) AS total_rows,
                SUM({bc}_qual) * 1.0 / COUNT(*) AS proportion_ones
            FROM map1
            """
            for bc in self.cols
        ]
        
        query = " UNION ALL ".join(union_queries)
        
        bc_qc_loss_table = self.con.execute(query).fetchdf()
        return bc_qc_loss_table 