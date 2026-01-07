    def check_designed(self, map_df):
        """
        Marks sequences as 'Designed' based on a provided design file (optional).
        If no design file is provided, returns the input DataFrame unchanged.

        Args:
            map_df (dask.DataFrame): DataFrame with an 'AD' column of sequences to check.

        Returns:
            dask.DataFrame: DataFrame with an additional 'Designed' column (1 if in design file, else 0).

        Example:
            >>> mapper = BarcodeMapper(...initialize...)
            >>> df_checked = mapper.check_designed(mapper.seq_df)
        """
        if self.design_file_path is None:
            map_df["Designed"] = 1
            return map_df
            
        design_file = dd.read_csv(self.design_file_path, header=None, names=["AD"], dtype=str)
        design_file["Designed"] = 1

        map_df["AD"] = map_df["AD"].str.strip()
        design_file["AD"] = design_file["AD"].str.strip()

        merged = dd.merge(map_df, design_file, on="AD", how="left")
        merged["Designed"] = merged["Designed"].fillna(0)
        
        return merged

    def create_map(self):
        """
        Processes the sequence file, extracts barcodes, and maps them to designed sequences.

        Returns:
            dask.DataFrame: DataFrame with barcodes mapped to design file.

        Example:
            >>> mapper = BarcodeMapper(...initialize...)
            >>> mapped_df = mapper.create_map()
        """
        dfs = []

        # 1. Load seq files as dask dataframe, translate =and use txt intermediate if needed
        all_seq_df = load_and_shorten_files(self.seq_file, self.reverse_complement)

        # 2. Extract all barcodes
        all_seq_df = add_multiple_barcodes(self.bc_objects, all_seq_df)

        if self.umi_length > 0:
            all_seq_df = add_umi(all_seq_df, self.umi_length)
        
        # 3. Apply design check
        self.mapped_df = self.check_designed(all_seq_df).drop(columns="sequence")
    
        return self.mapped_df

    

