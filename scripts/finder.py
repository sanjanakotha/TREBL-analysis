import dask.dataframe as dd
import pandas as pd

def add_barcode(seq_df, bc_name, preceder, post, bc_length=120):
    """
    Extracts a barcode from sequences between a preceding and following pattern.

    Args:
        seq_df (pd.DataFrame or dask.DataFrame): DataFrame containing a 'sequence' column.
        bc_name (str): Name of the barcode to create.
        preceder (str): Sequence preceding the barcode.
        post (str): Sequence following the barcode.
        bc_length (int, optional): Maximum length of barcode to extract. Defaults to 120.

    Returns:
        pd.DataFrame or dask.DataFrame: Updated DataFrame with barcode and quality columns.

    Example:
        >>> df = pd.DataFrame({"sequence": ["AAAXXXCCC", "AAAYYYCCC"]})
        >>> BarcodeMapper.add_barcode(df, "AD", "AAA", "CCC", 3)
           sequence   AD  AD_qual
        0  AAAXXXCCC  XXX     True
        1  AAAYYYCCC  YYY     True
    """
    regex = f"{preceder}(.*){post}"
    subseq_series = seq_df['sequence'].str.extract(regex)[0].str.slice(0, bc_length)
    seq_df[bc_name] = subseq_series
    seq_df[bc_name + "_qual"] = True #seq_df[bc_name].str.len() == bc_length
    seq_df[bc_name + "_qual"] = seq_df[bc_name + "_qual"].fillna(False)
    return seq_df
    
def add_multiple_barcodes(bc_names, preceders, posts, lengths, seq_df):
    """
    Extracts all defined barcodes from sequences in a DataFrame.

    Args:
        seq_df (pd.DataFrame or dask.DataFrame): DataFrame containing sequences.

    Returns:
        pd.DataFrame or dask.DataFrame: Updated DataFrame with all barcodes added.

    Example:
        >>> df = pd.DataFrame({"sequence": ["AAAXXXCCC", "AAAYYYCCC"]})
        >>> mapper = BarcodeMapper(...initialize...)
        >>> mapper.add_multiple_barcodes(df)
    """
    df = seq_df.copy()
    for bc_name, preceder, post, length in zip(bc_names, preceders, posts, lengths):
        df = add_barcode(df, bc_name, preceder, post, length)
    return df