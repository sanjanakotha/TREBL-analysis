import dask.dataframe as dd
import pandas as pd

class Barcode:
    def __init__(
        self,
        name,
        preceder,
        post, 
        length):
    
        self.name = name
        self.preceder = preceder
        self.post = post
        self.length = length
        
def add_barcode(seq_df, bc_object):
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
    regex = f"{bc_object.preceder}(.*){bc_object.post}"
    subseq_series = seq_df['sequence'].str.extract(regex)[0].str.slice(0, bc_object.length)
    seq_df[bc_object.name] = subseq_series
    seq_df[bc_object.name + "_qual"] = seq_df[bc_object.name].str.len() == bc_object.length
    seq_df[bc_object.name + "_qual"] = seq_df[bc_object.name + "_qual"].fillna(False)
    return seq_df

def add_umi(seq_df, umi_length=12):
    """
    Extracts the first `umi_length` bases from the 'sequence' column as the UMI.

    Args:
        seq_df (pd.DataFrame or dask.DataFrame): DataFrame containing a 'sequence' column.
        umi_length (int): Length of the UMI to extract. Default is 12.

    Returns:
        pd.DataFrame or dask.DataFrame: Updated DataFrame with 'UMI' and 'UMI_qual' columns.
    """
    seq_df['UMI'] = seq_df['sequence'].str.slice(-12,)
    return seq_df
    
def add_multiple_barcodes(bc_objects, seq_df):
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
    for bc_object in bc_objects:
        df = add_barcode(df, bc_object)
    return df