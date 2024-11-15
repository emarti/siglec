import numpy as np
import pandas as pd

def weighted_mean(x, sigma, axis=0):
    """
    Calculate weighted mean given known errors. Errors are given as standard
    deviations. If the errors are estimates, representing the sample standard
    deviation, then this is an approximation.

    The behavior of numpy arrays and pandas dataframes is different when some
    entries are NaN or infinite. Pandas is "smarter" and will calculate the
    weighted mean and number of degrees of freedoms for the non-NaN entries,
    while numpy will return nan as the result. See examples.

    Paramaters
    ----------
    x : numpy.array or pandas.DataFrame
        Array or dataframe of measured values.

    sigma : numpy.array or pandas.DataFrame
        Array or dataframe of standard deviations.

    axis : int
        Axis to sum over. Default is zero.

    Results
    -------
    x_mean : numpy.array or pandas.DataFrame
        Weighted mean, given by sum( x / sigma**2 ) / sum( 1 / sigma**2 )

    sigma_mean : numpy.array or pandas.DataFrame
        Weighted errors, given by sqrt( 1/(sum(1 / sigma**2) )

    chi2 : numpy.array or pandas.DataFrame

    dof : int
        Number of degrees of freedom. This is equal to the number of samples
        minus one.

    Examples
    --------
    >>> import numpy as np
    >>> import pandas as pd
    >>> import noisymice as nm
    >>> nm.statistics.weighted_mean([1, 2, 3], [1, 1, 0.5])
    (2.5, 0.408248290463863, 3.5, 2)
    >>> nm.statistics.weighted_mean(np.array([[1, 2, 3], [1, 2, np.nan]]),
                                    np.array([[1, 1, .5], [1, 1, np.nan]]),
                                    axis=1)
    (array([2.5, nan]), array([0.40824829,        nan]), array([3.5, nan]), array([2, 1]))
    >>> val = pd.DataFrame({'a': [1., 2., 3.], 'b': [1., 2., np.nan]})
    >>> err = pd.DataFrame({'a': [1., 1., .5], 'b': [1., 1., np.nan]})
    >>> mean, err, chi2, dof = nm.statistics.weighted_mean(val, err)
    >>> print(mean.to_string())
    a    2.5
    b    1.5
    >>> print(err.to_string())
    a    0.408248
    b    0.707107
    >>> print(chi2.to_string())
    a    3.5
    b    0.5
    >>> print(dof.to_string())
    a    2
    b    1
    """
    if isinstance(x, list):
        x = np.array(x)
    if isinstance(sigma, list):
        sigma = np.array(sigma)

    x_mean = (x * (sigma**-2)).sum(axis=axis) \
           / (sigma**-2).sum(axis=axis)
    sigma_mean =  ( 1/((sigma**-2).sum(axis=axis)) )**(0.5)
    try:
        chi2 =  ( (x - x_mean)**2 / sigma**2 ).sum(axis=axis)
    except ValueError:
        chi2 =  ( (x - np.expand_dims(x_mean, axis=axis))**2 / sigma**2 ).sum(axis=axis)

    dof = (np.isfinite(x) & np.isfinite(sigma)).sum(axis=axis) - 1

    return x_mean, sigma_mean, chi2, dof

def weighted_mean_from_panda_records(df,
                                     val_name,
                                     std_name,
                                     index_name='mouse',
                                     column_name='gene',
                                     fill_value=None):
    """
    Calculate weighted mean,  given known errors. This expects a DataFrame
    in records orientation, i.e. of the format

    mouse1 gene1 value1 error1
    mouse2 gene2 value2 error2
    ...

    Paramaters
    ----------
    df : pandas.DataFrame
        Dataframe in record orientation. Each row should include a gene, value,
        and error.

    val_column : str
        Name of column in df containing the values.

    err_column : str
        Name of column in df containing the errors.
        Errors are equal to the standard deviation.

    fill_value : float or None
        Value to fill missing entires in dataframe, default is None (fills with nan).

    Results
    -------
    x_mean : numpy.array or pandas.DataFrame
        Weighted mean, given by sum( x / sigma**2 ) / sum( 1 / sigma**2 )

    sigma_mean : numpy.array or pandas.DataFrame
        Weighted errors, given by sqrt( 1/(sum(1 / sigma**2) )

    chi2 : numpy.array or pandas.DataFrame

    dof : int
        Number of degrees of freedom. This is equal to the number of samples
        minus one.

    Examples
    --------

    """
    val = pd.pivot_table(df, values=val_name, index=index_name, columns=column_name, fill_value=fill_value)
    err = pd.pivot_table(df, values=std_name, index=index_name, columns=column_name, fill_value=fill_value)

    mean, weighted_error, chi2, dof = weighted_mean(val, err)
    weighted_error_corrected = weighted_error * np.sqrt( chi2/dof )

    results = pd.DataFrame({'weighted_mean': mean,
                            'weighted_error': weighted_error,
                            'weighted_error_corrected': weighted_error_corrected,
                            'chi2': chi2,
                            'dof': dof})

    return results

def Benjamini_Hochberg(p, alpha=0.1):
    """
    Calculate multiple hypothesis testing cutoff with the Benjamini-Hochberg
    procedure. Returns cutoff and index of p values that successfully
    passed the cutoff.

    Paramaters
    ----------
    p : list, numpy.ndarray, or pandas.Series
        List of p values. May be unsorted. 

    alpha : float
        The desired false discovery rate.

    Results
    -------

    p_cutoff : float
        p-value cutoff for significance

    index_pass : numpy.ndarray or pandas.Series
        Index of allegedly significant p-values 

    Examples
    --------
    >>> noisymice.statistics.Benjamini_Hochberg([0.1, 0.5, 1.0, 0.01])
    (0.04, array([False, False, False, True]))
    
    """
    
    p_sort = np.sort(p)
    n = len(p_sort)
    k = 1 + np.arange(n)
    p_k = alpha * k / n
    if np.any(p_sort < p_k):
        p_cutoff = p_k[np.argmax(p_sort[p_sort < p_k])]
    else:
        p_cutoff = p_k[0]

    index_pass = p < p_cutoff

    return p_cutoff, index_pass
